#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_euk.py
----------------
Fases 1 y 2 sobre los bins EUCARIOTAS (anexo de la tesis de la Estela de
Raimondi). Repite la logica de las fases bacterianas con herramientas para
eucariotas. Solo biblioteca estandar.

Subcomandos:
  seleccion     Copia los bins con mayoria eucariota segun Tiara (Fase 1).
  fase1         Cruza Tiara + EukCC + BUSCO, aplica >70 % / <5 % y escribe el
                reporte. Copia los seleccionados.
  refs          Tabla de genomas de referencia (accesion, especie, linaje NCBI).
  its-ventanas  Ventanas alrededor del ADNr detectado por barrnap, para ITSx.
  busco-genes   Un FASTA por gen BUSCO de copia unica compartido (arbol).
  concat        Concatena los alineamientos recortados (matriz del arbol).
  genomeinfo    genomeInfo.csv para dRep desde EukCC.
  fase2         Taxonomia en cascada + linajes dRep + prevalencia + reporte.

Taxonomia en cascada (Fase 2):
  1. skani contra genomas de referencia de hongos de NCBI:
     ANI >= 95 % y AF (fraccion del bin alineada) >= 65 %  -> especie.
  2. ITS (ITSx) contra UNITE, con los umbrales de identidad de Vu et al. 2019:
     >= 98.41 especie (hipotesis), >= 94.3 genero, >= 88.5 familia,
     >= 81.2 orden; cobertura de la consulta >= 90 %.
  3. Vecino mas cercano de skani (ANI < 95 %) -> genero, con advertencia.
  4. Linaje de la colocacion de EukCC (orientativo).
"""

import argparse
import csv
import glob
import math
import os
import re
import shutil
import sys
from datetime import datetime

RANGOS = ["dominio", "reino", "filo", "clase", "orden", "familia", "genero", "especie"]
# Rango NCBI -> columna. NCBI renombro 'superkingdom' a 'domain' en 2025.
NCBI_A_RANGO = {"superkingdom": "dominio", "domain": "dominio", "kingdom": "reino",
                "phylum": "filo", "class": "clase", "order": "orden",
                "family": "familia", "genus": "genero", "species": "especie"}
UNITE_PREFIJOS = {"k": "reino", "p": "filo", "c": "clase", "o": "orden",
                  "f": "familia", "g": "genero", "s": "especie"}
RE_ACC = re.compile(r"GC[AF]_\d+\.\d+")
RE_BUSCO = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


# Utilidades
def bin_name(path):
    """bin-1-55.fasta -> bin-1-55"""
    b = os.path.basename(path)
    for ext in (".gz", ".fasta", ".fa", ".fna"):
        if b.endswith(ext):
            b = b[: -len(ext)]
    return b


def sample_of(name):
    return name.split("-")[-1] if "-" in name else "NA"


def leer_tsv(path, delim="\t"):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def escribir_tsv(path, filas, cols):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for f in filas:
            w.writerow({c: ("" if f.get(c) is None else f.get(c)) for c in cols})


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def leer_fasta(path, solo=None):
    """Itera (id, secuencia). Si 'solo' es un conjunto, devuelve solo esos ids."""
    sid, partes = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None and (solo is None or sid in solo):
                    yield sid, "".join(partes)
                sid, partes = line[1:].split()[0], []
            else:
                partes.append(line.strip())
    if sid is not None and (solo is None or sid in solo):
        yield sid, "".join(partes)


def stats_fasta(path):
    """n contigs, longitud, N50 y GC (%) de un genoma."""
    largos, gc, total = [], 0, 0
    for _, s in leer_fasta(path):
        largos.append(len(s))
        su = s.upper()
        gc += su.count("G") + su.count("C")
        total += len(s) - su.count("N")
    largos.sort(reverse=True)
    acum, n50 = 0, 0
    for l in largos:
        acum += l
        if acum >= sum(largos) / 2:
            n50 = l
            break
    return {"n_contigs": len(largos), "longitud_bp": sum(largos), "n50": n50,
            "gc": round(100.0 * gc / total, 2) if total else ""}


# Taxonomia NCBI (taxdump)
class Taxdump:
    def __init__(self, d):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        with open(os.path.join(d, "nodes.dmp"), encoding="utf-8") as fh:
            for line in fh:
                p = line.split("\t|\t")
                self.parent[p[0].strip()] = p[1].strip()
                self.rank[p[0].strip()] = p[2].strip()
        with open(os.path.join(d, "names.dmp"), encoding="utf-8") as fh:
            for line in fh:
                p = line.split("\t|\t")
                if p[3].strip().rstrip("|").strip() == "scientific name":
                    self.name[p[0].strip()] = p[1].strip()
        mpath = os.path.join(d, "merged.dmp")
        if os.path.isfile(mpath):
            with open(mpath, encoding="utf-8") as fh:
                for line in fh:
                    p = line.split("\t|\t")
                    self.merged[p[0].strip()] = p[1].strip().rstrip("|").strip()

    def linaje(self, taxid):
        """taxid -> {rango: nombre} para los RANGOS."""
        t = self.merged.get(str(taxid), str(taxid))
        out, vistos = {}, set()
        while t and t not in vistos and t in self.parent:
            vistos.add(t)
            r = NCBI_A_RANGO.get(self.rank.get(t, ""))
            if r and r not in out:
                out[r] = self.name.get(t, "")
            if t == "1":
                break
            t = self.parent[t]
        return out


def linaje_eukcc(ncbi_lng, tax):
    """EukCC escribe el linaje de su colocacion como taxids unidos por '-'.
    Se toma el taxid con el linaje mas profundo."""
    mejor = {}
    for t in (ncbi_lng or "").split("-"):
        if t.strip().isdigit():
            lin = tax.linaje(t.strip())
            if len(lin) > len(mejor):
                mejor = lin
    return mejor


def texto_linaje(lin):
    return ";".join(lin.get(r, "") for r in RANGOS[1:]).strip(";") if lin else ""


# Subcomando: seleccion
def cmd_seleccion(args):
    filas = leer_tsv(args.tiara_summary)
    os.makedirs(args.dest, exist_ok=True)
    n, faltan = 0, []
    for f in filas:
        if f.get("majority") != "euk":
            continue
        src = os.path.join(args.bins_dir, f["bin"] + "." + args.bin_ext)
        if not os.path.isfile(src):
            faltan.append(f["bin"])
            continue
        shutil.copyfile(src, os.path.join(args.dest, os.path.basename(src)))
        n += 1
    print("    Bins con mayoria eucariota (Tiara): %d copiados a %s" % (n, args.dest))
    if faltan:
        sys.exit("ERROR: no se encontraron en %s: %s" % (args.bins_dir, ", ".join(faltan)))


# Subcomando: fase1
def leer_eukcc(eukcc_dir, name):
    """-> (completitud, contaminacion, ncbi_lng, estado)."""
    p = os.path.join(eukcc_dir, name, "eukcc.tsv")
    if not os.path.isfile(p):
        return None, None, "", "sin salida de EukCC"
    filas = leer_tsv(p)
    if not filas:
        return None, None, "", "EukCC sin filas"
    r = filas[0]
    comp, cont = num(r.get("completeness")), num(r.get("contamination"))
    if comp is None or cont is None:
        return None, None, r.get("ncbi_lng", ""), "EukCC sin estimacion (sin set de marcadores)"
    return comp, cont, r.get("ncbi_lng", ""), "ok"


def leer_busco(busco_dir, linaje, name):
    """Lee la linea C:..%[S:..,D:..],F:..,M:..,n:.. del short_summary."""
    pat = os.path.join(busco_dir, linaje, name, "short_summary*.txt")
    for p in sorted(glob.glob(pat)):
        with open(p, encoding="utf-8") as fh:
            m = RE_BUSCO.search(fh.read())
        if m:
            c, s, d, f, mi, n = m.groups()
            return {"C": float(c), "S": float(s), "D": float(d), "F": float(f),
                    "M": float(mi), "n": int(n)}
    return None


def cmd_fase1(args):
    tax = Taxdump(args.taxdump)
    tiara = {f["bin"]: f for f in leer_tsv(args.tiara_summary)}
    genomas = sorted(glob.glob(os.path.join(args.bins_dir, "*." + args.bin_ext)))
    linajes_busco = args.busco_lineages.split()
    os.makedirs(args.selected_dir, exist_ok=True)
    filas = []
    for g in genomas:
        name = bin_name(g)
        t = tiara.get(name, {})
        st = stats_fasta(g)
        comp, cont, lng, estado = leer_eukcc(args.eukcc_dir, name)
        lin = linaje_eukcc(lng, tax)
        fila = {"bin": name, "muestra": sample_of(name), **st,
                "pct_euk_tiara": t.get("pct_euk", ""),
                "pct_unknown_tiara": t.get("pct_unknown", ""),
                "eukcc_completitud": "" if comp is None else comp,
                "eukcc_contaminacion": "" if cont is None else cont,
                "eukcc_estado": estado,
                "eukcc_ncbi_lng": lng,
                "eukcc_linaje": texto_linaje(lin)}
        for r in RANGOS:
            fila["eukcc_" + r] = lin.get(r, "")
        for lb in linajes_busco:
            b = leer_busco(args.busco_dir, lb, name)
            clave = lb.split("_odb")[0]
            for k in ("C", "S", "D", "F", "M", "n"):
                fila["busco_%s_%s" % (clave, k)] = "" if b is None else b[k]
        if comp is None:
            fila["decision_final"], fila["motivo"] = "RECHAZADO", estado
        else:
            motivos = []
            if not comp > args.min_completeness:
                motivos.append("completitud %.2f<=%g" % (comp, args.min_completeness))
            if not cont < args.max_contamination:
                motivos.append("contaminacion %.2f>=%g" % (cont, args.max_contamination))
            fila["decision_final"] = "RECHAZADO" if motivos else "SELECCIONADO"
            fila["motivo"] = "; ".join(motivos) if motivos else "pasa"
        if fila["decision_final"] == "SELECCIONADO":
            shutil.copyfile(g, os.path.join(args.selected_dir, os.path.basename(g)))
        filas.append(fila)

    cols = (["bin", "muestra", "n_contigs", "longitud_bp", "n50", "gc",
             "pct_euk_tiara", "pct_unknown_tiara", "eukcc_completitud",
             "eukcc_contaminacion", "eukcc_estado"]
            + ["busco_%s_%s" % (lb.split("_odb")[0], k) for lb in linajes_busco
               for k in ("C", "S", "D", "F", "M", "n")]
            + ["eukcc_" + r for r in RANGOS] + ["eukcc_linaje", "eukcc_ncbi_lng",
                                                 "decision_final", "motivo"])
    escribir_tsv(args.out_tsv, filas, cols)
    reporte_fase1(args, filas, linajes_busco)
    sel = [f for f in filas if f["decision_final"] == "SELECCIONADO"]
    print("    Bins eucariotas evaluados: %d | seleccionados: %d" % (len(filas), len(sel)))
    print("    Tabla: %s\n    Reporte: %s" % (args.out_tsv, args.out_report))


def reporte_fase1(args, filas, linajes_busco):
    sel = [f for f in filas if f["decision_final"] == "SELECCIONADO"]
    sin = [f for f in filas if f["eukcc_completitud"] == ""]
    est = [f for f in filas if f["eukcc_completitud"] != ""]
    solo_comp = [f for f in est if "completitud" in f["motivo"] and "contaminacion" not in f["motivo"]]
    solo_cont = [f for f in est if "contaminacion" in f["motivo"] and "completitud" not in f["motivo"]]
    ambos = [f for f in est if "completitud" in f["motivo"] and "contaminacion" in f["motivo"]]
    lf = linajes_busco[0].split("_odb")[0]
    with open(args.out_report, "w", encoding="utf-8") as fh:
        fh.write("# Fase 1 eucariota — Control de calidad de los bins eucariotas\n\n")
        fh.write("**Job SLURM:** `%s` — **Fecha:** %s\n\n" % (args.job_id, datetime.now().strftime("%Y-%m-%d %H:%M")))
        fh.write("Entrada: bins con mayoria de longitud eucariota segun Tiara (Fase 1). "
                 "Calidad con EukCC (completitud > %g %% y contaminacion < %g %%). "
                 "BUSCO como medida complementaria (no decide).\n\n"
                 % (args.min_completeness, args.max_contamination))
        fh.write("| Metrica | Valor |\n|---|---:|\n")
        fh.write("| Bins eucariotas (Tiara) | %d |\n" % len(filas))
        fh.write("| Sin estimacion de EukCC | %d |\n" % len(sin))
        fh.write("| Rechazo solo por completitud | %d |\n" % len(solo_comp))
        fh.write("| Rechazo solo por contaminacion | %d |\n" % len(solo_cont))
        fh.write("| Rechazo por ambos | %d |\n" % len(ambos))
        fh.write("| **Seleccionados** | **%d** |\n\n" % len(sel))
        fh.write("## Detalle por bin\n\n")
        fh.write("| Bin | Mb | Contigs | N50 | EukCC compl. | EukCC cont. | BUSCO %s C/D | Linaje EukCC | Decision |\n" % lf)
        fh.write("|---|---:|---:|---:|---:|---:|---|---|---|\n")
        for f in sorted(filas, key=lambda x: (x["decision_final"] != "SELECCIONADO", x["bin"])):
            c, d = f.get("busco_%s_C" % lf, ""), f.get("busco_%s_D" % lf, "")
            fh.write("| %s | %.2f | %d | %d | %s | %s | %s | %s | %s |\n" % (
                f["bin"], f["longitud_bp"] / 1e6, f["n_contigs"], f["n50"],
                f["eukcc_completitud"], f["eukcc_contaminacion"],
                ("%s/%s" % (c, d)) if c != "" else "-",
                f["eukcc_linaje"] or "-",
                f["decision_final"] if f["decision_final"] == "SELECCIONADO" else f["motivo"]))
        fh.write("\nTabla completa: `%s`\n" % os.path.basename(args.out_tsv))


# Subcomando: refs
def cmd_refs(args):
    tax = Taxdump(args.taxdump)
    rutas = {}
    with open(args.fasta_list, encoding="utf-8") as fh:
        for line in fh:
            p = line.strip()
            m = RE_ACC.search(os.path.basename(p))
            if p and m:
                rutas[m.group(0)] = p
    filas, sin_fasta = [], 0
    for r in leer_tsv(args.datasets_tsv):
        acc = r.get("Assembly Accession") or r.get("accession") or ""
        taxid = r.get("Organism Taxonomic ID") or r.get("organism-tax-id") or ""
        org = r.get("Organism Name") or r.get("organism-name") or ""
        if acc not in rutas:
            sin_fasta += 1
            continue
        lin = tax.linaje(taxid)
        filas.append({"accesion": acc, "organismo": org, "taxid": taxid,
                      **{r_: lin.get(r_, "") for r_ in RANGOS}, "fasta": rutas[acc]})
    escribir_tsv(args.out, filas, ["accesion", "organismo", "taxid"] + RANGOS + ["fasta"])
    print("    Referencias con FASTA: %d (sin FASTA: %d) -> %s" % (len(filas), sin_fasta, args.out))


# Subcomando: its-ventanas
def cmd_its_ventanas(args):
    """barrnap marca 18S/5.8S/28S. Se toma cada hit +- flanco y se unen las
    ventanas solapadas; ITSx ubica despues los limites exactos del ITS."""
    n_total = 0
    with open(args.out, "w", encoding="utf-8") as out:
        for gff in sorted(glob.glob(os.path.join(args.gff_dir, "*.gff"))):
            name = bin_name(gff[:-4])
            hits = {}
            with open(gff, encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 9:
                        continue
                    hits.setdefault(p[0], []).append((int(p[3]), int(p[4])))
            if not hits:
                continue
            fasta = os.path.join(args.bins_dir, name + "." + args.bin_ext)
            for contig, seq in leer_fasta(fasta, solo=set(hits)):
                ivs = sorted((max(1, s - args.flank), min(len(seq), e + args.flank))
                             for s, e in hits[contig])
                unidas = [list(ivs[0])]
                for s, e in ivs[1:]:
                    if s <= unidas[-1][1]:
                        unidas[-1][1] = max(unidas[-1][1], e)
                    else:
                        unidas.append([s, e])
                for s, e in unidas:
                    out.write(">%s__%s__%d_%d\n%s\n" % (name, contig, s, e, seq[s - 1:e]))
                    n_total += 1
    print("    Ventanas de ADNr para ITSx: %d -> %s" % (n_total, args.out))


# Subcomandos: busco-genes y concat (matriz del arbol)
def cmd_busco_genes(args):
    nombres = [bin_name(g) for g in sorted(glob.glob(os.path.join(args.genomes_dir, "*." + args.bin_ext)))]
    genes = {}
    for n in nombres:
        d = glob.glob(os.path.join(args.busco_dir, n, "run_*", "busco_sequences",
                                   "single_copy_busco_sequences"))
        if not d:
            print("    AVISO: sin genes de copia unica de BUSCO para %s" % n)
            continue
        for faa in glob.glob(os.path.join(d[0], "*.faa")):
            seq = next(leer_fasta(faa), (None, ""))[1].replace("*", "")
            if seq:
                genes.setdefault(os.path.basename(faa)[:-4], {})[n] = seq
    minimo = max(2, math.ceil(args.min_frac * len(nombres)))
    os.makedirs(args.out_dir, exist_ok=True)
    k = 0
    for gid, seqs in sorted(genes.items()):
        if len(seqs) >= minimo:
            with open(os.path.join(args.out_dir, gid + ".faa"), "w", encoding="utf-8") as fh:
                for n, s in sorted(seqs.items()):
                    fh.write(">%s\n%s\n" % (n, s))
            k += 1
    print("    Genes BUSCO en >= %d de %d genomas: %d -> %s" % (minimo, len(nombres), k, args.out_dir))


def cmd_concat(args):
    nombres = [bin_name(g) for g in sorted(glob.glob(os.path.join(args.genomes_dir, "*." + args.bin_ext)))]
    matriz = {n: [] for n in nombres}
    particiones, pos = [], 1
    for aln in sorted(glob.glob(os.path.join(args.aln_dir, "*.aln"))):
        seqs = dict(leer_fasta(aln))
        if not seqs:
            continue
        largo = len(next(iter(seqs.values())))
        if largo == 0:
            continue
        for n in nombres:
            matriz[n].append(seqs.get(n, "-" * largo))
        particiones.append("LG, %s = %d-%d" % (bin_name(aln[:-4]), pos, pos + largo - 1))
        pos += largo
    with open(args.out_faa, "w", encoding="utf-8") as fh:
        for n in nombres:
            fh.write(">%s\n%s\n" % (n, "".join(matriz[n])))
    with open(args.out_partition, "w", encoding="utf-8") as fh:
        fh.write("\n".join(particiones) + "\n")
    total = pos - 1
    print("    Matriz concatenada: %d genes, %d posiciones -> %s" % (len(particiones), total, args.out_faa))
    for n in nombres:
        ocup = sum(1 for c in "".join(matriz[n]) if c != "-")
        print("      %s: %.1f %% de posiciones con aminoacido" % (n, 100.0 * ocup / total if total else 0))


# Subcomando: genomeinfo
def cmd_genomeinfo(args):
    sel = {f["bin"]: f for f in leer_tsv(args.phase1_tsv) if f["decision_final"] == "SELECCIONADO"}
    genomas = sorted(glob.glob(os.path.join(args.genomes_dir, "*." + args.bin_ext)))
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["genome", "completeness", "contamination"])
        for g in genomas:
            f = sel[bin_name(g)]
            # dRep espera el nombre de archivo exacto (con extension).
            w.writerow([os.path.basename(g), f["eukcc_completitud"], f["eukcc_contaminacion"]])
    print("    genomeInfo.csv (EukCC): %d genomas -> %s" % (len(genomas), args.out))


# Subcomando: fase2
def leer_skani(path, refs):
    """-> {bin: [hits ordenados por ANI y AF]}. AF en porcentaje."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    for r in leer_tsv(path):
        m = RE_ACC.search(os.path.basename(r.get("Ref_file", "")))
        acc = m.group(0) if m else ""
        hit = {"accesion": acc, "ani": num(r.get("ANI")),
               "af_query": num(r.get("Align_fraction_query")),
               "af_ref": num(r.get("Align_fraction_ref")),
               "ref": refs.get(acc, {})}
        out.setdefault(bin_name(r.get("Query_file", "")), []).append(hit)
    for v in out.values():
        v.sort(key=lambda h: (-(h["ani"] or 0), -(h["af_query"] or 0)))
    return out


def leer_blast_its(path, fuente):
    """blastn outfmt '6 qseqid pident length qlen qcovs evalue bitscore stitle'
    -> {bin: mejor hit por bitscore}."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            b = p[0].split("__")[0]
            hit = {"pid": float(p[1]), "qcov": float(p[4]), "bits": float(p[6]),
                   "titulo": p[7], "fuente": fuente}
            if b not in out or hit["bits"] > out[b]["bits"]:
                out[b] = hit
    return out


def taxo_unite(titulo):
    """'Fusarium_oxysporum|MH..|SH..|reps|k__Fungi;p__Ascomycota;...;s__Fusarium_oxysporum'"""
    lin = {}
    for seg in titulo.split("|"):
        if "k__" in seg:
            for t in seg.split(";"):
                if "__" in t:
                    pre, val = t.split("__", 1)
                    val = val.replace("_", " ").strip()
                    if pre in UNITE_PREFIJOS and val and not val.lower().startswith(("unidentified", "incertae")):
                        lin[UNITE_PREFIJOS[pre]] = val
    if lin:
        lin["dominio"] = "Eukaryota"
    return lin


def especie_refseq(titulo):
    """Titulo de RefSeq ITS -> 'Genero especie'. BLAST lo entrega sin accesion
    ('Papiliotrema fonsecae EX F-4087 ITS region; ...'); si viniera con ella
    ('NR_xxx.1 Fusarium ...'), se salta."""
    w = titulo.split()
    if w and re.match(r"^[A-Z]{1,2}_?\d+(\.\d+)?$", w[0]):
        w = w[1:]
    return " ".join(w[:2]) if len(w) >= 2 else ""


def recortar(lin, hasta):
    """Deja solo los rangos hasta 'hasta' (inclusive)."""
    i = RANGOS.index(hasta)
    return {r: lin.get(r, "") for r in RANGOS[: i + 1] if lin.get(r)}


def asignar(name, skani, unite, refseq, eukcc_lin, a):
    """Cascada de asignacion. Devuelve (linaje, nivel, metodo, conflicto)."""
    best = skani[0] if skani else None
    if best and best["ani"] is not None and best["af_query"] is not None \
            and best["ani"] >= a.min_ani and best["af_query"] >= a.min_af:
        return recortar(best["ref"], "especie"), "especie", "ANI (skani)", ""
    # ITS: primero UNITE (trae el linaje completo), luego RefSeq ITS (cepas tipo).
    for umbral, nivel in ((a.its_especie, "especie"), (a.its_genero, "genero")):
        for its in (unite, refseq):
            if not its or its["qcov"] < a.its_min_qcov or its["pid"] < umbral:
                continue
            if its["fuente"] == "UNITE":
                lin_its = taxo_unite(its["titulo"])
            else:
                sp = especie_refseq(its["titulo"])
                lin_its = {"especie": sp, "genero": sp.split()[0]} if sp else {}
            if lin_its.get(nivel):
                etiqueta = "especie (hipotesis ITS)" if nivel == "especie" else "genero"
                return recortar(lin_its, nivel), etiqueta, "ITS %s (%.2f %%)" % (its["fuente"], its["pid"]), ""
    its = unite
    if best and best["ref"].get("genero"):
        return (recortar(best["ref"], "genero"), "genero (vecino skani)",
                "skani ANI %.2f %% < %g" % (best["ani"], a.min_ani), "")
    if its and its["qcov"] >= a.its_min_qcov and its["fuente"] == "UNITE":
        lin_its = taxo_unite(its["titulo"])
        for umbral, nivel in ((a.its_familia, "familia"), (a.its_orden, "orden")):
            if its["pid"] >= umbral and lin_its.get(nivel):
                return recortar(lin_its, nivel), nivel, "ITS UNITE (%.2f %%)" % its["pid"], ""
    if eukcc_lin:
        prof = [r for r in RANGOS if eukcc_lin.get(r)]
        return eukcc_lin, prof[-1] if prof else "NA", "colocacion EukCC", ""
    return {}, "sin asignacion", "-", ""


def cmd_fase2(args):
    genomas = sorted(glob.glob(os.path.join(args.genomes_dir, "*." + args.bin_ext)))
    nombres = [bin_name(g) for g in genomas]
    p1 = {f["bin"]: f for f in leer_tsv(args.phase1_tsv)}
    refs = {r["accesion"]: r for r in leer_tsv(args.refs_tsv)} if args.refs_tsv and os.path.isfile(args.refs_tsv) else {}
    # genero -> rangos superiores (NCBI), para completar asignaciones hechas por ITS
    sobre_genero = {}
    for r in refs.values():
        if r.get("genero"):
            sobre_genero.setdefault(r["genero"], {k: r[k] for k in RANGOS[:RANGOS.index("genero")] if r.get(k)})
    skani = leer_skani(args.skani_tsv, refs)
    unite = leer_blast_its(args.its_unite, "UNITE")
    refseq = leer_blast_its(args.its_refseq, "RefSeq-ITS")
    copias = {}
    if args.its_fasta and os.path.isfile(args.its_fasta):
        for sid, _ in leer_fasta(args.its_fasta):
            b = sid.split("__")[0]
            copias[b] = copias.get(b, 0) + 1

    clusters, reps = {}, {}
    if args.cdb and os.path.isfile(args.cdb):
        for r in leer_tsv(args.cdb, ","):
            clusters[bin_name(r["genome"])] = r.get("secondary_cluster", "")
    if args.wdb and os.path.isfile(args.wdb):
        for r in leer_tsv(args.wdb, ","):
            reps[r["cluster"]] = bin_name(r["genome"])
    for i, n in enumerate(nombres):            # sin dRep (< 2 genomas): un linaje por bin
        clusters.setdefault(n, "solo_%d" % (i + 1))

    filas = []
    for n in nombres:
        f1 = p1.get(n, {})
        eukcc_lin = {r: f1.get("eukcc_" + r, "") for r in RANGOS if f1.get("eukcc_" + r)}
        sk = skani.get(n, [])
        lin, nivel, metodo, _ = asignar(n, sk, unite.get(n), refseq.get(n), eukcc_lin, args)
        for k, v in sobre_genero.get(lin.get("genero", ""), {}).items():
            lin.setdefault(k, v)
        best = sk[0] if sk else None
        u, rs = unite.get(n), refseq.get(n)
        # Conflicto: el genero del vecino skani y el del ITS no coinciden.
        g_sk = best["ref"].get("genero", "") if best else ""
        g_its = ""
        if u and u["pid"] >= args.its_genero:
            g_its = taxo_unite(u["titulo"]).get("genero", "")
        elif rs and rs["pid"] >= args.its_genero:
            g_its = especie_refseq(rs["titulo"]).split(" ")[0]
        conflicto = "genero skani=%s vs ITS=%s" % (g_sk, g_its) if g_sk and g_its and g_sk != g_its else ""
        # Otras especies que tambien pasan ANI/AF: la especie no queda resuelta
        # (p. ej. Alternaria seccion Alternaria). Se reporta, no se decide.
        otras = sorted({h["ref"].get("especie", "") for h in sk[1:]
                        if h["ani"] is not None and h["af_query"] is not None
                        and h["ani"] >= args.min_ani and h["af_query"] >= args.min_af
                        and h["ref"].get("especie") and h["ref"].get("especie") != lin.get("especie")})
        otras = "; ".join(otras)
        fila = {"genome": n, "muestra": sample_of(n), "linaje_cluster": clusters[n],
                **{r: lin.get(r, "") for r in RANGOS},
                "nivel_asignacion": nivel, "metodo": metodo,
                "skani_ref": best["accesion"] if best else "",
                "skani_organismo": best["ref"].get("organismo", "") if best else "",
                "ani": best["ani"] if best else "", "af_query": best["af_query"] if best else "",
                "af_ref": best["af_ref"] if best else "",
                "especie_confirmada_ani": ("si" if best and best["ani"] >= args.min_ani
                                           and best["af_query"] >= args.min_af else "no") if best else "NA",
                "its_copias": copias.get(n, 0),
                "its_unite": u["titulo"].split("|")[0].replace("_", " ") if u else "",
                "its_unite_pid": u["pid"] if u else "", "its_unite_qcov": u["qcov"] if u else "",
                "its_refseq": especie_refseq(rs["titulo"]) if rs else "",
                "its_refseq_pid": rs["pid"] if rs else "",
                "eukcc_linaje": f1.get("eukcc_linaje", ""),
                "conflicto": conflicto,
                "otras_especies_ani": otras}
        filas.append(fila)

    # Linajes: clusteres de dRep ordenados por prevalencia (prefijo E de eucariota)
    por_cl = {}
    for f in filas:
        por_cl.setdefault(f["linaje_cluster"], []).append(f)
    linajes = []
    for cl, miembros in por_cl.items():
        rep = reps.get(cl) or miembros[0]["genome"]
        frep = next(m for m in miembros if m["genome"] == rep)
        muestras = sorted(set(m["muestra"] for m in miembros), key=lambda x: int(x) if x.isdigit() else 0)
        linajes.append({"linaje_cluster": cl, "representante": rep, "n_genomas": len(miembros),
                        "prevalencia": len(muestras), "muestras": ",".join(muestras),
                        "taxonomia": ";".join(frep[r] for r in RANGOS[1:] if frep[r]),
                        "nivel_asignacion": frep["nivel_asignacion"],
                        "miembros": ",".join(m["genome"] for m in miembros)})
    linajes.sort(key=lambda x: (-x["prevalencia"], -x["n_genomas"], x["taxonomia"]))
    cl_a_id = {}
    for i, l in enumerate(linajes):
        l["linaje_id"] = "E%d" % (i + 1)
        cl_a_id[l["linaje_cluster"]] = l["linaje_id"]
    for f in filas:
        f["linaje_id"] = cl_a_id[f["linaje_cluster"]]

    cols = (["genome", "muestra", "linaje_id", "linaje_cluster"] + RANGOS
            + ["nivel_asignacion", "metodo", "skani_ref", "skani_organismo", "ani",
               "af_query", "af_ref", "especie_confirmada_ani", "its_copias", "its_unite",
               "its_unite_pid", "its_unite_qcov", "its_refseq", "its_refseq_pid",
               "eukcc_linaje", "conflicto", "otras_especies_ani"])
    escribir_tsv(args.out_tsv, filas, cols)
    escribir_tsv(args.out_linajes, linajes, ["linaje_id", "linaje_cluster", "representante",
                                             "n_genomas", "prevalencia", "muestras",
                                             "taxonomia", "nivel_asignacion", "miembros"])
    reporte_fase2(args, filas, linajes)
    print("    Genomas: %d | linajes: %d | especie por ANI: %d" % (
        len(filas), len(linajes), sum(1 for f in filas if f["nivel_asignacion"] == "especie")))
    print("    Tablas: %s, %s\n    Reporte: %s" % (args.out_tsv, args.out_linajes, args.out_report))


def reporte_fase2(args, filas, linajes):
    niveles = {}
    for f in filas:
        niveles[f["nivel_asignacion"]] = niveles.get(f["nivel_asignacion"], 0) + 1
    with open(args.out_report, "w", encoding="utf-8") as fh:
        fh.write("# Fase 2 eucariota — Taxonomia y linajes de los bins eucariotas\n\n")
        fh.write("**Job SLURM:** `%s` — **Fecha:** %s\n\n" % (args.job_id, datetime.now().strftime("%Y-%m-%d %H:%M")))
        fh.write("| Metrica | Valor |\n|---|---:|\n")
        fh.write("| Genomas (Fase 1 eucariota) | %d |\n" % len(filas))
        fh.write("| Linajes (dRep, %s ANI) | %d |\n" % (args.drep_sa, len(linajes)))
        for k, v in sorted(niveles.items()):
            fh.write("| Asignados a nivel: %s | %d |\n" % (k, v))
        fh.write("| Con ITS recuperado | %d |\n" % sum(1 for f in filas if f["its_copias"]))
        fh.write("| Conflictos de genero skani vs ITS | %d |\n" % sum(1 for f in filas if f["conflicto"]))
        fh.write("| Especie por ANI con otras especies sobre el umbral | %d |\n\n" % sum(1 for f in filas if f["otras_especies_ani"]))
        fh.write("## Linajes\n\n| Linaje | Prev. | Genomas | Muestras | Representante | Taxonomia | Nivel |\n")
        fh.write("|---|---:|---:|---|---|---|---|\n")
        for l in linajes:
            fh.write("| %s | %d | %d | %s | %s | %s | %s |\n" % (
                l["linaje_id"], l["prevalencia"], l["n_genomas"], l["muestras"],
                l["representante"], l["taxonomia"] or "-", l["nivel_asignacion"]))
        fh.write("\n## Detalle por genoma\n\n| Genoma | Linaje | Asignacion | Nivel | Metodo | "
                 "Vecino skani | ANI | AF | ITS UNITE (%) | ITS RefSeq (%) | Conflicto | Otras especies >= umbral |\n")
        fh.write("|---|---|---|---|---|---|---:|---:|---|---|---|---|\n")
        for f in sorted(filas, key=lambda x: x["linaje_id"]):
            nombre = f["especie"] or f["genero"] or f["familia"] or f["orden"] or f["clase"] or "-"
            fh.write("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" % (
                f["genome"], f["linaje_id"], nombre, f["nivel_asignacion"], f["metodo"],
                f["skani_organismo"] or "-", f["ani"] or "-", f["af_query"] or "-",
                ("%s (%s)" % (f["its_unite"], f["its_unite_pid"])) if f["its_unite"] else "-",
                ("%s (%s)" % (f["its_refseq"], f["its_refseq_pid"])) if f["its_refseq"] else "-",
                f["conflicto"] or "-", f["otras_especies_ani"] or "-"))
        fh.write("\n## Reglas\n\n")
        fh.write("1. skani: ANI >= %g %% y AF del bin >= %g %% -> especie.\n" % (args.min_ani, args.min_af))
        fh.write("2. ITS (cobertura >= %g %%), umbrales de Vu et al. 2019: especie >= %g, genero >= %g, "
                 "familia >= %g, orden >= %g %% de identidad.\n" % (
                     args.its_min_qcov, args.its_especie, args.its_genero, args.its_familia, args.its_orden))
        fh.write("3. Sin lo anterior: genero del vecino mas cercano de skani (ANI < %g %%).\n" % args.min_ani)
        fh.write("4. Ultimo recurso: linaje de la colocacion de EukCC (orientativo).\n")
        fh.write("\nEl 95 % de ANI es un umbral operativo: no hay un umbral de especie validado para hongos.\n")


# CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    s = sp.add_parser("seleccion")
    s.add_argument("--tiara-summary", required=True)
    s.add_argument("--bins-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--dest", required=True)
    s.set_defaults(func=cmd_seleccion)

    s = sp.add_parser("fase1")
    s.add_argument("--tiara-summary", required=True)
    s.add_argument("--bins-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--eukcc-dir", required=True)
    s.add_argument("--busco-dir", required=True)
    s.add_argument("--busco-lineages", required=True, help="lista separada por espacios")
    s.add_argument("--taxdump", required=True)
    s.add_argument("--min-completeness", type=float, default=70)
    s.add_argument("--max-contamination", type=float, default=5)
    s.add_argument("--selected-dir", required=True)
    s.add_argument("--out-tsv", required=True)
    s.add_argument("--out-report", required=True)
    s.add_argument("--job-id", default="NA")
    s.set_defaults(func=cmd_fase1)

    s = sp.add_parser("refs")
    s.add_argument("--datasets-tsv", required=True)
    s.add_argument("--fasta-list", required=True)
    s.add_argument("--taxdump", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_refs)

    s = sp.add_parser("its-ventanas")
    s.add_argument("--gff-dir", required=True)
    s.add_argument("--bins-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--flank", type=int, default=3000)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_its_ventanas)

    s = sp.add_parser("busco-genes")
    s.add_argument("--busco-dir", required=True, help="carpeta del linaje BUSCO (una subcarpeta por bin)")
    s.add_argument("--genomes-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--min-frac", type=float, default=0.5)
    s.add_argument("--out-dir", required=True)
    s.set_defaults(func=cmd_busco_genes)

    s = sp.add_parser("concat")
    s.add_argument("--aln-dir", required=True)
    s.add_argument("--genomes-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--out-faa", required=True)
    s.add_argument("--out-partition", required=True)
    s.set_defaults(func=cmd_concat)

    s = sp.add_parser("genomeinfo")
    s.add_argument("--phase1-tsv", required=True)
    s.add_argument("--genomes-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_genomeinfo)

    s = sp.add_parser("fase2")
    s.add_argument("--genomes-dir", required=True)
    s.add_argument("--bin-ext", default="fasta")
    s.add_argument("--phase1-tsv", required=True)
    s.add_argument("--refs-tsv")
    s.add_argument("--skani-tsv")
    s.add_argument("--its-fasta")
    s.add_argument("--its-unite")
    s.add_argument("--its-refseq")
    s.add_argument("--cdb")
    s.add_argument("--wdb")
    s.add_argument("--drep-sa", default="0.95")
    s.add_argument("--min-ani", type=float, default=95)
    s.add_argument("--min-af", type=float, default=65)
    s.add_argument("--its-min-qcov", type=float, default=90)
    s.add_argument("--its-especie", type=float, default=98.41)
    s.add_argument("--its-genero", type=float, default=94.3)
    s.add_argument("--its-familia", type=float, default=88.5)
    s.add_argument("--its-orden", type=float, default=81.2)
    s.add_argument("--out-tsv", required=True)
    s.add_argument("--out-linajes", required=True)
    s.add_argument("--out-report", required=True)
    s.add_argument("--job-id", default="NA")
    s.set_defaults(func=cmd_fase2)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
