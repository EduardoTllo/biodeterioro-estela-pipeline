#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07_particion.py
---------------
FASE 3 (OE3), Etapas 3.3-3.5: control de la anotacion y particion del pangenoma.

Solo biblioteca estandar.

Subcomandos:
  bakta      Checkpoint C5. Resume cada anotacion de Bakta (contigs, longitud,
             CDS, sORF, pseudogenes) y marca las referencias cuyo numero de CDS
             se aleja mas de +/- 20 % de la mediana de las referencias de su
             especie.
  particion  Particion core/shell/cloud calculada SOLO con las referencias
             (decision D22) a partir de gene_presence_absence.csv de Panaroo, y
             ubicacion de cada gen del bin de la Estela:
               core       f_ref >= 0.95  (sensibilidad: 0.90)
               shell      0.15 <= f_ref < 0.95
               cloud      0 < f_ref < 0.15
               exclusivo  f_ref = 0 y presente en >= 1 bin de la Estela
             Accesorio = shell + cloud; "solo cloud" queda como columna para la
             sensibilidad de la Fase 4 (decision D27). Calcula ademas la
             recuperacion del core de referencias en cada genoma: en las
             referencias debe ser >= 97 %, y en cada bin debe quedar a +/- 10
             puntos de su completitud CheckM2 (checkpoint C6).

Formato de Panaroo (v1.8.0, generate_output.py): gene_presence_absence.csv tiene
las columnas Gene, Non-unique Gene name, Annotation y una por genoma. Cada celda
lista los IDs de los genes (ID del GFF = locus tag de Bakta), separados por ';'
si hay paralogos; Panaroo puede anexar '_len' o '_pseudo' al ID, y los genes
re-encontrados aparecen como '<n>_refound_<m>'.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime


DEF_CORE = 0.95
DEF_CORE_SENS = 0.90
DEF_CLOUD = 0.15
DEF_MIN_RECUP_REF = 97.0     # % del core de referencias presente en cada referencia
DEF_MAX_DIF_BIN = 10.0       # puntos entre recuperacion del core y completitud CheckM2
DEF_TOL_CDS = 0.20           # +/- 20 % de la mediana de CDS de las referencias


def warn(msg):
    sys.stderr.write("[aviso] %s\n" % msg)


def die(msg):
    sys.stderr.write("[error] %s\n" % msg)
    sys.exit(1)


def leer_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def escribir_tsv(path, filas, campos):
    destino = os.path.dirname(os.path.abspath(path))
    if destino:
        os.makedirs(destino, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for f in filas:
            w.writerow(f)


def a_float(s, defecto=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return defecto


def mediana(valores):
    v = sorted(valores)
    if not v:
        return None
    m = len(v) // 2
    return v[m] if len(v) % 2 else (v[m - 1] + v[m]) / 2.0


def leer_bakta_tsv(path):
    """Filas de features del .tsv de Bakta (se saltan las lineas '#')."""
    campos = ["contig", "tipo", "inicio", "fin", "hebra", "locus", "gen", "producto",
              "dbxrefs"]
    filas = []
    with open(path, encoding="utf-8") as fh:
        for linea in fh:
            if linea.startswith("#") or not linea.strip():
                continue
            partes = linea.rstrip("\n").split("\t")
            partes += [""] * (len(campos) - len(partes))
            filas.append(dict(zip(campos, partes)))
    return filas


def leer_fna(path):
    """Devuelve {contig: longitud} y {contig: gc} del .fna de Bakta."""
    largos, gc, actual = {}, {}, None
    with open(path, encoding="ascii", errors="replace") as fh:
        for linea in fh:
            if linea.startswith(">"):
                actual = linea[1:].split()[0]
                largos[actual] = 0
                gc[actual] = 0
            elif actual is not None:
                s = linea.strip().upper()
                largos[actual] += len(s)
                gc[actual] += s.count("G") + s.count("C")
    return largos, {c: (gc[c] / float(largos[c]) if largos[c] else 0.0) for c in largos}


def completitud_checkm2(path):
    comp = {}
    if path and os.path.exists(path):
        for f in leer_tsv(path):
            comp[f.get("Name", "")] = a_float(f.get("Completeness"))
    elif path:
        warn("no se encontro el reporte de CheckM2: %s" % path)
    return comp


# ------------------------------------------------------------- subcomando bakta

def cmd_bakta(args):
    manif = leer_tsv(args.manifest)
    filas = []
    for m in manif:
        gid = m["id"]
        base = os.path.join(args.bakta_dir, gid, gid)
        fila = {"id": gid, "tipo": m["tipo"], "slug": m["slug"]}
        if not (os.path.exists(base + ".tsv") and os.path.exists(base + ".fna")):
            fila["estado"] = "sin_anotacion"
            filas.append(fila)
            continue
        feats = leer_bakta_tsv(base + ".tsv")
        largos, _ = leer_fna(base + ".fna")
        tipos = Counter(f["tipo"] for f in feats)
        fila.update({
            "estado": "ok", "n_contigs": len(largos), "longitud": sum(largos.values()),
            "n_cds": tipos.get("cds", 0), "n_sorf": tipos.get("sorf", 0),
            "n_pseudo": sum(1 for f in feats if f["producto"].startswith("(pseudo)")),
            "n_trna": tipos.get("tRNA", 0), "n_rrna": tipos.get("rRNA", 0),
        })
        filas.append(fila)

    por_especie = defaultdict(list)
    for f in filas:
        if f["tipo"] == "ref" and f.get("estado") == "ok":
            por_especie[f["slug"]].append(f["n_cds"])
    n_alerta = 0
    for f in filas:
        med = mediana(por_especie.get(f["slug"], []))
        f["mediana_cds_refs"] = med if med is not None else "NA"
        if f.get("estado") != "ok" or not med:
            f["alerta"] = "sin_anotacion" if f.get("estado") != "ok" else ""
            continue
        dif = (f["n_cds"] - med) / float(med)
        f["dif_rel_mediana"] = "%.3f" % dif
        if f["tipo"] == "ref" and abs(dif) > args.tolerancia:
            f["alerta"] = "cds_fuera_de_rango"
        elif f["tipo"] == "bin" and dif > args.tolerancia:
            f["alerta"] = "bin_con_exceso_de_cds"
        else:
            f["alerta"] = ""
        n_alerta += bool(f["alerta"])

    escribir_tsv(args.out, filas,
                 ["id", "tipo", "slug", "estado", "n_contigs", "longitud", "n_cds",
                  "n_sorf", "n_pseudo", "n_trna", "n_rrna", "mediana_cds_refs",
                  "dif_rel_mediana", "alerta"])
    faltan = [f["id"] for f in filas if f.get("estado") != "ok"]
    print("[ok] %s  (%d genomas, %d sin anotacion, %d con alerta)"
          % (args.out, len(filas), len(faltan), n_alerta))
    if faltan:
        warn("sin anotacion: %s" % ", ".join(faltan[:20]))
        sys.exit(1)


# --------------------------------------------------------- subcomando particion

def limpiar_token(t):
    t = t.strip()
    cambio = True
    while cambio:
        cambio = False
        for suf in ("_len", "_pseudo", "_stop"):
            if t.endswith(suf):
                t = t[: -len(suf)]
                cambio = True
    return t


def leer_panaroo(path, n_esperado=None):
    with open(path, encoding="utf-8") as fh:
        cab = fh.readline().rstrip("\n").split(",")
        genomas = cab[3:]
        n = len(genomas)
        familias = []
        for i, linea in enumerate(fh, 2):
            partes = linea.rstrip("\n").split(",")
            if len(partes) < 3 + n:
                die("linea %d de %s con %d campos (se esperaban >= %d)"
                    % (i, path, len(partes), 3 + n))
            celdas = partes[-n:]
            familias.append({
                "familia": partes[0],
                "nombre_no_unico": partes[1],
                # La anotacion puede contener comas: se reconstruye desde el centro.
                "anotacion": ",".join(partes[2:len(partes) - n]),
                "celdas": {g: [limpiar_token(t) for t in c.split(";") if t.strip()]
                           for g, c in zip(genomas, celdas)},
            })
    return genomas, familias


def categoria(n_ref, n_total, core, cloud, en_bin):
    if n_ref == 0:
        return "exclusivo" if en_bin else "ninguno"
    f = n_ref / float(n_total)
    if f >= core - 1e-9:
        return "core"
    if f >= cloud - 1e-9:
        return "shell"
    return "cloud"


def cmd_particion(args):
    p_gpa = os.path.join(args.panaroo_dir, "gene_presence_absence.csv")
    if not os.path.exists(p_gpa):
        die("falta %s" % p_gpa)
    manif = {m["id"]: m for m in leer_tsv(args.manifest) if m["slug"] == args.especie}
    genomas, familias = leer_panaroo(p_gpa)
    ajenos = [g for g in genomas if g not in manif]
    if ajenos:
        die("genomas en Panaroo que no estan en el manifiesto de %s: %s"
            % (args.especie, ", ".join(ajenos)))
    refs = [g for g in genomas if manif[g]["tipo"] == "ref"]
    bins = [g for g in genomas if manif[g]["tipo"] == "bin"]
    if len(refs) < 3 or not bins:
        die("se necesitan >= 3 referencias y >= 1 bin (hay %d y %d)" % (len(refs), len(bins)))
    nref = len(refs)
    comp = completitud_checkm2(args.checkm2)

    # --- particion por familia (frecuencias solo de referencias)
    filas_fam = []
    token_a_familia = {}
    for fam in familias:
        n_ref = sum(1 for g in refs if fam["celdas"][g])
        en_bin = [b for b in bins if fam["celdas"][b]]
        c95 = categoria(n_ref, nref, args.core, args.cloud, en_bin)
        c90 = categoria(n_ref, nref, args.core_sens, args.cloud, en_bin)
        fila = {
            "familia": fam["familia"], "nombre_no_unico": fam["nombre_no_unico"],
            "anotacion": fam["anotacion"], "n_ref": n_ref,
            "f_ref": "%.4f" % (n_ref / float(nref)),
            "categoria_95": c95, "categoria_90": c90,
            "accesorio": "si" if c95 in ("shell", "cloud") else "no",
            "solo_cloud": "si" if c95 == "cloud" else "no",
            "n_bins_presentes": len(en_bin),
        }
        for b in bins:
            fila["genes_" + b] = ";".join(fam["celdas"][b])
            for t in fam["celdas"][b]:
                token_a_familia[(b, t)] = fila
        fam["_fila"] = fila
        filas_fam.append(fila)

    os.makedirs(args.outdir, exist_ok=True)
    escribir_tsv(os.path.join(args.outdir, "particion_familias.tsv"), filas_fam,
                 ["familia", "nombre_no_unico", "anotacion", "n_ref", "f_ref",
                  "categoria_95", "categoria_90", "accesorio", "solo_cloud",
                  "n_bins_presentes"] + ["genes_" + b for b in bins])

    # --- matriz de presencia solo de referencias (para curvas y Heaps en R)
    with open(os.path.join(args.outdir, "pangenoma_refs.Rtab"), "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("\t".join(["Gene"] + refs) + "\n")
        for fam in familias:
            if fam["_fila"]["n_ref"] > 0:
                fh.write("\t".join([fam["familia"]] +
                                   ["1" if fam["celdas"][g] else "0" for g in refs]) + "\n")

    # --- recuperacion del core de referencias en cada genoma
    core_fams = [fam for fam in familias if fam["_fila"]["categoria_95"] == "core"]
    recup = []
    for g in genomas:
        pres = sum(1 for fam in core_fams if fam["celdas"][g])
        pct = 100.0 * pres / len(core_fams) if core_fams else 0.0
        fila = {"id": g, "tipo": manif[g]["tipo"],
                "genoma_original": manif[g]["genoma_original"],
                "familias_core_refs": len(core_fams), "presentes": pres,
                "recuperacion_core_pct": "%.2f" % pct}
        if manif[g]["tipo"] == "bin":
            c = comp.get(manif[g]["genoma_original"])
            fila["completitud_checkm2"] = "%.2f" % c if c is not None else "NA"
            if c is not None:
                fila["diferencia_pts"] = "%.2f" % (pct - c)
                fila["alerta"] = "revisar" if abs(pct - c) > args.max_dif_bin else ""
            else:
                fila["alerta"] = "sin_checkm2"
        else:
            fila["alerta"] = "revisar" if pct < args.min_recup_ref else ""
        recup.append(fila)
    escribir_tsv(os.path.join(args.outdir, "recuperacion_core.tsv"), recup,
                 ["id", "tipo", "genoma_original", "familias_core_refs", "presentes",
                  "recuperacion_core_pct", "completitud_checkm2", "diferencia_pts", "alerta"])

    # --- genes de cada bin con su categoria (interfaz con la Fase 4)
    genes, resumen_bins = [], []
    for b in bins:
        base = os.path.join(args.bakta_dir, b, b)
        if not os.path.exists(base + ".tsv"):
            die("falta %s.tsv" % base)
        largos, gc = leer_fna(base + ".fna")
        conteo = Counter()
        for f in leer_bakta_tsv(base + ".tsv"):
            if f["tipo"] not in ("cds", "sorf"):
                continue
            loc = f["locus"]
            fam = token_a_familia.get((b, loc)) or token_a_familia.get(
                (b, re.sub(r"\.\d+$", "", loc)))
            ini, fin = int(f["inicio"]), int(f["fin"])
            largo_contig = largos.get(f["contig"], 0)
            fila = {
                "bin": b, "genoma_original": manif[b]["genoma_original"],
                "locus_tag": loc, "tipo_feature": f["tipo"], "contig": f["contig"],
                "largo_contig": largo_contig, "gc_contig": "%.4f" % gc.get(f["contig"], 0),
                "inicio": ini, "fin": fin, "hebra": f["hebra"],
                "dist_borde_pb": min(ini - 1, largo_contig - fin) if largo_contig else "NA",
                "largo_aa": max((fin - ini + 1) // 3 - 1, 0),
                "gen": f["gen"], "producto": f["producto"],
                "pseudo": "si" if f["producto"].startswith("(pseudo)") else "no",
            }
            if fam is None:
                fila.update({"familia": "", "estado_panaroo": "no_en_pangenoma",
                             "n_ref": "", "f_ref": "", "categoria_95": "no_en_pangenoma",
                             "categoria_90": "no_en_pangenoma", "accesorio": "",
                             "solo_cloud": "", "exclusivo_candidato": "no"})
            else:
                fila.update({"familia": fam["familia"], "estado_panaroo": "en_familia",
                             "n_ref": fam["n_ref"], "f_ref": fam["f_ref"],
                             "categoria_95": fam["categoria_95"],
                             "categoria_90": fam["categoria_90"],
                             "accesorio": fam["accesorio"], "solo_cloud": fam["solo_cloud"],
                             "exclusivo_candidato": "si" if fam["categoria_95"] == "exclusivo" else "no"})
            conteo[fila["categoria_95"]] += 1
            genes.append(fila)
        r = next(x for x in recup if x["id"] == b)
        resumen_bins.append({
            "bin": b, "genoma_original": manif[b]["genoma_original"],
            "completitud_checkm2": r.get("completitud_checkm2", "NA"),
            "recuperacion_core_pct": r["recuperacion_core_pct"],
            "n_genes": sum(conteo.values()), "core": conteo["core"],
            "shell": conteo["shell"], "cloud": conteo["cloud"],
            "exclusivo_candidato": conteo["exclusivo"],
            "no_en_pangenoma": conteo["no_en_pangenoma"],
        })
    escribir_tsv(os.path.join(args.outdir, "genes_bin.tsv"), genes,
                 ["bin", "genoma_original", "locus_tag", "tipo_feature", "contig",
                  "largo_contig", "gc_contig", "inicio", "fin", "hebra", "dist_borde_pb",
                  "largo_aa", "gen", "producto", "pseudo", "familia", "estado_panaroo",
                  "n_ref", "f_ref", "categoria_95", "categoria_90", "accesorio",
                  "solo_cloud", "exclusivo_candidato"])
    escribir_tsv(os.path.join(args.outdir, "resumen_bins.tsv"), resumen_bins,
                 ["bin", "genoma_original", "completitud_checkm2", "recuperacion_core_pct",
                  "n_genes", "core", "shell", "cloud", "exclusivo_candidato",
                  "no_en_pangenoma"])

    # --- resumen de la especie
    c95 = Counter(f["categoria_95"] for f in filas_fam)
    c90 = Counter(f["categoria_90"] for f in filas_fam)
    resumen = {
        "slug": args.especie, "n_refs": nref, "n_bins": len(bins),
        "familias_refs": sum(1 for f in filas_fam if f["n_ref"] > 0),
        "core_95": c95["core"], "shell_95": c95["shell"], "cloud_95": c95["cloud"],
        "core_90": c90["core"], "shell_90": c90["shell"], "cloud_90": c90["cloud"],
        "familias_exclusivas": c95["exclusivo"],
    }
    escribir_tsv(os.path.join(args.outdir, "resumen_particion.tsv"), [resumen],
                 list(resumen.keys()))

    alertas = [r for r in recup if r.get("alerta")]
    with open(os.path.join(args.outdir, "particion_report.md"), "w", encoding="utf-8",
              newline="\n") as fh:
        w = fh.write
        w("# Particion del pangenoma - %s\n\n" % args.especie)
        w("Generado: %s. Frecuencias calculadas solo con las %d referencias (D22).\n\n"
          % (datetime.now().strftime("%Y-%m-%d %H:%M"), nref))
        w("| Categoria | Corte 95 % | Corte 90 % |\n|---|---|---|\n")
        for cat in ("core", "shell", "cloud"):
            w("| %s | %d | %d |\n" % (cat, c95[cat], c90[cat]))
        w("| familias solo en bins (exclusivas candidatas) | %d | %d |\n\n"
          % (c95["exclusivo"], c90["exclusivo"]))
        w("## Bins de la Estela\n\n| Bin | Completitud CheckM2 | Recuperacion del core | "
          "Genes | Core | Shell | Cloud | Exclusivos cand. | Fuera del pangenoma |\n"
          "|---|---|---|---|---|---|---|---|---|\n")
        for r in resumen_bins:
            w("| %s | %s | %s | %d | %d | %d | %d | %d | %d |\n" % (
                r["genoma_original"], r["completitud_checkm2"], r["recuperacion_core_pct"],
                r["n_genes"], r["core"], r["shell"], r["cloud"], r["exclusivo_candidato"],
                r["no_en_pangenoma"]))
        w("\n## Checkpoint C6\n\n")
        if alertas:
            for r in alertas:
                w("- REVISAR %s (%s): recuperacion del core %s %%%s\n" % (
                    r["id"], r["tipo"], r["recuperacion_core_pct"],
                    (", completitud %s" % r.get("completitud_checkm2")) if r["tipo"] == "bin" else ""))
        else:
            w("- Sin alertas: las referencias recuperan >= %.0f %% del core y los bins "
              "quedan a <= %.0f puntos de su completitud.\n" % (args.min_recup_ref, args.max_dif_bin))

    print("[ok] particion de %s: %d refs, %d bins | core=%d shell=%d cloud=%d exclusivas=%d"
          % (args.especie, nref, len(bins), c95["core"], c95["shell"], c95["cloud"],
             c95["exclusivo"]))
    for r in alertas:
        warn("checkpoint C6: revisar %s (%s), recuperacion del core %s %%"
             % (r["id"], r["tipo"], r["recuperacion_core_pct"]))


# ------------------------------------------------------------------------ main


# ------------------------------------------------------------ subcomando cog

def cogs_de(dbxrefs):
    """Letras de categoria COG de un feature de Bakta ('COG:L', no 'COG:COG0593')."""
    letras = []
    for x in dbxrefs.split(","):
        x = x.strip()
        if x.startswith("COG:") and re.fullmatch(r"[A-Z]+", x[4:]):
            letras.extend(x[4:])
    return letras


def cmd_cog(args):
    """Categoria COG de cada familia (la mas frecuente entre sus genes) y de cada
    gen del bin, para comparar la funcion de core, shell, cloud y exclusivos."""
    genomas, familias = leer_panaroo(os.path.join(args.panaroo_dir, "gene_presence_absence.csv"))
    cog_gen = {}
    for g in genomas:
        base = os.path.join(args.bakta_dir, g, g)
        if not os.path.exists(base + ".tsv"):
            die("falta %s.tsv" % base)
        for f in leer_bakta_tsv(base + ".tsv"):
            if f["tipo"] == "cds":
                cog_gen[(g, f["locus"])] = cogs_de(f["dbxrefs"])
    cat = {r["familia"]: r for r in leer_tsv(os.path.join(args.particion_dir, "particion_familias.tsv"))}
    filas = []
    for fam in familias:
        c = Counter()
        n = 0
        for g, genes in fam["celdas"].items():
            for loc in genes:
                letras = cog_gen.get((g, loc))
                if letras is None:
                    continue
                n += 1
                c.update(set(letras))
        top = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))
        info = cat.get(fam["familia"], {})
        filas.append({"familia": fam["familia"], "categoria_95": info.get("categoria_95", ""),
                      "n_ref": info.get("n_ref", ""), "cog": top[0][0] if top else "sin_COG",
                      "genes_con_cog": top[0][1] if top else 0,
                      "genes": n})
    escribir_tsv(args.out, filas, ["familia", "categoria_95", "n_ref", "cog", "genes_con_cog", "genes"])
    resumen = Counter((f["categoria_95"], f["cog"]) for f in filas)
    print("[ok] %d familias -> %s | sin COG: %d" % (len(filas), args.out,
          sum(v for (k, c), v in resumen.items() if c == "sin_COG")))

def main():
    ap = argparse.ArgumentParser(description="Fase 3: control de Bakta y particion del pangenoma.")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("bakta", help="resumen de las anotaciones (checkpoint C5)")
    p.add_argument("--manifest", required=True)
    p.add_argument("--bakta-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tolerancia", type=float, default=DEF_TOL_CDS)
    p.set_defaults(func=cmd_bakta)

    p = sub.add_parser("particion", help="particion solo con referencias + posicion del bin")
    p.add_argument("--panaroo-dir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--especie", required=True, help="slug de la especie (columna slug)")
    p.add_argument("--bakta-dir", required=True)
    p.add_argument("--checkm2", default="", help="quality_report.tsv de CheckM2 (Fase 1)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--core", type=float, default=DEF_CORE)
    p.add_argument("--core-sens", type=float, default=DEF_CORE_SENS)
    p.add_argument("--cloud", type=float, default=DEF_CLOUD)
    p.add_argument("--min-recup-ref", type=float, default=DEF_MIN_RECUP_REF)
    p.add_argument("--max-dif-bin", type=float, default=DEF_MAX_DIF_BIN)
    p.set_defaults(func=cmd_particion)

    p = sub.add_parser("cog", help="categoria COG de cada familia (desde Bakta)")
    p.add_argument("--panaroo-dir", required=True)
    p.add_argument("--particion-dir", required=True)
    p.add_argument("--bakta-dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_cog)

    args = ap.parse_args()
    if not getattr(args, "func", None):
        ap.print_help()
        sys.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
