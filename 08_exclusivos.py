#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_exclusivos.py
----------------
FASE 3 (OE3), Etapa 3.5: verificacion de los genes candidatos a exclusivos de
la Estela (decision D26). Un candidato es un gen del bin cuya familia de
Panaroo no aparece en ninguna referencia (07_particion.py).

Solo biblioteca estandar.

Filtros (de los baratos a los caros):
  F1 estructura      pseudogen de Bakta, < 100 aa, o CDS a < 100 pb de un
                     extremo del contig -> descartado.
  F4 especie amplia  tblastn contra TODOS los genomas post-QC descargados de la
                     especie (no solo los <= 50 del pangenoma); identidad >= 80 %
                     y cobertura de la consulta >= 80 % -> descartado (presente
                     en la especie).
  F2 anclaje         el contig lleva >= 1 gen de una familia presente en
                     referencias -> anclado; si no, contig huerfano (va a F6).
  F3 misma libreria  blastn contra los demas bins de la misma muestra; identidad
                     >= 95 % y cobertura >= 80 % -> bandera (elemento compartido
                     con otro organismo del cultivo o error de binning).
  F6 contig huerfano blastn remoto del contig contra core_nt; mejor hit de otro
                     genero -> descartado (contaminacion); del mismo genero ->
                     se trata como anclado; sin hit -> bandera.
  F5 origen          blastp remoto contra nr; el mejor hit clasifica el origen
                     probable (mismo genero / misma familia / mismo filo / otro
                     filo / sin hit = ORFan). Un hit de la MISMA especie con
                     identidad >= 80 % y cobertura >= 80 % -> descartado.

Clases finales: verificado | con_bandera | descartado (con motivo).

Subcomandos:
  preparar  F1 + F2; escribe las consultas para F4 y F3.       (SLURM, local)
  evaluar   lee F4 y F3; escribe las consultas para F5 y F6.    (SLURM, local)
  integrar  lee F5 y F6 (+ linajes de taxonkit); clase final.   (login/laptop)
  ebi       F5 por el REST del EBI (blastp contra UniProtKB), alternativa
            cuando la cola de NCBI no avanza. Salida en el mismo tabular.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime


DEF_MIN_AA = 100
DEF_MIN_BORDE = 100
DEF_F4_ID, DEF_F4_COV = 80.0, 80.0
DEF_F3_ID, DEF_F3_COV = 95.0, 80.0
DEF_F5_ID, DEF_F5_COV = 80.0, 80.0

OUTFMT_CAMPOS = ["qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore",
                 "staxids"]


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


def leer_fasta(path):
    seqs, actual = OrderedDict(), None
    with open(path, encoding="ascii", errors="replace") as fh:
        for linea in fh:
            if linea.startswith(">"):
                actual = linea[1:].split()[0]
                seqs[actual] = []
            elif actual is not None:
                seqs[actual].append(linea.strip())
    return OrderedDict((k, "".join(v)) for k, v in seqs.items())


def escribir_fasta(path, registros):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        for nombre, seq in registros:
            fh.write(">%s\n" % nombre)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def leer_blast(path):
    """Mejor hit por consulta (mayor bitscore) y lista completa de hits."""
    hits = defaultdict(list)
    if not path or not os.path.exists(path):
        return hits
    with open(path, encoding="utf-8") as fh:
        for linea in fh:
            if not linea.strip() or linea.startswith("#"):
                continue
            p = linea.rstrip("\n").split("\t")
            d = dict(zip(OUTFMT_CAMPOS, p))
            d["pident"] = a_float(d.get("pident"), 0.0)
            d["qcovs"] = a_float(d.get("qcovs"), 0.0)
            d["bitscore"] = a_float(d.get("bitscore"), 0.0)
            hits[d["qseqid"]].append(d)
    for q in hits:
        hits[q].sort(key=lambda d: -d["bitscore"])
    return hits


def primer_hit_que_cumple(hits, min_id, min_cov):
    for h in hits:
        if h["pident"] >= min_id and h["qcovs"] >= min_cov:
            return h
    return None


# ---------------------------------------------------------- subcomando preparar

def cmd_preparar(args):
    genes = leer_tsv(os.path.join(args.particion_dir, "genes_bin.tsv"))
    anclados = set((g["bin"], g["contig"]) for g in genes
                   if g["estado_panaroo"] == "en_familia" and a_float(g["n_ref"], 0) > 0)
    cand = [g for g in genes if g["exclusivo_candidato"] == "si"]
    os.makedirs(args.outdir, exist_ok=True)

    faa, ffn, fna = {}, {}, {}
    for b in sorted(set(g["bin"] for g in cand)):
        base = os.path.join(args.bakta_dir, b, b)
        faa[b], ffn[b], fna[b] = (leer_fasta(base + ext) for ext in (".faa", ".ffn", ".fna"))

    filas, q_faa, q_ffn, huerfanos = [], [], defaultdict(list), OrderedDict()
    for g in cand:
        b, loc = g["bin"], g["locus_tag"]
        dist = a_float(g["dist_borde_pb"], 0)
        if g["pseudo"] == "si":
            f1 = "descartado_F1_pseudogen"
        elif a_float(g["largo_aa"], 0) < args.min_aa:
            f1 = "descartado_F1_corto"
        elif dist < args.min_borde:
            f1 = "descartado_F1_borde_contig"
        else:
            f1 = "ok"
        anclado = (b, g["contig"]) in anclados
        fila = dict(g)
        fila.update({"consulta": loc, "f1": f1, "f2_anclaje": "anclado" if anclado else "huerfano"})
        filas.append(fila)
        if f1 != "ok":
            continue
        if loc not in faa[b] or loc not in ffn[b]:
            die("no se encontro %s en %s.faa/.ffn" % (loc, b))
        q_faa.append((loc, faa[b][loc]))
        q_ffn[b].append((loc, ffn[b][loc]))
        if not anclado:
            huerfanos["%s__%s" % (b, g["contig"])] = fna[b][g["contig"]]

    campos = ["consulta", "bin", "genoma_original", "locus_tag", "contig", "largo_contig",
              "gc_contig", "inicio", "fin", "hebra", "dist_borde_pb", "largo_aa", "gen",
              "producto", "familia", "f1", "f2_anclaje"]
    escribir_tsv(os.path.join(args.outdir, "candidatos.tsv"), filas, campos)
    escribir_fasta(os.path.join(args.outdir, "f1_ok.faa"), q_faa)
    for b, regs in q_ffn.items():
        escribir_fasta(os.path.join(args.outdir, "f1_ok_%s.ffn" % b), regs)
    escribir_fasta(os.path.join(args.outdir, "contigs_huerfanos.fna"), list(huerfanos.items()))
    c = Counter(f["f1"] for f in filas)
    print("[ok] %d candidatos | F1: %s | huerfanos (con F1 ok): %d"
          % (len(filas), ", ".join("%s=%d" % kv for kv in c.most_common()), len(huerfanos)))


# ----------------------------------------------------------- subcomando evaluar

def cmd_evaluar(args):
    filas = leer_tsv(os.path.join(args.outdir, "candidatos.tsv"))
    f4 = leer_blast(args.f4)
    f3 = leer_blast(args.f3)
    huerfanos = leer_fasta(os.path.join(args.outdir, "contigs_huerfanos.fna"))
    faa = leer_fasta(os.path.join(args.outdir, "f1_ok.faa"))

    q5, q6 = [], OrderedDict()
    for f in filas:
        f["f4"], f["f4_hit"], f["f3"], f["f3_hit"] = "", "", "", ""
        if f["f1"] != "ok":
            continue
        h4 = primer_hit_que_cumple(f4.get(f["consulta"], []), args.f4_id, args.f4_cov)
        if h4:
            f["f4"] = "descartado_F4_presente_en_especie"
            f["f4_hit"] = "%s (%.1f %% id, %.0f %% cov)" % (h4["sseqid"].split("__")[0],
                                                          h4["pident"], h4["qcovs"])
            continue
        f["f4"] = "ok"
        h3 = primer_hit_que_cumple(f3.get(f["consulta"], []), args.f3_id, args.f3_cov)
        if h3:
            f["f3"] = "bandera_F3_otro_bin_misma_libreria"
            f["f3_hit"] = "%s (%.1f %% id, %.0f %% cov)" % (h3["sseqid"].split("__")[0],
                                                          h3["pident"], h3["qcovs"])
        else:
            f["f3"] = "ok"
        q5.append((f["consulta"], faa[f["consulta"]]))
        if f["f2_anclaje"] == "huerfano":
            clave = "%s__%s" % (f["bin"], f["contig"])
            q6[clave] = huerfanos[clave]

    campos = list(filas[0].keys()) if filas else []
    escribir_tsv(os.path.join(args.outdir, "candidatos_local.tsv"), filas, campos)
    escribir_fasta(os.path.join(args.outdir, "f5_query.faa"), q5)
    escribir_fasta(os.path.join(args.outdir, "f6_query.fna"), list(q6.items()))
    c4 = Counter(f["f4"] for f in filas if f["f1"] == "ok")
    print("[ok] F4: %s | pasan a F5 (nr): %d | contigs huerfanos para F6: %d"
          % (", ".join("%s=%d" % kv for kv in c4.most_common()), len(q5), len(q6)))


# ---------------------------------------------------------- subcomando integrar

def leer_linajes(path):
    lin = {}
    with open(path, encoding="utf-8") as fh:
        for linea in fh:
            p = linea.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip():
                partes = (p[1].split(";") + ["", "", "", ""])[:4]
                lin[p[0].strip()] = dict(zip(("filo", "familia", "genero", "especie"),
                                             [x.strip() for x in partes]))
    return lin


def taxid_de(hit):
    return (hit.get("staxids") or "").split(";")[0].strip()


def origen(lin_hit, lin_ref):
    if not lin_hit:
        return "taxonomia_desconocida"
    for rango, etiqueta in (("especie", "misma_especie"), ("genero", "mismo_genero"),
                            ("familia", "misma_familia"), ("filo", "mismo_filo")):
        if lin_hit.get(rango) and lin_hit.get(rango) == lin_ref.get(rango):
            return etiqueta
    return "otro_filo"


def cmd_integrar(args):
    filas = leer_tsv(os.path.join(args.outdir, "candidatos_local.tsv"))
    for p in (args.f5, args.f6):
        if not os.path.exists(p):
            die("falta %s: corre 08_exclusivos_remoto.sh (deja el archivo vacio si no "
                "habia consultas)." % p)
    f5 = leer_blast(args.f5)
    f6 = leer_blast(args.f6)
    lin = leer_linajes(args.linajes)
    ref = lin.get(str(args.taxid_especie))
    if not ref:
        die("el taxid de la especie (%s) no esta en %s" % (args.taxid_especie, args.linajes))

    embudo = Counter()
    for f in filas:
        embudo["candidatos"] += 1
        f.update({"f6": "", "f6_hit": "", "f5_origen": "", "f5_hit": "", "f5_pident": "",
                  "f5_qcovs": "", "f5_taxon": "", "clase": "", "motivo": ""})
        if f["f1"] != "ok":
            f["clase"], f["motivo"] = "descartado", f["f1"]
            continue
        embudo["pasa_F1"] += 1
        if f["f4"] != "ok":
            f["clase"], f["motivo"] = "descartado", f["f4"]
            continue
        embudo["pasa_F4"] += 1
        banderas = []
        if f["f3"] != "ok":
            banderas.append(f["f3"])
        if f["f2_anclaje"] == "huerfano":
            h6 = (f6.get("%s__%s" % (f["bin"], f["contig"])) or [None])[0]
            if h6 is None:
                f["f6"] = "bandera_F6_huerfano_sin_hit"
                banderas.append(f["f6"])
            else:
                l6 = lin.get(taxid_de(h6), {})
                f["f6_hit"] = "%s %s (%.1f %% id)" % (taxid_de(h6), l6.get("especie", ""),
                                                    h6["pident"])
                if l6.get("genero") and l6.get("genero") == ref.get("genero"):
                    f["f6"] = "ok_mismo_genero"
                else:
                    f["f6"] = "descartado_F6_contaminacion"
                    f["clase"], f["motivo"] = "descartado", f["f6"]
                    continue
        embudo["pasa_F6"] += 1
        hits5 = f5.get(f["consulta"], [])
        if not hits5:
            f["f5_origen"] = "sin_hit_ORFan"
        else:
            h5 = hits5[0]
            l5 = lin.get(taxid_de(h5), {})
            f["f5_origen"] = origen(l5, ref)
            f["f5_hit"] = h5["sseqid"]
            f["f5_pident"] = "%.1f" % h5["pident"]
            f["f5_qcovs"] = "%.0f" % h5["qcovs"]
            f["f5_taxon"] = "%s | %s" % (taxid_de(h5), ";".join(
                l5.get(k, "") for k in ("filo", "familia", "genero", "especie")))
            misma = [h for h in hits5 if lin.get(taxid_de(h), {}).get("especie")
                     and lin[taxid_de(h)]["especie"] == ref.get("especie")
                     and h["pident"] >= args.f5_id and h["qcovs"] >= args.f5_cov]
            if misma:
                f["clase"], f["motivo"] = "descartado", "descartado_F5_misma_especie_en_nr"
                continue
        embudo["pasa_F5"] += 1
        if banderas:
            f["clase"], f["motivo"] = "con_bandera", ";".join(banderas)
            embudo["con_bandera"] += 1
        else:
            f["clase"], f["motivo"] = "verificado", ""
            embudo["verificado"] += 1

    campos = ["clase", "motivo", "bin", "genoma_original", "locus_tag", "familia", "gen",
              "producto", "largo_aa", "contig", "largo_contig", "gc_contig", "inicio", "fin",
              "hebra", "dist_borde_pb", "f1", "f4", "f4_hit", "f2_anclaje", "f3", "f3_hit",
              "f6", "f6_hit", "f5_origen", "f5_hit", "f5_pident", "f5_qcovs", "f5_taxon"]
    escribir_tsv(os.path.join(args.outdir, "exclusivos_verificados.tsv"), filas, campos)
    pasos = ["candidatos", "pasa_F1", "pasa_F4", "pasa_F6", "pasa_F5", "con_bandera",
             "verificado"]
    escribir_tsv(os.path.join(args.outdir, "embudo_exclusivos.tsv"),
                 [{"paso": p, "n": embudo[p]} for p in pasos], ["paso", "n"])
    origenes = Counter(f["f5_origen"] for f in filas if f["clase"] == "verificado")
    motivos = Counter(f["motivo"] for f in filas if f["clase"] == "descartado")
    with open(os.path.join(args.outdir, "exclusivos_report.md"), "w", encoding="utf-8",
              newline="\n") as fh:
        w = fh.write
        w("# Genes exclusivos - %s\n\n" % os.path.basename(os.path.normpath(args.outdir)))
        w("Generado: %s\n\n## Embudo\n\n| Paso | N |\n|---|---|\n"
          % datetime.now().strftime("%Y-%m-%d %H:%M"))
        for p in pasos:
            w("| %s | %d |\n" % (p, embudo[p]))
        w("\n## Motivos de descarte\n\n")
        for m, n in motivos.most_common():
            w("- %s: %d\n" % (m, n))
        w("\n## Origen probable de los exclusivos verificados (mejor hit en nr)\n\n")
        for o, n in origenes.most_common():
            w("- %s: %d\n" % (o, n))
        w("\nCheckpoint C7: todo candidato tiene clase; revisar a mano una muestra de 10 "
          "verificados (anotacion, contig, mejor hit).\n")
    print("[ok] embudo: " + " -> ".join("%s=%d" % (p, embudo[p]) for p in pasos))


# --------------------------------------------------------------- subcomando ebi

EBI_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"


def ebi_http(url, datos=None, intentos=5):
    """GET (datos=None) o POST al REST del EBI, con reintentos."""
    import time
    import urllib.parse
    import urllib.request
    cuerpo = urllib.parse.urlencode(datos).encode() if datos is not None else None
    for i in range(1, intentos + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=cuerpo),
                                        timeout=120) as r:
                return r.read().decode("utf-8")
        except Exception as e:
            if i == intentos:
                raise
            warn("EBI %s (intento %d): %s; reintento en 60 s" % (url.split("/")[-2], i, e))
            time.sleep(60)


def ebi_a_tabular(qseqid, qlen, res):
    """Hits del JSON del EBI -> filas outfmt 6 (OUTFMT_CAMPOS)."""
    filas = []
    for h in res.get("hits", []):
        hsps = h.get("hit_hsps") or []
        if not hsps:
            continue
        mejor = max(hsps, key=lambda s: s.get("hsp_bit_score", 0))
        # cobertura de la consulta: union de los tramos de todas las HSP (como qcovs)
        cubiertas = set()
        for s in hsps:
            a, b = sorted((s["hsp_query_from"], s["hsp_query_to"]))
            cubiertas.update(range(a, b + 1))
        qcov = 100.0 * len(cubiertas) / qlen if qlen else 0.0
        filas.append([qseqid, "%s|%s" % (h.get("hit_db", ""), h.get("hit_acc", "")),
                      "%.3f" % mejor["hsp_identity"], str(mejor["hsp_align_len"]),
                      "%.0f" % min(qcov, 100.0), "%g" % mejor["hsp_expect"],
                      "%.1f" % mejor["hsp_bit_score"], str(h.get("hit_uni_ox", ""))])
    return filas


def cmd_ebi(args):
    """blastp en el EBI (una secuencia por trabajo, hasta --paralelo a la vez).
    Cada resultado queda en <estado>/<consulta>.json: si se corta, al relanzar
    solo se envia lo que falta."""
    import json
    import time
    seqs = leer_fasta(args.query)
    os.makedirs(args.estado, exist_ok=True)
    def ruta(q):
        return os.path.join(args.estado, re.sub(r"[^A-Za-z0-9_.-]", "_", q))
    pendientes = [q for q in seqs if not os.path.exists(ruta(q) + ".json")]
    print("[ebi] %d consultas, %d pendientes (db %s)" % (len(seqs), len(pendientes), args.db))
    en_curso = {}
    for q in list(pendientes):
        if os.path.exists(ruta(q) + ".job"):
            en_curso[q] = open(ruta(q) + ".job").read().strip()
            pendientes.remove(q)
    t0 = time.time()
    while pendientes or en_curso:
        while pendientes and len(en_curso) < args.paralelo:
            q = pendientes.pop(0)
            job = ebi_http(EBI_URL + "/run", {
                "email": args.email, "program": "blastp", "stype": "protein",
                "database": args.db, "exp": args.evalue, "alignments": args.hits,
                "scores": args.hits, "sequence": ">%s\n%s\n" % (q, seqs[q])}).strip()
            if not job.startswith("ncbiblast-"):
                die("el EBI no acepto %s: %s" % (q, job[:300]))
            open(ruta(q) + ".job", "w").write(job)
            en_curso[q] = job
        time.sleep(15)
        for q, job in list(en_curso.items()):
            st = ebi_http("%s/status/%s" % (EBI_URL, job)).strip()
            if st == "FINISHED":
                with open(ruta(q) + ".json", "w", encoding="utf-8") as fh:
                    fh.write(ebi_http("%s/result/%s/json" % (EBI_URL, job)))
                os.remove(ruta(q) + ".job")
                del en_curso[q]
                print("[ebi] %s listo (%d min)" % (q, (time.time() - t0) // 60))
            elif st in ("ERROR", "FAILURE", "NOT_FOUND"):
                os.remove(ruta(q) + ".job")
                del en_curso[q]
                pendientes.append(q)
                warn("%s termino en %s; se reenvia" % (q, st))
    filas = []
    for q, seq in seqs.items():
        with open(ruta(q) + ".json", encoding="utf-8") as fh:
            filas += ebi_a_tabular(q, len(seq), json.load(fh))
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        for f in filas:
            fh.write("\t".join(f) + "\n")
    print("[ok] %d hits de %d consultas -> %s" % (len(filas), len(seqs), args.out))


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description="Fase 3: verificacion de genes exclusivos (F1-F6).")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("preparar", help="F1 + F2 y consultas de F3/F4")
    p.add_argument("--particion-dir", required=True)
    p.add_argument("--bakta-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-aa", type=int, default=DEF_MIN_AA)
    p.add_argument("--min-borde", type=int, default=DEF_MIN_BORDE)
    p.set_defaults(func=cmd_preparar)

    p = sub.add_parser("evaluar", help="aplica F4 y F3; consultas de F5/F6")
    p.add_argument("--outdir", required=True)
    p.add_argument("--f4", required=True, help="salida de tblastn (outfmt 6)")
    p.add_argument("--f3", required=True, help="salida de blastn contra otros bins")
    p.add_argument("--f4-id", type=float, default=DEF_F4_ID)
    p.add_argument("--f4-cov", type=float, default=DEF_F4_COV)
    p.add_argument("--f3-id", type=float, default=DEF_F3_ID)
    p.add_argument("--f3-cov", type=float, default=DEF_F3_COV)
    p.set_defaults(func=cmd_evaluar)

    p = sub.add_parser("integrar", help="aplica F5 y F6; clase final")
    p.add_argument("--outdir", required=True)
    p.add_argument("--f5", required=True, help="salida de blastp remoto contra nr")
    p.add_argument("--f6", required=True, help="salida de blastn remoto contra core_nt")
    p.add_argument("--linajes", required=True, help="taxid<TAB>filo;familia;genero;especie")
    p.add_argument("--taxid-especie", required=True, help="taxid NCBI de la especie del bin")
    p.add_argument("--f5-id", type=float, default=DEF_F5_ID)
    p.add_argument("--f5-cov", type=float, default=DEF_F5_COV)
    p.set_defaults(func=cmd_integrar)

    p = sub.add_parser("ebi", help="F5 por el BLAST del EBI (alternativa a NCBI)")
    p.add_argument("--query", required=True, help="FASTA de proteinas (f5_query.faa)")
    p.add_argument("--out", required=True, help="tabular outfmt 6 de salida")
    p.add_argument("--estado", required=True, help="carpeta con un .json por consulta")
    p.add_argument("--email", required=True, help="correo exigido por el EBI")
    p.add_argument("--db", default="uniprotkb")
    p.add_argument("--evalue", default="1e-5")
    p.add_argument("--hits", type=int, default=10)
    p.add_argument("--paralelo", type=int, default=5, help="trabajos simultaneos (EBI: <= 30)")
    p.set_defaults(func=cmd_ebi)

    args = ap.parse_args()
    if not getattr(args, "func", None):
        ap.print_help()
        sys.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
