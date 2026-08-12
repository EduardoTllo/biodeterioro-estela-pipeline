#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_phase2.py
-------------------
Utilidades de resumen para la FASE 2 (OE2): asignacion taxonomica y definicion
de linajes de la tesis de la Estela de Raimondi.

Solo biblioteca estandar (sin pandas), para correr con el python de cualquier
entorno conda del cluster.

Subcomandos:
  genomeinfo  Genera genomeInfo.csv (genome,completeness,contamination) para dRep
              a partir del quality_report.tsv de CheckM2 (Fase 1), de modo que
              dRep NO tenga que re-ejecutar CheckM.
  report      Cruza la clasificacion de GTDB-Tk con los clusteres de dRep, calcula
              la prevalencia espacial por linaje, confirma especie por ANI/AF,
              selecciona los top-N linajes y escribe phase2_report.md + TSV.

Reglas (tesis, seccion 4.1.4):
  - Confirmacion de especie: ANI >= 95 % y AF >= 65 %.
  - Linaje = cluster de especie a 95 % de ANI (dRep -sa 0.95).
  - Prevalencia espacial = nº de muestras distintas (libreria 49-72) que aportan
    al menos un genoma a cada cluster. Es el criterio primario de importancia.
  - Se priorizan los 3 linajes de mayor prevalencia para la Fase 3.
"""

import argparse
import csv
import glob
import os
import shutil
import sys
from datetime import datetime


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def bin_name(fname):
    """bin-1-49.fasta -> bin-1-49"""
    return os.path.splitext(os.path.basename(fname))[0]


def sample_of(name):
    """Extrae la muestra (libreria) del nombre bin-<N>-<muestra>."""
    core = bin_name(name)
    return core.split("-")[-1] if "-" in core else "NA"


def read_checkm2(path):
    """quality_report.tsv -> {name: (completeness, contamination)}."""
    out = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = [c.strip().lower() for c in header]
        ni, ci, ti = cols.index("name"), cols.index("completeness"), cols.index("contamination")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ni, ci, ti):
                continue
            try:
                out[p[ni]] = (float(p[ci]), float(p[ti]))
            except ValueError:
                continue
    return out


# ---------------------------------------------------------------------------
# Subcomando: genomeinfo
# ---------------------------------------------------------------------------
def cmd_genomeinfo(args):
    ext = args.bin_ext.lstrip(".")
    qc = read_checkm2(args.checkm2)
    genomes = sorted(glob.glob(os.path.join(args.genomes_dir, "*." + ext)))
    if not genomes:
        sys.exit("ERROR: no hay genomas *.%s en %s" % (ext, args.genomes_dir))

    missing = []
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genome", "completeness", "contamination"])
        for g in genomes:
            name = bin_name(g)
            if name not in qc:
                missing.append(name)
                continue
            comp, cont = qc[name]
            # dRep espera el nombre de archivo EXACTO (con extension).
            w.writerow([os.path.basename(g), comp, cont])

    print("    genomeInfo.csv escrito: %s (%d genomas)" % (args.out, len(genomes) - len(missing)))
    if missing:
        print("    ADVERTENCIA: sin metricas CheckM2 para: %s" % ", ".join(missing))


# ---------------------------------------------------------------------------
# Parseo de GTDB-Tk
# ---------------------------------------------------------------------------
def read_gtdbtk(gtdbtk_dir):
    """
    Lee los summary de CLASIFICACION de GTDB-Tk (bac120 y ar53) y devuelve
    {user_genome: {classification, ani, af, reference, method}}.

    IMPORTANTE: hay que excluir los archivos de identify/
    (gtdbtk.*.markers_summary.tsv y gtdbtk.translation_table_summary.tsv).
    Tambien terminan en "summary.tsv" y traen una columna 'name' con los
    nombres de genoma, pero NO tienen taxonomia. Si se leen, sobrescriben las
    entradas buenas de classify/ con valores vacios (se procesan despues por
    orden alfabetico) y el reporte sale sin taxonomia.
    """
    out = {}
    files = sorted(glob.glob(os.path.join(gtdbtk_dir, "**", "*summary.tsv"),
                             recursive=True))
    for f in files:
        base = os.path.basename(f)
        # Solo los summary de clasificacion.
        if "markers_summary" in base or "translation_table" in base:
            continue
        try:
            fh = open(f, encoding="utf-8")
        except OSError:
            # p.ej. enlace simbolico roto al copiar los resultados a otro SO.
            continue
        with fh:
            reader = csv.DictReader(fh, delimiter="\t")
            # Sin columna 'user_genome' no es un summary de clasificacion.
            if not reader.fieldnames or "user_genome" not in reader.fieldnames:
                continue
            for row in reader:
                g = row.get("user_genome")
                if not g:
                    continue
                # GTDB-Tk 2.x usa closest_genome_* (la asignacion de especie).
                # Se mantienen los nombres antiguos (fastani_*) por compatibilidad.
                ani = (_val(row, "closest_genome_ani")
                       or _val(row, "fastani_ani")
                       or _val(row, "closest_placement_ani"))
                af = (_val(row, "closest_genome_af")
                      or _val(row, "fastani_af")
                      or _val(row, "closest_placement_af"))
                ref = (_val(row, "closest_genome_reference")
                       or _val(row, "fastani_reference")
                       or _val(row, "closest_placement_reference"))
                out[g] = {
                    "classification": row.get("classification", ""),
                    "ani": ani,
                    "af": af,
                    "reference": ref,
                    "method": row.get("classification_method", ""),
                }
    return out


def _val(row, key):
    """Devuelve el valor de una columna, tratando 'N/A' como vacio."""
    v = (row.get(key) or "").strip()
    return "" if v.upper() in ("N/A", "NA", "NONE") else v


def read_drep_clusters(cdb_path):
    """Cdb.csv -> {genome_name_sin_ext: secondary_cluster}."""
    clusters = {}
    with open(cdb_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            g = row.get("genome", "")
            sc = row.get("secondary_cluster") or row.get("cluster") or ""
            if g:
                clusters[bin_name(g)] = sc
    return clusters


def read_drep_winners(wdb_path):
    """Wdb.csv -> {secondary_cluster: genome_representante_sin_ext}."""
    reps = {}
    if not os.path.isfile(wdb_path):
        return reps
    with open(wdb_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            g = row.get("genome", "")
            cl = row.get("cluster", "")
            if g and cl:
                reps[cl] = bin_name(g)
    return reps


# ---------------------------------------------------------------------------
# Subcomando: report
# ---------------------------------------------------------------------------
def cmd_report(args):
    ext = args.bin_ext.lstrip(".")
    genomes = sorted(glob.glob(os.path.join(args.genomes_dir, "*." + ext)))
    names = [bin_name(g) for g in genomes]

    gt = read_gtdbtk(args.gtdbtk_dir)
    clusters = read_drep_clusters(args.drep_cdb)
    winners = read_drep_winners(args.drep_wdb)

    min_ani, min_af = float(args.min_ani), float(args.min_af)

    # --- Fila por genoma ----------------------------------------------------
    per_genome = []
    for name in names:
        g = gt.get(name, {})
        cl = clusters.get(name, "NA")
        ani = g.get("ani", "")
        af_raw = g.get("af", "")
        # GTDB-Tk reporta ANI en porcentaje (98.5) y AF como fraccion (0.90).
        # Normalizamos AF a porcentaje para comparar con el umbral (65 %).
        af_pct = ""
        try:
            if af_raw != "":
                v = float(af_raw)
                if v <= 1.0:
                    v *= 100.0
                af_pct = round(v, 2)
        except ValueError:
            af_pct = ""
        # Confirmacion de especie por ANI/AF.
        conf = "NA"
        try:
            if ani != "" and af_pct != "":
                conf = "si" if (float(ani) >= min_ani and float(af_pct) >= min_af) else "no"
        except ValueError:
            conf = "NA"
        per_genome.append({
            "genome": name,
            "muestra": sample_of(name),
            "clasificacion": g.get("classification", ""),
            "ani": ani,
            "af": af_pct,
            "ref": g.get("reference", ""),
            "especie_confirmada": conf,
            "linaje_cluster": cl,
        })

    # --- Agrupar por cluster (linaje) --------------------------------------
    by_cluster = {}
    for r in per_genome:
        by_cluster.setdefault(r["linaje_cluster"], []).append(r)

    lineages = []
    for cl, members in by_cluster.items():
        samples = sorted(set(m["muestra"] for m in members))
        rep = winners.get(cl) or members[0]["genome"]
        # Taxonomia consenso: la del representante si existe, si no la del primero.
        rep_tax = next((m["clasificacion"] for m in members if m["genome"] == rep), "")
        if not rep_tax:
            rep_tax = members[0]["clasificacion"]
        lineages.append({
            "cluster": cl,
            "representante": rep,
            "n_genomas": len(members),
            "prevalencia": len(samples),
            "muestras": ",".join(samples),
            "taxonomia": rep_tax,
            "miembros": [m["genome"] for m in members],
        })

    # Ranking: prevalencia desc, luego nº de genomas desc.
    lineages.sort(key=lambda x: (-x["prevalencia"], -x["n_genomas"], x["cluster"]))
    top = int(args.top)
    for i, lin in enumerate(lineages):
        lin["linaje_id"] = "L%d" % (i + 1)
        lin["seleccionado_fase3"] = "SI" if i < top else "no"

    # --- Copiar genomas de los top-N linajes para la Fase 3 ----------------
    if args.selected_lineages_dir:
        os.makedirs(args.selected_lineages_dir, exist_ok=True)
        name_to_path = {bin_name(g): g for g in genomes}
        for lin in lineages[:top]:
            dst = os.path.join(args.selected_lineages_dir, lin["linaje_id"])
            os.makedirs(dst, exist_ok=True)
            for m in lin["miembros"]:
                src = name_to_path.get(m)
                if src:
                    shutil.copyfile(src, os.path.join(dst, os.path.basename(src)))

    # --- TSV por genoma -----------------------------------------------------
    cols = ["genome", "muestra", "clasificacion", "ani", "af", "ref",
            "especie_confirmada", "linaje_cluster"]
    # Anexar el linaje_id legible.
    cl_to_lid = {lin["cluster"]: lin["linaje_id"] for lin in lineages}
    with open(args.out_tsv, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols + ["linaje_id"]) + "\n")
        for r in per_genome:
            fh.write("\t".join(str(r[c]) for c in cols) + "\t"
                     + cl_to_lid.get(r["linaje_cluster"], "NA") + "\n")

    # --- Reporte legible ----------------------------------------------------
    write_report(args, per_genome, lineages, top, min_ani, min_af)

    print("    Genomas procesados: %d" % len(per_genome))
    print("    Linajes (clusteres) definidos: %d" % len(lineages))
    print("    Top-%d por prevalencia -> Fase 3" % top)
    print("    Reporte: %s" % args.out_report)


def write_report(args, per_genome, lineages, top, min_ani, min_af):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(args.out_report, "w", encoding="utf-8") as fh:
        fh.write("# Fase 2 — Asignacion taxonomica y definicion de linajes (OE2)\n\n")
        fh.write("**Proyecto:** Modelo predictivo de riesgo de deterioro de la "
                 "Estela de Raimondi\n\n")
        fh.write("**Cluster:** Khipu (UTEC) — job SLURM `%s`\n\n" % args.job_id)
        fh.write("**Fecha de ejecucion:** %s\n\n" % now)
        fh.write("---\n\n")

        fh.write("## 1. Resumen\n\n")
        fh.write("Partiendo de los %d genomas seleccionados en la Fase 1, se "
                 "asigno taxonomia con GTDB-Tk (datos R220), se confirmo la "
                 "especie por ANI y se agruparon los genomas redundantes en "
                 "**linajes** (clusteres de especie a 95 %% de ANI con dRep). "
                 "Cada linaje se prioriza por su **prevalencia espacial** "
                 "(nº de muestras distintas que lo aportan).\n\n"
                 % len(per_genome))
        fh.write("| Metrica | Valor |\n|---|---:|\n")
        fh.write("| Genomas de entrada (Fase 1) | %d |\n" % len(per_genome))
        fh.write("| Linajes definidos (clusteres 95%% ANI) | %d |\n" % len(lineages))
        n_conf = sum(1 for r in per_genome if r["especie_confirmada"] == "si")
        fh.write("| Genomas con especie confirmada (ANI>=%.0f%%, AF>=%.0f%%) | %d |\n"
                 % (min_ani, min_af, n_conf))
        fh.write("| Linajes priorizados para la Fase 3 | %d |\n\n" % min(top, len(lineages)))

        fh.write("## 2. Linajes priorizados (top-%d por prevalencia)\n\n" % top)
        fh.write("| Linaje | Prevalencia (nº muestras) | Genomas | Muestras | "
                 "Representante | Taxonomia (GTDB) | Fase 3 |\n")
        fh.write("|---|---:|---:|---|---|---|---|\n")
        for lin in lineages:
            fh.write("| %s | %d | %d | %s | %s | %s | %s |\n" % (
                lin["linaje_id"], lin["prevalencia"], lin["n_genomas"],
                lin["muestras"], lin["representante"],
                short_tax(lin["taxonomia"]), lin["seleccionado_fase3"]))
        fh.write("\n")

        fh.write("## 3. Detalle por genoma\n\n")
        fh.write("| Genoma | Muestra | Linaje | Clasificacion GTDB | ANI (%) | "
                 "AF (%) | Especie confirmada |\n")
        fh.write("|---|---|---|---|---:|---:|---|\n")
        cl_order = {lin["cluster"]: lin["linaje_id"] for lin in lineages}
        for r in sorted(per_genome, key=lambda x: cl_order.get(x["linaje_cluster"], "z")):
            fh.write("| %s | %s | %s | %s | %s | %s | %s |\n" % (
                r["genome"], r["muestra"], cl_order.get(r["linaje_cluster"], "NA"),
                short_tax(r["clasificacion"]),
                r["ani"] or "-", r["af"] or "-", r["especie_confirmada"]))
        fh.write("\n")

        fh.write("## 4. Parametros y trazabilidad\n\n")
        fh.write("- Clasificacion: GTDB-Tk `classify_wf` (marcadores bac120/ar53, datos R220).\n")
        fh.write("- Confirmacion de especie: ANI >= %.0f %% y AF >= %.0f %% "
                 "(FastANI interno de GTDB-Tk).\n" % (min_ani, min_af))
        fh.write("- Definicion de linaje: cluster de especie a 95 % de ANI "
                 "(dRep, `-pa 0.90 -sa 0.95`), reutilizando la calidad de CheckM2.\n")
        fh.write("- Prevalencia espacial: nº de muestras distintas por cluster.\n")
        fh.write("- Arbol filogenomico: **arbol de colocacion de `classify_wf`** "
                 "(`classify/*.classify.tree`). El **arbol de novo focalizado** "
                 "(query + referencias del taxon asignado + outgroup) se realizara "
                 "en la **Fase 4**, solo con los linajes seleccionados.\n")
        fh.write("- Tabla por genoma (TSV): `%s`\n" % os.path.basename(args.out_tsv))
        fh.write("- Genomas de los linajes top-%d (para Fase 3): "
                 "`phase2_selected_lineages/`\n\n" % top)
        fh.write("> Nota: solo se analizan los genomas que superaron la Fase 1. "
                 "La prevalencia no depende de las lecturas crudas (no disponibles) "
                 "y por eso es el criterio primario de importancia de un linaje.\n")


def short_tax(tax):
    """Acorta la cadena GTDB a genero/especie legible."""
    if not tax:
        return "-"
    parts = [p for p in tax.split(";") if len(p) > 3]
    # Tomar los dos ultimos rangos informativos (genero/especie).
    tail = parts[-2:] if len(parts) >= 2 else parts
    return " ".join(tail) if tail else tax


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Resumenes de la Fase 2 (OE2).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("genomeinfo", help="Generar genomeInfo.csv para dRep.")
    p1.add_argument("--genomes-dir", required=True)
    p1.add_argument("--bin-ext", default="fasta")
    p1.add_argument("--checkm2", required=True)
    p1.add_argument("--out", required=True)
    p1.set_defaults(func=cmd_genomeinfo)

    p2 = sub.add_parser("report", help="Reporte legible de la Fase 2.")
    p2.add_argument("--genomes-dir", required=True)
    p2.add_argument("--bin-ext", default="fasta")
    p2.add_argument("--gtdbtk-dir", required=True)
    p2.add_argument("--drep-cdb", required=True)
    p2.add_argument("--drep-wdb", required=True)
    p2.add_argument("--min-ani", default="95")
    p2.add_argument("--min-af", default="65")
    p2.add_argument("--top", default="3")
    p2.add_argument("--selected-lineages-dir", default="")
    p2.add_argument("--out-tsv", required=True)
    p2.add_argument("--out-report", required=True)
    p2.add_argument("--job-id", default="NA")
    p2.set_defaults(func=cmd_report)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
