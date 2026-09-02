#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_phase1.py
-------------------
Utilidades de resumen para la FASE 1 (OE1) de la tesis de la Estela de Raimondi.

Solo usa la biblioteca estandar de Python 3 (no requiere pandas), de modo que
corre con el python de cualquier entorno conda del cluster.

Subcomandos:
  tiara       Lee las salidas de Tiara + los FASTA de los bins y produce una
              tabla por bin con la composicion de dominio (ponderada por
              longitud) y la decision de filtrado.
  copy-prok   Copia al directorio de entrada de CheckM2 solo los bins que
              pasaron el filtro por dominio (procariotas).
  report      Cruza el resumen de Tiara con el reporte de CheckM2, aplica los
              umbrales de completitud/contaminacion, copia los genomas
              seleccionados y escribe un reporte legible (Markdown) + TSV.

Reglas (tesis, seccion 4.1.3):
  - Longitud minima de contig para Tiara: 3000 bp.
  - Se descarta un bin cuya fraccion MAYORITARIA de longitud sea eucariota
    u organelar.
  - Se conservan genomas con completitud > 70 % y contaminacion < 5 %.
"""

import argparse
import os
import shutil
import sys
from datetime import datetime

# Clases de Tiara (primera etapa) agrupadas por dominio de interes.
PROK_CLASSES = {"archaea", "bacteria", "prokarya"}
EUK_CLASSES = {"eukarya"}
ORGANELLE_CLASSES = {"organelle", "mitochondrion", "plastid"}


# Utilidades FASTA
def read_fasta_lengths(path):
    """Devuelve dict {contig_id: longitud} leyendo un FASTA (id = primer token)."""
    lengths = {}
    cur = None
    n = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = n
                cur = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
        if cur is not None:
            lengths[cur] = n
    return lengths


def list_bins(bins_dir, ext):
    ext = ext.lstrip(".")
    out = []
    for f in sorted(os.listdir(bins_dir)):
        if f.endswith("." + ext):
            out.append(os.path.join(bins_dir, f))
    return out


def bin_name(path):
    """Nombre del bin sin extension (p.ej. muestra49_bin1.fa -> muestra49_bin1)."""
    return os.path.splitext(os.path.basename(path))[0]


# Parseo de la salida de Tiara
def read_tiara(path):
    """
    Lee un archivo de salida de Tiara.
    Devuelve dict {contig_id: clase_primera_etapa (minusculas)}.
    Tiara escribe columnas separadas por tab; la 1a es el id de la secuencia y
    la 2a es la clasificacion de la primera etapa.
    """
    labels = {}
    if not os.path.isfile(path):
        return labels
    with open(path, encoding="utf-8") as fh:
        header = fh.readline()  # descartar encabezado
        # Detectar columnas por nombre si el encabezado las trae.
        cols = [c.strip().lower() for c in header.rstrip("\n").split("\t")]
        try:
            id_idx = cols.index("sequence_id")
        except ValueError:
            id_idx = 0
        try:
            cls_idx = cols.index("class_fst_stage")
        except ValueError:
            cls_idx = 1 if len(cols) > 1 else 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(id_idx, cls_idx):
                continue
            cid = parts[id_idx].split()[0]
            cls = parts[cls_idx].strip().lower()
            labels[cid] = cls
    return labels


def classify_bin(bin_path, tiara_path):
    """
    Calcula la composicion de dominio de un bin, ponderada por longitud.
    Devuelve un dict con longitudes por dominio, dominio mayoritario y decision.
    """
    lengths = read_fasta_lengths(bin_path)
    labels = read_tiara(tiara_path)

    total = sum(lengths.values())
    buckets = {"prok": 0, "euk": 0, "organelle": 0, "unknown": 0}
    for cid, L in lengths.items():
        cls = labels.get(cid, "unknown")
        if cls in PROK_CLASSES:
            buckets["prok"] += L
        elif cls in EUK_CLASSES:
            buckets["euk"] += L
        elif cls in ORGANELLE_CLASSES:
            buckets["organelle"] += L
        else:
            buckets["unknown"] += L

    # Dominio mayoritario por longitud.
    majority = max(buckets, key=lambda k: buckets[k])
    # Regla de la tesis: descartar si la mayoria es euk u organelar.
    if majority in ("euk", "organelle"):
        decision = "DESCARTADO_dominio"
    else:
        decision = "PASA_dominio"

    def pct(x):
        return round(100.0 * x / total, 2) if total > 0 else 0.0

    return {
        "n_contigs": len(lengths),
        "length_bp": total,
        "pct_prok": pct(buckets["prok"]),
        "pct_euk": pct(buckets["euk"]),
        "pct_organelle": pct(buckets["organelle"]),
        "pct_unknown": pct(buckets["unknown"]),
        "majority": majority,
        "decision": decision,
    }


# Subcomando: tiara
def cmd_tiara(args):
    bins = list_bins(args.bins_dir, args.bin_ext)
    if not bins:
        sys.exit("ERROR: no se hallaron bins en %s" % args.bins_dir)

    rows = []
    for b in bins:
        name = bin_name(b)
        tpath = os.path.join(args.tiara_dir, name + ".tiara.txt")
        info = classify_bin(b, tpath)
        info["bin"] = name
        rows.append(info)

    cols = ["bin", "n_contigs", "length_bp", "pct_prok", "pct_euk",
            "pct_organelle", "pct_unknown", "majority", "decision"]
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_pass = sum(1 for r in rows if r["decision"] == "PASA_dominio")
    print("    Bins evaluados: %d | pasan dominio: %d | descartados: %d"
          % (len(rows), n_pass, len(rows) - n_pass))
    print("    Resumen escrito en: %s" % args.out)


# Subcomando: copy-prok
def cmd_copy_prok(args):
    ext = args.bin_ext.lstrip(".")
    os.makedirs(args.dest, exist_ok=True)
    n = 0
    with open(args.summary, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        bi = header.index("bin")
        di = header.index("decision")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[di] == "PASA_dominio":
                src = os.path.join(args.bins_dir, parts[bi] + "." + ext)
                dst = os.path.join(args.dest, parts[bi] + "." + ext)
                if os.path.isfile(src):
                    shutil.copyfile(src, dst)
                    n += 1
    print("    Copiados %d bins procariotas a %s" % (n, args.dest))


# Parseo del reporte de CheckM2
def read_checkm2(path):
    """
    Lee quality_report.tsv de CheckM2.
    Devuelve dict {name: (completeness_float, contamination_float)}.
    """
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = [c.strip().lower() for c in header]
        ni = cols.index("name")
        ci = cols.index("completeness")
        ti = cols.index("contamination")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ni, ci, ti):
                continue
            name = parts[ni]
            try:
                comp = float(parts[ci])
                cont = float(parts[ti])
            except ValueError:
                continue
            out[name] = (comp, cont)
    return out


# Subcomando: report
def cmd_report(args):
    ext = args.bin_ext.lstrip(".")

    # Cargar resumen de Tiara.
    tiara = {}
    with open(args.tiara_summary, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            tiara[p[idx["bin"]]] = {
                "n_contigs": int(p[idx["n_contigs"]]),
                "length_bp": int(p[idx["length_bp"]]),
                "pct_euk": float(p[idx["pct_euk"]]),
                "pct_organelle": float(p[idx["pct_organelle"]]),
                "majority": p[idx["majority"]],
                "decision": p[idx["decision"]],
            }

    checkm2 = read_checkm2(args.checkm2)

    min_comp = float(args.min_completeness)
    max_cont = float(args.max_contamination)

    os.makedirs(args.selected_dir, exist_ok=True)

    rows = []
    n_total = len(tiara)
    n_pass_dom = 0
    n_selected = 0

    for name in sorted(tiara.keys()):
        t = tiara[name]
        dom_ok = (t["decision"] == "PASA_dominio")
        if dom_ok:
            n_pass_dom += 1
        comp, cont = checkm2.get(name, ("", ""))

        final = "RECHAZADO"
        reason = ""
        if not dom_ok:
            reason = "dominio mayoritario: %s" % t["majority"]
        elif name not in checkm2:
            reason = "sin resultado CheckM2"
        else:
            qc_ok = (comp > min_comp) and (cont < max_cont)
            if qc_ok:
                final = "SELECCIONADO"
                n_selected += 1
                src = os.path.join(args.bins_dir, name + "." + ext)
                if os.path.isfile(src):
                    shutil.copyfile(src, os.path.join(args.selected_dir,
                                                      name + "." + ext))
            else:
                fails = []
                if not comp > min_comp:
                    fails.append("completitud %.2f<=%.0f" % (comp, min_comp))
                if not cont < max_cont:
                    fails.append("contaminacion %.2f>=%.0f" % (cont, max_cont))
                reason = "; ".join(fails)

        rows.append({
            "bin": name,
            "n_contigs": t["n_contigs"],
            "length_bp": t["length_bp"],
            "dominio_majority": t["majority"],
            "dominio_decision": t["decision"],
            "completeness": comp,
            "contamination": cont,
            "decision_final": final,
            "motivo": reason,
        })

    # --- TSV de seleccion ----------------------------------------------------
    cols = ["bin", "n_contigs", "length_bp", "dominio_majority",
            "dominio_decision", "completeness", "contamination",
            "decision_final", "motivo"]
    with open(args.out_tsv, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # --- Reporte Markdown legible -------------------------------------------
    selected = [r for r in rows if r["decision_final"] == "SELECCIONADO"]
    write_markdown_report(args, rows, selected, n_total, n_pass_dom,
                          n_selected, min_comp, max_cont)

    print("    Genomas de entrada:        %d" % n_total)
    print("    Pasan filtro por dominio:  %d" % n_pass_dom)
    print("    SELECCIONADOS (OE1):       %d" % n_selected)
    print("    Reporte:  %s" % args.out_report)
    print("    Tabla:    %s" % args.out_tsv)


def write_markdown_report(args, rows, selected, n_total, n_pass_dom,
                          n_selected, min_comp, max_cont):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def md_table(rs):
        head = ("| Bin | Contigs | Long (bp) | Dom. mayoritario | "
                "Completitud (%) | Contaminacion (%) | Decision | Motivo |\n"
                "|---|---:|---:|---|---:|---:|---|---|\n")
        body = ""
        for r in rs:
            comp = ("%.2f" % r["completeness"]) if r["completeness"] != "" else "-"
            cont = ("%.2f" % r["contamination"]) if r["contamination"] != "" else "-"
            body += "| %s | %d | %d | %s | %s | %s | %s | %s |\n" % (
                r["bin"], r["n_contigs"], r["length_bp"], r["dominio_majority"],
                comp, cont, r["decision_final"], r["motivo"] or "-")
        return head + body

    with open(args.out_report, "w", encoding="utf-8") as fh:
        fh.write("# Fase 1 — Control de calidad y seleccion de genomas (OE1)\n\n")
        fh.write("**Proyecto:** Modelo predictivo de riesgo de deterioro de la "
                 "Estela de Raimondi\n\n")
        fh.write("**Cluster:** Khipu (UTEC) — job SLURM `%s`\n\n" % args.job_id)
        fh.write("**Fecha de ejecucion:** %s\n\n" % now)
        fh.write("---\n\n")

        fh.write("## 1. Resumen del proceso\n\n")
        fh.write("La Fase 1 aplica dos filtros secuenciales sobre los bins "
                 "genomicos crudos entregados por el centro de secuenciacion, "
                 "para quedarnos con genomas procariotas de calidad suficiente "
                 "para el analisis pangenomico posterior.\n\n")
        fh.write("| Etapa | Herramienta | Criterio | Entran | Pasan |\n")
        fh.write("|---|---|---|---:|---:|\n")
        fh.write("| 1. Filtrado por dominio | Tiara | contig >= 3000 bp; se "
                 "descarta el bin si su longitud mayoritaria es eucariota u "
                 "organelar | %d | %d |\n" % (n_total, n_pass_dom))
        fh.write("| 2. Calidad genomica | CheckM2 | completitud > %.0f%% y "
                 "contaminacion < %.0f%% | %d | %d |\n"
                 % (min_comp, max_cont, n_pass_dom, n_selected))
        fh.write("\n")
        fh.write("**Resultado final: %d de %d genomas seleccionados** para la "
                 "Fase 2.\n\n" % (n_selected, n_total))

        fh.write("## 2. Genomas seleccionados (OE1)\n\n")
        if selected:
            fh.write(md_table(sorted(
                selected, key=lambda r: (-float(r["completeness"] or 0)))))
        else:
            fh.write("_Ningun genoma supero ambos filtros._\n")
        fh.write("\n")

        fh.write("## 3. Detalle completo (todos los bins)\n\n")
        fh.write(md_table(rows))
        fh.write("\n")

        fh.write("## 4. Parametros y trazabilidad\n\n")
        fh.write("- Longitud minima de contig (Tiara): **%s bp**\n"
                 % os.environ.get("MIN_CONTIG_LEN", "3000"))
        fh.write("- Umbral de completitud: **> %.0f %%**\n" % min_comp)
        fh.write("- Umbral de contaminacion: **< %.0f %%**\n" % max_cont)
        fh.write("- Tabla completa en formato TSV: `%s`\n"
                 % os.path.basename(args.out_tsv))
        fh.write("- Clasificacion por contig (Tiara): `01_tiara/`\n")
        fh.write("- Reporte de calidad (CheckM2): `02_checkm2/quality_report.tsv`\n")
        fh.write("- Genomas seleccionados (FASTA): `phase1_selected_genomes/`\n\n")
        fh.write("> Nota: el umbral (>70%%, <5%%) no corresponde a una categoria "
                 "MIMAG estandar; es un compromiso propio del proyecto que "
                 "flexibiliza la completitud pero mantiene la contaminacion "
                 "estricta para no sesgar la reconstruccion del genoma "
                 "accesorio (tesis, seccion 4.1.3.2).\n")


# CLI
def main():
    ap = argparse.ArgumentParser(description="Resumenes de la Fase 1 (OE1).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("tiara", help="Resumen de dominio por bin.")
    p1.add_argument("--bins-dir", required=True)
    p1.add_argument("--bin-ext", default="fa")
    p1.add_argument("--tiara-dir", required=True)
    p1.add_argument("--out", required=True)
    p1.set_defaults(func=cmd_tiara)

    p2 = sub.add_parser("copy-prok", help="Copiar bins procariotas.")
    p2.add_argument("--summary", required=True)
    p2.add_argument("--bins-dir", required=True)
    p2.add_argument("--bin-ext", default="fa")
    p2.add_argument("--dest", required=True)
    p2.set_defaults(func=cmd_copy_prok)

    p3 = sub.add_parser("report", help="Seleccion final + reporte legible.")
    p3.add_argument("--bins-dir", required=True)
    p3.add_argument("--bin-ext", default="fa")
    p3.add_argument("--tiara-summary", required=True)
    p3.add_argument("--checkm2", required=True)
    p3.add_argument("--min-completeness", default="70")
    p3.add_argument("--max-contamination", default="5")
    p3.add_argument("--selected-dir", required=True)
    p3.add_argument("--out-tsv", required=True)
    p3.add_argument("--out-report", required=True)
    p3.add_argument("--job-id", default="NA")
    p3.set_defaults(func=cmd_report)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
