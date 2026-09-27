#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_censo_genomas.py
-------------------
FASE 3 (OE3), Etapas 3.0-3.2: censo de disponibilidad genomica publica para los
linajes definidos en la Fase 2, y seleccion de las especies diana para el
analisis pangenomico.

Solo biblioteca estandar (sin pandas), para correr con el python de cualquier
entorno conda del cluster Khipu.

Criterio de seleccion (acordado con asesoria; reemplaza al de prevalencia
espacial de la Fase 2): DISPONIBILIDAD GENOMICA. El umbral es de >= 15 genomas,
que es el minimo que PPanGGOLiN recomienda para que su particionado estadistico
en core/shell/cloud sea confiable (Gautreau et al. 2020, PLoS Comput Biol 16(3):
e1007732). Se aplica al conteo DEPURADO (post-QC y post-desreplicacion) y no al
bruto de NCBI: una especie con 17 MAGs fragmentados no es viable, y los genomas
redundantes inflan el core de forma artificial (Guerra 2026, Bioinform Adv 6(1):
vbag069).

Subcomandos:
  censo   Cruza phase2_selection.tsv con el metadata de GTDB (bac120_metadata_r220.tsv)
          y produce, por linaje: mapeo a NCBI Taxonomy, tamano del cluster de
          especie, conteo RefSeq/GenBank, conteo post-QC y habitats representados.
          Emite ademas la lista de accesiones candidatas que consume 04_fetch_refs.sh.
          NO requiere internet.
  ncbi    Consulta el CLI de NCBI Datasets para los conteos GenBank/RefSeq vigentes.
          REQUIERE INTERNET -> correr en el nodo de login de Khipu.
  tabla   Integra censo + ncbi, aplica el umbral y emite la tabla de viabilidad
          y el reporte en markdown para la tesis.

Notas taxonomicas que el script resuelve automaticamente:
  - El nombre NCBI se deriva por voto mayoritario de ncbi_species_taxid dentro del
    cluster de especie de GTDB. Esto resuelve solo casos como GTDB "Telluria timonae"
    -> NCBI "Massilia timonae", donde una busqueda por el nombre GTDB daria 0.
  - Los sufijos GTDB (_A, _B, _AB) marcan especies escindidas que NCBI no separa.
    En esos linajes el conteo NCBI esta INFLADO y el numero comparable es el del
    cluster de GTDB. El script lo marca en la columna 'conteo_ncbi_comparable'.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime


# ---------------------------------------------------------------- utilidades

# Umbrales de admision de genomas publicos al set de referencia (Etapa 3.1).
DEF_MIN_COMPLETITUD = 95.0
DEF_MAX_CONTAMINACION = 5.0
DEF_MAX_CONTIGS = 300
DEF_UMBRAL = 15          # minimo de genomas post-QC para declarar viabilidad
                         # (minimo recomendado por PPanGGOLiN para particionar)
DEF_TOPE = 50            # tope de genomas por especie (computo + sesgo clonal)

# Categorias de NCBI que indican que el genoma NO proviene de un aislado.
NO_AISLADO = ("derived from metagenome", "derived from environmental sample",
              "derived from single cell")

RE_PLACEHOLDER = re.compile(r"^sp\d+$")
RE_SUFIJO = re.compile(r"_[A-Z]{1,2}$")


def warn(msg):
    sys.stderr.write("[aviso] %s\n" % msg)


def die(msg):
    sys.stderr.write("[error] %s\n" % msg)
    sys.exit(1)


def leer_tsv(path):
    """Lee un TSV con cabecera y devuelve lista de dicts."""
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


def a_int(s, defecto=None):
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return defecto


def rango_taxon(taxonomia, rango):
    """De 'd__Bacteria;...;s__Bacillus licheniformis' extrae el rango pedido."""
    for parte in (taxonomia or "").split(";"):
        parte = parte.strip()
        if parte.startswith(rango + "__"):
            return parte[len(rango) + 2:].strip()
    return ""


def clasificar_nivel(especie):
    """Devuelve 'especie' | 'placeholder' | 'sin_especie'.

    'placeholder' son los nombres alfanumericos de GTDB (p. ej. sp021023135):
    son clusteres validos a nivel de especie, pero sin binomio publicado ni
    equivalente estable en NCBI, asi que no admiten censo por nombre.
    """
    if not especie:
        return "sin_especie"
    partes = especie.split()
    if len(partes) >= 2 and RE_PLACEHOLDER.match(partes[-1]):
        return "placeholder"
    return "especie"


def tiene_sufijo_gtdb(especie):
    """True si el genero o el epiteto llevan sufijo GTDB (_A, _B, _AB...)."""
    return any(RE_SUFIJO.search(p) for p in (especie or "").split())


def limpiar_sufijos(especie):
    """'Bacillus_AB infantis' -> 'Bacillus infantis' (nombre NCBI probable)."""
    return " ".join(RE_SUFIJO.sub("", p) for p in (especie or "").split())


def orden_linaje(lid):
    """Ordena L1, L2, ... L19 numericamente y no lexicograficamente."""
    m = re.match(r"^L(\d+)$", lid or "")
    return (0, int(m.group(1))) if m else (1, 0)


# ------------------------------------------------------- lectura de la Fase 2

def linajes_desde_fase2(path):
    """Agrupa phase2_selection.tsv por linaje_id.

    Devuelve dict linaje_id -> info del linaje (clasificacion, especie GTDB,
    bins que lo componen y prevalencia espacial).
    """
    por_linaje = defaultdict(lambda: {"bins": [], "muestras": set()})
    for fila in leer_tsv(path):
        lid = (fila.get("linaje_id") or "").strip()
        if not lid:
            continue
        d = por_linaje[lid]
        d["bins"].append((fila.get("genome") or "").strip())
        d["muestras"].add((fila.get("muestra") or "").strip())
        d["clasificacion"] = (fila.get("clasificacion") or "").strip()
        d["ani"] = fila.get("ani") or ""
        d["af"] = fila.get("af") or ""
        d["ref_gtdb"] = fila.get("ref") or ""

    linajes = {}
    for lid, d in por_linaje.items():
        especie = rango_taxon(d["clasificacion"], "s")
        linajes[lid] = {
            "linaje_id": lid,
            "especie_gtdb": especie,
            "genero_gtdb": rango_taxon(d["clasificacion"], "g"),
            "filo_gtdb": rango_taxon(d["clasificacion"], "p"),
            "clasificacion_gtdb": d["clasificacion"],
            "nivel": clasificar_nivel(especie),
            "sufijo_gtdb": "si" if tiene_sufijo_gtdb(especie) else "no",
            "n_bins": len(d["bins"]),
            "bins": ",".join(sorted(d["bins"])),
            "prevalencia_espacial": len([m for m in d["muestras"] if m]),
            "muestras": ",".join(sorted(m for m in d["muestras"] if m)),
            "ani_gtdb": d["ani"],
            "af_gtdb": d["af"],
            "ref_gtdb": d["ref_gtdb"],
        }
    return linajes


# ------------------------------------------------ lectura del metadata de GTDB

def resolver_columnas(cabecera):
    """Mapea nombres logicos -> columna real del metadata, con alternativas.

    GTDB cambio nombres entre releases (checkm_ -> checkm2_), asi que se resuelve
    en tiempo de ejecucion en vez de asumir un esquema fijo.
    """
    presentes = set(c.strip() for c in cabecera)
    alternativas = {
        "accession":          ["accession"],
        "gtdb_taxonomy":      ["gtdb_taxonomy"],
        "ncbi_taxonomy":      ["ncbi_taxonomy"],
        "ncbi_species_taxid": ["ncbi_species_taxid"],
        "ncbi_taxid":         ["ncbi_taxid"],
        "ncbi_organism_name": ["ncbi_organism_name"],
        "completitud":        ["checkm2_completeness", "checkm_completeness"],
        "contaminacion":      ["checkm2_contamination", "checkm_contamination"],
        "contig_count":       ["contig_count"],
        "n50":                ["n50_contigs", "n50_scaffolds"],
        "genome_size":        ["genome_size"],
        "assembly_level":     ["ncbi_assembly_level"],
        "genome_category":    ["ncbi_genome_category"],
        "isolation_source":   ["ncbi_isolation_source"],
        "representative":     ["gtdb_representative"],
        "genbank_acc":        ["ncbi_genbank_assembly_accession"],
    }
    resueltas, faltantes = {}, []
    for logico, opciones in alternativas.items():
        for op in opciones:
            if op in presentes:
                resueltas[logico] = op
                break
        else:
            faltantes.append(logico)
    return resueltas, faltantes


def escanear_metadata(path, especies_objetivo):
    """Recorre el metadata de GTDB una sola vez y retiene solo las especies diana.

    El archivo pesa cientos de MB, asi que se procesa en streaming y se descarta
    todo lo que no pertenezca a un cluster de interes.
    """
    genomas = defaultdict(list)
    with open(path, encoding="utf-8", newline="") as fh:
        lector = csv.reader(fh, delimiter="\t")
        try:
            cabecera = next(lector)
        except StopIteration:
            die("metadata de GTDB vacio: %s" % path)
        cols, faltantes = resolver_columnas(cabecera)
        if "gtdb_taxonomy" in faltantes or "accession" in faltantes:
            die("el metadata no tiene 'accession'/'gtdb_taxonomy'; "
                "verifica que sea bac120_metadata_r*.tsv")
        if faltantes:
            warn("columnas ausentes en el metadata (se rellenan como NA): %s"
                 % ", ".join(faltantes))
        idx = {k: cabecera.index(v) for k, v in cols.items()}

        def val(fila, clave):
            i = idx.get(clave)
            return fila[i].strip() if i is not None and i < len(fila) else ""

        for fila in lector:
            if not fila:
                continue
            especie = rango_taxon(val(fila, "gtdb_taxonomy"), "s")
            if especie not in especies_objetivo:
                continue
            genomas[especie].append({
                "accession": val(fila, "accession"),
                "genbank_acc": val(fila, "genbank_acc"),
                "ncbi_species_taxid": val(fila, "ncbi_species_taxid"),
                "ncbi_taxid": val(fila, "ncbi_taxid"),
                "ncbi_organism_name": val(fila, "ncbi_organism_name"),
                "ncbi_taxonomy": val(fila, "ncbi_taxonomy"),
                "completitud": a_float(val(fila, "completitud")),
                "contaminacion": a_float(val(fila, "contaminacion")),
                "contigs": a_int(val(fila, "contig_count")),
                "n50": a_int(val(fila, "n50")),
                "genome_size": a_int(val(fila, "genome_size")),
                "assembly_level": val(fila, "assembly_level"),
                "genome_category": val(fila, "genome_category"),
                "isolation_source": val(fila, "isolation_source"),
                "representative": val(fila, "representative"),
            })
    return genomas


def es_refseq(g):
    """En GTDB el prefijo de la accesion codifica la fuente: RS_ = RefSeq."""
    return (g.get("accession") or "").startswith("RS_")


def es_aislado(g):
    cat = (g.get("genome_category") or "").strip().lower()
    if not cat or cat in ("none", "na"):
        return True
    return not any(cat.startswith(x) for x in NO_AISLADO)


def pasa_qc(g, min_comp, max_cont, max_contigs):
    """Criterios de la Etapa 3.1. Ausencia de metricas = no admitido (conservador)."""
    if not es_aislado(g):
        return False, "no_aislado"
    comp, cont, ctg = g["completitud"], g["contaminacion"], g["contigs"]
    if comp is None or cont is None:
        return False, "sin_metricas"
    if comp < min_comp:
        return False, "completitud_baja"
    if cont > max_cont:
        return False, "contaminacion_alta"
    if ctg is not None and ctg > max_contigs:
        return False, "fragmentado"
    return True, "ok"


def mapear_a_ncbi(lista_genomas):
    """Deriva especie NCBI + taxid por voto mayoritario dentro del cluster GTDB.

    Se vota sobre ncbi_species_taxid porque el nombre de organismo incluye cepa.
    Se reporta la fraccion de acuerdo: valores bajos indican un cluster GTDB que
    NCBI reparte entre varias especies (o al reves), y exigen revision manual.
    """
    votos = Counter(g["ncbi_species_taxid"] for g in lista_genomas
                    if g.get("ncbi_species_taxid") not in (None, "", "none", "NA"))
    if not votos:
        return {"ncbi_especie": "NA", "ncbi_taxid": "NA",
                "acuerdo_mapeo": "NA", "n_taxids_distintos": 0}
    taxid, n = votos.most_common(1)[0]
    nombres = Counter(rango_taxon(g["ncbi_taxonomy"], "s") for g in lista_genomas
                      if g.get("ncbi_species_taxid") == taxid
                      and rango_taxon(g["ncbi_taxonomy"], "s"))
    if nombres:
        nombre = nombres.most_common(1)[0][0]
    else:
        # Fallback: primeros dos tokens del nombre de organismo (genero + epiteto).
        crudos = [g["ncbi_organism_name"] for g in lista_genomas
                  if g.get("ncbi_species_taxid") == taxid and g.get("ncbi_organism_name")]
        nombre = " ".join(crudos[0].split()[:2]) if crudos else "NA"
    return {
        "ncbi_especie": nombre,
        "ncbi_taxid": taxid,
        "acuerdo_mapeo": "%.2f" % (float(n) / len(lista_genomas)),
        "n_taxids_distintos": len(votos),
    }


def resumir_habitats(lista_genomas, tope=5):
    """Agrega ncbi_isolation_source. Relevante para biodeterioro petreo: un
    pangenoma armado solo con aislados clinicos dice poco sobre roca."""
    vacios = ("none", "na", "n/a", "not applicable", "missing", "not collected",
              "unknown", "not provided", "-")
    c = Counter()
    for g in lista_genomas:
        src = (g.get("isolation_source") or "").strip().lower()
        if src and src not in vacios:
            c[src[:40]] += 1
    if not c:
        return "NA", 0
    return "; ".join("%s (%d)" % (k, v) for k, v in c.most_common(tope)), len(c)


# ------------------------------------------------------------ subcomando censo

def cmd_censo(args):
    linajes = linajes_desde_fase2(args.selection)
    if not linajes:
        die("no se leyeron linajes de %s" % args.selection)

    # Solo se censan los linajes con binomio real: los placeholder y los sin
    # especie no tienen conjunto de referencia publico contra el cual comparar.
    diana = {lid: d for lid, d in linajes.items() if d["nivel"] == "especie"}
    excluidos = {lid: d for lid, d in linajes.items() if d["nivel"] != "especie"}
    print("[info] linajes totales en Fase 2: %d" % len(linajes))
    print("[info] con especie asignada (a censar): %d -> %s"
          % (len(diana), ", ".join(sorted(diana, key=orden_linaje))))
    print("[info] excluidos (sin especie o placeholder GTDB): %d -> %s"
          % (len(excluidos), ", ".join(sorted(excluidos, key=orden_linaje))))

    especies = set(d["especie_gtdb"] for d in diana.values())
    print("[info] escaneando metadata de GTDB: %s" % args.gtdb_metadata)
    genomas = escanear_metadata(args.gtdb_metadata, especies)
    for e in sorted(especies):
        if e not in genomas:
            warn("sin genomas en el metadata para 's__%s' (revisa el release de GTDB)" % e)

    filas, candidatos = [], []
    for lid in sorted(diana, key=orden_linaje):
        d = dict(diana[lid])
        lista = genomas.get(d["especie_gtdb"], [])

        aprobados, motivos = [], Counter()
        for g in lista:
            ok, motivo = pasa_qc(g, args.min_completitud, args.max_contaminacion,
                                 args.max_contigs)
            motivos[motivo] += 1
            if ok:
                aprobados.append(g)

        d.update(mapear_a_ncbi(lista))
        d["ncbi_especie_esperada"] = limpiar_sufijos(d["especie_gtdb"])
        d["conteo_ncbi_comparable"] = "no" if d["sufijo_gtdb"] == "si" else "si"
        d["n_cluster_gtdb"] = len(lista)
        d["n_refseq_gtdb"] = sum(1 for g in lista if es_refseq(g))
        d["n_genbank_gtdb"] = len(lista) - d["n_refseq_gtdb"]
        d["n_aislados"] = sum(1 for g in lista if es_aislado(g))
        d["n_mags_excluidos"] = motivos["no_aislado"]
        d["n_postqc"] = len(aprobados)
        d["descartes_qc"] = "; ".join("%s=%d" % (k, v) for k, v in sorted(motivos.items())
                                      if k != "ok") or "ninguno"

        niveles = Counter(g["assembly_level"] or "NA" for g in aprobados)
        d["n_complete"] = niveles.get("Complete Genome", 0)
        d["n_chromosome"] = niveles.get("Chromosome", 0)
        d["n_scaffold"] = niveles.get("Scaffold", 0)
        d["n_contig"] = niveles.get("Contig", 0)

        comps = [g["completitud"] for g in aprobados if g["completitud"] is not None]
        conts = [g["contaminacion"] for g in aprobados if g["contaminacion"] is not None]
        ctgs = sorted(g["contigs"] for g in aprobados if g["contigs"] is not None)
        d["completitud_media"] = "%.2f" % (sum(comps) / len(comps)) if comps else "NA"
        d["contaminacion_media"] = "%.2f" % (sum(conts) / len(conts)) if conts else "NA"
        d["contigs_mediana"] = str(ctgs[len(ctgs) // 2]) if ctgs else "NA"
        d["habitats"], d["n_habitats_distintos"] = resumir_habitats(aprobados)
        filas.append(d)

        for g in aprobados:
            candidatos.append({
                "linaje_id": lid,
                "especie_gtdb": d["especie_gtdb"],
                "accession_gtdb": g["accession"],
                # 04_fetch_refs.sh consume esta columna: acc sin el prefijo RS_/GB_.
                "accession_ncbi": re.sub(r"^(RS_|GB_)", "", g["accession"]),
                "organismo": g["ncbi_organism_name"],
                "assembly_level": g["assembly_level"],
                "completitud": g["completitud"],
                "contaminacion": g["contaminacion"],
                "contigs": g["contigs"],
                "n50": g["n50"],
                "genome_size": g["genome_size"],
                "isolation_source": g["isolation_source"],
                "representante_gtdb": g["representative"],
                "fuente": "RefSeq" if es_refseq(g) else "GenBank",
            })

    os.makedirs(args.outdir, exist_ok=True)
    p_censo = os.path.join(args.outdir, "phase3_censo_genomas.tsv")
    p_cand = os.path.join(args.outdir, "phase3_genomas_candidatos.tsv")
    p_map = os.path.join(args.outdir, "phase3_mapeo_ncbi.tsv")

    escribir_tsv(p_censo, filas,
                 ["linaje_id", "especie_gtdb", "ncbi_especie", "ncbi_taxid",
                  "sufijo_gtdb", "conteo_ncbi_comparable", "acuerdo_mapeo",
                  "n_cluster_gtdb", "n_refseq_gtdb", "n_genbank_gtdb",
                  "n_aislados", "n_mags_excluidos", "n_postqc",
                  "n_complete", "n_chromosome", "n_scaffold", "n_contig",
                  "completitud_media", "contaminacion_media", "contigs_mediana",
                  "n_habitats_distintos", "habitats", "descartes_qc",
                  "prevalencia_espacial", "n_bins", "bins", "muestras",
                  "genero_gtdb", "filo_gtdb", "clasificacion_gtdb"])
    escribir_tsv(p_cand, candidatos,
                 ["linaje_id", "especie_gtdb", "accession_gtdb", "accession_ncbi",
                  "organismo", "assembly_level", "completitud", "contaminacion",
                  "contigs", "n50", "genome_size", "isolation_source",
                  "representante_gtdb", "fuente"])
    escribir_tsv(p_map, filas,
                 ["linaje_id", "especie_gtdb", "ncbi_especie_esperada",
                  "ncbi_especie", "ncbi_taxid", "acuerdo_mapeo",
                  "n_taxids_distintos", "sufijo_gtdb", "conteo_ncbi_comparable"])

    print("\n[ok] %s  (%d linajes)" % (p_censo, len(filas)))
    print("[ok] %s  (%d genomas candidatos)" % (p_cand, len(candidatos)))
    print("[ok] %s" % p_map)

    # Avisos de mapeo: nombre NCBI distinto al esperado, o cluster heterogeneo.
    for d in filas:
        if d["ncbi_especie"] not in ("NA", d["ncbi_especie_esperada"]):
            print("[mapeo] %s: GTDB '%s' -> NCBI '%s' (taxid %s) *** nombre distinto ***"
                  % (d["linaje_id"], d["especie_gtdb"], d["ncbi_especie"], d["ncbi_taxid"]))
        if a_float(d["acuerdo_mapeo"], 1.0) < 0.90:
            warn("%s: acuerdo de mapeo %s con %s taxids distintos -> revisar a mano"
                 % (d["linaje_id"], d["acuerdo_mapeo"], d["n_taxids_distintos"]))


# ------------------------------------------------------------- subcomando ncbi

def datasets_total(taxid, extra):
    """Devuelve total_count de 'datasets summary genome taxon', o None si falla."""
    cmd = ["datasets", "summary", "genome", "taxon", str(taxid), "--limit", "1"] + extra
    try:
        salida = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    except (OSError, subprocess.TimeoutExpired) as exc:
        warn("fallo '%s': %s" % (" ".join(cmd), exc))
        return None
    if salida.returncode != 0:
        warn("datasets devolvio %d para taxid %s: %s"
             % (salida.returncode, taxid, (salida.stderr or "").strip()[:200]))
        return None
    try:
        return json.loads(salida.stdout).get("total_count", 0)
    except (ValueError, AttributeError):
        warn("no se pudo parsear la salida de datasets para taxid %s" % taxid)
        return None


def cmd_ncbi(args):
    try:
        hay_cli = subprocess.run(["datasets", "--version"],
                                 capture_output=True, text=True).returncode == 0
    except OSError:
        hay_cli = False
    if not hay_cli:
        die("no se encontro el CLI 'datasets'. Instalalo en el nodo de LOGIN:\n"
            "  conda create -y -n ncbi-datasets --override-channels "
            "-c conda-forge -c bioconda python=3.11 ncbi-datasets-cli\n"
            "Recuerda que solo el nodo de login de Khipu tiene internet.")

    p_censo = os.path.join(args.outdir, "phase3_censo_genomas.tsv")
    if not os.path.exists(p_censo):
        die("falta %s; corre primero el subcomando 'censo'." % p_censo)

    filas = []
    for d in leer_tsv(p_censo):
        taxid = (d.get("ncbi_taxid") or "NA").strip()
        if taxid in ("NA", ""):
            warn("%s sin taxid; se omite" % d["linaje_id"])
            continue
        print("[ncbi] %s  taxid %s (%s)" % (d["linaje_id"], taxid, d["ncbi_especie"]))
        filas.append({
            "linaje_id": d["linaje_id"],
            "especie_gtdb": d["especie_gtdb"],
            "ncbi_especie": d["ncbi_especie"],
            "ncbi_taxid": taxid,
            "ncbi_genbank_total": datasets_total(taxid, ["--assembly-source", "genbank"]),
            "ncbi_refseq_total": datasets_total(taxid, ["--assembly-source", "refseq"]),
            "ncbi_refseq_completos": datasets_total(
                taxid, ["--assembly-source", "refseq",
                        "--assembly-level", "complete,chromosome"]),
            "conteo_ncbi_comparable": d.get("conteo_ncbi_comparable", "si"),
            "consultado": datetime.now().strftime("%Y-%m-%d"),
        })

    p = os.path.join(args.outdir, "phase3_ncbi_counts.tsv")
    escribir_tsv(p, filas, ["linaje_id", "especie_gtdb", "ncbi_especie", "ncbi_taxid",
                            "ncbi_genbank_total", "ncbi_refseq_total",
                            "ncbi_refseq_completos", "conteo_ncbi_comparable",
                            "consultado"])
    print("\n[ok] %s  (%d linajes consultados)" % (p, len(filas)))
    print("[nota] en linajes con sufijo GTDB estos totales estan INFLADOS: "
          "NCBI no separa las especies que GTDB escinde.")


# ------------------------------------------------------------ subcomando tabla

def cmd_tabla(args):
    p_censo = os.path.join(args.outdir, "phase3_censo_genomas.tsv")
    if not os.path.exists(p_censo):
        die("falta %s; corre primero el subcomando 'censo'." % p_censo)
    censo = leer_tsv(p_censo)

    p_ncbi = os.path.join(args.outdir, "phase3_ncbi_counts.tsv")
    ncbi = {d["linaje_id"]: d for d in leer_tsv(p_ncbi)} if os.path.exists(p_ncbi) else {}
    if not ncbi:
        warn("sin phase3_ncbi_counts.tsv; la tabla usa solo conteos de GTDB. "
             "Corre el subcomando 'ncbi' en el nodo de login para completarla.")

    filas = []
    for d in censo:
        n = a_int(d["n_postqc"], 0)
        e = ncbi.get(d["linaje_id"], {})
        filas.append({
            "linaje_id": d["linaje_id"],
            "especie_gtdb": d["especie_gtdb"],
            "ncbi_especie": d["ncbi_especie"],
            "ncbi_taxid": d["ncbi_taxid"],
            "sufijo_gtdb": d["sufijo_gtdb"],
            "conteo_ncbi_comparable": d["conteo_ncbi_comparable"],
            "n_genbank": e.get("ncbi_genbank_total", d["n_genbank_gtdb"]),
            "n_refseq": e.get("ncbi_refseq_total", d["n_refseq_gtdb"]),
            "n_cluster_gtdb": d["n_cluster_gtdb"],
            "n_postqc": n,
            "n_objetivo_pangenoma": min(n, args.tope),
            "completitud_media": d["completitud_media"],
            "contaminacion_media": d["contaminacion_media"],
            "n_habitats_distintos": d["n_habitats_distintos"],
            "habitats": d["habitats"],
            "prevalencia_espacial": d["prevalencia_espacial"],
            "n_bins_propios": d["n_bins"],
            "viabilidad_pangenoma": "Si" if n >= args.umbral else "No",
            "_n": n,
        })

    # Ranking por disponibilidad depurada. El orden es PRELIMINAR: el definitivo
    # se fija tras la desreplicacion a 99 % ANI en 04_fetch_refs.sh, porque el
    # conteo bruto premia a las especies clinicas/industriales sobresecuenciadas.
    viables = sorted([f for f in filas if f["viabilidad_pangenoma"] == "Si"],
                     key=lambda f: -f["_n"])
    for i, f in enumerate(viables, 1):
        f["ranking_disponibilidad"] = i
        f["seleccion_fase3"] = "Si" if i <= args.top else "No"
    for f in filas:
        f.setdefault("ranking_disponibilidad", "NA")
        f.setdefault("seleccion_fase3", "No")
        f.pop("_n", None)

    filas.sort(key=lambda f: (1, 0) + orden_linaje(f["linaje_id"])
               if f["ranking_disponibilidad"] == "NA"
               else (0, f["ranking_disponibilidad"], 0, 0))

    campos = ["ranking_disponibilidad", "linaje_id", "especie_gtdb", "ncbi_especie",
              "ncbi_taxid", "sufijo_gtdb", "conteo_ncbi_comparable",
              "n_genbank", "n_refseq", "n_cluster_gtdb", "n_postqc",
              "n_objetivo_pangenoma", "completitud_media", "contaminacion_media",
              "n_habitats_distintos", "habitats", "prevalencia_espacial",
              "n_bins_propios", "viabilidad_pangenoma", "seleccion_fase3"]
    p_tabla = os.path.join(args.outdir, "phase3_viabilidad.tsv")
    escribir_tsv(p_tabla, filas, campos)

    seleccionadas = [f for f in filas if f["seleccion_fase3"] == "Si"]
    p_md = os.path.join(args.outdir, "phase3_censo_report.md")
    with open(p_md, "w", encoding="utf-8", newline="") as fh:
        w = fh.write
        w("# Fase 3 (OE3) - Censo de disponibilidad genomica y seleccion de especies diana\n\n")
        w("Generado: %s\n\n" % datetime.now().strftime("%Y-%m-%d %H:%M"))
        w("## Criterio de seleccion\n\n")
        w("Disponibilidad genomica publica. Umbral de viabilidad: **>= %d genomas** de\n"
          "calidad adecuada, aplicado sobre el conteo **post-QC** y no sobre el conteo\n"
          "bruto de NCBI.\n\n" % args.umbral)
        w("Admision de cada genoma publico al set de referencia:\n\n")
        w("- completitud >= %.1f %%, contaminacion <= %.1f %%\n"
          % (args.min_completitud, args.max_contaminacion))
        w("- <= %d contigs (la fragmentacion parte genes e infla el genoma accesorio)\n"
          % args.max_contigs)
        w("- aislados unicamente; se excluyen MAGs y genomas de muestra ambiental\n")
        w("- tope de %d genomas por especie (computo y sesgo por clonalidad)\n\n" % args.tope)
        w("## Tabla de viabilidad\n\n")
        enc = ["#", "Linaje", "Especie GTDB", "Especie NCBI / TaxID", "GenBank",
               "RefSeq", "Cluster GTDB", "Post-QC", "Prev. esp.", "Viable"]
        w("| " + " | ".join(enc) + " |\n")
        w("|" + "---|" * len(enc) + "\n")
        for f in filas:
            w("| %s | %s | *%s* | %s / %s | %s | %s | %s | %s | %s | **%s** |\n" % (
                f["ranking_disponibilidad"], f["linaje_id"], f["especie_gtdb"],
                f["ncbi_especie"], f["ncbi_taxid"], f["n_genbank"], f["n_refseq"],
                f["n_cluster_gtdb"], f["n_postqc"], f["prevalencia_espacial"],
                f["viabilidad_pangenoma"]))
        w("\n## Especies seleccionadas para el analisis pangenomico\n\n")
        if seleccionadas:
            for f in seleccionadas:
                w("%s. **%s** (%s) - %s genomas post-QC, %s para el pangenoma; "
                  "%s habitat(s) distinto(s); prevalencia espacial %s; %s bin(s) propio(s).\n"
                  % (f["ranking_disponibilidad"], f["especie_gtdb"], f["linaje_id"],
                     f["n_postqc"], f["n_objetivo_pangenoma"],
                     f["n_habitats_distintos"], f["prevalencia_espacial"],
                     f["n_bins_propios"]))
        else:
            w("Ningun linaje alcanzo el umbral de %d genomas.\n" % args.umbral)
        w("\n## Advertencias para la redaccion\n\n")
        conflicto = [f for f in filas if f["viabilidad_pangenoma"] == "No"
                     and a_int(f["prevalencia_espacial"], 0) > 1]
        if conflicto:
            w("- Linajes con prevalencia espacial > 1 que NO alcanzan el umbral: %s.\n"
              "  Se excluyen del pangenoma por factibilidad analitica, no por falta de\n"
              "  relevancia ecologica; deben reportarse igualmente como resultado de OE2.\n"
              % ", ".join("%s (%s)" % (f["linaje_id"], f["especie_gtdb"]) for f in conflicto))
        incomp = [f for f in filas if f["conteo_ncbi_comparable"] == "no"]
        if incomp:
            w("- Conteos NCBI NO comparables por sufijo GTDB (NCBI no separa la especie\n"
              "  escindida, el total esta inflado): %s. Usar el tamano del cluster GTDB.\n"
              % ", ".join("%s (%s)" % (f["linaje_id"], f["especie_gtdb"]) for f in incomp))
        w("- El ranking es PRELIMINAR: el orden definitivo se fija tras la desreplicacion\n"
          "  a 99 % ANI (04_fetch_refs.sh), que mide diversidad de cepas no redundante.\n")
        w("- Los bins propios son MAGs (70-100 %% de completitud) frente a aislados al\n"
          "  >= %.0f %%. Las conclusiones sobre PRESENCIA de genes son validas; las de\n"
          "  AUSENCIA no lo son sin controlar por la completitud del MAG.\n"
          % args.min_completitud)

    print("[ok] %s" % p_tabla)
    print("[ok] %s" % p_md)
    print("\nViables (>= %d genomas post-QC): %d de %d"
          % (args.umbral, len(viables), len(filas)))
    for f in seleccionadas:
        print("  %s. %-42s %s post-QC" % (f["ranking_disponibilidad"],
                                          f["especie_gtdb"], f["n_postqc"]))


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Fase 3: censo de disponibilidad genomica y seleccion de especies diana.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Uso completo en Khipu:\n"
               "  python 03_censo_genomas.py censo \\\n"
               "      --selection results/fase2/phase2_selection.tsv \\\n"
               "      --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv \\\n"
               "      --outdir results/fase3\n"
               "  python 03_censo_genomas.py ncbi  --outdir results/fase3   # nodo de LOGIN\n"
               "  python 03_censo_genomas.py tabla --outdir results/fase3\n")
    sub = ap.add_subparsers(dest="cmd")

    def comunes(p):
        p.add_argument("--outdir", default="results/fase3",
                       help="directorio de salida (def: results/fase3)")

    def qc(p):
        p.add_argument("--min-completitud", type=float, default=DEF_MIN_COMPLETITUD)
        p.add_argument("--max-contaminacion", type=float, default=DEF_MAX_CONTAMINACION)
        p.add_argument("--max-contigs", type=int, default=DEF_MAX_CONTIGS)

    p1 = sub.add_parser("censo", help="cruza Fase 2 con el metadata de GTDB (sin internet)")
    comunes(p1)
    qc(p1)
    p1.add_argument("--selection", default="results/fase2/phase2_selection.tsv",
                    help="phase2_selection.tsv de la Fase 2")
    p1.add_argument("--gtdb-metadata", required=True,
                    help="bac120_metadata_r220.tsv del release de GTDB usado en Fase 2")
    p1.set_defaults(func=cmd_censo)

    p2 = sub.add_parser("ncbi", help="conteos vigentes via NCBI Datasets (requiere internet)")
    comunes(p2)
    p2.set_defaults(func=cmd_ncbi)

    p3 = sub.add_parser("tabla", help="tabla de viabilidad + reporte markdown")
    comunes(p3)
    qc(p3)
    p3.add_argument("--umbral", type=int, default=DEF_UMBRAL,
                    help="minimo de genomas post-QC para viabilidad (def: %d)" % DEF_UMBRAL)
    p3.add_argument("--tope", type=int, default=DEF_TOPE,
                    help="maximo de genomas por pangenoma (def: %d)" % DEF_TOPE)
    p3.add_argument("--top", type=int, default=3,
                    help="cuantas especies seleccionar para la Fase 3 (def: 3)")
    p3.set_defaults(func=cmd_tabla)

    args = ap.parse_args()
    if not getattr(args, "func", None):
        ap.print_help()
        sys.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
