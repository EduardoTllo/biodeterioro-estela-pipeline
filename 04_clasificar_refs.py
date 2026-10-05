#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_clasificar_refs.py
---------------------
FASE 3 (OE3), Etapa 3.2b: genomas de referencia. Toma el censo de la etapa
3.0-3.2 (03_censo_genomas.py), verifica la descarga de NCBI, clasifica cada
genoma por prioridad, habitat y continente, prepara la desreplicacion con dRep
al 99 % ANI, fija el ranking definitivo de especies (top-3 por genomas no
redundantes) y selecciona hasta 50 referencias por especie.

Solo biblioteca estandar (sin pandas).

Subcomandos, en el orden en que se usan:
  listas       Lista de especies viables y archivos acc_<especie>.txt con las
               accesiones a descargar (entrada de 04_fetch_refs.sh).
  verificar    Tras la descarga: controla que cada FASTA exista y que su
               longitud total coincida con genome_size de GTDB (+/- 1 %).
               Aplica las exclusiones manuales (--exclusiones) y recalcula la
               viabilidad con los genomas que quedan.
  clasificar   Nivel de prioridad (T0-T2), categoria de habitat y continente de
               cada genoma. Prepara la entrada de dRep (lista de genomas,
               genomeInfo con CheckM2 del metadata y tabla de pesos extra).
  seleccionar  Tras dRep: numero de genomas no redundantes por especie,
               ranking definitivo, top-N y seleccion de hasta 50 referencias.
               Emite el manifiesto de Bakta (referencias + bins propios).

Niveles de prioridad (decisiones D23 y D31):
  T0  representante de especie de GTDB o cepa tipo de la especie (siempre entra)
  T1  aislado de sustrato petreo o ambiente arido (ncbi_isolation_source)
  T2  el resto; completa los cupos por rotacion entre habitats y continentes
No hay preferencia por region: el pais solo se usa para la rotacion entre
continentes, que recorre los continentes en orden alfabetico.
Dentro de cada cluster de dRep al 99 % gana el genoma prioritario gracias a la
tabla de pesos extra (-extraW): T0 +1000, T1 +500, T2 0.
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict, deque
from datetime import datetime


# ---------------------------------------------------------------- parametros

DEF_UMBRAL = 15          # minimo de genomas no redundantes (D11)
DEF_TOPE = 50            # maximo de referencias por especie (D12, D23)
DEF_TOP = 3              # especies que pasan al pangenoma (D10, D25)
DEF_TOLERANCIA = 0.01    # diferencia relativa maxima de longitud vs GTDB (D29)

PESO_T0 = 1000
PESO_T1 = 500

# Orden fijo de la rotacion de T2 entre categorias de habitat: primero los
# ambientes mas cercanos al nicho del monumento.
ORDEN_HABITAT = ["suelo", "agua_sedimento", "planta", "otro_ambiental",
                 "petreo_arido", "animal", "alimento_industrial",
                 "clinico_humano", "sin_clasificar", "desconocido"]
# Continentes en orden alfabetico (sin preferencia regional, D31); los genomas
# sin pais van al final.
ORDEN_CONTINENTE = ["Africa", "America_Latina_Caribe", "America_del_Norte",
                    "Antartida", "Asia", "Europa", "Oceania", "Oceano",
                    "sin_mapear", "desconocido"]
ORDEN_ENSAMBLAJE = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}

VACIOS = ("", "none", "na", "n/a", "nan", "not applicable", "missing",
          "not collected", "unknown", "not provided", "not available", "-",
          "not determined", "unspecified", "null")


# ---------------------------------------------------------------- utilidades

def warn(msg):
    sys.stderr.write("[aviso] %s\n" % msg)


def die(msg):
    sys.stderr.write("[error] %s\n" % msg)
    sys.exit(1)


def leer_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def leer_csv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


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


def orden_linaje(lid):
    m = re.match(r"^L(\d+)$", lid or "")
    return int(m.group(1)) if m else 10 ** 6


def slug_especie(linaje_id, especie):
    """'L14', 'Bacillus licheniformis' -> 'L14_Bacillus_licheniformis'."""
    return re.sub(r"[^A-Za-z0-9]+", "_", "%s_%s" % (linaje_id, especie)).strip("_")


def id_genoma(accesion):
    """'GCF_000011645.1' -> 'GCF_000011645_1'. Solo [A-Za-z0-9_]: el ID se usa
    como nombre de archivo, prefijo y locus tag de Bakta y nombre de genoma en
    Panaroo."""
    return re.sub(r"[^A-Za-z0-9]", "_", accesion)


def id_bin(nombre):
    """'bin-6-50' -> 'bin_6_50'."""
    return re.sub(r"[^A-Za-z0-9]", "_", nombre)


def vacio(valor):
    v = (valor or "").strip().lower()
    return v in VACIOS or v.startswith("missing") or v.startswith("not ")


def longitud_fasta(path):
    """Longitud total (pb) de un FASTA, sin cargarlo en memoria."""
    total = 0
    with open(path, encoding="ascii", errors="replace") as fh:
        for linea in fh:
            if not linea.startswith(">"):
                total += len(linea.strip())
    return total


def leer_candidatos(censo_dir):
    p = os.path.join(censo_dir, "phase3_genomas_candidatos.tsv")
    if not os.path.exists(p):
        die("falta %s; corre primero 03_censo_genomas.py censo." % p)
    por_linaje = defaultdict(list)
    for c in leer_tsv(p):
        por_linaje[c["linaje_id"]].append(c)
    return por_linaje


def leer_especies(path):
    if not os.path.exists(path):
        die("falta %s" % path)
    return leer_tsv(path)


# ------------------------------------------------------------ subcomando listas

def cmd_listas(args):
    p_via = os.path.join(args.censo_dir, "phase3_viabilidad.tsv")
    p_censo = os.path.join(args.censo_dir, "phase3_censo_genomas.tsv")
    for p in (p_via, p_censo):
        if not os.path.exists(p):
            die("falta %s; corre 03_censo_genomas.py censo y tabla." % p)
    censo = {d["linaje_id"]: d for d in leer_tsv(p_censo)}
    candidatos = leer_candidatos(args.censo_dir)

    especies = []
    for v in leer_tsv(p_via):
        if v.get("viabilidad_pangenoma") != "Si" or v.get("seleccion_fase3") != "Si":
            continue
        lid = v["linaje_id"]
        c = censo.get(lid, {})
        slug = slug_especie(lid, v["especie_gtdb"])
        accs = [x["accession_ncbi"] for x in candidatos.get(lid, [])]
        if not accs:
            warn("%s sin accesiones candidatas; se omite" % lid)
            continue
        p_acc = os.path.join(args.outdir, "acc_%s.txt" % slug)
        os.makedirs(args.outdir, exist_ok=True)
        with open(p_acc, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(accs) + "\n")
        especies.append({
            "slug": slug,
            "linaje_id": lid,
            "especie_gtdb": v["especie_gtdb"],
            "ncbi_especie": v.get("ncbi_especie", "NA"),
            "ncbi_taxid": v.get("ncbi_taxid", "NA"),
            "n_postqc": len(accs),
            "bins": c.get("bins", ""),
            "n_bins": c.get("n_bins", ""),
        })
        print("[ok] %-45s %4d accesiones -> %s" % (slug, len(accs), p_acc))

    if not especies:
        die("ninguna especie viable en %s: revisa el umbral con el asesor (D25)." % p_via)
    escribir_tsv(os.path.join(args.outdir, "especies_viables.tsv"), especies,
                 ["slug", "linaje_id", "especie_gtdb", "ncbi_especie", "ncbi_taxid",
                  "n_postqc", "bins", "n_bins"])
    print("[ok] %d especies viables -> %s" % (len(especies),
                                              os.path.join(args.outdir, "especies_viables.tsv")))


# --------------------------------------------------------- subcomando verificar

def leer_exclusiones(path):
    """accesion -> motivo. TSV con columnas accession, linaje_id, motivo."""
    if not path:
        return {}
    if not os.path.exists(path):
        die("no existe el archivo de exclusiones %s" % path)
    excl = {}
    for f in leer_tsv(path):
        acc = (f.get("accession") or "").strip()
        if acc and not acc.startswith("#"):
            excl[acc] = (f.get("motivo") or "").strip() or "sin motivo"
    return excl


def cmd_verificar(args):
    especies = leer_especies(os.path.join(args.descarga_dir, "especies_viables.tsv"))
    candidatos = leer_candidatos(args.censo_dir)
    exclusiones = leer_exclusiones(args.exclusiones)
    usadas = set()

    estado_filas, resumen, para_drep = [], [], []
    for e in especies:
        slug, lid = e["slug"], e["linaje_id"]
        carpeta = os.path.join(args.refs_dir, slug)
        conteo = Counter()
        ok_paths = []
        for c in candidatos.get(lid, []):
            acc = c["accession_ncbi"]
            gid = id_genoma(acc)
            path = os.path.join(carpeta, gid + ".fna")
            esperado = a_int(c.get("genome_size"))
            fila = {"slug": slug, "linaje_id": lid, "accession": acc, "id": gid,
                    "longitud_fasta": "NA", "genome_size_gtdb": esperado or "NA",
                    "dif_relativa": "NA", "motivo_exclusion": ""}
            if acc in exclusiones:
                # Exclusion manual documentada (metadata/exclusiones.tsv): el
                # genoma no entra al dRep, al pangenoma ni a la BD del filtro F4.
                fila["estado"] = "excluido"
                fila["motivo_exclusion"] = exclusiones[acc]
                usadas.add(acc)
            elif not os.path.exists(path):
                base = id_genoma(acc.rsplit(".", 1)[0])
                otras = glob.glob(os.path.join(carpeta, base + "_*.fna"))
                fila["estado"] = "version_distinta" if otras else "no_descargado"
            else:
                n = longitud_fasta(path)
                fila["longitud_fasta"] = n
                if n == 0:
                    fila["estado"] = "vacio"
                elif not esperado:
                    fila["estado"] = "sin_tamano_gtdb"
                else:
                    dif = abs(n - esperado) / float(esperado)
                    fila["dif_relativa"] = "%.4f" % dif
                    fila["estado"] = "ok" if dif <= args.tolerancia else "tamano_discrepante"
            conteo[fila["estado"]] += 1
            if fila["estado"] == "ok":
                ok_paths.append(os.path.abspath(path))
            estado_filas.append(fila)

        viable = len(ok_paths) >= args.umbral
        resumen.append({
            "slug": slug, "linaje_id": lid, "especie_gtdb": e["especie_gtdb"],
            "solicitados": sum(conteo.values()), "ok": conteo["ok"],
            "no_descargado": conteo["no_descargado"],
            "version_distinta": conteo["version_distinta"],
            "tamano_discrepante": conteo["tamano_discrepante"],
            "vacio": conteo["vacio"], "sin_tamano_gtdb": conteo["sin_tamano_gtdb"],
            "excluido": conteo["excluido"],
            "viable_tras_descarga": "si" if viable else "no",
        })
        with open(os.path.join(args.descarga_dir, "genomas_ok_%s.txt" % slug), "w",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(ok_paths) + ("\n" if ok_paths else ""))
        if viable:
            para_drep.append({"slug": slug, "linaje_id": lid,
                              "especie_gtdb": e["especie_gtdb"], "n_ok": len(ok_paths)})
        print("[%s] %-45s ok=%d de %d%s" % ("ok" if viable else "--", slug, len(ok_paths),
                                            sum(conteo.values()),
                                            "" if viable else "  (bajo el umbral: deja de ser viable)"))

    escribir_tsv(os.path.join(args.descarga_dir, "descarga_estado.tsv"), estado_filas,
                 ["slug", "linaje_id", "accession", "id", "estado", "longitud_fasta",
                  "genome_size_gtdb", "dif_relativa", "motivo_exclusion"])
    escribir_tsv(os.path.join(args.descarga_dir, "descarga_resumen.tsv"), resumen,
                 ["slug", "linaje_id", "especie_gtdb", "solicitados", "ok",
                  "no_descargado", "version_distinta", "tamano_discrepante", "vacio",
                  "sin_tamano_gtdb", "excluido", "viable_tras_descarga"])
    for acc in sorted(set(exclusiones) - usadas):
        warn("la exclusion %s no corresponde a ningun candidato de las especies viables" % acc)
    if exclusiones:
        print("[info] genomas excluidos a mano: %d" % len(usadas))
    escribir_tsv(os.path.join(args.descarga_dir, "especies_para_drep.tsv"), para_drep,
                 ["slug", "linaje_id", "especie_gtdb", "n_ok"])
    print("\n[ok] %d especies pasan al dRep -> %s"
          % (len(para_drep), os.path.join(args.descarga_dir, "especies_para_drep.tsv")))
    if not para_drep:
        die("ninguna especie conserva >= %d genomas tras la descarga (D25)." % args.umbral)


# -------------------------------------------------------- subcomando clasificar

def cargar_keywords(path):
    reglas = []
    for fila in leer_tsv(path):
        cat = (fila.get("categoria") or "").strip()
        pat = (fila.get("patron") or "").strip()
        if not cat or not pat or cat.startswith("#"):
            continue
        try:
            reglas.append((cat, pat, re.compile(pat, re.IGNORECASE)))
        except re.error as exc:
            die("patron invalido en %s: %r (%s)" % (path, pat, exc))
    if not reglas:
        die("sin reglas en %s" % path)
    return reglas


def cargar_paises(path):
    mapa = {}
    for fila in leer_tsv(path):
        pais = (fila.get("pais") or "").strip().lower()
        if pais:
            mapa[pais] = ((fila.get("continente") or "").strip(),
                          (fila.get("latam") or "no").strip().lower() == "si")
    return mapa


def clasificar_habitat(fuente, reglas):
    if vacio(fuente):
        return "desconocido", ""
    for cat, pat, rx in reglas:
        if rx.search(fuente):
            return cat, pat
    return "sin_clasificar", ""


def clasificar_pais(pais_raw, mapa):
    if vacio(pais_raw):
        return "desconocido", "desconocido", False
    pais = pais_raw.split(":")[0].strip()
    info = mapa.get(pais.lower())
    if not info:
        return pais, "sin_mapear", False
    return pais, info[0], info[1]


def cmd_clasificar(args):
    para_drep = leer_especies(os.path.join(args.descarga_dir, "especies_para_drep.tsv"))
    candidatos = leer_candidatos(args.censo_dir)
    estado = {(f["slug"], f["id"]): f["estado"]
              for f in leer_tsv(os.path.join(args.descarga_dir, "descarga_estado.tsv"))}
    reglas = cargar_keywords(args.keywords)
    paises = cargar_paises(args.paises)

    filas, fuentes, sin_mapear = [], defaultdict(Counter), Counter()
    for e in para_drep:
        slug, lid = e["slug"], e["linaje_id"]
        dir_in = os.path.join(args.outdir, slug, "input")
        os.makedirs(dir_in, exist_ok=True)
        genomas, ginfo, extraw = [], [], []
        for c in candidatos.get(lid, []):
            gid = id_genoma(c["accession_ncbi"])
            if estado.get((slug, gid)) != "ok":
                continue
            fuente = (c.get("isolation_source") or "").strip()
            habitat, patron = clasificar_habitat(fuente, reglas)
            pais, continente, latam = clasificar_pais(c.get("ncbi_country"), paises)
            if continente == "sin_mapear":
                sin_mapear[pais] += 1
            es_rep = (c.get("representante_gtdb") or "").strip().lower() in ("t", "true")
            es_tipo = (c.get("tipo_designacion") or "").strip().lower() == "type strain of species"
            if es_rep or es_tipo:
                nivel, peso = "T0", PESO_T0
            elif habitat == "petreo_arido":
                nivel, peso = "T1", PESO_T1
            else:
                nivel, peso = "T2", 0
            fuentes[fuente or "(vacio)"][(habitat, patron)] += 1
            filas.append({
                "slug": slug, "linaje_id": lid, "id": gid,
                "accession": c["accession_ncbi"], "nivel": nivel, "peso_drep": peso,
                "representante_gtdb": "si" if es_rep else "no",
                "cepa_tipo": "si" if es_tipo else "no",
                "habitat": habitat, "patron_habitat": patron,
                "isolation_source": fuente, "pais": pais, "continente": continente,
                "latam": "si" if latam else "no",
                "assembly_level": c.get("assembly_level", ""),
                "contigs": c.get("contigs", ""), "n50": c.get("n50", ""),
                "completitud": c.get("completitud", ""),
                "contaminacion": c.get("contaminacion", ""),
                "genome_size": c.get("genome_size", ""), "fuente": c.get("fuente", ""),
            })
            fasta = os.path.abspath(os.path.join(args.refs_dir, slug, gid + ".fna"))
            genomas.append(fasta)
            ginfo.append({"genome": gid + ".fna", "completeness": c.get("completitud", ""),
                          "contamination": c.get("contaminacion", "")})
            extraw.append((gid + ".fna", peso))

        with open(os.path.join(dir_in, "genomas.txt"), "w", encoding="utf-8",
                  newline="\n") as fh:
            fh.write("\n".join(genomas) + "\n")
        with open(os.path.join(dir_in, "genomeInfo.csv"), "w", encoding="utf-8",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["genome", "completeness", "contamination"],
                               lineterminator="\n")
            w.writeheader()
            for g in ginfo:
                w.writerow(g)
        # Todos los genomas van en la tabla (T2 con peso 0): si dRep no pudiera
        # leerla, seguiria sin pesos y sin avisar; seleccionar lo verifica en Sdb.
        with open(os.path.join(dir_in, "extraW.tsv"), "w", encoding="utf-8",
                  newline="\n") as fh:
            for g, p in extraw:
                fh.write("%s\t%d\n" % (g, p))
        n_niv = Counter(f["nivel"] for f in filas if f["slug"] == slug)
        print("[ok] %-45s %4d genomas  (%s)" % (slug, len(genomas),
              ", ".join("%s=%d" % kv for kv in sorted(n_niv.items()))))

    campos = ["slug", "linaje_id", "id", "accession", "nivel", "peso_drep",
              "representante_gtdb", "cepa_tipo", "habitat", "patron_habitat",
              "isolation_source", "pais", "continente", "latam", "assembly_level",
              "contigs", "n50", "completitud", "contaminacion", "genome_size", "fuente"]
    escribir_tsv(os.path.join(args.outdir, "refs_clasificadas.tsv"), filas, campos)

    fu = []
    for fuente, c in sorted(fuentes.items(), key=lambda kv: -sum(kv[1].values())):
        for (hab, pat), n in c.items():
            fu.append({"isolation_source": fuente, "n_genomas": n,
                       "habitat": hab, "patron": pat})
    escribir_tsv(os.path.join(args.outdir, "fuentes_unicas.tsv"), fu,
                 ["isolation_source", "n_genomas", "habitat", "patron"])
    escribir_tsv(os.path.join(args.outdir, "paises_sin_mapear.tsv"),
                 [{"pais": p, "n_genomas": n} for p, n in sin_mapear.most_common()],
                 ["pais", "n_genomas"])

    hab = Counter(f["habitat"] for f in filas)
    print("\nHabitats: " + ", ".join("%s=%d" % kv for kv in hab.most_common()))
    if hab.get("sin_clasificar"):
        warn("%d genomas con fuente sin clasificar: revisa fuentes_unicas.tsv y amplia "
             "metadata/habitat_keywords.tsv (checkpoint C3)." % hab["sin_clasificar"])
    if sin_mapear:
        warn("%d paises sin mapear: agregalos a metadata/paises_continentes.tsv "
             "(ver paises_sin_mapear.tsv)." % len(sin_mapear))


# ------------------------------------------------------- subcomando seleccionar

def clave_calidad(r):
    return (ORDEN_ENSAMBLAJE.get(r.get("assembly_level"), 9),
            a_int(r.get("contigs"), 10 ** 6),
            -(a_float(r.get("completitud"), 0.0)),
            a_float(r.get("contaminacion"), 100.0),
            r["id"])


def seleccionar_tope(ganadores, tope):
    """Aplica D23/D31 sobre los ganadores de dRep de una especie."""
    if len(ganadores) <= tope:
        return [(r, "todos (N_nr <= tope)") for r in sorted(ganadores, key=clave_calidad)]

    elegidos = []
    t0 = sorted([r for r in ganadores if r["nivel"] == "T0"], key=clave_calidad)
    for r in t0:
        elegidos.append((r, "T0 representante/cepa tipo"))
    t1 = sorted([r for r in ganadores if r["nivel"] == "T1"], key=clave_calidad)
    for r in t1:
        if len(elegidos) >= tope:
            break
        elegidos.append((r, "T1 petreo/arido"))

    # T2: rotacion por habitat y, dentro de cada habitat, por continente.
    celdas = defaultdict(lambda: defaultdict(deque))
    for r in sorted([r for r in ganadores if r["nivel"] == "T2"], key=clave_calidad):
        celdas[r["habitat"]][r["continente"]].append(r)
    orden_hab = ORDEN_HABITAT + sorted(h for h in celdas if h not in ORDEN_HABITAT)
    puntero = Counter()
    while len(elegidos) < tope and any(q for h in celdas.values() for q in h.values()):
        for hab in orden_hab:
            if len(elegidos) >= tope:
                break
            conts = [c for c in ORDEN_CONTINENTE + sorted(celdas[hab]) if celdas[hab].get(c)]
            conts = list(OrderedDict.fromkeys(conts))
            if not conts:
                continue
            c = conts[puntero[hab] % len(conts)]
            puntero[hab] += 1
            elegidos.append((celdas[hab][c].popleft(), "T2 rotacion habitat/continente"))
    return elegidos


def completitud_bins(path):
    comp = {}
    if path and os.path.exists(path):
        for f in leer_tsv(path):
            comp[f.get("Name", "")] = a_float(f.get("Completeness"), 0.0)
    elif path:
        warn("no se encontro el reporte de CheckM2 %s; el desempate por "
             "completitud del bin queda en 0" % path)
    return comp


def cmd_seleccionar(args):
    para_drep = leer_especies(os.path.join(args.descarga_dir, "especies_para_drep.tsv"))
    viables_info = {e["slug"]: e for e in leer_especies(
        os.path.join(args.descarga_dir, "especies_viables.tsv"))}
    clasif = defaultdict(dict)
    for r in leer_tsv(os.path.join(args.drep_dir, "refs_clasificadas.tsv")):
        clasif[r["slug"]][r["id"]] = r
    comp_bins = completitud_bins(args.checkm2)

    ranking = []
    for e in para_drep:
        slug = e["slug"]
        tablas = os.path.join(args.drep_dir, slug, "data_tables")
        p_w, p_c, p_s = (os.path.join(tablas, x) for x in ("Wdb.csv", "Cdb.csv", "Sdb.csv"))
        if not os.path.exists(p_w):
            die("falta %s: el dRep de %s no termino (04_drep_refs.slurm)." % (p_w, slug))
        cluster = {os.path.splitext(r["genome"])[0]: r["secondary_cluster"]
                   for r in leer_csv(p_c)}
        score = {os.path.splitext(r["genome"])[0]: a_float(r["score"], 0.0)
                 for r in leer_csv(p_s)} if os.path.exists(p_s) else {}
        # Control: si dRep no leyo extraW, los T0 no tendrian el bono de +1000.
        t0 = [g for g, r in clasif[slug].items() if r["nivel"] == "T0" and g in score]
        if t0 and max(score[g] for g in t0) < PESO_T0 * 0.5:
            die("en %s los genomas T0 no recibieron el peso extra: dRep no leyo "
                "extraW.tsv. Revisa el log de dRep y vuelve a correrlo." % slug)
        ganadores = []
        for r in leer_csv(p_w):
            gid = os.path.splitext(r["genome"])[0]
            info = dict(clasif[slug].get(gid, {"id": gid, "nivel": "T2",
                                                "habitat": "desconocido",
                                                "continente": "desconocido"}))
            info["cluster_drep"] = cluster.get(gid, r.get("cluster", ""))
            info["score_drep"] = score.get(gid, "")
            ganadores.append(info)
        v = viables_info.get(slug, {})
        bins = [b for b in (v.get("bins") or "").split(",") if b]
        ranking.append({
            "slug": slug, "linaje_id": e["linaje_id"], "especie_gtdb": e["especie_gtdb"],
            "ncbi_especie": v.get("ncbi_especie", ""), "ncbi_taxid": v.get("ncbi_taxid", ""),
            "n_postqc": v.get("n_postqc", ""), "n_descargados_ok": e["n_ok"],
            "n_no_redundantes": len(ganadores),
            "viable": "si" if len(ganadores) >= args.umbral else "no",
            "bins": ",".join(bins), "n_bins": len(bins),
            "completitud_bin_max": max([comp_bins.get(b, 0.0) for b in bins] or [0.0]),
            "_ganadores": ganadores,
        })

    ranking.sort(key=lambda f: (0 if f["viable"] == "si" else 1, -f["n_no_redundantes"],
                                -f["n_bins"], -f["completitud_bin_max"],
                                orden_linaje(f["linaje_id"])))
    n_sel = 0
    for i, f in enumerate(ranking, 1):
        f["rango"] = i if f["viable"] == "si" else "NA"
        if f["viable"] == "si" and n_sel < args.top:
            f["seleccionada"] = "si"
            n_sel += 1
        else:
            f["seleccionada"] = "no"
    if n_sel == 0:
        die("ninguna especie alcanza %d genomas no redundantes: detener y "
            "decidir con el asesor (D25)." % args.umbral)
    if n_sel < args.top:
        warn("solo %d especie(s) viable(s); se continua con las que hay (D25)." % n_sel)

    os.makedirs(args.outdir, exist_ok=True)
    escribir_tsv(os.path.join(args.outdir, "phase3_ranking_especies.tsv"), ranking,
                 ["rango", "linaje_id", "especie_gtdb", "slug", "ncbi_especie",
                  "ncbi_taxid", "n_postqc", "n_descargados_ok", "n_no_redundantes",
                  "viable", "seleccionada", "n_bins", "bins", "completitud_bin_max"])

    finales, composicion, manifiesto, especies_sel = [], [], [], []
    for f in [r for r in ranking if r["seleccionada"] == "si"]:
        slug = f["slug"]
        elegidos = seleccionar_tope(f["_ganadores"], args.tope)
        for orden, (r, motivo) in enumerate(elegidos, 1):
            fila = dict(r)
            fila.update({"slug": slug, "orden_seleccion": orden, "motivo": motivo})
            finales.append(fila)
            manifiesto.append({"id": r["id"], "fasta": os.path.abspath(
                os.path.join(args.refs_dir, slug, r["id"] + ".fna")),
                "tipo": "ref", "slug": slug, "genoma_original": r.get("accession", r["id"])})
        for dim in ("nivel", "habitat", "continente"):
            for valor, n in Counter(r[dim] for r, _ in elegidos).most_common():
                composicion.append({"slug": slug, "dimension": dim, "valor": valor, "n": n})
        bins_ids = []
        for b in f["bins"].split(","):
            fasta = os.path.abspath(os.path.join(args.bins_dir, b + ".fasta"))
            if not os.path.exists(fasta):
                die("no existe el bin %s (esperado en %s)" % (b, fasta))
            manifiesto.append({"id": id_bin(b), "fasta": fasta, "tipo": "bin",
                               "slug": slug, "genoma_original": b})
            bins_ids.append(id_bin(b))
        especies_sel.append({
            "rango": f["rango"], "slug": slug, "linaje_id": f["linaje_id"],
            "especie_gtdb": f["especie_gtdb"], "ncbi_especie": f["ncbi_especie"],
            "ncbi_taxid": f["ncbi_taxid"], "n_refs": len(elegidos),
            "n_no_redundantes": f["n_no_redundantes"], "bins": f["bins"],
            "bins_ids": ",".join(bins_ids),
        })
        print("[ok] %s. %-45s N_nr=%d  referencias=%d  bins=%s"
              % (f["rango"], slug, f["n_no_redundantes"], len(elegidos), f["bins"]))

    escribir_tsv(os.path.join(args.outdir, "phase3_referencias_finales.tsv"), finales,
                 ["slug", "orden_seleccion", "id", "accession", "nivel", "motivo",
                  "habitat", "patron_habitat", "isolation_source", "pais", "continente",
                  "latam", "assembly_level", "contigs", "n50", "completitud",
                  "contaminacion", "cluster_drep", "score_drep"])
    escribir_tsv(os.path.join(args.outdir, "phase3_composicion_referencias.tsv"),
                 composicion, ["slug", "dimension", "valor", "n"])
    escribir_tsv(os.path.join(args.outdir, "especies_seleccionadas.tsv"), especies_sel,
                 ["rango", "slug", "linaje_id", "especie_gtdb", "ncbi_especie",
                  "ncbi_taxid", "n_refs", "n_no_redundantes", "bins", "bins_ids"])
    campos_man = ["id", "fasta", "tipo", "slug", "genoma_original"]
    escribir_tsv(os.path.join(args.outdir, "bakta_manifest.tsv"), manifiesto, campos_man)

    # Piloto P0: bins + 3 referencias (en orden de seleccion) de la especie n.o 1.
    primera = especies_sel[0]["slug"]
    refs1 = [m for m in manifiesto if m["slug"] == primera and m["tipo"] == "ref"][:3]
    bins1 = [m for m in manifiesto if m["slug"] == primera and m["tipo"] == "bin"]
    escribir_tsv(os.path.join(args.outdir, "piloto_manifest.tsv"), refs1 + bins1, campos_man)

    with open(os.path.join(args.outdir, "phase3_seleccion_report.md"), "w",
              encoding="utf-8", newline="\n") as fh:
        w = fh.write
        w("# Fase 3 - Ranking definitivo y seleccion de referencias\n\n")
        w("Generado: %s\n\n" % datetime.now().strftime("%Y-%m-%d %H:%M"))
        w("Umbral: >= %d genomas no redundantes (dRep 99 %% ANI). Tope: %d por especie. "
          "Top: %d.\n\n" % (args.umbral, args.tope, args.top))
        w("| Rango | Linaje | Especie GTDB | Post-QC | Descargados | No redundantes | "
          "Viable | Seleccionada |\n|---|---|---|---|---|---|---|---|\n")
        for f in ranking:
            w("| %s | %s | *%s* | %s | %s | %s | %s | %s |\n" % (
                f["rango"], f["linaje_id"], f["especie_gtdb"], f["n_postqc"],
                f["n_descargados_ok"], f["n_no_redundantes"], f["viable"],
                f["seleccionada"]))
        w("\n## Composicion de las referencias seleccionadas\n\n")
        for e in especies_sel:
            w("### %s (%s)\n\n" % (e["especie_gtdb"], e["linaje_id"]))
            for dim in ("nivel", "habitat", "continente"):
                partes = ["%s=%d" % (c["valor"], c["n"]) for c in composicion
                          if c["slug"] == e["slug"] and c["dimension"] == dim]
                w("- %s: %s\n" % (dim, ", ".join(partes)))
            w("- bins propios: %s\n\n" % e["bins"])

    print("\n[ok] manifiesto de Bakta: %d genomas -> %s"
          % (len(manifiesto), os.path.join(args.outdir, "bakta_manifest.tsv")))
    print("[ok] manifiesto del piloto: %d genomas -> %s"
          % (len(refs1) + len(bins1), os.path.join(args.outdir, "piloto_manifest.tsv")))


# ------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description="Fase 3, etapa 3.2b: verificacion, clasificacion y seleccion de referencias.")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("listas", help="especies viables y accesiones a descargar")
    p.add_argument("--censo-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_listas)

    p = sub.add_parser("verificar", help="controla la descarga contra GTDB")
    p.add_argument("--censo-dir", required=True)
    p.add_argument("--descarga-dir", required=True)
    p.add_argument("--refs-dir", required=True)
    p.add_argument("--umbral", type=int, default=DEF_UMBRAL)
    p.add_argument("--tolerancia", type=float, default=DEF_TOLERANCIA)
    p.add_argument("--exclusiones", default="",
                   help="TSV (accession, linaje_id, motivo) de genomas excluidos a mano")
    p.set_defaults(func=cmd_verificar)

    p = sub.add_parser("clasificar", help="prioridad, habitat, continente + entrada de dRep")
    p.add_argument("--censo-dir", required=True)
    p.add_argument("--descarga-dir", required=True)
    p.add_argument("--refs-dir", required=True)
    p.add_argument("--keywords", required=True, help="metadata/habitat_keywords.tsv")
    p.add_argument("--paises", required=True, help="metadata/paises_continentes.tsv")
    p.add_argument("--outdir", required=True, help="carpeta de dRep (results/02_drep)")
    p.set_defaults(func=cmd_clasificar)

    p = sub.add_parser("seleccionar", help="ranking definitivo + seleccion <= tope")
    p.add_argument("--descarga-dir", required=True)
    p.add_argument("--drep-dir", required=True)
    p.add_argument("--refs-dir", required=True)
    p.add_argument("--bins-dir", required=True, help="FASTA de los bins propios")
    p.add_argument("--checkm2", default="",
                   help="quality_report.tsv de CheckM2 (Fase 1) para el desempate")
    p.add_argument("--outdir", required=True)
    p.add_argument("--umbral", type=int, default=DEF_UMBRAL)
    p.add_argument("--tope", type=int, default=DEF_TOPE)
    p.add_argument("--top", type=int, default=DEF_TOP)
    p.set_defaults(func=cmd_seleccionar)

    args = ap.parse_args()
    if not getattr(args, "func", None):
        ap.print_help()
        sys.exit(2)
    args.func(args)


if __name__ == "__main__":
    main()
