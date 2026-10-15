#!/bin/bash
# eucariotas/00_setup_euk_khipu.sh
# Preparacion de las Fases 1 y 2 EUCARIOTAS (anexo) en el cluster Khipu (UTEC).
#
# Dos etapas, para poder lanzar la Fase 1 mientras bajan las referencias:
#   bash 00_setup_euk_khipu.sh base   -> entornos + BD de EukCC + BUSCO + ITS
#                                        (lo necesario para run_fase1_euk.slurm)
#   bash 00_setup_euk_khipu.sh refs   -> genomas de referencia de hongos de NCBI
#                                        (un genoma de referencia por especie,
#                                        ~5 600 genomas; lo usa la Fase 2)
#
# Entornos nuevos (versiones fijadas, verificadas en bioconda el 2026-10-14):
#   eukcc   -> EukCC 2.2.0 (completitud/contaminacion de eucariotas)
#   busco   -> BUSCO 6.1.0 (datasets odb12.2)
#   euktax  -> skani 0.3.2, barrnap 0.9, ITSx 1.1.3, MAFFT 7.525, trimAl 1.5.1
# Se reutilizan los entornos ya existentes: drep (Fase 2), blast y
# ncbi-datasets (Fase 3) y el taxdump de ~/dbs/taxdump.
#
# Donde se ejecuta: NODO DE LOGIN (unico con internet), dentro de tmux.
# Es idempotente: omite entornos y descargas que ya estan completos.

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

ETAPA="${1:-}"
if [[ "$ETAPA" != "base" && "$ETAPA" != "refs" ]]; then
  echo "Uso: bash 00_setup_euk_khipu.sh base|refs" >&2; exit 1
fi

# Rutas configurables
SCRIPTS_DIR="${SCRIPTS_DIR:-$HOME/estela/fase3/scripts}"   # clon del repo
EUK_DIR="$SCRIPTS_DIR/eucariotas"
DBS_DIR="${DBS_DIR:-$HOME/dbs}"
EUKCC_DIR="$DBS_DIR/eukcc"
BUSCO_DB="$DBS_DIR/busco"
ITS_DIR="$DBS_DIR/its"
REFS_DIR="$DBS_DIR/euk_refs"
TAXDUMP_DIR="$DBS_DIR/taxdump"
# UNITE no tiene enlace de descarga estable: se baja a mano (ver eucariotas/README.md)
# y se deja aqui. Si no existe, la Fase 2 usa solo RefSeq ITS.
UNITE_FASTA="${UNITE_FASTA:-$ITS_DIR/unite_fungi.fasta}"

EUKCC_DB_URL="http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz"
BUSCO_LINAJES="fungi_odb12.2 eukaryota_odb12.2"
REFS_TAXON="${REFS_TAXON:-4751}"     # 4751 = Fungi

PYTHON_VERSION="3.11"
EUKCC_VERSION="2.2.0"
BUSCO_VERSION="6.1.0"
SKANI_VERSION="0.3.2"
BARRNAP_VERSION="0.9"
ITSX_VERSION="1.1.3"
MAFFT_VERSION="7.525"
TRIMAL_VERSION="1.5.1"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
mkdir -p "$DBS_DIR"

SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; fi

# crear_env <nombre> <paquete_a_verificar> <version_esperada> <specs...>
crear_env() {
  local nombre="$1" pkg="$2" version="$3"; shift 3
  local instalada
  instalada="$(conda list -n "$nombre" --export 2>/dev/null | grep -E "^${pkg}=" | cut -d= -f2 || true)"
  if [[ "$instalada" == "$version" ]]; then
    echo "    '$nombre' ya existe con $pkg=$instalada. Se omite."
    return 0
  fi
  if conda env list | grep -qE "^\s*${nombre}\s"; then
    echo "    '$nombre' tiene $pkg='${instalada:-ninguna}' (se requiere $version). Se recrea."
    conda env remove -n "$nombre" -y >/dev/null 2>&1 || true
  fi
  echo "    Creando '$nombre': $*"
  $SOLVER create -y -n "$nombre" --override-channels -c conda-forge -c bioconda "$@"
}

if [[ "$ETAPA" == "base" ]]; then
  echo "==> [1/5] Entornos conda"
  crear_env eukcc  eukcc "$EUKCC_VERSION" "eukcc=$EUKCC_VERSION" "python=$PYTHON_VERSION"
  crear_env busco  busco "$BUSCO_VERSION" "busco=$BUSCO_VERSION" "python=$PYTHON_VERSION"
  crear_env euktax skani "$SKANI_VERSION" "skani=$SKANI_VERSION" "barrnap=$BARRNAP_VERSION" \
      "itsx=$ITSX_VERSION" "mafft=$MAFFT_VERSION" "trimal=$TRIMAL_VERSION" "python=$PYTHON_VERSION"
  for env in drep blast ncbi-datasets; do
    conda env list | grep -qE "^\s*${env}\s" || {
      echo "ERROR: falta el entorno '$env' (se crea en 02_setup_fase2 / 03_setup_fase3)." >&2; exit 1; }
  done
  [[ -s "$TAXDUMP_DIR/nodes.dmp" ]] || {
    echo "ERROR: falta $TAXDUMP_DIR/nodes.dmp (se baja en 03_setup_fase3_khipu.sh)." >&2; exit 1; }

  echo "==> [2/5] Base de datos de EukCC 2 (v1.2)"
  if [[ -d "$EUKCC_DIR/eukcc2_db_ver_1.2" ]]; then
    echo "    Ya existe. Se omite."
  else
    mkdir -p "$EUKCC_DIR"
    curl -fL --retry 5 -C - -o "$EUKCC_DIR/eukcc2_db_ver_1.2.tar.gz" "$EUKCC_DB_URL"
    tar xzf "$EUKCC_DIR/eukcc2_db_ver_1.2.tar.gz" -C "$EUKCC_DIR"
    rm -f "$EUKCC_DIR/eukcc2_db_ver_1.2.tar.gz"
  fi

  echo "==> [3/5] Datasets de BUSCO: $BUSCO_LINAJES"
  conda activate busco
  for lin in $BUSCO_LINAJES; do
    if [[ -d "$BUSCO_DB/lineages/$lin" ]]; then
      echo "    $lin ya existe. Se omite."
    else
      busco --download_path "$BUSCO_DB" --download "$lin"
    fi
  done
  conda deactivate

  echo "==> [4/5] Bases de ITS"
  mkdir -p "$ITS_DIR"
  conda activate blast
  if [[ -f "$ITS_DIR/ITS_RefSeq_Fungi.nsq" || -f "$ITS_DIR/ITS_RefSeq_Fungi.00.nsq" ]]; then
    echo "    RefSeq ITS ya existe. Se omite."
  else
    (cd "$ITS_DIR" && update_blastdb.pl --decompress ITS_RefSeq_Fungi)
    date +%F > "$ITS_DIR/ITS_RefSeq_Fungi.fecha_descarga.txt"
  fi
  if [[ -f "$UNITE_FASTA" ]]; then
    if [[ -f "$ITS_DIR/unite.nsq" || -f "$ITS_DIR/unite.00.nsq" ]]; then
      echo "    BD BLAST de UNITE ya existe. Se omite."
    else
      makeblastdb -in "$UNITE_FASTA" -dbtype nucl -out "$ITS_DIR/unite" -title "UNITE $(basename "$UNITE_FASTA")"
    fi
  else
    echo "    AVISO: no esta $UNITE_FASTA. Bajalo a mano (eucariotas/README.md, paso 2) y vuelve a correr 'base'."
  fi
  conda deactivate

  echo "==> [5/5] Lock files y versiones"
  ENVS_DIR="$EUK_DIR/envs"
  mkdir -p "$ENVS_DIR"
  VERS="$ENVS_DIR/versiones_herramientas_euk.txt"
  { echo "Versiones (analisis eucariota) - $(date)"; echo; } > "$VERS"
  declare -A PKGS=( [eukcc]="eukcc|metaeuk|pplacer|epa-ng|hmmer" [busco]="busco|miniprot|hmmer" \
                    [euktax]="skani|barrnap|itsx|mafft|trimal" [blast]="blast" [drep]="drep|fastani|mash" )
  for env in eukcc busco euktax; do
    conda env export -n "$env" > "$ENVS_DIR/${env}.lock.yml" 2>/dev/null || true
  done
  for env in eukcc busco euktax blast drep; do
    { echo "[$env]"; conda list -n "$env" --export 2>/dev/null | grep -E "^(${PKGS[$env]}|python)=" || true; echo; } >> "$VERS"
  done
  { echo "EukCC DB: eukcc2_db_ver_1.2"; echo "BUSCO: $BUSCO_LINAJES";
    echo "RefSeq ITS: $(cat "$ITS_DIR/ITS_RefSeq_Fungi.fecha_descarga.txt" 2>/dev/null || echo NA)";
    echo "UNITE: $(basename "$UNITE_FASTA") $( [[ -f "$UNITE_FASTA" ]] && echo OK || echo FALTA)"; } >> "$VERS"
  cat "$VERS"
  echo
  echo "Etapa 'base' lista. Ya puedes enviar run_fase1_euk.slurm."
  echo "Despues corre (en tmux):  bash $EUK_DIR/00_setup_euk_khipu.sh refs"
fi

if [[ "$ETAPA" == "refs" ]]; then
  echo "==> Genomas de referencia (taxon $REFS_TAXON, solo genomas de referencia de NCBI)"
  mkdir -p "$REFS_DIR"
  conda activate ncbi-datasets
  cd "$REFS_DIR"
  if [[ ! -s refs_datasets.tsv ]]; then
    datasets summary genome taxon "$REFS_TAXON" --reference --as-json-lines \
      | dataformat tsv genome --fields accession,organism-name,organism-tax-id > refs_datasets.tsv
    date +%F > fecha_descarga.txt
  fi
  echo "    Genomas en el resumen: $(($(wc -l < refs_datasets.tsv) - 1))"
  if [[ ! -d genomas/ncbi_dataset/data ]]; then
    datasets download genome taxon "$REFS_TAXON" --reference --include genome \
      --dehydrated --filename refs_dehydrated.zip
    unzip -q -o refs_dehydrated.zip -d genomas
  fi
  # rehydrate es reanudable: si se corta, se vuelve a correr y sigue.
  datasets rehydrate --directory genomas --gzip --max-workers 20
  find "$REFS_DIR/genomas/ncbi_dataset/data" -name '*_genomic.fna.gz' | sort > refs_fasta.txt
  echo "    FASTA descargados: $(wc -l < refs_fasta.txt)"
  conda deactivate

  conda activate euktax
  python "$EUK_DIR/bin/summarize_euk.py" refs \
      --datasets-tsv refs_datasets.tsv --fasta-list refs_fasta.txt \
      --taxdump "$TAXDUMP_DIR" --out refs_fungi.tsv
  conda deactivate
  echo
  echo "Etapa 'refs' lista: $REFS_DIR/refs_fungi.tsv"
  echo "El sketch de skani se arma dentro de run_fase2_euk.slurm (usa 32 CPU)."
fi
