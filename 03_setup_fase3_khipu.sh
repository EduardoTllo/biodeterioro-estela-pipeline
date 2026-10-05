#!/bin/bash
# 03_setup_fase3_khipu.sh
# Preparacion del entorno para la FASE 3 (analisis pangenomico comparativo, OE3)
# en el cluster Khipu (UTEC).
#
# Que hace:
#   1. Verifica que ~/.condarc use solo conda-forge + bioconda (el canal
#      'defaults' de repo.anaconda.com no resuelve en Khipu).
#   2. Crea los entornos con versiones FIJADAS (verificadas en bioconda el
#      2026-09-27) y python=3.11:
#        ncbi-datasets  -> descarga de genomas (CLI de NCBI Datasets)
#        bakta          -> anotacion (bakta 1.12.1, requiere BD esquema 6)
#        panaroo        -> pangenoma (panaroo 1.8.0; trae mafft, gffutils, biopython)
#        iqtree         -> arbol del core (iqtree 3.1.3 + snp-sites 2.5.1)
#        blast          -> verificacion de genes exclusivos (blast 2.17.0 + taxonkit 0.20.0)
#      El dRep es el MISMO entorno 'drep' de la Fase 2 (no se recrea).
#   3. Descarga la BD de Bakta v6.0 full (30 GB comprimida, 84 GB descomprimida).
#   4. Descarga el conversor Bakta -> Prokka de Panaroo desde el tag v1.8.0.
#      Panaroo no lee de forma completa los GFF3 de Bakta; el conversor si. No se
#      instala como comando en el paquete de bioconda, por eso se baja aparte.
#   5. Descarga el NCBI taxdump (para taxonkit) y el metadata de GTDB R220
#      (bac120_metadata_r220.tsv, que no viene en el tarball de GTDB-Tk).
#   6. Exporta lock files de TODOS los entornos del proyecto (tiara, checkm2,
#      gtdbtk, drep y los de esta fase) y un resumen de versiones.
#
# Donde se ejecuta: NODO DE LOGIN de Khipu (unico con internet), dentro de tmux
# (la BD de Bakta tarda en bajar y el SSH puede cortarse).
#   Uso:   bash 03_setup_fase3_khipu.sh
#
# Es idempotente: omite entornos que ya existen con la version correcta y
# descargas que ya estan completas.

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Rutas configurables
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SCRIPTS_DIR="$WORKDIR/scripts"
DBS_DIR="${DBS_DIR:-$HOME/dbs}"
BAKTA_DB_PARENT="$DBS_DIR/bakta"
TAXDUMP_DIR="$DBS_DIR/taxdump"
GTDB_RELEASE_DIR="${GTDB_RELEASE_DIR:-$HOME/gtdbtk_data/release220}"
GTDB_METADATA_URL="https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
PANAROO_CONVERTER_URL="https://raw.githubusercontent.com/gtonkinhill/panaroo/v1.8.0/scripts/convert_bakta_to_prokka_gff.py"
TAXDUMP_URL="https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Versiones fijadas (verificadas en bioconda/conda-forge el 2026-09-27)
PYTHON_VERSION="3.11"
BAKTA_VERSION="1.12.1"
BAKTA_DB_MAJOR="6"
PANAROO_VERSION="1.8.0"
IQTREE_VERSION="3.1.3"
SNPSITES_VERSION="2.5.1"
BLAST_VERSION="2.17.0"
TAXONKIT_VERSION="0.20.0"
# NCBI Datasets es solo el cliente de descarga (no afecta el contenido de los
# genomas, que se piden por accesion con version exacta); se instala la
# version vigente y queda registrada en el lock file.

MIN_FREE_GB_BAKTA=130

echo "==> [1/6] Cargando modulo miniconda"
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
mkdir -p "$WORKDIR" "$SCRIPTS_DIR" "$DBS_DIR"

CONDARC="$HOME/.condarc"
if [[ ! -f "$CONDARC" ]] || grep -qE 'defaults|repo\.anaconda\.com' "$CONDARC"; then
  [[ -f "$CONDARC" ]] && cp -f "$CONDARC" "${CONDARC}.bak.$(date +%Y%m%d%H%M%S)"
  cat > "$CONDARC" << 'EOF'
channels:
  - conda-forge
  - bioconda
channel_priority: strict
remote_connect_timeout_secs: 30
remote_read_timeout_secs: 120
remote_max_retries: 10
EOF
  echo "    ~/.condarc reescrito (solo conda-forge + bioconda)"
else
  echo "    ~/.condarc ya usa solo conda-forge + bioconda"
fi

SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; fi
echo "    Usando solver: $SOLVER"

# crear_env <nombre> <paquete_a_verificar> <version_esperada|""> <specs...>
crear_env() {
  local nombre="$1" pkg="$2" version="$3"; shift 3
  local instalada
  instalada="$(conda list -n "$nombre" --export 2>/dev/null | grep -E "^${pkg}=" | cut -d= -f2 || true)"
  if [[ -n "$instalada" && ( -z "$version" || "$instalada" == "$version" ) ]]; then
    echo "    '$nombre' ya existe con $pkg=$instalada. Se omite."
    return 0
  fi
  if conda env list | grep -qE "^\s*${nombre}\s"; then
    echo "    '$nombre' existe con $pkg='${instalada:-ninguna}' (se requiere '${version:-cualquiera}'). Se recrea."
    conda env remove -n "$nombre" -y >/dev/null 2>&1 || true
  fi
  echo "    Creando '$nombre': $*"
  $SOLVER create -y -n "$nombre" --override-channels -c conda-forge -c bioconda "$@"
}

echo "==> [2/6] Entornos conda"
crear_env ncbi-datasets ncbi-datasets-cli "" ncbi-datasets-cli "python=$PYTHON_VERSION"
crear_env bakta   bakta   "$BAKTA_VERSION"   "bakta=$BAKTA_VERSION" "python=$PYTHON_VERSION"
crear_env panaroo panaroo "$PANAROO_VERSION" "panaroo=$PANAROO_VERSION" "python=$PYTHON_VERSION"
crear_env iqtree  iqtree  "$IQTREE_VERSION"  "iqtree=$IQTREE_VERSION" "snp-sites=$SNPSITES_VERSION" "python=$PYTHON_VERSION"
crear_env blast   blast   "$BLAST_VERSION"   "blast=$BLAST_VERSION" "taxonkit=$TAXONKIT_VERSION" "python=$PYTHON_VERSION"
if ! conda run -n drep dRep -h >/dev/null 2>&1; then
  echo "    AVISO: no se encontro el entorno 'drep' de la Fase 2. Corre 02_setup_fase2_khipu.sh." >&2
fi

echo "==> [3/6] Base de datos de Bakta (v${BAKTA_DB_MAJOR}.x, tipo full)"
BAKTA_DB=""
VJSON="$(find "$BAKTA_DB_PARENT" -maxdepth 2 -name version.json 2>/dev/null | head -n1 || true)"
if [[ -n "$VJSON" ]] && grep -qE "\"major\": *${BAKTA_DB_MAJOR}\b" "$VJSON"; then
  BAKTA_DB="$(dirname "$VJSON")"
  echo "    Ya existe: $BAKTA_DB. Se omite la descarga."
else
  mkdir -p "$BAKTA_DB_PARENT"
  FREE_GB=$(df -BG --output=avail "$BAKTA_DB_PARENT" | tail -1 | tr -dc '0-9')
  echo "    Espacio libre: ${FREE_GB} GB (se requieren >= ${MIN_FREE_GB_BAKTA} GB)"
  if [[ "${FREE_GB:-0}" -lt "$MIN_FREE_GB_BAKTA" ]]; then
    echo "ERROR: espacio insuficiente para la BD de Bakta en $BAKTA_DB_PARENT" >&2; exit 1
  fi
  conda activate bakta
  # bakta_db verifica el MD5 del paquete antes de descomprimirlo.
  bakta_db download --output "$BAKTA_DB_PARENT" --type full
  conda deactivate
  VJSON="$(find "$BAKTA_DB_PARENT" -maxdepth 2 -name version.json | head -n1)"
  BAKTA_DB="$(dirname "$VJSON")"
fi
conda activate bakta
# Recomendado por la documentacion de Bakta: AMRFinderPlus prepara su propia BD.
amrfinder_update --force_update --database "$BAKTA_DB/amrfinderplus-db/" || \
  echo "    AVISO: amrfinder_update fallo; revisa antes de correr Bakta." >&2
conda deactivate
echo "$BAKTA_DB" > "$WORKDIR/bakta_db_path.txt"
echo "    BD de Bakta: $BAKTA_DB (ruta guardada en $WORKDIR/bakta_db_path.txt)"

echo "==> [4/6] Conversor Bakta -> Prokka de Panaroo (tag v${PANAROO_VERSION})"
CONV="$SCRIPTS_DIR/convert_bakta_to_prokka_gff.py"
if [[ ! -s "$CONV" ]]; then
  curl -fsSL "$PANAROO_CONVERTER_URL" -o "$CONV"
fi
conda activate panaroo
python "$CONV" -h >/dev/null && echo "    OK: $CONV"
conda deactivate

echo "==> [5/6] NCBI taxdump y metadata de GTDB R220"
mkdir -p "$TAXDUMP_DIR"
if [[ ! -s "$TAXDUMP_DIR/nodes.dmp" ]]; then
  curl -fsSL "$TAXDUMP_URL" -o "$TAXDUMP_DIR/taxdump.tar.gz"
  tar xzf "$TAXDUMP_DIR/taxdump.tar.gz" -C "$TAXDUMP_DIR" names.dmp nodes.dmp delnodes.dmp merged.dmp
  echo "$(date +%Y-%m-%d)" > "$TAXDUMP_DIR/fecha_descarga.txt"
fi
echo "    taxdump: $TAXDUMP_DIR (descargado el $(cat "$TAXDUMP_DIR/fecha_descarga.txt" 2>/dev/null || echo NA))"
mkdir -p "$GTDB_RELEASE_DIR"
if [[ ! -s "$GTDB_RELEASE_DIR/bac120_metadata_r220.tsv" ]]; then
  curl -fsSL "$GTDB_METADATA_URL" -o "$GTDB_RELEASE_DIR/bac120_metadata_r220.tsv.gz"
  gunzip -f "$GTDB_RELEASE_DIR/bac120_metadata_r220.tsv.gz"
fi
for col in ncbi_country gtdb_type_designation_ncbi_taxa checkm2_completeness; do
  if ! head -n1 "$GTDB_RELEASE_DIR/bac120_metadata_r220.tsv" | tr '\t' '\n' | grep -qx "$col"; then
    echo "ERROR: el metadata no tiene la columna $col" >&2; exit 1
  fi
done
echo "    metadata: $GTDB_RELEASE_DIR/bac120_metadata_r220.tsv"

echo "==> [6/6] Lock files y versiones de todos los entornos del proyecto"
ENVS_DIR="$SCRIPTS_DIR/envs"
mkdir -p "$ENVS_DIR"
VERS="$ENVS_DIR/versiones_herramientas.txt"
{
  echo "Versiones de herramientas - generado $(date)"
  echo "Fuente: conda list --export (paquete=version=build)"
  echo
} > "$VERS"
declare -A PKGS=( [tiara]="tiara" [checkm2]="checkm2" [gtdbtk]="gtdbtk" \
                  [drep]="drep|fastani|mash" [ncbi-datasets]="ncbi-datasets-cli" \
                  [bakta]="bakta|pyrodigal|diamond|amrfinder" [panaroo]="panaroo|mafft|cd-hit" \
                  [iqtree]="iqtree|snp-sites" [blast]="blast|taxonkit" )
for env in tiara checkm2 gtdbtk drep ncbi-datasets bakta panaroo iqtree blast; do
  if conda env list | grep -qE "^\s*${env}\s"; then
    conda env export -n "$env" > "$ENVS_DIR/${env}.lock.yml" 2>/dev/null || true
    {
      echo "[$env]"
      conda list -n "$env" --export 2>/dev/null | grep -E "^(${PKGS[$env]}|python)=" || true
      echo
    } >> "$VERS"
  else
    echo "    AVISO: entorno '$env' no existe; sin lock file." >&2
  fi
done
echo "Bakta DB: $(cat "$BAKTA_DB/version.json" | tr -d '\n ')" >> "$VERS"
cat "$VERS"

echo
echo "======================================================================"
echo " SETUP DE LA FASE 3 COMPLETO"
echo " Checkpoint C0:"
echo "   - Lock files en $ENVS_DIR  (copiarlos al repo: envs/)"
echo "   - BD de Bakta:  $BAKTA_DB"
echo "   - Conversor:    $CONV"
echo "   - taxdump:      $TAXDUMP_DIR"
echo "   - metadata:     $GTDB_RELEASE_DIR/bac120_metadata_r220.tsv"
echo "======================================================================"
