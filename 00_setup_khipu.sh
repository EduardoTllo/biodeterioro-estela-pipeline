#!/bin/bash
# 00_setup_khipu.sh
# Preparacion del entorno para la FASE 1 (Control de calidad y seleccion de
# genomas, OE1) de la tesis de la Estela de Raimondi, en el cluster Khipu (UTEC).
#
# Que hace:
#   - Crea dos entornos conda separados (para evitar conflictos de dependencias):
#       * tiara    -> filtrado por dominio (paso 1)
#       * checkm2  -> calidad genomica: completitud/contaminacion (paso 2)
#   - Descarga la base de datos de CheckM2 (~3 GB, DIAMOND).
#
# Donde se ejecuta:
#   SOLO en el NODO DE LOGIN de Khipu (es el unico con acceso a internet).
#   NO lo mandes con sbatch: correlo directamente en la terminal del login.
#
#   Uso:   bash 00_setup_khipu.sh
#
# Es idempotente: si un entorno ya existe, no lo vuelve a crear.

# Nota: NO usar 'set -u' (nounset): Lmod (module) referencia LD_PRELOAD sin
# definir y aborta la carga de modulos en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Rutas configurables
# Carpeta donde vivira la base de datos de CheckM2 (en tu HOME por defecto).
CHECKM2_DB_DIR="${CHECKM2_DB_DIR:-$HOME/dbs/checkm2}"

echo "==> [1/4] Cargando modulo miniconda"
module purge 2>/dev/null || true
module load miniconda/3.0
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: 'conda' no quedo en el PATH tras 'module load miniconda/3.0'." >&2
  echo "       Verifica el nombre exacto del modulo con:" >&2
  echo "         module avail 2>&1 | grep -i conda" >&2
  echo "       y ajusta la linea 'module load ...' en este script." >&2
  exit 1
fi
eval "$(conda shell.bash hook)"

# Usar mamba si esta disponible (mas rapido); si no, conda.
SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; fi
echo "    Usando solver: $SOLVER"

# Entorno 1: Tiara
echo "==> [2/4] Creando entorno 'tiara'"
if conda env list | grep -qE '^\s*tiara\s'; then
  echo "    El entorno 'tiara' ya existe. Se omite."
else
  $SOLVER create -y -n tiara -c conda-forge -c bioconda tiara
fi

# Entorno 2: CheckM2
echo "==> [3/4] Creando entorno 'checkm2'"
if conda env list | grep -qE '^\s*checkm2\s'; then
  echo "    El entorno 'checkm2' ya existe. Se omite."
else
  $SOLVER create -y -n checkm2 -c conda-forge -c bioconda checkm2
fi

# Base de datos de CheckM2
echo "==> [4/4] Descargando base de datos de CheckM2 en: $CHECKM2_DB_DIR"
mkdir -p "$CHECKM2_DB_DIR"
conda activate checkm2

DB_FILE="$CHECKM2_DB_DIR/CheckM2_database/uniref100.KO.1.dmnd"
if [[ -f "$DB_FILE" ]]; then
  echo "    La base de datos ya existe en $DB_FILE. Se omite la descarga."
else
  checkm2 database --download --path "$CHECKM2_DB_DIR"
fi
conda deactivate

# Exportar entornos exactos (reproducibilidad)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="$SCRIPT_DIR/envs"
mkdir -p "$ENVS_DIR"
echo "==> Exportando lock files de conda en: $ENVS_DIR"
conda env export -n tiara   > "$ENVS_DIR/tiara.lock.yml"   2>/dev/null || true
conda env export -n checkm2 > "$ENVS_DIR/checkm2.lock.yml" 2>/dev/null || true
echo "    (commitea envs/*.lock.yml al repo para fijar las versiones exactas)"

echo
echo "============================================================"
echo " SETUP COMPLETO"
echo "------------------------------------------------------------"
echo " Entornos conda creados:  tiara, checkm2"
echo " Base de datos CheckM2:   $DB_FILE"
echo
echo " Verificacion de versiones:"
conda activate tiara   && echo -n "   Tiara:   " && tiara   --version 2>&1 | head -n1 ; conda deactivate
conda activate checkm2 && echo -n "   CheckM2: " && checkm2 --version 2>&1 | head -n1 ; conda deactivate
echo
echo " SIGUIENTE PASO:"
echo "   1) Sube tus 97 bins (FASTA) a la carpeta de trabajo en Khipu."
echo "   2) Edita las rutas en run_fase1.slurm si hace falta."
echo "   3) Envia el job:   sbatch run_fase1.slurm"
echo "============================================================"
