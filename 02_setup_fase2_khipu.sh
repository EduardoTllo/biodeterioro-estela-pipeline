#!/bin/bash
# =============================================================================
# 02_setup_fase2_khipu.sh
# -----------------------------------------------------------------------------
# Preparacion del entorno para la FASE 2 (Asignacion taxonomica y definicion de
# linajes, OE2) en el cluster Khipu (UTEC).
#
# QUE HACE:
#   - Crea dos entornos conda:
#       * gtdbtk  -> clasificacion taxonomica + ANI + arbol (GTDB-Tk, datos R220)
#       * drep    -> desreplicacion y definicion de linajes (dRep + fastANI + Mash)
#   - Verifica que haya espacio suficiente y descomprime la base GTDB R220
#     (~110 GB descomprimidos) desde gtdbtk_r220_data.tar.gz.
#
# DONDE SE EJECUTA:
#   En el NODO DE LOGIN de Khipu (unico con internet, para crear los entornos).
#   La descompresion es I/O de archivos (permitida en login).
#
#   Uso:   bash 02_setup_fase2_khipu.sh
#
# Es idempotente: no recrea entornos ni re-descomprime si ya existen.
# =============================================================================
# NOTA: no usamos 'set -u' por compatibilidad con Lmod en Khipu (igual que Fase 1).
set -eo pipefail

# --- Rutas configurables -----------------------------------------------------
# Carpeta que contiene gtdbtk_r220_data.tar.gz:
GTDB_TARDIR="${GTDB_TARDIR:-$HOME/gtdbtk_data}"
GTDB_TAR="$GTDB_TARDIR/gtdbtk_r220_data.tar.gz"
# Carpeta destino donde se descomprimira la base (necesita ~110 GB libres):
GTDB_DEST="${GTDB_DEST:-$HOME/gtdbtk_data}"
# Espacio minimo requerido en el destino (GB):
MIN_FREE_GB=120

echo "==> [1/4] Cargando modulo miniconda"
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; fi
echo "    Usando solver: $SOLVER"

# --- Entorno GTDB-Tk ---------------------------------------------------------
echo "==> [2/4] Creando entorno 'gtdbtk'"
if conda env list | grep -qE '^\s*gtdbtk\s'; then
  echo "    El entorno 'gtdbtk' ya existe. Se omite."
else
  # gtdbtk >=2.4 es necesario para los datos R220.
  $SOLVER create -y -n gtdbtk -c conda-forge -c bioconda 'gtdbtk>=2.4.0'
fi

# --- Entorno dRep ------------------------------------------------------------
echo "==> [3/4] Creando entorno 'drep'"
if conda env list | grep -qE '^\s*drep\s'; then
  echo "    El entorno 'drep' ya existe. Se omite."
else
  $SOLVER create -y -n drep -c conda-forge -c bioconda drep fastani mash
fi

# --- Descompresion de la base GTDB R220 --------------------------------------
echo "==> [4/4] Preparando base GTDB R220"
RELEASE_DIR="$GTDB_DEST/release220"
if [[ -d "$RELEASE_DIR" && -n "$(ls -A "$RELEASE_DIR" 2>/dev/null)" ]]; then
  echo "    Ya existe $RELEASE_DIR (base descomprimida). Se omite."
else
  if [[ ! -f "$GTDB_TAR" ]]; then
    echo "ERROR: no se encuentra $GTDB_TAR" >&2
    echo "       Ajusta GTDB_TARDIR o coloca ahi gtdbtk_r220_data.tar.gz." >&2
    exit 1
  fi
  mkdir -p "$GTDB_DEST"
  # Verificar espacio libre en el destino.
  FREE_GB=$(df -BG --output=avail "$GTDB_DEST" | tail -1 | tr -dc '0-9')
  echo "    Espacio libre en $GTDB_DEST: ${FREE_GB} GB (se requieren >= ${MIN_FREE_GB} GB)"
  if [[ "${FREE_GB:-0}" -lt "$MIN_FREE_GB" ]]; then
    echo "ERROR: espacio insuficiente. Cambia GTDB_DEST a un area con mas espacio" >&2
    echo "       (p.ej. scratch o proyecto) exportando GTDB_DEST=/ruta antes de correr." >&2
    exit 1
  fi
  echo "    Descomprimiendo $GTDB_TAR (puede tardar)..."
  tar xzf "$GTDB_TAR" -C "$GTDB_DEST"
  echo "    Descompresion lista."
fi

# Detectar la carpeta real de la base (algunos paquetes extraen a release220/).
if [[ -d "$RELEASE_DIR" ]]; then
  GTDBTK_DATA_PATH_DETECTED="$RELEASE_DIR"
else
  # Buscar una carpeta que contenga los marcadores esperados.
  GTDBTK_DATA_PATH_DETECTED="$(find "$GTDB_DEST" -maxdepth 2 -type d -name 'markers' -printf '%h\n' 2>/dev/null | head -1)"
  [[ -z "$GTDBTK_DATA_PATH_DETECTED" ]] && GTDBTK_DATA_PATH_DETECTED="$GTDB_DEST"
fi

# Fijar la variable en el entorno gtdbtk para futuras corridas.
conda activate gtdbtk
conda env config vars set GTDBTK_DATA_PATH="$GTDBTK_DATA_PATH_DETECTED" >/dev/null 2>&1 || true
conda deactivate

# --- Exportar lock files (reproducibilidad) ----------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="$SCRIPT_DIR/envs"
mkdir -p "$ENVS_DIR"
conda env export -n gtdbtk > "$ENVS_DIR/gtdbtk.lock.yml" 2>/dev/null || true
conda env export -n drep   > "$ENVS_DIR/drep.lock.yml"   2>/dev/null || true

echo
echo "============================================================"
echo " SETUP FASE 2 COMPLETO"
echo "------------------------------------------------------------"
echo " Entornos conda:        gtdbtk, drep"
echo " GTDBTK_DATA_PATH:      $GTDBTK_DATA_PATH_DETECTED"
echo " Verificacion rapida:"
conda activate gtdbtk && echo -n "   GTDB-Tk: " && gtdbtk --version 2>&1 | head -n1 ; conda deactivate || true
conda activate drep   && echo -n "   dRep:    " && dRep --version 2>&1 | head -n1 ; conda deactivate || true
echo
echo " Si GTDBTK_DATA_PATH quedo mal, verifica la ruta que contiene"
echo " las carpetas 'markers/', 'msa/', 'taxonomy/', etc. y ajustala"
echo " en run_fase2.slurm."
echo "============================================================"
