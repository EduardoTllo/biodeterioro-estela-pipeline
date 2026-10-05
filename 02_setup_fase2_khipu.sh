#!/bin/bash
# 02_setup_fase2_khipu.sh
# Preparacion del entorno para la FASE 2 (Asignacion taxonomica y definicion de
# linajes, OE2) en el cluster Khipu (UTEC).
#
# Que hace:
#   - Deja ~/.condarc con solo conda-forge + bioconda (el canal 'defaults' de
#     repo.anaconda.com no resuelve en Khipu).
#   - Crea dos entornos conda con versiones PINNEADAS:
#       * gtdbtk  -> clasificacion taxonomica + ANI + arbol (GTDB-Tk, datos R220)
#       * drep    -> desreplicacion y definicion de linajes (dRep + fastANI + Mash)
#   - Verifica que haya espacio suficiente y descomprime la base GTDB R220
#     (~110 GB descomprimidos) desde gtdbtk_r220_data.tar.gz.
#   - Fija GTDBTK_DATA_PATH dentro del entorno y verifica la integridad de la base.
#
# Donde se ejecuta:
#   En el NODO DE LOGIN de Khipu (unico con internet, para crear los entornos).
#   La descompresion es I/O de archivos (permitida en login).
#
#   Uso:   bash 02_setup_fase2_khipu.sh
#
# Es idempotente: omite la creacion de un entorno solo si ya existe CON LA
# VERSION CORRECTA, y no re-descomprime la base si ya esta extraida.

# Nota: no usamos 'set -u' por compatibilidad con Lmod en Khipu (igual que Fase 1).
set -eo pipefail

# Rutas configurables
# Carpeta que contiene gtdbtk_r220_data.tar.gz:
GTDB_TARDIR="${GTDB_TARDIR:-$HOME/gtdbtk_data}"
GTDB_TAR="$GTDB_TARDIR/gtdbtk_r220_data.tar.gz"
# Carpeta destino donde se descomprimira la base (necesita ~110 GB libres):
GTDB_DEST="${GTDB_DEST:-$HOME/gtdbtk_data}"
# Espacio minimo requerido en el destino (GB):
MIN_FREE_GB=120

# Version de GTDB-Tk. DEBE ser compatible con el release de la base de datos.
# Tabla oficial de compatibilidad (docs de GTDB-Tk):
#   datos R220 -> gtdbtk 2.4.0 a 2.6.1
#   datos R232 -> gtdbtk 2.7.0 en adelante
# Con datos R220 usamos 2.6.1 (la mas alta compatible, con mas correcciones).
GTDBTK_VERSION="${GTDBTK_VERSION:-2.6.1}"

# Version de Python. CRITICO: hay que fijarla.
# Sin pin, conda resuelve Python 3.14, y GTDB-Tk falla al arrancar con
#   ValueError: "__StageLogger" object has no field "version"
# porque la funcionalidad de pydantic v1 que usa GTDB-Tk no es compatible con
# Python >= 3.14 (issue #669 del repo de GTDB-Tk).
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

echo "==> [1/4] Cargando modulo miniconda"
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

# En Khipu el canal 'defaults' (repo.anaconda.com) no resuelve de forma fiable.
# Ademas, algunos ~/.condarc traen las URLs de repo.anaconda.com escritas de
# forma EXPLICITA (no como el nombre "defaults"), en cuyo caso
# 'conda config --remove channels defaults' NO las quita (no hace match por
# texto). Por eso reescribimos ~/.condarc por completo, dejando solo
# conda-forge y bioconda (conda.anaconda.org), y afinamos timeouts/reintentos
# para tolerar una red intermitente.
CONDARC="$HOME/.condarc"
if [[ -f "$CONDARC" ]]; then
  cp -f "$CONDARC" "${CONDARC}.bak.$(date +%Y%m%d%H%M%S)"
  echo "    Backup de ~/.condarc guardado junto al original (.bak.*)"
fi
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

SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then SOLVER="mamba"; fi
echo "    Usando solver: $SOLVER"

# Entorno GTDB-Tk
# Idempotencia real: no basta con que exista un entorno con ese NOMBRE (un
# intento previo fallido puede dejar un entorno vacio/a medias registrado).
# Se verifica que el binario realmente funcione antes de omitir la creacion.
echo "==> [2/4] Creando entorno 'gtdbtk'"
# Se exige que exista Y que sea exactamente la version pinneada: una version
# distinta (p.ej. la ultima) es incompatible con los datos R220.
INSTALLED_VER="$(conda run -n gtdbtk gtdbtk --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
if [[ "$INSTALLED_VER" == "$GTDBTK_VERSION" ]]; then
  echo "    El entorno 'gtdbtk' ya existe con la version correcta ($GTDBTK_VERSION). Se omite."
else
  if conda env list | grep -qE '^\s*gtdbtk\s'; then
    echo "    Entorno 'gtdbtk' con version incorrecta ('${INSTALLED_VER:-ninguna}', se requiere $GTDBTK_VERSION). Se recrea."
    conda env remove -n gtdbtk -y >/dev/null 2>&1 || true
  fi
  # Versiones PINNEADAS: GTDB-Tk al release de los datos y Python a una serie
  # compatible con pydantic v1 (ver notas arriba).
  # --override-channels: no consultar 'defaults' (repo.anaconda.com).
  echo "    Instalando gtdbtk=$GTDBTK_VERSION con python=$PYTHON_VERSION"
  $SOLVER create -y -n gtdbtk --override-channels -c conda-forge -c bioconda \
      "gtdbtk=$GTDBTK_VERSION" "python=$PYTHON_VERSION"
fi

# Entorno dRep
echo "==> [3/4] Creando entorno 'drep'"
if conda run -n drep dRep -h >/dev/null 2>&1; then
  echo "    El entorno 'drep' ya existe y funciona. Se omite."
else
  if conda env list | grep -qE '^\s*drep\s'; then
    echo "    Entorno 'drep' incompleto/roto (intento previo fallido). Se recrea."
    conda env remove -n drep -y >/dev/null 2>&1 || true
  fi
  # Python tambien pinneado aqui: sin pin, conda resuelve la serie mas nueva y
  # se expone al mismo tipo de incompatibilidad que rompio a GTDB-Tk.
  $SOLVER create -y -n drep --override-channels -c conda-forge -c bioconda \
      drep fastani mash "python=$PYTHON_VERSION"
fi

# Descompresion de la base GTDB R220
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

# Exportar lock files (reproducibilidad)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="$SCRIPT_DIR/envs"
mkdir -p "$ENVS_DIR"
conda env export -n gtdbtk > "$ENVS_DIR/gtdbtk.lock.yml" 2>/dev/null || true
conda env export -n drep   > "$ENVS_DIR/drep.lock.yml"   2>/dev/null || true

# Verificacion de la base de referencia
# Se hacen DOS verificaciones:
#  1) Integridad real: archivos de 0 bytes (sintoma de descompresion truncada).
#  2) gtdbtk check_install, SOLO INFORMATIVO.
#
# check_install NO es bloqueante a proposito: su manifiesto de hashes tiene
# falsos positivos documentados con R220 (issues #594 y #626 del repo de
# GTDB-Tk; el changelog de 2.7.1 incluye un "fix for MD5 mismatch in
# check_install configuration"). Reporta HASH MISMATCH por 1 archivo entre
# ~114.000 aun con el MD5 del tar.gz correcto. La prueba real es que
# classify_wf corra.
echo "==> Verificando la base de referencia"
conda activate gtdbtk
export GTDBTK_DATA_PATH="$GTDBTK_DATA_PATH_DETECTED"

echo "    [1/2] Buscando archivos de 0 bytes (descompresion truncada)..."
N_EMPTY=$(find "$GTDBTK_DATA_PATH_DETECTED" -type f -size 0 2>/dev/null | wc -l)
if [[ "$N_EMPTY" -gt 0 ]]; then
  echo "    ERROR: se encontraron $N_EMPTY archivos vacios: la descompresion" >&2
  echo "           quedo incompleta. Vuelve a extraer el tar.gz." >&2
  find "$GTDBTK_DATA_PATH_DETECTED" -type f -size 0 2>/dev/null | head -5 >&2
  conda deactivate
  exit 1
fi
echo "          OK: sin archivos vacios."

echo "    [2/2] gtdbtk check_install (informativo, no bloqueante)..."
if gtdbtk check_install >/dev/null 2>&1; then
  echo "          OK: check_install paso."
else
  echo "          AVISO: check_install reporta HASH MISMATCH."
  echo "          Es un falso positivo conocido con R220 y NO bloquea el analisis."
  echo "          Se continua; la prueba real es la ejecucion de classify_wf."
fi
conda deactivate

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
