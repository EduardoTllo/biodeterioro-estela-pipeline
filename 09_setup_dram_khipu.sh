#!/bin/bash
# 09_setup_dram_khipu.sh
# Preparacion de DRAM para la FASE 4 (potencial metabolico y riesgo, OE4).
#
# Que hace:
#   1. Crea el entorno 'dram' con el environment.yaml oficial de DRAM 1.5.0
#      (python 3.10 y versiones fijadas por los autores; DRAM-bio 1.5.0 por pip).
#   2. Prepara SOLO las bases de datos que usa la Fase 4, con --select_db:
#        kofam_hmm + kofam_ko_list      -> identificadores KEGG (KO) sin licencia
#        dbcan + actividades + EC        -> CAZymes (polisacaridos, EPS)
#        peptidase (MEROPS)              -> peptidasas
#        formularios de DRAM             -> resumen ('distill') de modulos y funciones
#      Se omiten UniRef90 (requiere ~220-512 GB de RAM), Pfam, VOGDB y RefSeq
#      viral: no hacen falta para la tabla de marcadores y son las que mas
#      memoria piden al procesarse. KEGG Genes no se usa (requiere licencia);
#      los KO salen de KOfam.
#
# Donde se ejecuta: NODO DE LOGIN (necesita internet), dentro de tmux. Sin
# UniRef ni Pfam el procesamiento es liviano (descompresion y hmmpress).
#   Uso:   bash 09_setup_dram_khipu.sh
# Es idempotente: no recrea el entorno si ya existe con DRAM 1.5.0.

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

WORKDIR="${WORKDIR:-$HOME/estela/fase4}"
DBS_DIR="${DBS_DIR:-$HOME/dbs}"
DRAM_DB="$DBS_DIR/dram"
DRAM_VERSION="1.5.0"
ENV_YAML_URL="https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v${DRAM_VERSION}/environment.yaml"
THREADS=8

mkdir -p "$WORKDIR/logs" "$DBS_DIR"
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

echo "==> [1/3] Entorno 'dram' (DRAM $DRAM_VERSION)"
INSTALADA="$(conda run -n dram DRAM.py --help >/dev/null 2>&1 && conda list -n dram --export 2>/dev/null | grep -iE '^dram-bio=' | cut -d= -f2 || true)"
if [[ "$INSTALADA" == "$DRAM_VERSION" ]]; then
  echo "    Ya existe con DRAM $INSTALADA. Se omite."
else
  conda env remove -n dram -y >/dev/null 2>&1 || true
  curl -fsSL "$ENV_YAML_URL" -o "$WORKDIR/dram_environment.yaml"
  SOLVER="conda"; command -v mamba >/dev/null 2>&1 && SOLVER="mamba"
  $SOLVER env create -n dram -f "$WORKDIR/dram_environment.yaml"
fi
conda activate dram
DRAM.py --help >/dev/null && echo "    OK: $(conda list --export | grep -iE '^dram-bio=')"

echo "==> [2/3] Bases de datos de DRAM en $DRAM_DB (solo las de la Fase 4)"
if DRAM-setup.py print_config 2>/dev/null | grep -qE "KOfam db: .*kofam"; then
  echo "    Ya configuradas (ver print_config abajo). Se omite."
else
  DRAM-setup.py prepare_databases --output_dir "$DRAM_DB" --skip_uniref --threads "$THREADS" \
      --select_db kofam_hmm --select_db kofam_ko_list \
      --select_db dbcan --select_db dbcan_fam_activities --select_db dbcan_subfam_ec \
      --select_db peptidase \
      --select_db genome_summary_form --select_db module_step_form \
      --select_db function_heatmap_form --select_db amg_database --select_db etc_module_database \
      --clear_config --verbose
fi

echo "==> [3/3] Configuracion final"
DRAM-setup.py print_config | tee "$WORKDIR/dram_config.txt"
conda env export -n dram > "$WORKDIR/dram.lock.yml" 2>/dev/null || true
echo
echo "======================================================================"
echo " SETUP DE DRAM COMPLETO"
echo "   Bases:     $DRAM_DB"
echo "   Config:    $WORKDIR/dram_config.txt (deben figurar KOfam, dbCAN y MEROPS)"
echo "   Lock file: $WORKDIR/dram.lock.yml (copiarlo al repo: envs/)"
echo "======================================================================"
