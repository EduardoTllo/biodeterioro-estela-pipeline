#!/bin/bash
# run_fase1.slurm
# FASE 1 - Control de calidad y seleccion de genomas (OE1)
# Tesis: Estela de Raimondi | Cluster: Khipu (UTEC) | SLURM 25.11
#
# Flujo:
#   Paso 1 (Tiara)   -> filtrado por dominio; descarta bins cuya fraccion
#                       mayoritaria de longitud sea eucariota u organelar.
#                       Longitud minima de contig = 3000 bp.
#   Paso 2 (CheckM2) -> completitud y contaminacion de los bins procariotas.
#                       Se conservan: completitud > 70 % y contaminacion < 5 %.
#   Paso 3 (resumen) -> tablas y reporte legible del proceso.
#
# Envio:   sbatch run_fase1.slurm

# Recursos dentro de los limites de la cuenta (max: 32 cores, 98 GB RAM, 24 h):
#SBATCH --job-name=fase1_estela
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=fase1_%j.out
#SBATCH --error=fase1_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=diegotello2000@hotmail.com

# Nota: NO usar 'set -u' (nounset): Lmod (module) referencia LD_PRELOAD sin
# definir y aborta la carga de modulos en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion: edita estas rutas segun tu cuenta en Khipu
# Carpeta raiz del proyecto de la Fase 1 en Khipu:
WORKDIR="$HOME/estela/fase1"

# Carpeta con los 97 bins de entrada (un archivo FASTA por bin):
BINS_DIR="$WORKDIR/bins"

# Extension de los archivos de bins (fa | fna | fasta):
# Tras 01_prep_bins.sh los bins quedan como bin_0001.fasta ...
BIN_EXT="fasta"

# Base de datos de CheckM2 descargada en el setup:
CHECKM2_DB="$HOME/dbs/checkm2/CheckM2_database/uniref100.KO.1.dmnd"

# Umbrales de la Fase 1 (segun tesis, seccion 4.1.3):
MIN_CONTIG_LEN=3000     # longitud minima de contig para Tiara
MIN_COMPLETENESS=70     # completitud > 70 %
MAX_CONTAMINATION=5     # contaminacion < 5 %

# Directorio donde estan los scripts (esta carpeta):
SCRIPTS_DIR="$WORKDIR/scripts"

# No edites debajo de esta linea salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-16}"
RESULTS="$WORKDIR/results"
TIARA_DIR="$RESULTS/01_tiara"
CHECKM2_INPUT="$RESULTS/checkm2_input"
CHECKM2_OUT="$RESULTS/02_checkm2"
SELECTED_DIR="$RESULTS/phase1_selected_genomes"

echo "======================================================================"
echo " FASE 1 - Estela de Raimondi"
echo " Job SLURM: ${SLURM_JOB_ID:-NA}  |  Nodo: $(hostname)  |  Hilos: $THREADS"
echo " Inicio: $(date)"
echo "======================================================================"

# Comprobaciones previas
if [[ ! -d "$BINS_DIR" ]]; then
  echo "ERROR: no existe la carpeta de bins: $BINS_DIR" >&2; exit 1
fi
shopt -s nullglob
BIN_FILES=("$BINS_DIR"/*."$BIN_EXT")
N_BINS=${#BIN_FILES[@]}
if [[ "$N_BINS" -eq 0 ]]; then
  echo "ERROR: no se encontraron archivos *.$BIN_EXT en $BINS_DIR" >&2; exit 1
fi
if [[ ! -f "$CHECKM2_DB" ]]; then
  echo "ERROR: no existe la base de datos de CheckM2: $CHECKM2_DB" >&2
  echo "       Corre primero 00_setup_khipu.sh en el nodo de login." >&2; exit 1
fi
echo "Bins de entrada encontrados: $N_BINS"

mkdir -p "$TIARA_DIR" "$CHECKM2_INPUT" "$CHECKM2_OUT" "$SELECTED_DIR"

# Entorno conda
module purge 2>/dev/null || true
module load miniconda/3.0
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: 'conda' no quedo en el PATH tras 'module load miniconda/3.0'." >&2
  echo "       Verifica con: module avail 2>&1 | grep -i conda" >&2
  exit 1
fi
eval "$(conda shell.bash hook)"

# Registro de versiones (reproducibilidad)
VERS="$RESULTS/versions.txt"
{
  echo "Fase 1 - registro de versiones"
  echo "Fecha:   $(date)"
  echo "Job:     ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  echo "Umbrales: min_contig=$MIN_CONTIG_LEN, completitud>$MIN_COMPLETENESS, contaminacion<$MAX_CONTAMINATION"
  echo "N bins entrada: $N_BINS"
  echo "----"
  conda activate tiara   2>/dev/null && echo "Tiara:   $(tiara --version 2>&1 | head -n1)"   && conda deactivate || true
  conda activate checkm2 2>/dev/null && echo "CheckM2: $(checkm2 --version 2>&1 | head -n1)" && conda deactivate || true
  echo "CheckM2 DB: $CHECKM2_DB"
} > "$VERS"
echo "Versiones registradas en: $VERS"

# Paso 1: filtrado por dominio (Tiara)
echo
echo "---> PASO 1: Filtrado por dominio con Tiara (min_len=$MIN_CONTIG_LEN)"
conda activate tiara
echo "     Tiara version: $(tiara --version 2>&1 | head -n1)"

cd "$TIARA_DIR"
i=0
for f in "${BIN_FILES[@]}"; do
  i=$((i+1))
  base="$(basename "$f")"
  name="${base%.*}"
  echo "     [$i/$N_BINS] Tiara -> $name"
  # -o = archivo de clasificacion por contig; --min_len filtra contigs cortos.
  tiara -i "$f" -o "$TIARA_DIR/${name}.tiara.txt" \
        --min_len "$MIN_CONTIG_LEN" -t "$THREADS" >/dev/null 2>&1 || {
          echo "     ADVERTENCIA: Tiara fallo en $name; se marcara como 'error'." >&2
        }
done
conda deactivate

# Resumen de Tiara + decision por bin (script python)
echo "     Generando resumen de dominio por bin..."
conda activate checkm2   # usamos su python (solo stdlib) para el resumen
python "$SCRIPTS_DIR/summarize_phase1.py" tiara \
    --bins-dir "$BINS_DIR" \
    --bin-ext "$BIN_EXT" \
    --tiara-dir "$TIARA_DIR" \
    --out "$TIARA_DIR/tiara_bin_summary.tsv"

# Copiar (enlazar) al input de CheckM2 solo los bins procariotas que pasaron.
python "$SCRIPTS_DIR/summarize_phase1.py" copy-prok \
    --summary "$TIARA_DIR/tiara_bin_summary.tsv" \
    --bins-dir "$BINS_DIR" \
    --bin-ext "$BIN_EXT" \
    --dest "$CHECKM2_INPUT"

N_PROK=$(find "$CHECKM2_INPUT" -maxdepth 1 -name "*.$BIN_EXT" | wc -l)
echo "     Bins que pasan el filtro por dominio (procariotas): $N_PROK"

if [[ "$N_PROK" -eq 0 ]]; then
  echo "ERROR: ningun bin paso el filtro por dominio. Revisa tiara_bin_summary.tsv" >&2
  conda deactivate; exit 1
fi
conda deactivate

# Paso 2: calidad genomica (CheckM2)
echo
echo "---> PASO 2: Calidad genomica con CheckM2"
conda activate checkm2
echo "     CheckM2 version: $(checkm2 --version 2>&1 | head -n1)"
export CHECKM2DB="$CHECKM2_DB"

checkm2 predict \
    --threads "$THREADS" \
    --input "$CHECKM2_INPUT" \
    --extension ".$BIN_EXT" \
    --output-directory "$CHECKM2_OUT" \
    --database_path "$CHECKM2_DB" \
    --force

echo "     Reporte de calidad: $CHECKM2_OUT/quality_report.tsv"

# Paso 3: seleccion final y reporte
echo
echo "---> PASO 3: Seleccion final y reporte"
python "$SCRIPTS_DIR/summarize_phase1.py" report \
    --bins-dir "$BINS_DIR" \
    --bin-ext "$BIN_EXT" \
    --tiara-summary "$TIARA_DIR/tiara_bin_summary.tsv" \
    --checkm2 "$CHECKM2_OUT/quality_report.tsv" \
    --min-completeness "$MIN_COMPLETENESS" \
    --max-contamination "$MAX_CONTAMINATION" \
    --selected-dir "$SELECTED_DIR" \
    --out-tsv "$RESULTS/phase1_selection.tsv" \
    --out-report "$RESULTS/phase1_report.md" \
    --job-id "${SLURM_JOB_ID:-NA}"
conda deactivate

echo
echo "======================================================================"
echo " FASE 1 COMPLETA - $(date)"
echo "----------------------------------------------------------------------"
echo " Resultados en: $RESULTS"
echo "   - Resumen por bin (dominio):   01_tiara/tiara_bin_summary.tsv"
echo "   - Calidad CheckM2:             02_checkm2/quality_report.tsv"
echo "   - Tabla de seleccion final:    phase1_selection.tsv"
echo "   - Reporte legible:             phase1_report.md"
echo "   - Genomas seleccionados:       phase1_selected_genomes/"
echo "======================================================================"
