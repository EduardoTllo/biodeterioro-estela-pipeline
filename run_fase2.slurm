#!/bin/bash
# =============================================================================
# run_fase2.slurm
# -----------------------------------------------------------------------------
# FASE 2 - Asignacion taxonomica y definicion de linajes (OE2)
# Tesis: Estela de Raimondi | Cluster: Khipu (UTEC) | SLURM 25.11
#
# Flujo (solo sobre los genomas SELECCIONADOS en la Fase 1):
#   Paso 1 (GTDB-Tk classify_wf) -> taxonomia + ANI/AF + arbol de colocacion.
#   Paso 2 (dRep)                -> desreplicacion; linaje = cluster de especie
#                                   a 95 % de ANI. Reusa la calidad de CheckM2
#                                   (Fase 1) via genomeInfo (no re-corre CheckM).
#   Paso 3 (resumen)             -> prevalencia espacial, confirmacion de especie,
#                                   top-3 linajes para la Fase 3 y reporte legible.
#
# NOTA: el arbol de novo FOCALIZADO se hara en la Fase 4 (no aqui), solo con los
#       linajes seleccionados. En Fase 2 se usa el arbol de colocacion.
#
# Envio:   sbatch run_fase2.slurm
# =============================================================================
# Recursos dentro de los limites de la cuenta (max: 32 cores, 98 GB RAM, 24 h):
#SBATCH --job-name=fase2_estela
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --output=fase2_%j.out
#SBATCH --error=fase2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=diegotello2000@hotmail.com

# No usamos 'set -u' por compatibilidad con Lmod en Khipu (igual que Fase 1).
set -eo pipefail

# =============================================================================
# CONFIGURACION  (edita segun tu cuenta en Khipu)
# =============================================================================
# Carpeta de la Fase 1 (para tomar los genomas seleccionados y la calidad):
FASE1_DIR="$HOME/estela/fase1"
SELECTED_DIR="$FASE1_DIR/results/phase1_selected_genomes"
CHECKM2_REPORT="$FASE1_DIR/results/02_checkm2/quality_report.tsv"

# Carpeta de trabajo de la Fase 2:
WORKDIR="$HOME/estela/fase2"

# Carpeta con los scripts (summarize_phase2.py):
SCRIPTS_DIR="$WORKDIR/scripts"

# Base de datos de GTDB-Tk R220 (carpeta con markers/, msa/, taxonomy/, ...):
GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-$HOME/gtdbtk_data/release220}"

# Extension de los genomas:
BIN_EXT="fasta"

# Umbrales (tesis, seccion 4.1.4):
DREP_PA=0.90        # umbral primario (Mash)
DREP_SA=0.95        # umbral secundario (ANI) -> linaje = especie
ANI_ESPECIE=95      # confirmacion de especie: ANI >=
AF_ESPECIE=65       # confirmacion de especie: AF  >=
TOP_LINAJES=3       # linajes priorizados para la Fase 3

# =============================================================================
# No editar debajo salvo que sepas lo que haces
# =============================================================================
THREADS="${SLURM_CPUS_PER_TASK:-32}"
RESULTS="$WORKDIR/results"
GTDBTK_OUT="$RESULTS/03_gtdbtk"
DREP_OUT="$RESULTS/04_drep"
DREP_INPUT="$RESULTS/drep_input"
SEL_LINEAGES="$RESULTS/phase2_selected_lineages"

echo "======================================================================"
echo " FASE 2 - Estela de Raimondi"
echo " Job SLURM: ${SLURM_JOB_ID:-NA}  |  Nodo: $(hostname)  |  Hilos: $THREADS"
echo " Inicio: $(date)"
echo "======================================================================"

# --- Comprobaciones previas --------------------------------------------------
if [[ ! -d "$SELECTED_DIR" ]]; then
  echo "ERROR: no existe la carpeta de genomas seleccionados: $SELECTED_DIR" >&2; exit 1
fi
shopt -s nullglob
GENOMES=("$SELECTED_DIR"/*."$BIN_EXT")
N_GEN=${#GENOMES[@]}
if [[ "$N_GEN" -eq 0 ]]; then
  echo "ERROR: no hay genomas *.$BIN_EXT en $SELECTED_DIR" >&2; exit 1
fi
if [[ ! -f "$CHECKM2_REPORT" ]]; then
  echo "ERROR: no existe el reporte de CheckM2: $CHECKM2_REPORT" >&2; exit 1
fi
if [[ ! -d "$GTDBTK_DATA_PATH" ]]; then
  echo "ERROR: no existe GTDBTK_DATA_PATH: $GTDBTK_DATA_PATH" >&2
  echo "       Corre 02_setup_fase2_khipu.sh o ajusta la ruta." >&2; exit 1
fi
echo "Genomas de entrada (seleccionados en Fase 1): $N_GEN"

mkdir -p "$GTDBTK_OUT" "$DREP_OUT" "$DREP_INPUT" "$SEL_LINEAGES"

# --- Entorno conda -----------------------------------------------------------
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

# --- Registro de versiones ---------------------------------------------------
VERS="$RESULTS/versions.txt"
{
  echo "Fase 2 - registro de versiones"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  echo "N genomas entrada: $N_GEN"
  echo "GTDBTK_DATA_PATH: $GTDBTK_DATA_PATH"
  echo "Umbrales: dRep -pa $DREP_PA -sa $DREP_SA | especie ANI>=$ANI_ESPECIE AF>=$AF_ESPECIE | top=$TOP_LINAJES"
  echo "----"
  conda activate gtdbtk 2>/dev/null && echo "GTDB-Tk: $(gtdbtk --version 2>&1 | head -n1)" && conda deactivate || true
  conda activate drep   2>/dev/null && echo "dRep:    $(dRep --version 2>&1 | head -n1)"   && conda deactivate || true
} > "$VERS"
echo "Versiones registradas en: $VERS"

# =============================================================================
# PASO 1 - CLASIFICACION TAXONOMICA (GTDB-Tk classify_wf)
# =============================================================================
echo
echo "---> PASO 1: GTDB-Tk classify_wf"
conda activate gtdbtk
export GTDBTK_DATA_PATH="$GTDBTK_DATA_PATH"
echo "     GTDB-Tk: $(gtdbtk --version 2>&1 | head -n1)"

# --skip_ani_screen: omite el pre-filtro Mash (no necesita mash_db); la especie
#   igual se confirma con FastANI interno en el paso de clasificacion.
# --pplacer_cpus 1: limita la memoria de pplacer (bac120 ~55 GB con 1 cpu).
gtdbtk classify_wf \
    --genome_dir "$SELECTED_DIR" \
    --extension "$BIN_EXT" \
    --out_dir "$GTDBTK_OUT" \
    --cpus "$THREADS" \
    --pplacer_cpus 1 \
    --skip_ani_screen

echo "     Resumen(es): $(ls "$GTDBTK_OUT"/classify/*summary.tsv 2>/dev/null | tr '\n' ' ')"
conda deactivate

# =============================================================================
# PASO 2 - DESREPLICACION Y DEFINICION DE LINAJES (dRep)
# =============================================================================
echo
echo "---> PASO 2: dRep (linaje = especie a ${DREP_SA} ANI)"
conda activate drep
echo "     dRep: $(dRep --version 2>&1 | head -n1)"

# Enlazar/copiar genomas al input de dRep.
for g in "${GENOMES[@]}"; do cp -f "$g" "$DREP_INPUT/"; done

# genomeInfo.csv desde CheckM2 (Fase 1) para que dRep NO re-ejecute CheckM.
python "$SCRIPTS_DIR/summarize_phase2.py" genomeinfo \
    --genomes-dir "$DREP_INPUT" \
    --bin-ext "$BIN_EXT" \
    --checkm2 "$CHECKM2_REPORT" \
    --out "$RESULTS/genomeInfo.csv"

# -comp 0 -con 100: se desactiva el filtro de calidad de dRep porque ya se aplico
#   en la Fase 1 (>70 % completitud, <5 % contaminacion). Asi no se pierde ningun
#   genoma seleccionado.
dRep dereplicate "$DREP_OUT" \
    -g "$DREP_INPUT"/*."$BIN_EXT" \
    -p "$THREADS" \
    -pa "$DREP_PA" -sa "$DREP_SA" \
    --S_algorithm fastANI \
    -comp 0 -con 100 \
    --genomeInfo "$RESULTS/genomeInfo.csv"

echo "     Clusteres (linajes): $DREP_OUT/data_tables/Cdb.csv"
conda deactivate

# =============================================================================
# PASO 3 - PREVALENCIA, SELECCION Y REPORTE
# =============================================================================
echo
echo "---> PASO 3: prevalencia espacial, top-${TOP_LINAJES} y reporte"
conda activate drep   # su python (stdlib) sirve para el resumen
python "$SCRIPTS_DIR/summarize_phase2.py" report \
    --genomes-dir "$SELECTED_DIR" \
    --bin-ext "$BIN_EXT" \
    --gtdbtk-dir "$GTDBTK_OUT" \
    --drep-cdb "$DREP_OUT/data_tables/Cdb.csv" \
    --drep-wdb "$DREP_OUT/data_tables/Wdb.csv" \
    --min-ani "$ANI_ESPECIE" \
    --min-af "$AF_ESPECIE" \
    --top "$TOP_LINAJES" \
    --selected-lineages-dir "$SEL_LINEAGES" \
    --out-tsv "$RESULTS/phase2_selection.tsv" \
    --out-report "$RESULTS/phase2_report.md" \
    --job-id "${SLURM_JOB_ID:-NA}"
conda deactivate

echo
echo "======================================================================"
echo " FASE 2 COMPLETA - $(date)"
echo "----------------------------------------------------------------------"
echo " Resultados en: $RESULTS"
echo "   - Taxonomia GTDB-Tk:        03_gtdbtk/classify/*summary.tsv"
echo "   - Arbol de colocacion:      03_gtdbtk/classify/*.classify.tree"
echo "   - Clusteres/linajes dRep:   04_drep/data_tables/Cdb.csv"
echo "   - Tabla por genoma:         phase2_selection.tsv"
echo "   - Reporte legible:          phase2_report.md"
echo "   - Genomas top-3 (Fase 3):   phase2_selected_lineages/"
echo "======================================================================"
