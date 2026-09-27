#!/bin/bash
# 04_drep_refs.slurm
# FASE 3 (OE3), Etapa 3.2b: desreplicacion de las referencias al 99 % ANI.
#
# Un elemento del array por especie (fila de especies_para_drep.tsv). Solo entran
# las referencias descargadas y verificadas; los bins propios NO se desreplican.
#   - genomeInfo: completitud/contaminacion de CheckM2 tomadas del metadata de
#     GTDB (no se corre CheckM2).
#   - extraW: T0 +1000, T1/T2 +500, T3 0, para que dentro de cada cluster al 99 %
#     gane el genoma prioritario (decision D23).
#   - Mismo algoritmo secundario que la Fase 2 (fastANI) y mismo entorno 'drep'.
# El numero de clusters secundarios es el N de genomas no redundantes que decide
# la viabilidad (>= 15) y el ranking de especies.
#
# Antes de enviar, en el login:
#   python scripts/04_clasificar_refs.py clasificar ...   (ver README_fase3.md)
# Envio (N = numero de especies en especies_para_drep.tsv):
#   N=$(tail -n +2 results/01_descarga/especies_para_drep.tsv | wc -l)
#   sbatch --array=1-$N scripts/04_drep_refs.slurm

#SBATCH --job-name=f3_drep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fase3_drep_%A_%a.out
#SBATCH --error=logs/fase3_drep_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
DESCARGA_DIR="$WORKDIR/results/01_descarga"
DREP_DIR="$WORKDIR/results/02_drep"
DREP_PA=0.90        # umbral primario (Mash), igual que la Fase 2
DREP_SA=0.99        # umbral secundario (ANI): genomas casi clonales

# No editar debajo salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-16}"
TASK="${SLURM_ARRAY_TASK_ID:?Este script se envia como array: sbatch --array=1-N}"
LINEA="$(tail -n +2 "$DESCARGA_DIR/especies_para_drep.tsv" | sed -n "${TASK}p")"
SLUG="$(echo "$LINEA" | cut -f1)"
if [[ -z "$SLUG" ]]; then
  echo "ERROR: no hay especie en la fila $TASK de especies_para_drep.tsv" >&2; exit 1
fi
IN="$DREP_DIR/$SLUG/input"
OUT="$DREP_DIR/$SLUG"
for f in genomas.txt genomeInfo.csv extraW.tsv; do
  if [[ ! -s "$IN/$f" ]]; then
    echo "ERROR: falta $IN/$f (corre 04_clasificar_refs.py clasificar)." >&2; exit 1
  fi
done

echo "======================================================================"
echo " FASE 3 - dRep ${DREP_SA} ANI | especie: $SLUG"
echo " Job: ${SLURM_JOB_ID:-NA} (tarea $TASK) | Nodo: $(hostname) | Hilos: $THREADS"
echo " Genomas: $(grep -c . "$IN/genomas.txt") | Inicio: $(date)"
echo "======================================================================"

if [[ -s "$OUT/data_tables/Wdb.csv" ]]; then
  echo "Ya existe $OUT/data_tables/Wdb.csv. Se omite (borra la carpeta para rehacer)."
  exit 0
fi

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate drep

{
  echo "Fase 3 - dRep referencias - $SLUG"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}_${TASK}   Nodo: $(hostname)"
  echo "Parametros: -pa $DREP_PA -sa $DREP_SA --S_algorithm fastANI -comp 0 -con 100"
  conda list --export | grep -E '^(drep|fastani|mash|python)=' || true
} > "$OUT/versions.txt"

# -comp 0 -con 100: el QC ya se aplico en el censo (>= 95 % / <= 5 %).
dRep dereplicate "$OUT" \
    -g "$IN/genomas.txt" \
    -p "$THREADS" \
    -pa "$DREP_PA" -sa "$DREP_SA" \
    --S_algorithm fastANI \
    -comp 0 -con 100 \
    --genomeInfo "$IN/genomeInfo.csv" \
    -extraW "$IN/extraW.tsv"

N_NR=$(( $(wc -l < "$OUT/data_tables/Wdb.csv") - 1 ))
echo
echo "======================================================================"
echo " dRep COMPLETO - $SLUG - $(date)"
echo " Genomas no redundantes (clusters al ${DREP_SA}): $N_NR"
echo "======================================================================"
