#!/bin/bash
# 06_panaroo.slurm
# FASE 3 (OE3), Etapa 3.4: pangenoma por especie con Panaroo + particion
# calculada solo con referencias (07_particion.py).
#
# Un elemento del array por especie (fila de especies_seleccionadas.tsv).
#   1. Convierte cada GFF3 de Bakta al formato Prokka que Panaroo lee, con el
#      conversor oficial de Panaroo 1.8.0 (convert_bakta_to_prokka_gff.py). El
#      conversor descarta CDS con stop interno, largo no multiplo de 3 o < 34 pb
#      (incluye los CDS con selenocisteina, issue #373 de Panaroo). Lo descartado
#      se registra en conversion_descartes.tsv.
#   2. Panaroo --clean-mode moderate (decision D15): strict exige que un gen
#      posiblemente contaminante este en >= 5 % de los genomas (3 de 51) y
#      eliminaria los genes presentes solo en el bin; moderate exige 1 %.
#      Parametros de agrupamiento en sus valores por defecto (-c 0.98, -f 0.7,
#      --len_dif_percent 0.98, --refind-mode default).
#      --core_threshold 0.95 aqui solo define el ALINEAMIENTO del core para el
#      arbol (06b_iqtree.slurm); la particion que se reporta la recalcula
#      07_particion.py solo con referencias (decision D22).
#   3. 07_particion.py particion: core/shell/cloud/exclusivo con frecuencias de
#      referencias, posicion de cada bin y recuperacion del core.
#
# Envio desde ~/estela/fase3:
#   Piloto P0:  sbatch --export=ALL,PILOTO=1 scripts/06_panaroo.slurm
#   Completo:   N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l)
#               sbatch --array=1-$N%1 scripts/06_panaroo.slurm
#
# Analisis de sensibilidad del modo de limpieza (08-10-2026):
#   sbatch --array=1-3%1 --export=ALL,CLEAN_MODE=sensitive scripts/06_panaroo.slurm
#   Con un modo distinto de moderate escribe en results/05_panaroo_<modo>/, reutiliza
#   los GFF ya convertidos de results/05_panaroo/ y no alinea el core (ALINEAR=0 por
#   defecto en ese caso; el arbol se hace solo con moderate). En sensitive Panaroo
#   no poda los genes de extremo de contig (trailing_recursive = 0).

#SBATCH --job-name=f3_panaroo
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fase3_panaroo_%A_%a.out
#SBATCH --error=logs/fase3_panaroo_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SCRIPTS_DIR="$WORKDIR/scripts"
SEL_DIR="$WORKDIR/results/03_seleccion"
BAKTA_OUT="$WORKDIR/results/04_bakta"
PANAROO_BASE="$WORKDIR/results/05_panaroo"
CHECKM2_REPORT="${CHECKM2_REPORT:-$HOME/estela/fase1/results/02_checkm2/quality_report.tsv}"
CLEAN_MODE="${CLEAN_MODE:-moderate}"
CORE_ALN_THRESHOLD=0.95
if [[ "$CLEAN_MODE" == "moderate" ]]; then
  PANAROO_OUT="$PANAROO_BASE"
  ALINEAR="${ALINEAR:-1}"
else
  PANAROO_OUT="$WORKDIR/results/05_panaroo_$CLEAN_MODE"
  ALINEAR="${ALINEAR:-0}"
fi

# No editar debajo salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-32}"
if [[ "${PILOTO:-0}" == "1" ]]; then
  MANIFEST="$SEL_DIR/piloto_manifest.tsv"
  SLUG="$(tail -n +2 "$MANIFEST" | head -n1 | cut -f4)"
  OUT="$PANAROO_OUT/_piloto_$SLUG"
else
  MANIFEST="$SEL_DIR/bakta_manifest.tsv"
  TASK="${SLURM_ARRAY_TASK_ID:?Envia como array: sbatch --array=1-N%1 (o PILOTO=1)}"
  SLUG="$(tail -n +2 "$SEL_DIR/especies_seleccionadas.tsv" | sed -n "${TASK}p" | cut -f2)"
  OUT="$PANAROO_OUT/$SLUG"
fi
if [[ -z "$SLUG" ]]; then
  echo "ERROR: no se pudo determinar la especie" >&2; exit 1
fi
mkdir -p "$OUT/gff"

echo "======================================================================"
echo " FASE 3 - Panaroo ($CLEAN_MODE) | especie: $SLUG | piloto: ${PILOTO:-0}"
echo " Job: ${SLURM_JOB_ID:-NA} | Nodo: $(hostname) | Hilos: $THREADS | Inicio: $(date)"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate panaroo

# 1) Conversion de GFF3 de Bakta -> formato Prokka
LISTA="$OUT/gff_list.txt"
DESC="$OUT/conversion_descartes.tsv"
: > "$LISTA"
printf "id\ttipo\tcds_bakta\tcds_convertidos\tcds_descartados\n" > "$DESC"
N_GEN=0
while IFS=$'\t' read -r ID FASTA TIPO ESP ORIG; do
  [[ "$ESP" == "$SLUG" ]] || continue
  GFF3="$BAKTA_OUT/$ID/$ID.gff3"
  if [[ ! -s "$GFF3" ]]; then
    echo "ERROR: falta $GFF3 (Bakta no termino para $ID)" >&2; exit 1
  fi
  GFF="$OUT/gff/$ID.gff"
  GFF_BASE="$PANAROO_BASE/$SLUG/gff/$ID.gff"
  if [[ ! -s "$GFF" && "$OUT" != "$PANAROO_BASE/$SLUG" && -s "$GFF_BASE" ]]; then
    cp "$GFF_BASE" "$GFF"                  # misma conversion que la corrida moderate
  elif [[ ! -s "$GFF" ]]; then
    python "$SCRIPTS_DIR/convert_bakta_to_prokka_gff.py" -g "$GFF3" -o "$GFF"
  fi
  N_B=$(awk -F'\t' '/^##FASTA/{exit} $3=="CDS"' "$GFF3" | wc -l)
  N_C=$(awk -F'\t' '/^##FASTA/{exit} $3=="CDS"' "$GFF" | wc -l)
  printf "%s\t%s\t%s\t%s\t%s\n" "$ID" "$TIPO" "$N_B" "$N_C" "$((N_B - N_C))" >> "$DESC"
  echo "$GFF" >> "$LISTA"
  N_GEN=$((N_GEN + 1))
done < <(tail -n +2 "$MANIFEST")
echo "Genomas convertidos: $N_GEN (descartes por genoma en $DESC)"
if [[ "$N_GEN" -lt 3 ]]; then
  echo "ERROR: menos de 3 genomas para $SLUG" >&2; exit 1
fi

{
  echo "Fase 3 - Panaroo - $SLUG"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  if [[ "$ALINEAR" == "1" ]]; then
    echo "Parametros: --clean-mode $CLEAN_MODE --remove-invalid-genes -a core --aligner mafft --core_threshold $CORE_ALN_THRESHOLD (resto por defecto)"
  else
    echo "Parametros: --clean-mode $CLEAN_MODE --remove-invalid-genes (sin alineamiento; resto por defecto)"
  fi
  echo "Genomas: $N_GEN"
  conda list --export | grep -E '^(panaroo|mafft|cd-hit|python)=' || true
} > "$OUT/versions.txt"

# 2) Panaroo
if [[ -s "$OUT/summary_statistics.txt" && -s "$OUT/gene_presence_absence.csv" ]]; then
  echo "Panaroo ya corrio para $SLUG. Se omite (borra la carpeta para rehacer)."
elif [[ "$ALINEAR" == "1" ]]; then
  panaroo -i "$LISTA" -o "$OUT" \
      --clean-mode "$CLEAN_MODE" \
      --remove-invalid-genes \
      -a core --aligner mafft \
      --core_threshold "$CORE_ALN_THRESHOLD" \
      -t "$THREADS"
else
  panaroo -i "$LISTA" -o "$OUT" \
      --clean-mode "$CLEAN_MODE" \
      --remove-invalid-genes \
      -t "$THREADS"
fi
cat "$OUT/summary_statistics.txt"; echo

# 3) Particion solo con referencias + posicion del bin
python "$SCRIPTS_DIR/07_particion.py" particion \
    --panaroo-dir "$OUT" \
    --manifest "$MANIFEST" \
    --especie "$SLUG" \
    --bakta-dir "$BAKTA_OUT" \
    --checkm2 "$CHECKM2_REPORT" \
    --outdir "$OUT/particion"

echo
echo "======================================================================"
echo " PANAROO COMPLETO - $SLUG - $(date)"
echo "   Pangenoma:        $OUT/gene_presence_absence.csv"
echo "   Alineamiento core: $OUT/core_gene_alignment_filtered.aln"
echo "   Particion:        $OUT/particion/"
echo " Checkpoint C6: revisar particion/recuperacion_core.tsv"
echo "======================================================================"
