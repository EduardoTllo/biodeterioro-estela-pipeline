#!/bin/bash
# 10_dram.slurm
# FASE 4 (OE4): anotacion funcional con DRAM de las proteinas de los bins.
#
# DRAM se corre sobre las proteinas que predijo Bakta (annotate_genes), no con
# su propia prediccion de genes: asi cada gen conserva el locus tag de Bakta y
# se puede cruzar con la particion del pangenoma (genes_bin.tsv de la Fase 3).
# Entrada: results/04_bakta/<bin>/<bin>.faa de los bins de la Fase 1 (los 21
# seleccionados; los 3 del pangenoma ya estan anotados, los demas se anotan con
# 05_bakta_all.slurm y el manifiesto bins21_manifest.tsv).
#   1. DRAM.py annotate_genes -> KO (KOfam), CAZymes (dbCAN), peptidasas (MEROPS)
#   2. DRAM.py distill        -> resumen de modulos KEGG, ciclos C/N/S y CAZymes
#
# Envio desde ~/estela/fase3:
#   sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" scripts/10_dram.slurm

#SBATCH --job-name=f4_dram
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fase4_dram_%j.out
#SBATCH --error=logs/fase4_dram_%j.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
BAKTA_OUT="$WORKDIR/results/04_bakta"
BINS_DIR="$WORKDIR/data/bins"
OUT="${OUT:-$HOME/estela/fase4/results/dram}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Proteinas de Bakta de cada bin (un archivo por genoma; DRAM usa el nombre del
# archivo como nombre del genoma).
ENTRADA="$OUT/entrada_faa"
mkdir -p "$ENTRADA"
FALTAN=0
for f in "$BINS_DIR"/*.fasta; do
  B="$(basename "$f" .fasta)"
  ID="$(echo "$B" | sed 's/[^A-Za-z0-9]/_/g')"
  FAA="$BAKTA_OUT/$ID/$ID.faa"
  if [[ -s "$FAA" ]]; then
    cp -f "$FAA" "$ENTRADA/$B.faa"
  else
    echo "FALTA la anotacion de Bakta de $B ($FAA)" >&2; FALTAN=$((FALTAN + 1))
  fi
done
if [[ "$FALTAN" -gt 0 ]]; then
  echo "ERROR: faltan $FALTAN bins por anotar con Bakta (ver README_fase4.md)." >&2; exit 1
fi

echo "======================================================================"
echo " FASE 4 - DRAM | bins: $(ls "$ENTRADA"/*.faa | wc -l) | Hilos: $THREADS"
echo " Job: ${SLURM_JOB_ID:-NA} | Nodo: $(hostname) | Inicio: $(date)"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate dram

{
  echo "Fase 4 - DRAM"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  conda list --export | grep -iE '^(dram-bio|hmmer|mmseqs2|python)=' || true
  DRAM-setup.py print_config 2>/dev/null | grep -iE "kofam|dbcan|merops|peptidase" || true
} > "$OUT/versions.txt"

if [[ ! -s "$OUT/anotacion/annotations.tsv" ]]; then
  rm -rf "$OUT/anotacion"
  DRAM.py annotate_genes -i "$ENTRADA/*.faa" -o "$OUT/anotacion" --threads "$THREADS"
fi

rm -rf "$OUT/resumen"
DRAM.py distill -i "$OUT/anotacion/annotations.tsv" -o "$OUT/resumen"

echo
echo "======================================================================"
echo " DRAM COMPLETO - $(date)"
echo "   Anotaciones por gen: $OUT/anotacion/annotations.tsv"
echo "   Resumen:             $OUT/resumen/ (metabolism_summary.xlsx, product.tsv, genome_stats.tsv)"
echo "======================================================================"
