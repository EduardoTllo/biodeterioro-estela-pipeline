#!/bin/bash
# 05_bakta_all.slurm
# FASE 3 (OE3), Etapa 3.3: anotacion homogenea con Bakta.
#
# Un elemento del array por genoma del manifiesto (referencias seleccionadas +
# bins propios de las especies del top-3). Todos los genomas se anotan con la
# MISMA version de Bakta, la MISMA BD y las MISMAS opciones (decision D14):
#   - sin --meta, --complete, --partial, --genus ni --species;
#   - sin --partial, Bakta predice genes con extremos cerrados y no anota genes
#     que corren fuera del borde de un contig (lo tiene en cuenta el filtro F1 de
#     genes exclusivos);
#   - --locus-tag = ID del genoma, para que los identificadores no choquen en
#     Panaroo.
# Todas las coordenadas y secuencias posteriores se toman de las salidas de
# Bakta (.gff3, .tsv, .faa, .ffn, .fna), no de los FASTA originales.
#
# Envio desde ~/estela/fase3 (crear antes la carpeta logs/):
#   Piloto P0 (3 referencias + bins de la especie n.o 1):
#     N=$(tail -n +2 results/03_seleccion/piloto_manifest.tsv | wc -l)
#     sbatch --array=1-$N --export=ALL,MANIFEST=results/03_seleccion/piloto_manifest.tsv scripts/05_bakta_all.slurm
#   Corrida completa (maximo 4 tareas simultaneas = 32 cores, 96 GB):
#     N=$(tail -n +2 results/03_seleccion/bakta_manifest.tsv | wc -l)
#     sbatch --array=1-$N%4 scripts/05_bakta_all.slurm
# Los genomas ya anotados se omiten, asi que el piloto no se repite.

#SBATCH --job-name=f3_bakta
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fase3_bakta_%A_%a.out
#SBATCH --error=logs/fase3_bakta_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
MANIFEST="${MANIFEST:-$WORKDIR/results/03_seleccion/bakta_manifest.tsv}"
BAKTA_OUT="$WORKDIR/results/04_bakta"
BAKTA_DB="${BAKTA_DB:-$(cat "$WORKDIR/bakta_db_path.txt" 2>/dev/null)}"

# No editar debajo salvo que sepas lo que haces
[[ "$MANIFEST" = /* ]] || MANIFEST="$WORKDIR/$MANIFEST"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
TASK="${SLURM_ARRAY_TASK_ID:?Este script se envia como array: sbatch --array=1-N}"
LINEA="$(tail -n +2 "$MANIFEST" | sed -n "${TASK}p")"
ID="$(echo "$LINEA" | cut -f1)"
FASTA="$(echo "$LINEA" | cut -f2)"
TIPO="$(echo "$LINEA" | cut -f3)"
if [[ -z "$ID" || ! -s "$FASTA" ]]; then
  echo "ERROR: fila $TASK del manifiesto invalida o FASTA inexistente: '$ID' '$FASTA'" >&2; exit 1
fi
if [[ -z "$BAKTA_DB" || ! -s "$BAKTA_DB/version.json" ]]; then
  echo "ERROR: BD de Bakta no encontrada ('$BAKTA_DB'). Corre 03_setup_fase3_khipu.sh." >&2; exit 1
fi

OUT="$BAKTA_OUT/$ID"
if [[ -s "$OUT/$ID.gff3" && -s "$OUT/$ID.faa" && -s "$OUT/$ID.json" && -s "$OUT/$ID.tsv" ]]; then
  echo "$ID ya anotado. Se omite."
  exit 0
fi

echo "======================================================================"
echo " FASE 3 - Bakta | $ID ($TIPO)"
echo " Job: ${SLURM_JOB_ID:-NA} (tarea $TASK) | Nodo: $(hostname) | Hilos: $THREADS"
echo " FASTA: $FASTA | Inicio: $(date)"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate bakta

TMP="${TMPDIR:-$WORKDIR/tmp}/bakta_${ID}_${SLURM_JOB_ID:-0}"
mkdir -p "$TMP" "$BAKTA_OUT"
trap 'rm -rf "$TMP"' EXIT

bakta --db "$BAKTA_DB" \
      --output "$OUT" \
      --prefix "$ID" \
      --locus-tag "$ID" \
      --threads "$THREADS" \
      --tmp-dir "$TMP" \
      --skip-plot \
      --force \
      "$FASTA"

{
  echo "Fase 3 - Bakta - $ID"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}_${TASK}   Nodo: $(hostname)"
  conda list --export | grep -E '^(bakta|pyrodigal|diamond|python)=' || true
  echo "BD: $(tr -d '\n ' < "$BAKTA_DB/version.json")"
} > "$OUT/versions.txt"

echo "Listo: $(grep -c $'\tcds\t' "$OUT/$ID.tsv" || true) CDS | $(date)"
