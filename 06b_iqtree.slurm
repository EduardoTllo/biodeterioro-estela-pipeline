#!/bin/bash
# 06b_iqtree.slurm
# FASE 3 (OE3), Etapa 3.4: arbol de maxima verosimilitud del core por especie
# (decision D28). Ubica el bin de la Estela entre las cepas de su especie.
#
# Entrada: core_gene_alignment_filtered.aln de Panaroo (genes en >= 95 % de los
# genomas, bin incluido; se toleran huecos).
#   MODO=completo (defecto): IQ-TREE 3 con ModelFinder (-m MFP) y 1000 replicas
#                            de UFBoot sobre el alineamiento completo.
#   MODO=snps:               si el anterior no cabe en 24 h, snp-sites -c y
#                            IQ-TREE con -m MFP+ASC (correccion por sesgo de
#                            muestreo de sitios variables).
# El modo usado queda registrado en versions.txt. Limitacion declarada: no se
# corrige por recombinacion.
#
# Envio desde ~/estela/fase3:
#   N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l)
#   sbatch --array=1-$N%1 scripts/06b_iqtree.slurm
#   sbatch --array=2 --export=ALL,MODO=snps scripts/06b_iqtree.slurm   # respaldo

#SBATCH --job-name=f3_iqtree
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fase3_iqtree_%A_%a.out
#SBATCH --error=logs/fase3_iqtree_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SEL_DIR="$WORKDIR/results/03_seleccion"
PANAROO_OUT="$WORKDIR/results/05_panaroo"
IQTREE_OUT="$WORKDIR/results/06_iqtree"
MODO="${MODO:-completo}"
UFBOOT=1000
SEED=12345

# No editar debajo salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-32}"
TASK="${SLURM_ARRAY_TASK_ID:?Envia como array: sbatch --array=1-N%1}"
SLUG="$(tail -n +2 "$SEL_DIR/especies_seleccionadas.tsv" | sed -n "${TASK}p" | cut -f2)"
ALN="$PANAROO_OUT/$SLUG/core_gene_alignment_filtered.aln"
OUT="$IQTREE_OUT/$SLUG"
if [[ -z "$SLUG" || ! -s "$ALN" ]]; then
  echo "ERROR: falta el alineamiento del core de '$SLUG': $ALN" >&2; exit 1
fi
mkdir -p "$OUT"

echo "======================================================================"
echo " FASE 3 - IQ-TREE ($MODO) | especie: $SLUG"
echo " Job: ${SLURM_JOB_ID:-NA} | Nodo: $(hostname) | Hilos: $THREADS | Inicio: $(date)"
echo " Secuencias en el alineamiento: $(grep -c '^>' "$ALN")"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate iqtree

{
  echo "Fase 3 - IQ-TREE - $SLUG"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  echo "Modo: $MODO   UFBoot: $UFBOOT   Semilla: $SEED"
  conda list --export | grep -E '^(iqtree|snp-sites)=' || true
} > "$OUT/versions_${MODO}.txt"

if [[ "$MODO" == "snps" ]]; then
  snp-sites -c -o "$OUT/core_snps.aln" "$ALN"
  echo "Sitios variables: $(awk '/^>/{if(n++)exit;next}{l+=length($0)}END{print l}' "$OUT/core_snps.aln")"
  iqtree3 -s "$OUT/core_snps.aln" -m MFP+ASC -B "$UFBOOT" -T AUTO --threads-max "$THREADS" \
          --seed "$SEED" --prefix "$OUT/core_snps" -redo
  cp -f "$OUT/core_snps.treefile" "$OUT/core.treefile"
else
  iqtree3 -s "$ALN" -m MFP -B "$UFBOOT" -T AUTO --threads-max "$THREADS" \
          --seed "$SEED" --prefix "$OUT/core" -redo
fi
echo "$MODO" > "$OUT/modo_usado.txt"

echo
echo "======================================================================"
echo " IQ-TREE COMPLETO - $SLUG - $(date)"
echo "   Arbol:  $OUT/core.treefile"
echo "   Modelo: $(grep -m1 'Best-fit model' "$OUT"/core*.iqtree 2>/dev/null || echo 'ver .iqtree')"
echo "======================================================================"
