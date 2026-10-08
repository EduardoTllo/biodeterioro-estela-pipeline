#!/bin/bash
# 06c_ani_bins.slurm
# FASE 3 (OE3): ANI de cada bin seleccionado contra TODOS los genomas descargados
# de su especie (genomas_ok_<especie>.txt, no solo las referencias del
# pangenoma). Sirve para ubicar el bin en la especie (cepa nueva o miembro de un
# clon ya secuenciado) y comparar con el arbol del core.
#
# Un solo trabajo para las tres especies (fastANI es rapido).
# Envio desde ~/estela/fase3:
#   sbatch --account=tesis scripts/06c_ani_bins.slurm
#
# Salida: results/06c_ani/<especie>/ani_bin_genomas.tsv
#   consulta  referencia  ANI  fragmentos_compartidos  fragmentos_totales

#SBATCH --job-name=f3_ani
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fase3_ani_%j.out
#SBATCH --error=logs/fase3_ani_%j.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SEL_DIR="$WORKDIR/results/03_seleccion"
DESCARGA_DIR="$WORKDIR/results/01_descarga"
ANI_OUT="$WORKDIR/results/06c_ani"

# No editar debajo salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-16}"
module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate drep
echo "fastANI: $(fastANI --version 2>&1 | head -n1) | $(date) | $(hostname)"

tail -n +2 "$SEL_DIR/bakta_manifest.tsv" | awk -F'\t' '$3 == "bin" {print $4 "\t" $2}' |
while IFS=$'\t' read -r SLUG FASTA; do
  [[ "$FASTA" = /* ]] || FASTA="$WORKDIR/$FASTA"
  OUT="$ANI_OUT/$SLUG"
  mkdir -p "$OUT"
  LISTA="$DESCARGA_DIR/genomas_ok_${SLUG}.txt"
  [[ -s "$LISTA" ]] || { echo "ERROR: falta $LISTA" >&2; exit 1; }
  echo "==> $SLUG: $(basename "$FASTA") contra $(wc -l < "$LISTA") genomas"
  fastANI -q "$FASTA" --rl "$LISTA" -t "$THREADS" -o "$OUT/fastani_crudo.tsv"
  awk -F'\t' 'BEGIN{OFS="\t"; print "consulta","referencia","ani","fragmentos_compartidos","fragmentos_totales"}
       {q=$1; r=$2; sub(/.*\//,"",q); sub(/\.[^.]*$/,"",q); sub(/.*\//,"",r); sub(/\.[^.]*$/,"",r);
        print q, r, $3, $4, $5}' "$OUT/fastani_crudo.tsv" > "$OUT/ani_bin_genomas.tsv"
  echo "    mayor ANI: $(tail -n +2 "$OUT/ani_bin_genomas.tsv" | sort -t$'\t' -k3,3gr | head -n3 | cut -f2,3 | tr '\t\n' ' ;')"
done
echo "Listo: $(date)"
