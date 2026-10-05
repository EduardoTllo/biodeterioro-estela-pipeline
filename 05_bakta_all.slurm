#!/bin/bash
# 05_bakta_all.slurm
# FASE 3 (OE3), Etapa 3.3: anotacion homogenea con Bakta.
#
# Anota todos los genomas del manifiesto (referencias seleccionadas + bins
# propios de las especies del top-3). Todos los genomas se anotan con la MISMA
# version de Bakta, la MISMA BD y las MISMAS opciones (decision D14):
#   - sin --meta, --complete, --partial, --genus ni --species;
#   - sin --partial, Bakta predice genes con extremos cerrados y no anota genes
#     que corren fuera del borde de un contig (lo tiene en cuenta el filtro F1 de
#     genes exclusivos);
#   - --locus-tag = ID del genoma, para que los identificadores no choquen en
#     Panaroo.
# Todas las coordenadas y secuencias posteriores se toman de las salidas de
# Bakta (.gff3, .tsv, .faa, .ffn, .fna), no de los FASTA originales.
#
# Trabajo por LOTES: muchos clusters limitan cuantos trabajos puede tener un
# usuario en cola (en Khipu, cuenta tesis: 5 enviados, 3 corriendo, 32 CPU y
# 98 GB en total). Por eso cada tarea del array anota un lote de genomas, uno
# tras otro: la tarea k de un array 1-LOTES toma las filas k, k+LOTES,
# k+2*LOTES... del manifiesto. Los genomas ya anotados se saltan, asi que si
# una tarea llega al limite de tiempo basta con volver a enviarla. Si un genoma
# falla, el lote sigue con los demas y el fallo queda en el registro.
#
# Envio desde ~/estela/fase3 (crear antes la carpeta logs/). El numero de
# tareas del array debe ser igual a LOTES:
#   Piloto P0 (3 referencias + bins de la especie n.o 1, un solo lote):
#     sbatch --account=tesis --array=1 --export=ALL,LOTES=1,MANIFEST=results/03_seleccion/piloto_manifest.tsv scripts/05_bakta_all.slurm
#   Corrida completa (3 lotes en paralelo = 30 CPU, 90 GB):
#     sbatch --account=tesis --array=1-3 --export=ALL,LOTES=3 scripts/05_bakta_all.slurm

#SBATCH --job-name=f3_bakta
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fase3_bakta_%A_%a.out
#SBATCH --error=logs/fase3_bakta_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
MANIFEST="${MANIFEST:-$WORKDIR/results/03_seleccion/bakta_manifest.tsv}"
LOTES="${LOTES:-3}"
BAKTA_OUT="$WORKDIR/results/04_bakta"
BAKTA_DB="${BAKTA_DB:-$(cat "$WORKDIR/bakta_db_path.txt" 2>/dev/null)}"

# No editar debajo salvo que sepas lo que haces
[[ "$MANIFEST" = /* ]] || MANIFEST="$WORKDIR/$MANIFEST"
THREADS="${SLURM_CPUS_PER_TASK:-10}"
TASK="${SLURM_ARRAY_TASK_ID:?Este script se envia como array: sbatch --array=1-LOTES}"
if [[ ! -s "$MANIFEST" ]]; then
  echo "ERROR: no existe el manifiesto $MANIFEST" >&2; exit 1
fi
if [[ "$TASK" -gt "$LOTES" ]]; then
  echo "ERROR: tarea $TASK mayor que LOTES=$LOTES (usa --array=1-$LOTES)" >&2; exit 1
fi
if [[ -z "$BAKTA_DB" || ! -s "$BAKTA_DB/version.json" ]]; then
  echo "ERROR: BD de Bakta no encontrada ('$BAKTA_DB'). Corre 03_setup_fase3_khipu.sh." >&2; exit 1
fi

N_TOTAL=$(tail -n +2 "$MANIFEST" | grep -c . || true)
echo "======================================================================"
echo " FASE 3 - Bakta | lote $TASK de $LOTES | genomas en el manifiesto: $N_TOTAL"
echo " Job: ${SLURM_JOB_ID:-NA} | Nodo: $(hostname) | Hilos: $THREADS | Inicio: $(date)"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate bakta
mkdir -p "$BAKTA_OUT"
VERSIONES="$(conda list --export | grep -E '^(bakta|pyrodigal|diamond|python)=' || true)"

N_HECHOS=0; N_SALTADOS=0; FALLIDOS=()
FILA=0
while IFS=$'\t' read -r ID FASTA TIPO ESP ORIG; do
  FILA=$((FILA + 1))
  # La tarea k toma las filas k, k+LOTES, k+2*LOTES...
  [[ $(( (FILA - 1) % LOTES + 1 )) -eq "$TASK" ]] || continue
  OUT="$BAKTA_OUT/$ID"
  if [[ -s "$OUT/$ID.gff3" && -s "$OUT/$ID.faa" && -s "$OUT/$ID.json" && -s "$OUT/$ID.tsv" ]]; then
    N_SALTADOS=$((N_SALTADOS + 1))
    continue
  fi
  if [[ -z "$ID" || ! -s "$FASTA" ]]; then
    echo "FALLO: fila $FILA invalida o FASTA inexistente: '$ID' '$FASTA'" >&2
    FALLIDOS+=("$ID"); continue
  fi
  echo "---> [$(date +%H:%M)] $ID ($TIPO)"
  TMP="${TMPDIR:-$WORKDIR/tmp}/bakta_${ID}_${SLURM_JOB_ID:-0}"
  mkdir -p "$TMP"
  if bakta --db "$BAKTA_DB" \
           --output "$OUT" \
           --prefix "$ID" \
           --locus-tag "$ID" \
           --threads "$THREADS" \
           --tmp-dir "$TMP" \
           --skip-plot \
           --force \
           "$FASTA"; then
    {
      echo "Fase 3 - Bakta - $ID"
      echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}_${TASK}   Nodo: $(hostname)"
      echo "$VERSIONES"
      echo "BD: $(tr -d '\n ' < "$BAKTA_DB/version.json")"
    } > "$OUT/versions.txt"
    echo "     listo: $(grep -c $'\tcds\t' "$OUT/$ID.tsv" || true) CDS"
    N_HECHOS=$((N_HECHOS + 1))
  else
    echo "FALLO: Bakta termino con error en $ID" >&2
    FALLIDOS+=("$ID")
  fi
  rm -rf "$TMP"
done < <(tail -n +2 "$MANIFEST")

echo
echo "======================================================================"
echo " LOTE $TASK COMPLETO - $(date)"
echo "   Anotados ahora: $N_HECHOS | Ya estaban: $N_SALTADOS | Fallidos: ${#FALLIDOS[@]}"
[[ ${#FALLIDOS[@]} -gt 0 ]] && echo "   Fallidos: ${FALLIDOS[*]}"
echo "======================================================================"
[[ ${#FALLIDOS[@]} -eq 0 ]]
