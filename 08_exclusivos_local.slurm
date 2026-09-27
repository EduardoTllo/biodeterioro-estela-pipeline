#!/bin/bash
# 08_exclusivos_local.slurm
# FASE 3 (OE3), Etapa 3.5: verificacion local de genes exclusivos (F1-F4).
#
# Un elemento del array por especie seleccionada.
#   F1 + F2  08_exclusivos.py preparar (estructura del gen y anclaje del contig).
#   F4       tblastn de las proteinas candidatas contra TODOS los genomas post-QC
#            descargados de la especie (genomas_ok_<especie>.txt), no solo los
#            <= 50 del pangenoma.
#   F3       blastn de los genes candidatos contra los demas bins de la misma
#            libreria (muestra) de cada bin.
#   evaluar  08_exclusivos.py evaluar: deja listas las consultas de F5 (nr) y F6
#            (contigs huerfanos contra core_nt), que corre 08_exclusivos_remoto.sh
#            en el nodo de login porque necesitan internet.
#
# Requiere en data/bins_libreria/ los 97 bins crudos con su nombre original
# (bin-<N>-<muestra>.fasta, carpeta bins_clean/ de la Fase 1).
#
# Envio desde ~/estela/fase3:
#   N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l)
#   sbatch --array=1-$N scripts/08_exclusivos_local.slurm

#SBATCH --job-name=f3_exclus
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/fase3_exclusivos_%A_%a.out
#SBATCH --error=logs/fase3_exclusivos_%A_%a.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SCRIPTS_DIR="$WORKDIR/scripts"
SEL_DIR="$WORKDIR/results/03_seleccion"
DESCARGA_DIR="$WORKDIR/results/01_descarga"
BAKTA_OUT="$WORKDIR/results/04_bakta"
PANAROO_OUT="$WORKDIR/results/05_panaroo"
EXCL_OUT="$WORKDIR/results/07_exclusivos"
BINS_LIB="$WORKDIR/data/bins_libreria"
MIN_AA=100
MIN_BORDE=100
EVALUE=1e-10

# No editar debajo salvo que sepas lo que haces
THREADS="${SLURM_CPUS_PER_TASK:-16}"
TASK="${SLURM_ARRAY_TASK_ID:?Envia como array: sbatch --array=1-N}"
LINEA="$(tail -n +2 "$SEL_DIR/especies_seleccionadas.tsv" | sed -n "${TASK}p")"
SLUG="$(echo "$LINEA" | cut -f2)"
BINS_ORIG="$(echo "$LINEA" | cut -f9)"
OUT="$EXCL_OUT/$SLUG"
PART="$PANAROO_OUT/$SLUG/particion"
if [[ -z "$SLUG" || ! -s "$PART/genes_bin.tsv" ]]; then
  echo "ERROR: falta la particion de '$SLUG' ($PART/genes_bin.tsv). Corre 06_panaroo.slurm." >&2
  exit 1
fi
mkdir -p "$OUT/db"

echo "======================================================================"
echo " FASE 3 - Exclusivos (F1-F4) | especie: $SLUG | bins: $BINS_ORIG"
echo " Job: ${SLURM_JOB_ID:-NA} | Nodo: $(hostname) | Hilos: $THREADS | Inicio: $(date)"
echo "======================================================================"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate blast

{
  echo "Fase 3 - Exclusivos (local) - $SLUG"
  echo "Fecha: $(date)   Job: ${SLURM_JOB_ID:-NA}   Nodo: $(hostname)"
  echo "F1: >= $MIN_AA aa, >= $MIN_BORDE pb del borde | F3/F4 evalue <= $EVALUE"
  conda list --export | grep -E '^(blast|python)=' || true
} > "$OUT/versions_local.txt"

# F1 + F2
python "$SCRIPTS_DIR/08_exclusivos.py" preparar \
    --particion-dir "$PART" --bakta-dir "$BAKTA_OUT" --outdir "$OUT" \
    --min-aa "$MIN_AA" --min-borde "$MIN_BORDE"

# F4: BD nucleotidica con todos los genomas post-QC de la especie.
# Encabezados reescritos a <ID>__<contig> para saber de que genoma es cada hit.
: > "$OUT/f4_tblastn.tsv"
if [[ -s "$OUT/f1_ok.faa" ]]; then
  if [[ ! -s "$OUT/db/especie.nsq" && ! -s "$OUT/db/especie.00.nsq" ]]; then
    : > "$OUT/db/especie.fna"
    while read -r FNA; do
      [[ -s "$FNA" ]] || continue
      ID="$(basename "$FNA" .fna)"
      awk -v id="$ID" '/^>/{split(substr($0,2),a," "); print ">" id "__" a[1]; next} {print}' \
          "$FNA" >> "$OUT/db/especie.fna"
    done < "$DESCARGA_DIR/genomas_ok_${SLUG}.txt"
    makeblastdb -in "$OUT/db/especie.fna" -dbtype nucl -parse_seqids -out "$OUT/db/especie"
  fi
  tblastn -query "$OUT/f1_ok.faa" -db "$OUT/db/especie" -evalue "$EVALUE" \
          -max_target_seqs 50 -num_threads "$THREADS" \
          -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
          > "$OUT/f4_tblastn.tsv"
fi
echo "F4: $(cut -f1 "$OUT/f4_tblastn.tsv" | sort -u | wc -l) consultas con hit en la especie"

# F3: otros bins de la misma libreria (muestra) de cada bin de la especie
: > "$OUT/f3_blastn.tsv"
for B in ${BINS_ORIG//,/ }; do
  BID="$(echo "$B" | sed 's/[^A-Za-z0-9]/_/g')"
  Q="$OUT/f1_ok_${BID}.ffn"
  [[ -s "$Q" ]] || continue
  MUESTRA="${B##*-}"
  OTROS="$OUT/db/otros_bins_${BID}.fna"
  : > "$OTROS"
  for F in "$BINS_LIB"/bin-*-"${MUESTRA}".fasta; do
    [[ -s "$F" ]] || continue
    NOMBRE="$(basename "$F" .fasta)"
    [[ "$NOMBRE" == "$B" ]] && continue
    awk -v id="$NOMBRE" '/^>/{split(substr($0,2),a," "); print ">" id "__" a[1]; next} {print}' \
        "$F" >> "$OTROS"
  done
  if [[ -s "$OTROS" ]]; then
    blastn -query "$Q" -subject "$OTROS" -evalue "$EVALUE" \
           -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
           >> "$OUT/f3_blastn.tsv"
  else
    echo "F3: $B no tiene otros bins en la libreria $MUESTRA"
  fi
done

# Evaluacion de F4 y F3; consultas para F5 y F6
python "$SCRIPTS_DIR/08_exclusivos.py" evaluar \
    --outdir "$OUT" --f4 "$OUT/f4_tblastn.tsv" --f3 "$OUT/f3_blastn.tsv"

echo
echo "======================================================================"
echo " EXCLUSIVOS (LOCAL) COMPLETO - $SLUG - $(date)"
echo "   Consultas para nr (F5):      $OUT/f5_query.faa ($(grep -c '^>' "$OUT/f5_query.faa" || true))"
echo "   Contigs para core_nt (F6):   $OUT/f6_query.fna ($(grep -c '^>' "$OUT/f6_query.fna" || true))"
echo " Siguiente paso (login): bash scripts/08_exclusivos_remoto.sh"
echo "======================================================================"
