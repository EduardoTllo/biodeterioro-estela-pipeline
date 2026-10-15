#!/bin/bash
# eucariotas/run_fase1_euk.slurm
# FASE 1 EUCARIOTA - Control de calidad de los bins eucariotas (anexo)
# Tesis: Estela de Raimondi | Cluster: Khipu (UTEC) | SLURM 25.11
#
# Flujo:
#   Paso 1 (Tiara)  -> se reutiliza la clasificacion de la Fase 1 bacteriana
#                      (mismos 97 bins, min contig 3000 bp). Se toman los bins
#                      con mayoria de longitud eucariota (31). No se re-corre.
#   Paso 2 (EukCC)  -> completitud y contaminacion con marcadores eucariotas.
#                      Se conservan: completitud > 70 % y contaminacion < 5 %.
#   Paso 3 (BUSCO)  -> medida complementaria (fungi_odb12.2 y eukaryota_odb12.2);
#                      no decide la seleccion. Sus genes de copia unica se usan
#                      para el arbol de la Fase 2.
#   Paso 4 (resumen)-> tabla, reporte y copia de los seleccionados.
#
# El job es reanudable: si se corta, se reenvia y salta los bins ya procesados.
#
# Envio (desde ~/estela/euk):
#   sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" \
#          ~/estela/fase3/scripts/eucariotas/run_fase1_euk.slurm

#SBATCH --job-name=fase1_euk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=fase1_euk_%j.out
#SBATCH --error=fase1_euk_%j.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
SCRIPTS_DIR="$HOME/estela/fase3/scripts"
BINS_DIR="$HOME/estela/fase1/bins"          # los 97 bins crudos
BIN_EXT="fasta"
TIARA_SUMMARY="$SCRIPTS_DIR/results/fase1/phase1_tiara_bin_summary.tsv"
WORKDIR="$HOME/estela/euk/fase1"
EUKCC_DB="$HOME/dbs/eukcc/eukcc2_db_ver_1.2"
BUSCO_DB="$HOME/dbs/busco"
BUSCO_LINAJES="fungi_odb12.2 eukaryota_odb12.2"
TAXDUMP="$HOME/dbs/taxdump"
MIN_COMPLETENESS=70
MAX_CONTAMINATION=5

THREADS="${SLURM_CPUS_PER_TASK:-32}"
RESULTS="$WORKDIR/results"
EUK_BINS="$RESULTS/bins_euk"
EUKCC_OUT="$RESULTS/01_eukcc"
BUSCO_OUT="$RESULTS/02_busco"
SELECTED_DIR="$RESULTS/phase1_euk_selected_genomes"
PY="$SCRIPTS_DIR/eucariotas/bin/summarize_euk.py"

echo "======================================================================"
echo " FASE 1 EUCARIOTA - Job ${SLURM_JOB_ID:-NA} | $(hostname) | $THREADS hilos | $(date)"
echo "======================================================================"

for p in "$BINS_DIR" "$EUKCC_DB" "$BUSCO_DB" "$TAXDUMP"; do
  [[ -e "$p" ]] || { echo "ERROR: no existe $p (corre 00_setup_euk_khipu.sh base)" >&2; exit 1; }
done
[[ -f "$TIARA_SUMMARY" ]] || { echo "ERROR: no existe $TIARA_SUMMARY (haz git pull)" >&2; exit 1; }
mkdir -p "$EUKCC_OUT" "$BUSCO_OUT" "$SELECTED_DIR"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

VERS="$RESULTS/versions.txt"
{
  echo "Fase 1 eucariota - $(date) - job ${SLURM_JOB_ID:-NA}"
  echo "Umbrales: completitud>$MIN_COMPLETENESS, contaminacion<$MAX_CONTAMINATION (EukCC)"
  echo "EukCC DB: $EUKCC_DB | BUSCO: $BUSCO_LINAJES"
  conda list -n eukcc --export | grep -E '^(eukcc|metaeuk|python)=' || true
  conda list -n busco --export | grep -E '^(busco|miniprot|python)=' || true
} > "$VERS"

echo "---> PASO 1: bins con mayoria eucariota segun Tiara"
conda activate euktax
python "$PY" seleccion --tiara-summary "$TIARA_SUMMARY" --bins-dir "$BINS_DIR" \
    --bin-ext "$BIN_EXT" --dest "$EUK_BINS"
conda deactivate
shopt -s nullglob
BINS=("$EUK_BINS"/*."$BIN_EXT")
N=${#BINS[@]}

echo "---> PASO 2: EukCC ($N bins)"
conda activate eukcc
i=0
for f in "${BINS[@]}"; do
  i=$((i+1)); name="$(basename "$f" ."$BIN_EXT")"; out="$EUKCC_OUT/$name"
  if [[ -s "$out/eukcc.tsv" ]]; then echo "     [$i/$N] $name ya procesado"; continue; fi
  echo "     [$i/$N] EukCC -> $name ($(date +%H:%M))"
  rc=0
  eukcc single --out "$out" --threads "$THREADS" --db "$EUKCC_DB" "$f" \
      > "$EUKCC_OUT/$name.log" 2>&1 || rc=$?
  # 201 = sin set de marcadores, 204 = sin proteinas: esperable en bins pequenos.
  if [[ "$rc" -ne 0 ]]; then
    echo "     AVISO: EukCC salio con codigo $rc en $name (ver $name.log)"
    printf "%s\t%s\n" "$name" "$rc" >> "$EUKCC_OUT/codigos_salida.tsv"
  fi
done
conda deactivate

echo "---> PASO 3: BUSCO ($BUSCO_LINAJES)"
conda activate busco
for lin in $BUSCO_LINAJES; do
  mkdir -p "$BUSCO_OUT/$lin"
  i=0
  for f in "${BINS[@]}"; do
    i=$((i+1)); name="$(basename "$f" ."$BIN_EXT")"
    if compgen -G "$BUSCO_OUT/$lin/$name/short_summary*.txt" >/dev/null; then continue; fi
    echo "     [$i/$N] BUSCO $lin -> $name ($(date +%H:%M))"
    busco -i "$f" -m genome -l "$lin" --offline --download_path "$BUSCO_DB" \
        -o "$name" --out_path "$BUSCO_OUT/$lin" -c "$THREADS" -f \
        > "$BUSCO_OUT/$lin/$name.log" 2>&1 \
      || echo "     AVISO: BUSCO fallo en $name ($lin), ver $name.log"
  done
done
conda deactivate

echo "---> PASO 4: seleccion y reporte"
conda activate euktax
rm -f "$SELECTED_DIR"/*."$BIN_EXT"
python "$PY" fase1 --tiara-summary "$TIARA_SUMMARY" --bins-dir "$EUK_BINS" --bin-ext "$BIN_EXT" \
    --eukcc-dir "$EUKCC_OUT" --busco-dir "$BUSCO_OUT" --busco-lineages "$BUSCO_LINAJES" \
    --taxdump "$TAXDUMP" --min-completeness "$MIN_COMPLETENESS" \
    --max-contamination "$MAX_CONTAMINATION" --selected-dir "$SELECTED_DIR" \
    --out-tsv "$RESULTS/phase1_euk_selection.tsv" --out-report "$RESULTS/phase1_euk_report.md" \
    --job-id "${SLURM_JOB_ID:-NA}"
conda deactivate

echo "======================================================================"
echo " FASE 1 EUCARIOTA COMPLETA - $(date)"
echo "   Tabla:   $RESULTS/phase1_euk_selection.tsv"
echo "   Reporte: $RESULTS/phase1_euk_report.md"
echo "   Seleccionados: $(ls "$SELECTED_DIR" | wc -l) en $SELECTED_DIR"
echo "======================================================================"
