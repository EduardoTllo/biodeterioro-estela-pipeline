#!/bin/bash
# eucariotas/run_fase2_euk.slurm
# FASE 2 EUCARIOTA - Taxonomia y linajes de los bins eucariotas (anexo)
# Tesis: Estela de Raimondi | Cluster: Khipu (UTEC) | SLURM 25.11
#
# Flujo (solo sobre los genomas seleccionados en la Fase 1 eucariota):
#   Paso 1 (skani)  -> ANI/AF contra los genomas de referencia de hongos de NCBI
#                      (reemplaza al ANI de GTDB-Tk). Especie: ANI>=95, AF>=65.
#   Paso 2 (ITS)    -> barrnap ubica el ADNr, ITSx extrae el ITS y BLAST lo
#                      compara con UNITE y RefSeq ITS.
#   Paso 3 (arbol)  -> genes BUSCO de copia unica -> MAFFT -> trimAl -> matriz
#                      concatenada (equivalente al user_msa bac120). El NJ se
#                      calcula en R, igual que la Fig. 7.4.
#   Paso 4 (dRep)   -> linaje = cluster a 95 % de ANI (umbral operativo),
#                      con la calidad de EukCC via genomeInfo.
#   Paso 5 (resumen)-> taxonomia en cascada, linajes y reporte (sin prevalencia).
#
# Requiere: 00_setup_euk_khipu.sh base y refs, y run_fase1_euk.slurm terminado.
# Envio (desde ~/estela/euk):
#   sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" \
#          ~/estela/fase3/scripts/eucariotas/run_fase2_euk.slurm

#SBATCH --job-name=fase2_euk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=24:00:00
#SBATCH --output=fase2_euk_%j.out
#SBATCH --error=fase2_euk_%j.err

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
SCRIPTS_DIR="$HOME/estela/fase3/scripts"
FASE1="$HOME/estela/euk/fase1/results"
SELECTED="$FASE1/phase1_euk_selected_genomes"
PHASE1_TSV="$FASE1/phase1_euk_selection.tsv"
BIN_EXT="fasta"
WORKDIR="$HOME/estela/euk/fase2"
REFS_DIR="$HOME/dbs/euk_refs"
SKANI_DB="$REFS_DIR/skani_db"
ITS_DIR="$HOME/dbs/its"
BUSCO_ARBOL="fungi_odb12.2"   # cambiar a eukaryota_odb12.2 si algun seleccionado no es hongo
MIN_FRAC_BUSCO=0.5            # gen en >= 50 % de los genomas para entrar al arbol
DREP_PA=0.90
DREP_SA=0.95
ANI_ESPECIE=95
AF_ESPECIE=65

THREADS="${SLURM_CPUS_PER_TASK:-32}"
RESULTS="$WORKDIR/results"
SKANI_OUT="$RESULTS/01_skani"
ITS_OUT="$RESULTS/02_its"
ARBOL_OUT="$RESULTS/03_arbol_busco"
DREP_OUT="$RESULTS/04_drep"
PY="$SCRIPTS_DIR/eucariotas/bin/summarize_euk.py"

echo "======================================================================"
echo " FASE 2 EUCARIOTA - Job ${SLURM_JOB_ID:-NA} | $(hostname) | $THREADS hilos | $(date)"
echo "======================================================================"

shopt -s nullglob
GENOMES=("$SELECTED"/*."$BIN_EXT")
N=${#GENOMES[@]}
[[ "$N" -gt 0 ]] || { echo "ERROR: no hay genomas en $SELECTED" >&2; exit 1; }
[[ -s "$REFS_DIR/refs_fungi.tsv" ]] || { echo "ERROR: falta $REFS_DIR/refs_fungi.tsv (00_setup_euk_khipu.sh refs)" >&2; exit 1; }
echo "Genomas de entrada: $N"
mkdir -p "$SKANI_OUT" "$ITS_OUT/gff" "$ARBOL_OUT"/{genes,aln,trim} "$DREP_OUT"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"

UNITE_DB=""
if [[ -f "$ITS_DIR/unite.nsq" || -f "$ITS_DIR/unite.00.nsq" ]]; then UNITE_DB="$ITS_DIR/unite"; fi
{
  echo "Fase 2 eucariota - $(date) - job ${SLURM_JOB_ID:-NA}"
  echo "Referencias: $(($(wc -l < "$REFS_DIR/refs_fungi.tsv") - 1)) genomas (NCBI, descarga $(cat "$REFS_DIR/fecha_descarga.txt" 2>/dev/null))"
  echo "UNITE: ${UNITE_DB:-no disponible} | RefSeq ITS: $ITS_DIR/ITS_RefSeq_Fungi"
  echo "Umbrales: especie ANI>=$ANI_ESPECIE AF>=$AF_ESPECIE | dRep -pa $DREP_PA -sa $DREP_SA | BUSCO arbol $BUSCO_ARBOL (>= $MIN_FRAC_BUSCO)"
  conda list -n euktax --export | grep -E '^(skani|barrnap|itsx|mafft|trimal)=' || true
  conda list -n blast --export | grep -E '^blast=' || true
  conda list -n drep --export | grep -E '^(drep|fastani|mash)=' || true
} > "$RESULTS/versions.txt"
printf "%s\n" "${GENOMES[@]}" > "$RESULTS/genomas.txt"

echo "---> PASO 1: skani contra referencias"
conda activate euktax
if [[ ! -s "$SKANI_DB/markers.bin" ]]; then
  echo "     Armando el sketch de las referencias (una sola vez)..."
  rm -rf "$SKANI_DB"
  skani sketch -l "$REFS_DIR/refs_fasta.txt" -o "$SKANI_DB" -t "$THREADS"
fi
skani search -d "$SKANI_DB" --ql "$RESULTS/genomas.txt" -o "$SKANI_OUT/skani_bins_vs_refs.tsv" \
    -n 10 -t "$THREADS"

echo "---> PASO 2: ITS (barrnap -> ITSx -> BLAST)"
for f in "${GENOMES[@]}"; do
  name="$(basename "$f" ."$BIN_EXT")"
  [[ -s "$ITS_OUT/gff/$name.gff" ]] || barrnap --kingdom euk --threads "$THREADS" --quiet "$f" \
      > "$ITS_OUT/gff/$name.gff" || true
done
python "$PY" its-ventanas --gff-dir "$ITS_OUT/gff" --bins-dir "$SELECTED" --bin-ext "$BIN_EXT" \
    --out "$ITS_OUT/ventanas_adnr.fasta"
ITS_FULL="$ITS_OUT/itsx.full.fasta"
if [[ -s "$ITS_OUT/ventanas_adnr.fasta" ]]; then
  ITSx -i "$ITS_OUT/ventanas_adnr.fasta" -o "$ITS_OUT/itsx" -t all --cpu "$THREADS" \
      --save_regions all --graphical F --preserve T > "$ITS_OUT/itsx.log" 2>&1 || true
fi
conda deactivate
FMT="6 qseqid pident length qlen qcovs evalue bitscore stitle"
if [[ -s "$ITS_FULL" ]]; then
  conda activate blast
  [[ -n "$UNITE_DB" ]] && blastn -query "$ITS_FULL" -db "$UNITE_DB" -outfmt "$FMT" \
      -max_target_seqs 20 -evalue 1e-20 -num_threads "$THREADS" > "$ITS_OUT/its_vs_unite.tsv"
  blastn -query "$ITS_FULL" -db "$ITS_DIR/ITS_RefSeq_Fungi" -outfmt "$FMT" \
      -max_target_seqs 20 -evalue 1e-20 -num_threads "$THREADS" > "$ITS_OUT/its_vs_refseq.tsv"
  conda deactivate
  echo "     ITS completos: $(grep -c '>' "$ITS_FULL")"
else
  echo "     AVISO: no se recupero ningun ITS completo (ADNr ausente o fragmentado en los bins)."
fi

echo "---> PASO 3: matriz BUSCO para el arbol"
conda activate euktax
if [[ "$N" -ge 3 ]]; then
  rm -f "$ARBOL_OUT"/genes/*.faa "$ARBOL_OUT"/aln/*.aln "$ARBOL_OUT"/trim/*.aln
  python "$PY" busco-genes --busco-dir "$FASE1/02_busco/$BUSCO_ARBOL" --genomes-dir "$SELECTED" \
      --bin-ext "$BIN_EXT" --min-frac "$MIN_FRAC_BUSCO" --out-dir "$ARBOL_OUT/genes"
  for g in "$ARBOL_OUT"/genes/*.faa; do
    id="$(basename "$g" .faa)"
    mafft --auto --thread "$THREADS" --quiet "$g" > "$ARBOL_OUT/aln/$id.aln"
    trimal -in "$ARBOL_OUT/aln/$id.aln" -out "$ARBOL_OUT/trim/$id.aln" -automated1 >/dev/null 2>&1 || true
  done
  python "$PY" concat --aln-dir "$ARBOL_OUT/trim" --genomes-dir "$SELECTED" --bin-ext "$BIN_EXT" \
      --out-faa "$RESULTS/phase2_euk_busco_concat.faa" \
      --out-partition "$RESULTS/phase2_euk_busco_particiones.txt"
else
  echo "     Menos de 3 genomas: no se arma arbol."
fi
conda deactivate

echo "---> PASO 4: dRep (linaje = ${DREP_SA} ANI)"
conda activate drep
CDB=""; WDB=""
if [[ "$N" -ge 2 ]]; then
  python "$PY" genomeinfo --phase1-tsv "$PHASE1_TSV" --genomes-dir "$SELECTED" \
      --bin-ext "$BIN_EXT" --out "$RESULTS/genomeInfo_eukcc.csv"
  # -comp 0 -con 100: la calidad ya se filtro en la Fase 1 con EukCC.
  dRep dereplicate "$DREP_OUT" -g "${GENOMES[@]}" -p "$THREADS" \
      -pa "$DREP_PA" -sa "$DREP_SA" --S_algorithm fastANI \
      -comp 0 -con 100 --genomeInfo "$RESULTS/genomeInfo_eukcc.csv"
  CDB="$DREP_OUT/data_tables/Cdb.csv"; WDB="$DREP_OUT/data_tables/Wdb.csv"
else
  echo "     Un solo genoma: no se corre dRep."
fi
conda deactivate

echo "---> PASO 5: taxonomia en cascada, linajes y reporte"
conda activate euktax
python "$PY" fase2 --genomes-dir "$SELECTED" --bin-ext "$BIN_EXT" --phase1-tsv "$PHASE1_TSV" \
    --refs-tsv "$REFS_DIR/refs_fungi.tsv" --skani-tsv "$SKANI_OUT/skani_bins_vs_refs.tsv" \
    --its-fasta "$ITS_FULL" --its-unite "$ITS_OUT/its_vs_unite.tsv" \
    --its-refseq "$ITS_OUT/its_vs_refseq.tsv" --cdb "$CDB" --wdb "$WDB" --drep-sa "$DREP_SA" \
    --min-ani "$ANI_ESPECIE" --min-af "$AF_ESPECIE" \
    --out-tsv "$RESULTS/phase2_euk_selection.tsv" --out-linajes "$RESULTS/phase2_euk_linajes.tsv" \
    --out-report "$RESULTS/phase2_euk_report.md" --job-id "${SLURM_JOB_ID:-NA}"
conda deactivate

echo "======================================================================"
echo " FASE 2 EUCARIOTA COMPLETA - $(date)"
echo "   Por genoma: $RESULTS/phase2_euk_selection.tsv"
echo "   Linajes:    $RESULTS/phase2_euk_linajes.tsv"
echo "   Reporte:    $RESULTS/phase2_euk_report.md"
echo "======================================================================"
