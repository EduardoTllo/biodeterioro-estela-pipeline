#!/bin/bash
# 04_fetch_refs.sh
# FASE 3 (OE3), Etapa 3.2b: descarga de los genomas de referencia desde NCBI.
#
# Que hace, para cada especie viable del censo:
#   1. Arma la lista de accesiones candidatas (post-QC) con 04_clasificar_refs.py listas.
#   2. Descarga los FASTA con NCBI Datasets en modo deshidratado y los rehidrata
#      (flujo que NCBI recomienda para descargas grandes). Solo --include genome:
#      nunca se bajan los GFF de NCBI (todo se reanota con Bakta, decision D14).
#   3. Renombra cada FASTA a <ID>.fna, con ID = accesion con '.' -> '_'.
#   4. Verifica la descarga (04_clasificar_refs.py verificar): cada FASTA debe
#      existir y su longitud total debe coincidir con genome_size de GTDB (+/- 1 %).
#      No hay sustituciones: una accesion que no baja o que no coincide queda fuera
#      y se registra (decision D29). Tambien se aplican las exclusiones
#      manuales de metadata/exclusiones.tsv.
#
# Se descargan TODAS las candidatas post-QC de cada especie viable, no solo las
# <= 50 que entraran al pangenoma: el conjunto completo se usa luego en el filtro
# F4 de genes exclusivos (tblastn contra toda la especie).
#
# Donde se ejecuta: NODO DE LOGIN de Khipu (unico con internet), dentro de tmux.
#   Uso:   bash 04_fetch_refs.sh
#
# Es idempotente: una especie cuya carpeta ya tiene todos los FASTA se omite, y
# una descarga interrumpida se retoma al volver a correr el script.

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SCRIPTS_DIR="$WORKDIR/scripts"
CENSO_DIR="$WORKDIR/results/00_censo"
DESCARGA_DIR="$WORKDIR/results/01_descarga"
REFS_DIR="$WORKDIR/data/refs"
UMBRAL=15
TOLERANCIA=0.01
# Genomas excluidos a mano, con su motivo (accession, linaje_id, motivo).
EXCLUSIONES="$SCRIPTS_DIR/metadata/exclusiones.tsv"

# No editar debajo salvo que sepas lo que haces
TMP_DIR="$WORKDIR/tmp/descarga"
mkdir -p "$DESCARGA_DIR" "$REFS_DIR" "$TMP_DIR"
LOG="$DESCARGA_DIR/descarga.log"

module purge
module load miniconda/3.0
eval "$(conda shell.bash hook)"
conda activate ncbi-datasets

echo "======================================================================"
echo " FASE 3 - Descarga de referencias (NCBI Datasets)"
echo " Inicio: $(date)   datasets: $(datasets --version 2>&1 | head -n1)"
echo "======================================================================"
echo "Inicio: $(date)  datasets: $(datasets --version 2>&1 | head -n1)" >> "$LOG"

if [[ ! -f "$CENSO_DIR/phase3_viabilidad.tsv" ]]; then
  echo "ERROR: falta $CENSO_DIR/phase3_viabilidad.tsv (corre 03_censo_genomas.py censo/tabla)." >&2
  exit 1
fi

# 1) Listas de accesiones por especie viable
python "$SCRIPTS_DIR/04_clasificar_refs.py" listas \
    --censo-dir "$CENSO_DIR" --outdir "$DESCARGA_DIR"

# 2-3) Descarga especie por especie
tail -n +2 "$DESCARGA_DIR/especies_viables.tsv" | cut -f1 | while read -r SLUG; do
  ACC_FILE="$DESCARGA_DIR/acc_${SLUG}.txt"
  DEST="$REFS_DIR/$SLUG"
  mkdir -p "$DEST"
  N_ACC=$(grep -c . "$ACC_FILE" || true)
  N_TIENE=$(find "$DEST" -maxdepth 1 -name '*.fna' | wc -l)
  echo
  echo "---> $SLUG: $N_ACC accesiones ($N_TIENE ya descargadas)"
  if [[ "$N_TIENE" -ge "$N_ACC" ]]; then
    echo "     Completa. Se omite."
    continue
  fi

  PKG="$TMP_DIR/$SLUG"
  ZIP="$TMP_DIR/${SLUG}.zip"
  if [[ ! -d "$PKG/ncbi_dataset" ]]; then
    rm -f "$ZIP"
    datasets download genome accession --inputfile "$ACC_FILE" \
        --include genome --dehydrated --filename "$ZIP"
    unzip -q -o "$ZIP" -d "$PKG"
  fi
  # La rehidratacion se reintenta: baja solo lo que falta.
  for intento in 1 2 3; do
    if datasets rehydrate --directory "$PKG/"; then break; fi
    echo "     rehydrate fallo (intento $intento); se reintenta en 30 s" | tee -a "$LOG"
    sleep 30
  done

  # Cada genoma queda en ncbi_dataset/data/<ACCESION>/<...>_genomic.fna
  N_OK=0
  for d in "$PKG"/ncbi_dataset/data/GC[AF]_*; do
    [[ -d "$d" ]] || continue
    ACC="$(basename "$d")"
    ID="$(echo "$ACC" | sed 's/[^A-Za-z0-9]/_/g')"
    FNA="$(find "$d" -maxdepth 1 -name '*_genomic.fna' | head -n1)"
    if [[ -n "$FNA" && -s "$FNA" ]]; then
      cp -f "$FNA" "$DEST/$ID.fna"
      N_OK=$((N_OK + 1))
    else
      echo "$SLUG	$ACC	sin_fasta_tras_rehidratar" >> "$LOG"
    fi
  done
  echo "     FASTA copiados: $N_OK de $N_ACC" | tee -a "$LOG"
done

# 4) Verificacion contra GTDB y viabilidad tras la descarga
echo
echo "---> Verificacion de la descarga (existencia y longitud vs GTDB)"
python "$SCRIPTS_DIR/04_clasificar_refs.py" verificar \
    --censo-dir "$CENSO_DIR" --descarga-dir "$DESCARGA_DIR" --refs-dir "$REFS_DIR" \
    --umbral "$UMBRAL" --tolerancia "$TOLERANCIA" --exclusiones "$EXCLUSIONES"

echo
echo "======================================================================"
echo " DESCARGA COMPLETA - $(date)"
echo "   Estado por genoma:   $DESCARGA_DIR/descarga_estado.tsv"
echo "   Resumen por especie: $DESCARGA_DIR/descarga_resumen.tsv"
echo "   Pasan al dRep:       $DESCARGA_DIR/especies_para_drep.tsv"
echo " Checkpoint C2: descargados + fallidos = solicitados; revisar"
echo " 'tamano_discrepante' y 'version_distinta' en descarga_estado.tsv."
echo " Los temporales de $TMP_DIR se pueden borrar tras revisar."
echo "======================================================================"
