#!/bin/bash
# =============================================================================
# 01_prep_bins.sh
# -----------------------------------------------------------------------------
# Normaliza los nombres de los bins antes de la Fase 1.
#
# CONTEXTO:
#   Los 97 bins vienen de una carpeta (MAGs_all) donde, al consolidarlos, el
#   sistema operativo renombro los duplicados con sufijos tipo:
#       "bin-1.fasta 10", "bin-1.fasta 22"
#   Esto rompe la extension y genera colisiones. Este script copia cada bin a
#   un nombre limpio y unico ( bin_0001.fasta ... bin_0097.fasta ) y guarda una
#   tabla de mapeo para trazabilidad (bin_rename_map.tsv).
#
#   IMPORTANTE: este paso NO recupera la muestra de origen (libreria 49-72), que
#   se perdio al aplanar la carpeta. Eso no afecta la Fase 1 (el QC es por bin),
#   pero para la Fase 2 (prevalencia por muestra) hay que recuperar el mapeo
#   bin->muestra del binning original (UPCH).
#
# USO (local o en el cluster):
#   bash 01_prep_bins.sh <carpeta_bins_crudos> <carpeta_salida>
# Ejemplo:
#   bash 01_prep_bins.sh "/c/Tesis-EstelaRaimondi/Datos/MAGs_all" ./bins
# =============================================================================
set -euo pipefail

SRC="${1:?Falta la carpeta de bins crudos}"
DEST="${2:?Falta la carpeta de salida}"
OUT_EXT="fasta"

if [[ ! -d "$SRC" ]]; then
  echo "ERROR: no existe la carpeta de origen: $SRC" >&2; exit 1
fi

mkdir -p "$DEST"
MAP="$DEST/bin_rename_map.tsv"
printf "nuevo_nombre\tarchivo_original\tn_contigs\tlongitud_bp\n" > "$MAP"

# Recorremos TODOS los archivos regulares del origen (la carpeta solo trae
# bins). Orden determinista por nombre para que el resultado sea reproducible.
i=0
while IFS= read -r -d '' f; do
  i=$((i+1))
  new=$(printf "bin_%04d.%s" "$i" "$OUT_EXT")
  cp -f "$f" "$DEST/$new"
  # Metricas basicas para la tabla de mapeo.
  nctg=$(grep -c '^>' "$DEST/$new" || echo 0)
  len=$(grep -v '^>' "$DEST/$new" | tr -d '\n' | wc -c)
  printf "%s\t%s\t%s\t%s\n" "$new" "$(basename "$f")" "$nctg" "$len" >> "$MAP"
done < <(find "$SRC" -maxdepth 1 -type f -print0 | sort -z)

echo "Bins normalizados: $i"
echo "Carpeta de salida: $DEST"
echo "Tabla de mapeo:    $MAP"
if [[ "$i" -ne 97 ]]; then
  echo "ADVERTENCIA: se esperaban 97 bins y se encontraron $i. Revisa el origen." >&2
fi
