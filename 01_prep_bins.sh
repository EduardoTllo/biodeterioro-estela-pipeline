#!/bin/bash
# =============================================================================
# 01_prep_bins.sh
# -----------------------------------------------------------------------------
# Prepara la carpeta de bins limpia y trazable para la Fase 1.
#
# FUENTE RECOMENDADA:
#   C:\Tesis-EstelaRaimondi\Datos\Shotgun Analysis\bins completos\
#   -> 97 genomas con nombres "bin-<N>-<muestra>.fasta" (p.ej. bin-1-49.fasta).
#      El nombre YA codifica el bin y la muestra de origen (libreria 49-72), lo
#      que preserva la trazabilidad necesaria para la Fase 2 (prevalencia).
#
# QUE HACE:
#   - Copia los .fasta de genoma (patron bin-*-*.fasta) a la carpeta de salida,
#     conservando el nombre original (limpio y unico).
#   - Ignora los ensamblajes por muestra "<muestra>_mag.fasta" (no son bins).
#   - Escribe bin_rename_map.tsv con: archivo, muestra, bin, n_contigs, long_bp.
#
# USO:
#   bash 01_prep_bins.sh <carpeta_origen> <carpeta_salida>
# Ejemplo:
#   bash 01_prep_bins.sh "/c/Tesis-EstelaRaimondi/Datos/Shotgun Analysis/bins completos" \
#                        "/c/Tesis-EstelaRaimondi/Datos/bins_clean"
#
# Tambien acepta como origen el arbol por muestra ("Shotgun Analysis") y busca
# recursivamente los bin-*-*.fasta.
# =============================================================================
set -euo pipefail

SRC="${1:?Falta la carpeta de origen}"
DEST="${2:?Falta la carpeta de salida}"

if [[ ! -d "$SRC" ]]; then
  echo "ERROR: no existe la carpeta de origen: $SRC" >&2; exit 1
fi

mkdir -p "$DEST"
MAP="$DEST/bin_rename_map.tsv"
printf "archivo\tmuestra\tbin\tn_contigs\tlongitud_bp\n" > "$MAP"

i=0
# Buscamos genomas de bin (bin-<N>-<muestra>.fasta), excluyendo <muestra>_mag.fasta.
# Deduplicamos por nombre de archivo (por si el origen es el arbol por muestra
# y ademas contiene la carpeta "bins completos").
declare -A seen
while IFS= read -r -d '' f; do
  bn="$(basename "$f")"
  [[ -n "${seen[$bn]:-}" ]] && continue
  seen[$bn]=1
  i=$((i+1))
  cp -f "$f" "$DEST/$bn"
  # Parseo del nombre: bin-<N>-<muestra>.fasta
  core="${bn%.fasta}"                 # bin-1-49
  sample="${core##*-}"                # 49
  rest="${core%-*}"                   # bin-1
  binn="${rest##*-}"                  # 1
  nctg=$(grep -c '^>' "$DEST/$bn" || echo 0)
  len=$(grep -v '^>' "$DEST/$bn" | tr -d '\n' | wc -c)
  printf "%s\t%s\t%s\t%s\t%s\n" "$bn" "$sample" "$binn" "$nctg" "$len" >> "$MAP"
done < <(find "$SRC" -type f -name 'bin-*-*.fasta' ! -name '*_mag.fasta' -print0 | sort -z)

echo "Genomas preparados: $i"
echo "Carpeta de salida:  $DEST"
echo "Tabla de mapeo:     $MAP  (con muestra de origen para trazabilidad)"
if [[ "$i" -ne 97 ]]; then
  echo "ADVERTENCIA: se esperaban 97 bins y se encontraron $i. Revisa el origen." >&2
fi
