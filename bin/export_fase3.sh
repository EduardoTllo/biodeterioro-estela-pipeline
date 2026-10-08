#!/bin/bash
# bin/export_fase3.sh
# Reune en ~/estela/fase3/export_fase3/ las tablas, reportes y arboles de la
# Fase 3 que se versionan en el repo (results/fase3/). No copia FASTA, GFF ni
# alineamientos: esos quedan en el cluster.
#
# Uso (en Khipu, en cualquier momento; copia lo que exista):
#   bash scripts/bin/export_fase3.sh
# Luego, desde el laptop (WSL):
#   rsync -avh eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/export_fase3/ /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/results/fase3/

set -eo pipefail

WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
R="$WORKDIR/results"
E="$WORKDIR/export_fase3"
mkdir -p "$E"

cp_si() {  # cp_si <origen> <destino>: copia solo si existe
  if [[ -s "$1" ]]; then mkdir -p "$(dirname "$2")"; cp -f "$1" "$2"; fi
}

# Censo (etapas 3.0-3.2)
for f in phase3_censo_report.md phase3_viabilidad.tsv phase3_censo_genomas.tsv \
         phase3_genomas_candidatos.tsv phase3_mapeo_ncbi.tsv phase3_ncbi_counts.tsv; do
  cp_si "$R/00_censo/$f" "$E/$f"
done
# Descarga y clasificacion (etapa 3.2b)
for f in descarga_resumen descarga_estado especies_viables especies_para_drep; do
  cp_si "$R/01_descarga/$f.tsv" "$E/phase3_$f.tsv"
done
for f in refs_clasificadas fuentes_unicas paises_sin_mapear; do
  cp_si "$R/02_drep/$f.tsv" "$E/phase3_$f.tsv"
done
# Seleccion
for f in phase3_ranking_especies.tsv phase3_referencias_finales.tsv \
         phase3_composicion_referencias.tsv phase3_seleccion_report.md; do
  cp_si "$R/03_seleccion/$f" "$E/$f"
done
cp_si "$R/03_seleccion/especies_seleccionadas.tsv" "$E/phase3_especies_seleccionadas.tsv"
cp_si "$R/03_seleccion/bakta_manifest.tsv" "$E/phase3_bakta_manifest.tsv"
cp_si "$R/04_bakta_resumen.tsv" "$E/phase3_bakta_resumen.tsv"
cp_si "$WORKDIR/scripts/envs/versiones_herramientas.txt" "$E/phase3_versiones_herramientas.txt"

# Por especie
for d in "$R"/02_drep/*/; do
  s="$(basename "$d")"
  for t in Cdb Wdb Sdb; do cp_si "$d/data_tables/$t.csv" "$E/$s/drep_$t.csv"; done
  cp_si "$d/versions.txt" "$E/$s/drep_versions.txt"
done
for d in "$R"/05_panaroo/*/; do
  s="$(basename "$d")"
  [[ "$s" == _piloto_* ]] && continue
  cp_si "$d/summary_statistics.txt" "$E/$s/panaroo_summary_statistics.txt"
  cp_si "$d/conversion_descartes.tsv" "$E/$s/panaroo_conversion_descartes.tsv"
  cp_si "$d/versions.txt" "$E/$s/panaroo_versions.txt"
  for f in particion_familias.tsv genes_bin.tsv recuperacion_core.tsv resumen_bins.tsv \
           resumen_particion.tsv particion_report.md pangenoma_refs.Rtab; do
    cp_si "$d/particion/$f" "$E/$s/$f"
  done
done
for d in "$R"/06_iqtree/*/; do
  s="$(basename "$d")"
  cp_si "$d/core.treefile" "$E/$s/core.treefile"
  cp_si "$d/modo_usado.txt" "$E/$s/iqtree_modo.txt"
  for f in "$d"/core.iqtree "$d"/core_snps.iqtree; do cp_si "$f" "$E/$s/iqtree_$(basename "$f" .iqtree).txt"; done
  for f in "$d"/versions_*.txt; do cp_si "$f" "$E/$s/iqtree_$(basename "$f")"; done
done
for d in "$R"/07_exclusivos/*/; do
  s="$(basename "$d")"
  for f in exclusivos_verificados.tsv embudo_exclusivos.tsv exclusivos_report.md \
           candidatos_local.tsv versions_local.txt versions_remoto.txt \
           f5_blastp_nr.tsv f5b_ani.tsv linajes.tsv; do
    cp_si "$d/$f" "$E/$s/$f"
  done
done
# Anotacion de Bakta y secuencia de cada bin (analisis de genes especificos:
# KO/EC y GC de las islas). El .fna no se versiona (.gitignore), solo se usa local.
tail -n +2 "$R/03_seleccion/bakta_manifest.tsv" | awk -F'\t' '$3 == "bin" {print $1 "\t" $4}' |
while IFS=$'\t' read -r B SLUG; do
  for ext in tsv fna; do
    cp_si "$R/04_bakta/$B/$B.$ext" "$E/$SLUG/bakta_bin/$B.$ext"
  done
done
# ANI de cada bin contra los genomas de su especie (06c_ani_bins.slurm)
for d in "$R"/06c_ani/*/; do
  [[ -d "$d" ]] || continue
  cp_si "$d/ani_bin_genomas.tsv" "$E/$(basename "$d")/ani_bin_genomas.tsv"
done
# Categoria COG de cada familia (07_particion.py cog)
for d in "$R"/05_panaroo/*/; do
  cp_si "$d/particion/familias_cog.tsv" "$E/$(basename "$d")/familias_cog.tsv"
done

# Exclusivos con F5 por el EBI (MOTOR=ebi), si se corrio: en <especie>/ebi/
for d in "$R"/07_exclusivos_ebi/*/; do
  [[ -d "$d" ]] || continue
  s="$(basename "$d")"
  for f in exclusivos_verificados.tsv embudo_exclusivos.tsv exclusivos_report.md            versions_remoto.txt; do
    cp_si "$d/$f" "$E/$s/ebi/$f"
  done
done

echo "Exportado a $E:"
find "$E" -type f | sed "s|$E/||" | sort
