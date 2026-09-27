#!/bin/bash
# 08_exclusivos_remoto.sh
# FASE 3 (OE3), Etapa 3.5: verificacion remota de genes exclusivos (F5, F6) y
# clasificacion final.
#
#   F5  blastp -remote contra nr de las proteinas que sobrevivieron a F1 y F4.
#       Se envia por lotes de 50 secuencias, uno tras otro (BLAST+ respeta los
#       limites de uso de NCBI). Cada lote terminado queda marcado y no se
#       repite si el script se interrumpe.
#   F6  blastn -remote -task megablast contra core_nt de los contigs huerfanos.
#       Si core_nt no respondiera, exportar NT_DB=nt y registrarlo.
#   Taxonomia: taxonkit (NCBI taxdump local) resuelve filo;familia;genero;especie
#       de cada hit y de la especie del bin.
#   Integrar: 08_exclusivos.py integrar -> exclusivos_verificados.tsv.
#
# Donde se ejecuta: NODO DE LOGIN de Khipu, dentro de tmux (necesita internet).
# Si el login no admite procesos largos, corre igual desde WSL en el laptop con
# blast+, taxonkit y python instalados: el BLAST remoto solo necesita internet.
#   Uso:   bash 08_exclusivos_remoto.sh              (todas las especies)
#          bash 08_exclusivos_remoto.sh <slug>       (una especie)

# No usamos 'set -u' por compatibilidad con Lmod en Khipu.
set -eo pipefail
export LD_PRELOAD="${LD_PRELOAD:-}"

# Configuracion
WORKDIR="${WORKDIR:-$HOME/estela/fase3}"
SCRIPTS_DIR="$WORKDIR/scripts"
SEL_DIR="$WORKDIR/results/03_seleccion"
EXCL_OUT="$WORKDIR/results/07_exclusivos"
TAXDUMP_DIR="${TAXDUMP_DIR:-$HOME/dbs/taxdump}"
NR_DB="${NR_DB:-nr}"
NT_DB="${NT_DB:-core_nt}"
LOTE=50
PAUSA_S=15
OUTFMT="6 qseqid sseqid pident length qcovs evalue bitscore staxids"

# No editar debajo salvo que sepas lo que haces
if command -v module >/dev/null 2>&1; then
  module purge
  module load miniconda/3.0
fi
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate blast
fi
for b in blastp blastn taxonkit python; do
  command -v "$b" >/dev/null || { echo "ERROR: falta '$b' en el PATH" >&2; exit 1; }
done

SOLO="$1"
tail -n +2 "$SEL_DIR/especies_seleccionadas.tsv" | while IFS=$'\t' read -r RANGO SLUG LID ESP NCBI_ESP TAXID RESTO; do
  [[ -n "$SOLO" && "$SLUG" != "$SOLO" ]] && continue
  OUT="$EXCL_OUT/$SLUG"
  if [[ ! -s "$OUT/candidatos_local.tsv" ]]; then
    echo "ERROR: falta $OUT/candidatos_local.tsv (corre 08_exclusivos_local.slurm)." >&2; exit 1
  fi
  echo
  echo "===> $SLUG (taxid NCBI $TAXID)"
  {
    echo "Fase 3 - Exclusivos (remoto) - $SLUG"
    echo "Fecha: $(date)   Host: $(hostname)"
    echo "F5: blastp -remote -db $NR_DB -evalue 1e-5 -max_target_seqs 10 (lotes de $LOTE)"
    echo "F6: blastn -remote -db $NT_DB -task megablast -evalue 1e-10 -max_target_seqs 10"
    echo "blast: $(blastp -version | head -n1) | taxonkit: $(taxonkit version 2>&1 | head -n1)"
    echo "taxdump: $TAXDUMP_DIR ($(cat "$TAXDUMP_DIR/fecha_descarga.txt" 2>/dev/null || echo NA))"
  } > "$OUT/versions_remoto.txt"

  # F5: lotes de LOTE proteinas
  mkdir -p "$OUT/f5_lotes"
  if [[ -s "$OUT/f5_query.faa" ]]; then
    awk -v n="$LOTE" -v d="$OUT/f5_lotes" '/^>/{if(c%n==0){f=sprintf("%s/lote_%04d.faa",d,int(c/n)+1)} c++} {print > f}' \
        "$OUT/f5_query.faa"
    for Q in "$OUT"/f5_lotes/lote_*.faa; do
      R="${Q%.faa}.tsv"
      if [[ -f "${R}.ok" ]]; then continue; fi
      echo "     F5 $(basename "$Q"): $(grep -c '^>' "$Q") proteinas ($(date +%H:%M))"
      blastp -remote -db "$NR_DB" -query "$Q" -evalue 1e-5 -max_target_seqs 10 \
             -outfmt "$OUTFMT" > "$R"
      touch "${R}.ok"
      sleep "$PAUSA_S"
    done
    cat "$OUT"/f5_lotes/lote_*.tsv > "$OUT/f5_blastp_nr.tsv"
  else
    : > "$OUT/f5_blastp_nr.tsv"
  fi

  # F6: contigs huerfanos
  if [[ -s "$OUT/f6_query.fna" && ! -f "$OUT/f6_blastn_nt.tsv.ok" ]]; then
    echo "     F6: $(grep -c '^>' "$OUT/f6_query.fna") contigs huerfanos"
    blastn -remote -db "$NT_DB" -task megablast -query "$OUT/f6_query.fna" \
           -evalue 1e-10 -max_target_seqs 10 -outfmt "$OUTFMT" > "$OUT/f6_blastn_nt.tsv"
    touch "$OUT/f6_blastn_nt.tsv.ok"
  elif [[ ! -s "$OUT/f6_query.fna" ]]; then
    : > "$OUT/f6_blastn_nt.tsv"
  fi

  # Taxonomia de los hits y de la especie del bin
  { echo "$TAXID"; cut -f8 "$OUT/f5_blastp_nr.tsv" "$OUT/f6_blastn_nt.tsv" | tr ';' '\n'; } \
      | grep -E '^[0-9]+$' | sort -u > "$OUT/taxids.txt"
  taxonkit reformat --data-dir "$TAXDUMP_DIR" -I 1 -f "{p};{f};{g};{s}" "$OUT/taxids.txt" \
      > "$OUT/linajes.tsv"

  python "$SCRIPTS_DIR/08_exclusivos.py" integrar \
      --outdir "$OUT" --f5 "$OUT/f5_blastp_nr.tsv" --f6 "$OUT/f6_blastn_nt.tsv" \
      --linajes "$OUT/linajes.tsv" --taxid-especie "$TAXID"
done

echo
echo "======================================================================"
echo " EXCLUSIVOS (REMOTO) COMPLETO - $(date)"
echo "   Resultado por especie: $EXCL_OUT/<especie>/exclusivos_verificados.tsv"
echo "   Embudo:                $EXCL_OUT/<especie>/embudo_exclusivos.tsv"
echo " Checkpoint C7: revisar exclusivos_report.md de cada especie."
echo "======================================================================"
