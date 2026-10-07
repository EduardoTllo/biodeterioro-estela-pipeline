#!/bin/bash
# 08_exclusivos_remoto.sh
# FASE 3 (OE3), Etapa 3.5: verificacion remota de genes exclusivos (F5, F6) y
# clasificacion final.
#
#   F5  blastp contra nr de las proteinas que sobrevivieron a F1 y F4.
#       Se envia por lotes de 50 secuencias, uno tras otro. Cada lote terminado
#       queda marcado y no se repite si el script se interrumpe.
#   F6  blastn megablast contra core_nt de los contigs huerfanos.
#       Si core_nt no respondiera, exportar NT_DB=nt y registrarlo.
#
#   Las busquedas NO usan 'blastp -remote': ese cliente se cuelga o falla con
#   "Connection stream is in bad state (Blast4-request)". Se envian por la
#   URL API de NCBI con curl (CMD=Put), se consulta el estado cada 60 s
#   (minimo que pide NCBI) y el resultado se recoge con 'blast_formatter -rid',
#   que da el mismo formato tabular (qcovs, staxids). El RID queda guardado en
#   <lote>.tsv.rid: si el script se corta, al relanzarlo se reutiliza (NCBI lo
#   conserva ~36 h) en vez de enviar la busqueda de nuevo.
#   Taxonomia: taxonkit (NCBI taxdump local) resuelve filo;familia;genero;especie
#       de cada hit y de la especie del bin.
#   Integrar: 08_exclusivos.py integrar -> exclusivos_verificados.tsv.
#
# Donde se ejecuta: NODO DE LOGIN de Khipu, dentro de tmux (necesita internet).
# Si el login no admite procesos largos, corre igual desde WSL en el laptop con
# blast+, taxonkit y python instalados: el BLAST remoto solo necesita internet.
#   Uso:   bash 08_exclusivos_remoto.sh              (todas las especies)
#          bash 08_exclusivos_remoto.sh <slug>       (una especie)
#
# Alternativa si la cola de NCBI no avanza (MOTOR=ebi): F5 con el blastp del
# EBI contra UniProtKB (REST, un trabajo por proteina, 5 a la vez; 08_exclusivos.py
# ebi). Mismo tabular y misma integracion; los taxid de UniProt (OX) son los de
# NCBI. Escribe en results/07_exclusivos_ebi/<especie>/ para no mezclarse con la
# corrida de NCBI. El EBI exige un correo: EBI_EMAIL (por defecto $CORREO).
# F6 no tiene equivalente aqui: si hay contigs huerfanos, usar MOTOR=ncbi.
#   Uso:   MOTOR=ebi bash 08_exclusivos_remoto.sh

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
BLAST_URL="https://blast.ncbi.nlm.nih.gov/Blast.cgi"
ESPERA_MAX_H="${ESPERA_MAX_H:-24}"   # horas maximas de espera por busqueda
MOTOR="${MOTOR:-ncbi}"               # ncbi | ebi
EBI_EMAIL="${EBI_EMAIL:-$CORREO}"
EBI_DB="${EBI_DB:-uniprotkb}"
if [[ "$MOTOR" == "ebi" ]]; then
  [[ -n "$EBI_EMAIL" ]] || { echo "ERROR: MOTOR=ebi necesita EBI_EMAIL (o CORREO)" >&2; exit 1; }
  EXCL_LOCAL="$EXCL_OUT"
  EXCL_OUT="$WORKDIR/results/07_exclusivos_ebi"
elif [[ "$MOTOR" != "ncbi" ]]; then
  echo "ERROR: MOTOR debe ser ncbi o ebi" >&2; exit 1
fi

# No editar debajo salvo que sepas lo que haces
if command -v module >/dev/null 2>&1; then
  module purge
  module load miniconda/3.0
fi
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate blast
fi
for b in curl blast_formatter blastp taxonkit python; do
  command -v "$b" >/dev/null || { echo "ERROR: falta '$b' en el PATH" >&2; exit 1; }
done

# blast_ncbi <programa> <bd> <evalue> <consulta> <salida> [parametros URL extra...]
# Envia (o reutiliza el RID guardado), espera y recoge el tabular en <salida>.
blast_ncbi() {
  local prog="$1" db="$2" ev="$3" q="$4" out="$5"; shift 5
  local rid="" st="" i
  [[ -s "${out}.rid" ]] && rid="$(cat "${out}.rid")"
  if [[ -z "$rid" ]]; then
    for i in 1 2 3 4 5; do
      rid="$(curl -s --max-time 300 -X POST "$BLAST_URL" \
               --data-urlencode "CMD=Put" --data-urlencode "PROGRAM=$prog" \
               --data-urlencode "DATABASE=$db" --data-urlencode "EXPECT=$ev" \
               --data-urlencode "HITLIST_SIZE=10" --data-urlencode "TOOL=estela_fase3" \
               "$@" --data-urlencode "QUERY@$q" \
             | sed -n 's/^ *RID = *//p' | head -n1 | tr -d '[:space:]' || true)"
      [[ -n "$rid" ]] && break
      echo "       envio sin RID (intento $i); reintento en 2 min" >&2
      sleep 120
    done
    [[ -n "$rid" ]] || { echo "ERROR: NCBI no acepto la busqueda $(basename "$q")" >&2; return 1; }
    echo "$rid" > "${out}.rid"
  fi
  echo "       RID $rid ($(date +%H:%M))"
  for (( i = 1; i <= ESPERA_MAX_H * 60; i++ )); do
    sleep 60
    st="$(curl -s --max-time 120 "$BLAST_URL?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid" \
            | grep -oE 'Status=[A-Z]+' | head -n1 | cut -d= -f2 || true)"
    case "$st" in
      READY) break ;;
      UNKNOWN|FAILED)
        echo "ERROR: RID $rid termino en estado $st; borra ${out}.rid y relanza" >&2
        return 1 ;;
    esac
    (( i % 30 == 0 )) && echo "       ... sigue en cola ($(date +%H:%M), ${st:-sin respuesta})"
  done
  [[ "$st" == "READY" ]] || { echo "ERROR: RID $rid sin terminar tras $ESPERA_MAX_H h" >&2; return 1; }
  for i in 1 2 3 4 5; do
    if blast_formatter -rid "$rid" -outfmt "$OUTFMT" > "$out"; then
      echo "       listo ($(date +%H:%M)): $(wc -l < "$out") hits"
      return 0
    fi
    echo "       blast_formatter fallo (intento $i); reintento en 1 min" >&2
    sleep 60
  done
  echo "ERROR: no se pudo recuperar el RID $rid" >&2
  return 1
}

SOLO="$1"
tail -n +2 "$SEL_DIR/especies_seleccionadas.tsv" | while IFS=$'\t' read -r RANGO SLUG LID ESP NCBI_ESP TAXID RESTO; do
  [[ -n "$SOLO" && "$SLUG" != "$SOLO" ]] && continue
  OUT="$EXCL_OUT/$SLUG"
  if [[ "$MOTOR" == "ebi" ]]; then
    mkdir -p "$OUT"
    for f in candidatos_local.tsv f5_query.faa f6_query.fna versions_local.txt; do
      [[ -f "$EXCL_LOCAL/$SLUG/$f" ]] && cp -f "$EXCL_LOCAL/$SLUG/$f" "$OUT/$f"
    done
    if [[ -s "$OUT/f6_query.fna" ]]; then
      echo "ERROR: $SLUG tiene contigs huerfanos (F6); el EBI no cubre core_nt. Usa MOTOR=ncbi." >&2; exit 1
    fi
  fi
  if [[ ! -s "$OUT/candidatos_local.tsv" ]]; then
    echo "ERROR: falta $OUT/candidatos_local.tsv (corre 08_exclusivos_local.slurm)." >&2; exit 1
  fi
  echo
  echo "===> $SLUG (taxid NCBI $TAXID)"
  {
    echo "Fase 3 - Exclusivos (remoto) - $SLUG"
    echo "Fecha: $(date)   Host: $(hostname)"
    if [[ "$MOTOR" == "ebi" ]]; then
      echo "F5: blastp en el EBI (REST ncbiblast) -db $EBI_DB -evalue 1e-5, 10 hits, una proteina por trabajo"
    else
      echo "F5: blastp (URL API de NCBI + blast_formatter -rid) -db $NR_DB -evalue 1e-5, 10 hits (lotes de $LOTE)"
    fi
    echo "F6: blastn megablast (URL API de NCBI) -db $NT_DB -evalue 1e-10, 10 hits"
    echo "blast: $(blastp -version | head -n1) | taxonkit: $(taxonkit version 2>&1 | head -n1)"
    echo "taxdump: $TAXDUMP_DIR ($(cat "$TAXDUMP_DIR/fecha_descarga.txt" 2>/dev/null || echo NA))"
  } > "$OUT/versions_remoto.txt"

  # F5: lotes de LOTE proteinas (NCBI) o una por trabajo (EBI)
  mkdir -p "$OUT/f5_lotes"
  if [[ "$MOTOR" == "ebi" && -s "$OUT/f5_query.faa" ]]; then
    echo "     F5 EBI: $(grep -c '^>' "$OUT/f5_query.faa") proteinas ($(date +%H:%M))"
    python "$SCRIPTS_DIR/08_exclusivos.py" ebi --query "$OUT/f5_query.faa"         --out "$OUT/f5_blastp_nr.tsv" --estado "$OUT/f5_lotes" --email "$EBI_EMAIL" --db "$EBI_DB"
  elif [[ -s "$OUT/f5_query.faa" ]]; then
    awk -v n="$LOTE" -v d="$OUT/f5_lotes" '/^>/{if(c%n==0){f=sprintf("%s/lote_%04d.faa",d,int(c/n)+1)} c++} {print > f}' \
        "$OUT/f5_query.faa"
    for Q in "$OUT"/f5_lotes/lote_*.faa; do
      R="${Q%.faa}.tsv"
      if [[ -f "${R}.ok" ]]; then continue; fi
      echo "     F5 $(basename "$Q"): $(grep -c '^>' "$Q") proteinas ($(date +%H:%M))"
      blast_ncbi blastp "$NR_DB" 1e-5 "$Q" "$R"
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
    blast_ncbi blastn "$NT_DB" 1e-10 "$OUT/f6_query.fna" "$OUT/f6_blastn_nt.tsv" \
               --data-urlencode "MEGABLAST=on"
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
