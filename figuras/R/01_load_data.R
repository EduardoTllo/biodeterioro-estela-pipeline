# =============================================================================
# 01_load_data.R — Carga y armonizacion de los resultados de las Fases 1 y 2
# Produce el objeto `maestra` (una fila por genoma seleccionado) y tablas
# auxiliares usadas por el resto de los scripts.
# =============================================================================

# --- Fase 1: Tiara (dominio) -------------------------------------------------
tiara <- read_tsv(file.path(RES_F1, "phase1_tiara_bin_summary.tsv"),
                  show_col_types = FALSE)

# --- Fase 1: seleccion (embudo completo, 97 bins) ----------------------------
seleccion_f1 <- read_tsv(file.path(RES_F1, "phase1_selection.tsv"),
                         show_col_types = FALSE) %>%
  mutate(
    categoria = case_when(
      decision_final == "SELECCIONADO"                ~ "Seleccionado",
      dominio_decision == "DESCARTADO_dominio"        ~ "Descartado: dominio",
      str_detect(motivo, "completitud") &
        str_detect(motivo, "contaminacion")           ~ "Rechazado: ambos",
      str_detect(motivo, "completitud")               ~ "Rechazado: completitud",
      str_detect(motivo, "contaminacion")             ~ "Rechazado: contaminacion",
      TRUE                                            ~ "Otro"
    ),
    categoria = factor(categoria, levels = names(PAL_DECISION))
  )

# --- Fase 1: CheckM2 (66 bins procariotas) -----------------------------------
checkm2 <- read_tsv(file.path(RES_F1, "phase1_checkm2_quality_report.tsv"),
                    show_col_types = FALSE) %>%
  rename(genoma = Name) %>%
  mutate(
    gc_pct        = GC_Content * 100,
    tamano_mb     = Genome_Size / 1e6,
    n50_kb        = Contig_N50 / 1e3,
    densidad_cod  = Coding_Density * 100
  )

# --- Fase 2: GTDB-Tk summary -------------------------------------------------
gtdbtk <- read_tsv(file.path(RES_F2, "phase2_gtdbtk_bac120_summary.tsv"),
                   show_col_types = FALSE, na = c("", "NA", "N/A")) %>%
  rename(genoma = user_genome) %>%
  separar_taxonomia("classification") %>%
  mutate(
    ani = as.numeric(closest_genome_ani),
    af  = as.numeric(closest_genome_af) * 100,   # GTDB-Tk reporta AF en fraccion
    # Cuando classify_wf no hace ANI, usa la colocacion (placement)
    ani = coalesce(ani, as.numeric(closest_placement_ani)),
    af  = coalesce(af,  as.numeric(closest_placement_af) * 100),
    referencia = coalesce(closest_genome_reference, closest_placement_reference)
  )

# La AF en el resumen de la Fase 2 ya venia en porcentaje; comprobamos escala.
if (max(gtdbtk$af, na.rm = TRUE) <= 1.5) gtdbtk <- mutate(gtdbtk, af = af * 100)

# --- Fase 2: seleccion por linaje --------------------------------------------
seleccion_f2 <- read_tsv(file.path(RES_F2, "phase2_selection.tsv"),
                         show_col_types = FALSE,
                         col_types = cols(muestra = col_character())) %>%
  rename(genoma = genome) %>%
  mutate(linaje_num = as.integer(str_remove(linaje_id, "^L")))

# --- Fase 2: dRep ------------------------------------------------------------
lee_drep <- function(archivo) {
  read_csv(file.path(RES_F2, archivo), show_col_types = FALSE) %>%
    mutate(across(any_of(c("genome", "reference", "querry", "genome1", "genome2")),
                  ~ str_remove(.x, "\\.fasta$")))
}
drep_cdb <- lee_drep("phase2_drep_Cdb.csv")   # clusters secundarios
drep_ndb <- lee_drep("phase2_drep_Ndb.csv")   # fastANI (solo dentro de cluster primario)
drep_mdb <- lee_drep("phase2_drep_Mdb.csv")   # MASH todos contra todos
drep_wdb <- lee_drep("phase2_drep_Wdb.csv")   # representantes ganadores

# --- Fase 2: marcadores bac120 -----------------------------------------------
marcadores <- read_tsv(file.path(RES_F2, "phase2_gtdbtk_markers_summary.tsv"),
                       show_col_types = FALSE) %>%
  select(genoma = name,
         unicos = number_unique_genes,
         multiples = number_multiple_genes,
         faltantes = number_missing_genes) %>%
  mutate(pct_marcadores = 100 * unicos / 120)

# =============================================================================
# TABLA MAESTRA — una fila por genoma seleccionado (n = 21)
# =============================================================================
maestra <- seleccion_f2 %>%
  select(genoma, muestra, linaje_id, linaje_num, linaje_cluster,
         especie_confirmada) %>%
  left_join(select(gtdbtk, genoma, classification, all_of(RANGOS),
                   referencia, ani, af, classification_method, note,
                   msa_percent, red_value, warnings),
            by = "genoma") %>%
  left_join(select(checkm2, genoma, Completeness, Contamination, tamano_mb,
                   gc_pct, n50_kb, Total_Contigs, Max_Contig_Length,
                   densidad_cod, Total_Coding_Sequences),
            by = "genoma") %>%
  left_join(select(marcadores, genoma, pct_marcadores, unicos, faltantes),
            by = "genoma") %>%
  left_join(drep_wdb %>%
              select(genoma = genome, score) %>%
              mutate(representante = TRUE),
            by = "genoma") %>%
  mutate(
    representante = coalesce(representante, FALSE),
    filo          = factor(filo, levels = names(PAL_FILO)),
    estado_especie = if_else(especie_confirmada == "si",
                             "Especie confirmada", "Sin especie asignada"),
    taxon_corto   = if_else(!is.na(especie), especie, paste0(genero, " sp.")),
    etiqueta      = paste0(genoma, " | ", taxon_corto),
    # Calidad segun umbrales de CheckM2 (nota: MIMAG exige ademas rRNA/tRNA)
    calidad_checkm2 = case_when(
      Completeness >= 90 & Contamination < 5  ~ "Alta (>=90/<5)",
      Completeness >= 50 & Contamination < 10 ~ "Media (>=50/<10)",
      TRUE                                    ~ "Baja"
    )
  ) %>%
  arrange(linaje_num, genoma)

# --- Prevalencia por linaje --------------------------------------------------
linajes <- maestra %>%
  group_by(linaje_id, linaje_num) %>%
  summarise(
    n_genomas   = n(),
    muestras    = paste(sort(unique(muestra)), collapse = ", "),
    prevalencia = n_distinct(muestra),
    taxon       = first(taxon_corto),
    filo        = first(filo),
    genero      = first(genero),
    repres      = paste(genoma[representante], collapse = ", "),
    completitud_rep = Completeness[representante][1],
    contam_rep      = Contamination[representante][1],
    especie_conf    = first(especie_confirmada),
    .groups = "drop"
  ) %>%
  mutate(
    fase3 = linaje_id %in% c("L1", "L2", "L3"),
    repres = if_else(repres == "", first(repres), repres)
  ) %>%
  arrange(desc(prevalencia), linaje_num)

message("Datos cargados: ", nrow(maestra), " genomas, ", nrow(linajes), " linajes, ",
        nrow(seleccion_f1), " bins en el embudo de la Fase 1.")
