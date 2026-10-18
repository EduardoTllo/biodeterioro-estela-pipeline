# =============================================================================
# euk_01_load_data.R — Carga y une los resultados de las Fases 1 y 2 eucariotas
# =============================================================================

leer <- function(...) read_tsv(file.path(...), show_col_types = FALSE, progress = FALSE)

# Tiara de la Fase 1 (los 97 bins)
tiara <- leer(RES_F1, "phase1_tiara_bin_summary.tsv")

# Fase 1 eucariota: los 31 bins con mayoria eucariota
sel1 <- leer(RES_E1, "phase1_euk_selection.tsv") %>%
  mutate(
    muestra = as.character(muestra),
    tamano_mb = longitud_bp / 1e6,
    categoria = case_when(
      decision_final == "SELECCIONADO"        ~ "Seleccionado",
      is.na(eukcc_completitud)                ~ "Sin estimacion de EukCC",
      str_detect(motivo, "contaminacion")     ~ "Rechazado: contaminacion",
      str_detect(motivo, "completitud")       ~ "Rechazado: completitud"
    ),
    categoria = factor(categoria, levels = names(PAL_DECISION_EUK))
  )

# Fase 2 eucariota: los 13 genomas seleccionados
sel2 <- leer(RES_E2, "phase2_euk_selection.tsv") %>%
  mutate(muestra = as.character(muestra))
linajes <- leer(RES_E2, "phase2_euk_linajes.tsv")

maestra <- sel2 %>%
  left_join(select(sel1, genome = bin, n_contigs, tamano_mb, n50, gc,
                   eukcc_completitud, eukcc_contaminacion,
                   starts_with("busco_fungi_")), by = "genome") %>%
  left_join(select(linajes, linaje_id, representante, prevalencia), by = "linaje_id") %>%
  mutate(
    representante = genome == representante,
    taxon_corto = taxon_corto_euk(especie, genero, familia, orden, clase),
    clase = factor(clase, levels = names(PAL_CLASE)),
    n50_kb = n50 / 1e3,
    linaje_num = as.integer(str_remove(linaje_id, "E"))
  )

linajes <- linajes %>%
  left_join(distinct(select(maestra, representante = genome, clase, filo,
                            taxon_corto)), by = "representante") %>%
  mutate(linaje_num = as.integer(str_remove(linaje_id, "E")))

# dRep: Mash todos contra todos y fastANI dentro de los clusteres
drep_mdb <- read_csv(file.path(RES_E2, "04_drep", "data_tables", "Mdb.csv"),
                     show_col_types = FALSE) %>%
  mutate(across(c(genome1, genome2), ~ str_remove(.x, "\\.fasta$")))
drep_ndb <- read_csv(file.path(RES_E2, "04_drep", "data_tables", "Ndb.csv"),
                     show_col_types = FALSE) %>%
  mutate(across(c(reference, querry), ~ str_remove(.x, "\\.fasta$")))

# skani: todos los hits (para la figura de novedad)
skani <- leer(RES_E2, "01_skani", "skani_bins_vs_refs.tsv") %>%
  transmute(genome = str_remove(basename(Query_file), "\\.fasta$"),
            ani = ANI, af_query = Align_fraction_query,
            ref = str_extract(Ref_name, "[A-Z][a-z]+ [a-z]+"))

message("Fase 1: ", nrow(sel1), " bins eucariotas; Fase 2: ", nrow(maestra),
        " genomas en ", nrow(linajes), " linajes")
