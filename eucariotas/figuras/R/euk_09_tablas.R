# =============================================================================
# euk_09_tablas.R — Tablas del anexo eucariota (TSV)
#   TE1 genomas | TE2 linajes | TE3 contraste con Kraken2 | TE5 los 31 bins
#   (TE4, ANI intracluster, la escribe euk_05_fig4_mash.R)
# =============================================================================

te1 <- maestra %>%
  arrange(linaje_num, genome) %>%
  transmute(Genoma = genome, Muestra = muestra, Linaje = linaje_id,
            Representante = if_else(representante, "si", ""),
            Filo = filo, Clase = as.character(clase), Orden = orden, Familia = familia,
            Genero = genero, Especie = especie, `Nivel de asignacion` = nivel_asignacion,
            `Vecino skani` = skani_organismo, `ANI (%)` = ani, `AF (%)` = af_query,
            `Otras especies >= umbral` = otras_especies_ani,
            `Tamano (Mb)` = round(tamano_mb, 2), Contigs = n_contigs, `N50 (kb)` = round(n50_kb, 1),
            `GC (%)` = gc, `EukCC compl. (%)` = eukcc_completitud,
            `EukCC cont. (%)` = eukcc_contaminacion, `BUSCO C (%)` = busco_fungi_C,
            `BUSCO D (%)` = busco_fungi_D)
write_tsv(te1, file.path(DIR_TAB, "TE1_genomas.tsv"), na = "")

te2 <- linajes %>%
  arrange(linaje_num) %>%
  transmute(Linaje = linaje_id, Taxon = taxon_corto, Clase = as.character(clase),
            Genomas = n_genomas, Prevalencia = prevalencia, Muestras = muestras,
            Representante = representante, Miembros = miembros, Nivel = nivel_asignacion)
write_tsv(te2, file.path(DIR_TAB, "TE2_linajes.tsv"), na = "")

# --- TE3: Kraken2 del centro de secuenciacion frente a la asignacion final ----
kraken_top <- function(genoma, muestra) {
  k <- str_match(genoma, "^bin-(\\d+)-")[, 2]
  f <- file.path(DIR_KRAKEN, muestra, genoma, paste0("bin-", k, "_taxonomy.tsv"))
  if (!file.exists(f)) return(tibble(kraken_genero = NA, kraken_pct = NA, kraken_hongos = NA))
  t <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    setNames(c("pct", "n", "n_dir", "rango", "taxid", "nombre")) %>%
    mutate(pct = as.numeric(pct), nombre = str_trim(nombre))
  g <- t %>% filter(rango == "G") %>% slice_max(pct, n = 1, with_ties = FALSE)
  tibble(kraken_genero = g$nombre[1], kraken_pct = g$pct[1],
         kraken_hongos = t$pct[t$rango == "K" & t$nombre == "Fungi"][1])
}

te3 <- maestra %>%
  arrange(linaje_num) %>%
  mutate(k = map2(genome, muestra, kraken_top)) %>%
  unnest(k) %>%
  transmute(Genoma = genome, `Asignacion final` = taxon_corto, `Genero final` = genero,
            `Genero mas frecuente en Kraken2` = kraken_genero,
            `% de contigs en ese genero` = kraken_pct,
            `% de contigs en Fungi` = kraken_hongos,
            Coincide = case_when(is.na(genero) ~ "sin genero final",
                                 genero == kraken_genero ~ "si", TRUE ~ "no"))
write_tsv(te3, file.path(DIR_TAB, "TE3_kraken2_vs_final.tsv"), na = "")
message("Kraken2 coincide con el genero final en ", sum(te3$Coincide == "si"), " de ",
        sum(te3$Coincide != "sin genero final"), " genomas con genero")
print(as.data.frame(te3))

te5 <- sel1 %>%
  arrange(categoria, desc(tamano_mb)) %>%
  transmute(Bin = bin, Muestra = muestra, `Tamano (Mb)` = round(tamano_mb, 2),
            Contigs = n_contigs, `N50 (kb)` = round(n50 / 1e3, 1), `GC (%)` = gc,
            `Tiara eucariota (%)` = pct_euk_tiara, `EukCC compl. (%)` = eukcc_completitud,
            `EukCC cont. (%)` = eukcc_contaminacion, `BUSCO C (%)` = busco_fungi_C,
            `BUSCO D (%)` = busco_fungi_D, `Linaje EukCC` = eukcc_linaje,
            Decision = as.character(categoria))
write_tsv(te5, file.path(DIR_TAB, "TE5_bins_eucariotas.tsv"), na = "")
