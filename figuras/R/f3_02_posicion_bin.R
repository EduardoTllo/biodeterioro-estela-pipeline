# =============================================================================
# f3_02_posicion_bin.R — Posicion de cada bin de la Estela en el pangenoma
#
# (a) Genes del bin por categoria de su familia (frecuencias de referencias).
# (b) Recuperacion del core de referencias vs completitud CheckM2 del bin:
#     control de coherencia del checkpoint C6 (banda de +/- 10 puntos).
# Tabla resumen del pangenoma por especie (corte core 95 % y 90 %).
# =============================================================================

resumen_bins <- map_dfr(especies_f3$slug, function(s) {
  d <- leer_especie(s, "resumen_bins.tsv", col_types = cols(.default = "c"))
  if (is.null(d)) NULL else mutate(d, slug = s)
})
resumen_part <- map_dfr(especies_f3$slug, function(s) {
  d <- leer_especie(s, "resumen_particion.tsv")
  if (is.null(d)) NULL else d
})

if (nrow(resumen_bins) > 0) {
  rb <- resumen_bins |>
    mutate(across(c(n_genes, core, shell, cloud, exclusivo_candidato, no_en_pangenoma,
                    completitud_checkm2, recuperacion_core_pct), as.numeric)) |>
    left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
    mutate(bin_lab = paste0(genoma_original, "\n", etiqueta))

  df_a <- rb |>
    select(bin_lab, core, shell, cloud, exclusivo = exclusivo_candidato, no_en_pangenoma) |>
    pivot_longer(-bin_lab, names_to = "categoria", values_to = "n") |>
    group_by(bin_lab) |> mutate(prop = n / sum(n)) |> ungroup() |>
    mutate(categoria = factor(categoria, levels = names(PAL_CATEGORIA)))

  ga <- ggplot(df_a, aes(prop, bin_lab, fill = categoria)) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = PAL_CATEGORIA, drop = FALSE,
                      labels = c("Core", "Shell", "Cloud", "Exclusivo (candidato)",
                                 "Fuera del pangenoma")) +
    scale_x_continuous(labels = label_percent(), expand = c(0, 0)) +
    labs(x = "Genes del bin", y = NULL, fill = "Categoria de la familia",
         tag = "a") +
    guides(fill = guide_legend(nrow = 2, title.position = "top")) +
    theme(legend.position = "bottom")

  gb <- ggplot(rb, aes(completitud_checkm2, recuperacion_core_pct)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_ribbon(data = tibble(x = c(50, 100)), inherit.aes = FALSE,
                aes(x = x, ymin = x - 10, ymax = x + 10), fill = "grey85", alpha = 0.5) +
    geom_point(size = 2.5, colour = "#B2182B") +
    geom_text_repel(aes(label = genoma_original), size = 3, seed = 1) +
    coord_cartesian(xlim = c(min(60, min(rb$completitud_checkm2, na.rm = TRUE)), 100),
                    ylim = c(min(60, min(rb$recuperacion_core_pct, na.rm = TRUE)), 100)) +
    labs(x = "Completitud CheckM2 (%)", y = "Core de referencias recuperado (%)",
         tag = "b")

  guardar_fig(ga / gb + plot_layout(heights = c(1, 1)), "fig_f3_posicion_bin",
              ancho = 7, alto = 3 + 1.2 * nrow(rb))
  write_tsv(rb |> select(-bin_lab), file.path(DIR_TAB, "tab_f3_posicion_bin.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_posicion_bin.tsv"))
}

if (nrow(resumen_part) > 0) {
  tab <- resumen_part |>
    left_join(especies_f3 |> select(slug, linaje_id, especie_gtdb), by = "slug") |>
    select(linaje_id, especie_gtdb, n_refs, n_bins, familias_refs,
           core_95, shell_95, cloud_95, core_90, familias_exclusivas)
  write_tsv(tab, file.path(DIR_TAB, "tab_f3_resumen_pangenoma.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_resumen_pangenoma.tsv"))
}
