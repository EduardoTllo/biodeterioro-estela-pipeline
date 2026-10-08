# =============================================================================
# f3_04_exclusivos.R — Embudo de verificacion de genes exclusivos (F1-F6) y
# origen probable de los verificados segun su mejor hit en nr (decision D26).
# =============================================================================

PASOS <- c(candidatos = "Candidatos", pasa_F1 = "F1 estructura",
           pasa_F4 = "F4 especie ampliada", pasa_F6 = "F6 contigs huerfanos",
           pasa_F5 = "F5 nr (misma especie)", verificado = "Verificados")

embudo <- map_dfr(especies_f3$slug, function(s) {
  d <- leer_especie(s, "embudo_exclusivos.tsv")
  if (is.null(d)) NULL else mutate(d, slug = s)
})
exclus <- map_dfr(especies_f3$slug, function(s) {
  d <- leer_especie(s, "exclusivos_verificados.tsv", col_types = cols(.default = "c"))
  if (is.null(d)) NULL else mutate(d, slug = s)
})

if (nrow(embudo) > 0) {
  df_e <- embudo |>
    filter(paso %in% names(PASOS)) |>
    left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
    mutate(paso = factor(PASOS[paso], levels = rev(PASOS)))
  ga <- ggplot(df_e, aes(n, paso)) +
    geom_col(fill = "#B2182B", width = 0.7) +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    facet_wrap(~ etiqueta, scales = "free_x", labeller = label_wrap_gen(width = 20)) +
    theme(strip.text = element_text(size = 8)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Genes", y = NULL, tag = "a")

  graficos <- list(ga)
  if (nrow(exclus) > 0) {
    df_o <- exclus |>
      filter(clase == "verificado") |>
      count(slug, f5_origen) |>
      left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
      mutate(f5_origen = factor(f5_origen, levels = c(
        "misma_especie", "mismo_genero", "misma_familia", "mismo_filo", "otro_filo",
        "sin_hit_ORFan", "taxonomia_desconocida")))
    if (nrow(df_o) > 0) {
      graficos[[2]] <- ggplot(df_o, aes(n, etiqueta, fill = f5_origen)) +
        geom_col(width = 0.6, colour = "grey30", linewidth = 0.2) +
        scale_fill_brewer(palette = "RdYlBu", direction = -1, drop = TRUE,
                          name = "Mejor hit en nr") +
        labs(x = "Exclusivos verificados", y = NULL, tag = "b")
    }
  }
  guardar_fig(wrap_plots(graficos, ncol = 1), "fig_f3_exclusivos",
              ancho = 7.5, alto = 3 + 2.2 * length(graficos))
  write_tsv(embudo, file.path(DIR_TAB, "tab_f3_embudo_exclusivos.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_embudo_exclusivos.tsv"))
}
