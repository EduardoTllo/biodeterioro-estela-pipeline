# =============================================================================
# f3_01_curvas_heaps.R — Curvas de acumulacion del pangenoma y ley de Heaps
#
# Se usa la matriz de presencia/ausencia SOLO de referencias
# (pangenoma_refs.Rtab, decision D22): el bin incompleto no entra en la curva.
#   - Curvas del pan (familias acumuladas) y del core (familias presentes en
#     todos los genomas anadidos) con N_PERM_CURVA ordenes aleatorios.
#   - Ley de Heaps con heaps_rapido() (f3_00_setup.R; equivalente exacto de
#     micropan::heaps, modelo de Tettelin et al. 2008): alpha < 1
#     indica un pangenoma abierto; alpha > 1, cerrado.
# =============================================================================

set.seed(42)

curvas <- list()
heaps_tab <- list()

for (slug in especies_f3$slug) {
  rtab <- leer_especie(slug, "pangenoma_refs.Rtab")
  if (is.null(rtab)) { message("  sin pangenoma_refs.Rtab para ", slug); next }
  pm <- t(as.matrix(rtab[, -1]))          # genomas x familias
  colnames(pm) <- rtab$Gene
  n <- nrow(pm)

  pm <- (pm > 0) * 1L
  pm <- pm[, colSums(pm) > 0, drop = FALSE]
  res <- map_dfr(seq_len(N_PERM_CURVA), function(p) {
    cur <- curva_acumulacion(pm[sample(n), , drop = FALSE])
    tibble(perm = p, n_genomas = seq_len(n), pan = cur$pan, core = cur$core)
  })
  curvas[[slug]] <- res |>
    pivot_longer(c(pan, core), names_to = "curva", values_to = "familias") |>
    group_by(n_genomas, curva) |>
    summarise(mediana = median(familias), q1 = quantile(familias, 0.25),
              q3 = quantile(familias, 0.75), .groups = "drop") |>
    mutate(slug = slug)

  h <- heaps_rapido(pm, n.perm = N_PERM_HEAPS)   # equivalente a micropan::heaps
  heaps_tab[[slug]] <- tibble(
    slug = slug, n_referencias = n, familias_refs = ncol(pm),
    heaps_intercepto = unname(h["Intercept"]), heaps_alpha = unname(h["alpha"]),
    pangenoma = if_else(unname(h["alpha"]) < 1, "abierto", "cerrado")
  )
}

if (length(curvas) > 0) {
  df_curvas <- bind_rows(curvas) |>
    left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
    mutate(curva = factor(curva, levels = c("pan", "core"),
                          labels = c("Pangenoma", "Genoma core")))

  g <- ggplot(df_curvas, aes(n_genomas, mediana, colour = curva, fill = curva)) +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.25, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ etiqueta, scales = "free") +
    scale_colour_manual(values = c("Pangenoma" = "#B2182B", "Genoma core" = "#2166AC")) +
    scale_fill_manual(values = c("Pangenoma" = "#B2182B", "Genoma core" = "#2166AC")) +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "Genomas de referencia anadidos", y = "Familias de genes",
         colour = NULL, fill = NULL,
         caption = sprintf("Mediana y rango intercuartil de %d ordenes aleatorios; solo referencias.",
                           N_PERM_CURVA)) +
    theme(legend.position = "bottom")
  guardar_fig(g, "fig_f3_curvas_acumulacion", ancho = 7.5,
              alto = 2.8 * ceiling(length(curvas) / 3) + 0.6)

  tab_heaps <- bind_rows(heaps_tab) |>
    left_join(especies_f3 |> select(slug, linaje_id, especie_gtdb), by = "slug") |>
    select(linaje_id, especie_gtdb, n_referencias, familias_refs,
           heaps_intercepto, heaps_alpha, pangenoma)
  write_tsv(tab_heaps, file.path(DIR_TAB, "tab_f3_heaps.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_heaps.tsv"))
}
