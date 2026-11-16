# =============================================================================
# f4_03_cazy.R — Asignaciones CAZy por bin y clase (Fase 4, anexo)
#
# Usa orden_bins de f4_02_matriz.R para mantener el mismo orden y etiquetas.
# =============================================================================

ETIQ_CAZY <- c(GH = "GH (glucósido hidrolasas)", GT = "GT (glicosiltransferasas)",
               CE = "CE (carbohidrato esterasas)", CBM = "CBM (módulos de unión)",
               AA = "AA (actividades auxiliares)", PL = "PL (polisacárido liasas)")

cz <- read_tsv(file.path(DIR_TAB, "tab_f4_cazy_bin.tsv"), show_col_types = FALSE) |>
  count(genoma, clase_cazy, wt = genes, name = "n") |>
  left_join(orden_bins |> select(genoma, etiqueta), by = "genoma") |>
  mutate(etiqueta = factor(etiqueta, levels = rev(orden_bins$etiqueta)),
         clase_cazy = factor(clase_cazy, levels = names(ETIQ_CAZY), labels = ETIQ_CAZY))
stopifnot(!anyNA(cz$etiqueta), !anyNA(cz$clase_cazy))

g <- ggplot(cz, aes(n, etiqueta, fill = clase_cazy)) +
  geom_col(width = 0.75) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Asignaciones gen-familia CAZy", y = NULL) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  guides(fill = guide_legend(nrow = 3))

guardar_fig(g, "fig_f4_cazy_bin", ancho = 8, alto = 0.28 * nrow(orden_bins) + 2)
