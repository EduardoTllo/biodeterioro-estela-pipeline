# =============================================================================
# euk_07_fig6_prevalencia.R — Fig. E6: matriz linaje x muestra (como la Fig. 7.7)
# Muestras: las que aportan al menos un bin eucariota (31 bins, 17 muestras).
# =============================================================================

muestras <- sort(unique(as.integer(sel1$muestra)))

matriz <- maestra %>%
  count(linaje_id, muestra, name = "n") %>%
  mutate(muestra = as.integer(muestra)) %>%
  complete(linaje_id = linajes$linaje_id, muestra = muestras, fill = list(n = 0)) %>%
  left_join(select(linajes, linaje_id, linaje_num, prevalencia, taxon_corto, clase), by = "linaje_id") %>%
  mutate(etq = paste0(linaje_id, " | ", taxon_corto),
         etq = fct_reorder(etq, -linaje_num),
         muestra = factor(muestra, levels = muestras))

pa <- ggplot(matriz, aes(x = muestra, y = etq)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.7) +
  geom_point(data = filter(matriz, n > 0), aes(colour = clase, shape = clase), size = 3.4) +
  scale_colour_manual(values = PAL_CLASE, name = "Clase", drop = TRUE) +
  scale_shape_manual(values = SHP_CLASE, name = "Clase", drop = TRUE) +
  labs(x = "Muestra (hisopo B, fraccion viable)", y = NULL, tag = "a") +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 8))

pb <- matriz %>% distinct(etq, prevalencia, clase) %>%
  ggplot(aes(x = prevalencia, y = etq)) +
  geom_col(aes(fill = clase), width = 0.62, colour = "grey25", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = prevalencia), hjust = -0.45, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASE) +
  scale_x_continuous(breaks = 0:3, limits = c(0, 3.6), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Prevalencia\n(nº de muestras)", y = NULL, tag = "b") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())

pc <- maestra %>%
  mutate(muestra = factor(as.integer(muestra), levels = muestras)) %>%
  count(muestra, clase, .drop = FALSE) %>%
  ggplot(aes(x = muestra, y = n, fill = clase)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = PAL_CLASE) +
  scale_y_continuous(breaks = 0:3, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Genomas\npor muestra", tag = "c") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig6 <- pc + plot_spacer() + pa + pb +
  plot_layout(ncol = 2, widths = c(3.1, 1), heights = c(0.35, 2.2), guides = "collect")

guardar_fig(fig6, "figE6_prevalencia_linajes", ancho = 9.5, alto = 5.4)
