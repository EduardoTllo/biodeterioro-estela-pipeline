# =============================================================================
# 07_fig6_prevalencia.R — Figura 6: matriz linaje x muestra y priorizacion
# para la Fase 3. Criterio primario del OE2: prevalencia espacial.
# =============================================================================

# Todas las muestras con bins en el analisis original (para no ocultar ceros)
muestras_todas <- sort(unique(maestra$muestra))

matriz <- maestra %>%
  count(linaje_id, muestra, name = "n_genomas") %>%
  complete(linaje_id = linajes$linaje_id, muestra = muestras_todas,
           fill = list(n_genomas = 0)) %>%
  left_join(select(linajes, linaje_id, prevalencia, taxon, filo, fase3),
            by = "linaje_id") %>%
  mutate(
    etq = paste0(linaje_id, " | ", taxon),
    etq = fct_reorder(etq, prevalencia * 100 - as.integer(str_remove(linaje_id, "L"))),
    presente = n_genomas > 0
  )

# --- Panel (a): matriz de presencia ------------------------------------------
pa <- ggplot(matriz, aes(x = muestra, y = etq)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.7) +
  geom_point(data = filter(matriz, presente), aes(colour = filo), size = 3.4) +
  scale_colour_manual(values = PAL_FILO, name = "Filo") +
  labs(x = "Muestra (hisopo B, fraccion viable)", y = NULL, tag = "a") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 8))

# --- Panel (b): prevalencia y seleccion --------------------------------------
# Las barras se colorean por filo, igual que los paneles (a) y (c). No se
# resalta ningun linaje: la seleccion para la Fase 3 se decide en esa fase por
# disponibilidad genomica publica y no por prevalencia espacial.
pb <- matriz %>%
  distinct(etq, prevalencia, filo) %>%
  ggplot(aes(x = prevalencia, y = etq)) +
  geom_col(aes(fill = filo), width = 0.62, colour = "grey25", linewidth = 0.3,
           show.legend = FALSE) +
  geom_text(aes(label = prevalencia), hjust = -0.45, size = 2.8,
            fontface = "bold") +
  scale_fill_manual(values = PAL_FILO) +
  scale_x_continuous(breaks = 0:3, limits = c(0, 2.6),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Prevalencia\n(nº de muestras)", y = NULL, tag = "b") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())

# --- Panel (c): riqueza de linajes por muestra --------------------------------
pc <- maestra %>%
  count(muestra, filo) %>%
  ggplot(aes(x = muestra, y = n, fill = filo)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = PAL_FILO) +
  scale_y_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Genomas\npor muestra", tag = "c") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig6 <- pc + plot_spacer() + pa + pb +
  plot_layout(ncol = 2, widths = c(3.1, 1), heights = c(0.35, 2.4),
              guides = "collect")

guardar_fig(fig6, "fig6_prevalencia_linajes", ancho = 9.5, alto = 6.4)
