# =============================================================================
# euk_03_fig2_calidad.R — Fig. E2: calidad de los bins eucariotas
# (a) EukCC completitud vs contaminacion | (b) tamano vs GC de los 13
# (c) EukCC vs BUSCO completitud        | (d) EukCC contaminacion vs BUSCO D
# =============================================================================

estim <- sel1 %>% filter(!is.na(eukcc_completitud))

etq_rech <- function(d) filter(d, categoria != "Seleccionado")

pa <- ggplot(estim, aes(x = eukcc_completitud, y = eukcc_contaminacion)) +
  annotate("rect", xmin = UMBRAL_COMPLETITUD, xmax = 101,
           ymin = -0.5, ymax = UMBRAL_CONTAMINACION, fill = "#1B7837", alpha = 0.07) +
  geom_vline(xintercept = UMBRAL_COMPLETITUD, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = UMBRAL_CONTAMINACION, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_point(aes(colour = categoria, shape = categoria), size = 2.4, stroke = 0.6) +
  geom_text_repel(data = etq_rech(estim), aes(label = bin, colour = categoria),
                  size = 2.5, seed = 3, min.segment.length = 0.1, show.legend = FALSE) +
  scale_colour_manual(values = PAL_DECISION_EUK, name = NULL, drop = TRUE) +
  scale_shape_manual(values = c("Seleccionado" = 16, "Rechazado: contaminacion" = 17,
                                "Rechazado: completitud" = 15), name = NULL, drop = TRUE) +
  scale_x_continuous(limits = c(35, 101), breaks = seq(40, 100, 10)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 1, 2, 5, 10, 25, 50, 100), limits = c(0, 125)) +
  labs(x = "Completitud EukCC (%)", y = "Contaminacion EukCC (%, escala log)", tag = "a") +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

pb <- ggplot(maestra, aes(x = tamano_mb, y = gc)) +
  geom_point(aes(colour = clase, shape = clase, size = n50_kb), alpha = 0.9) +
  geom_text_repel(aes(label = genome, colour = clase), size = 2.5, seed = 42,
                  min.segment.length = 0.15, segment.linewidth = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = PAL_CLASE, name = "Clase", drop = TRUE) +
  scale_shape_manual(values = SHP_CLASE, name = "Clase", drop = TRUE) +
  scale_size_continuous(name = "N50 (kb)", range = c(1.8, 6.5), breaks = c(10, 100, 400)) +
  labs(x = "Tamano del genoma (Mb)", y = "Contenido GC (%)", tag = "b") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)),
         shape = guide_legend(order = 1), size = guide_legend(order = 2))

pc <- ggplot(estim, aes(x = eukcc_completitud, y = busco_fungi_C)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
  geom_point(aes(colour = categoria, shape = categoria), size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = PAL_DECISION_EUK) +
  scale_shape_manual(values = c("Seleccionado" = 16, "Rechazado: contaminacion" = 17,
                                "Rechazado: completitud" = 15)) +
  coord_equal(xlim = c(35, 101), ylim = c(35, 101)) +
  labs(x = "Completitud EukCC (%)", y = "Completitud BUSCO fungi (%)", tag = "c")

pd <- ggplot(estim, aes(x = eukcc_contaminacion, y = busco_fungi_D)) +
  geom_point(aes(colour = categoria, shape = categoria), size = 2.4, show.legend = FALSE) +
  geom_text_repel(data = etq_rech(estim) %>% filter(eukcc_contaminacion > 5),
                  aes(label = bin, colour = categoria), size = 2.5, seed = 5,
                  show.legend = FALSE) +
  scale_colour_manual(values = PAL_DECISION_EUK) +
  scale_shape_manual(values = c("Seleccionado" = 16, "Rechazado: contaminacion" = 17,
                                "Rechazado: completitud" = 15)) +
  scale_x_continuous(limits = c(0, 102)) +
  scale_y_continuous(limits = c(0, 80)) +
  labs(x = "Contaminacion EukCC (%)", y = "BUSCO duplicados, D (%)", tag = "d")

fig2 <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.2, 1))

guardar_fig(fig2, "figE2_calidad_genomica", ancho = 11, alto = 8.4)
