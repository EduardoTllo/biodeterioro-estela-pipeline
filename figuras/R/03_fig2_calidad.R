# =============================================================================
# 03_fig2_calidad.R — Figura 2: calidad y metricas de ensamblaje de los genomas
# Panel (a) completitud vs contaminacion | (b) tamano vs GC | (c) metricas por filo
# =============================================================================

# 66 bins procariotas con CheckM2, etiquetados con su destino en la Fase 1
calidad <- checkm2 %>%
  left_join(select(seleccion_f1, genoma = bin, categoria), by = "genoma") %>%
  left_join(select(maestra, genoma, filo, genero, taxon_corto, linaje_id),
            by = "genoma")

# --- Panel (a): completitud vs contaminacion ---------------------------------
pa <- ggplot(calidad, aes(x = Completeness, y = Contamination)) +
  # region de calidad alta segun umbrales del estudio
  annotate("rect", xmin = UMBRAL_COMPLETITUD, xmax = 101,
           ymin = -0.5, ymax = UMBRAL_CONTAMINACION,
           fill = "#1B7837", alpha = 0.07) +
  geom_vline(xintercept = UMBRAL_COMPLETITUD, linetype = "dashed",
             colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = UMBRAL_CONTAMINACION, linetype = "dashed",
             colour = "grey35", linewidth = 0.4) +
  geom_point(aes(colour = categoria, shape = categoria),
             size = 2.2, stroke = 0.6, alpha = 0.9) +
  scale_colour_manual(values = PAL_DECISION, name = NULL, drop = TRUE) +
  scale_shape_manual(values = c("Seleccionado" = 16, "Rechazado: contaminacion" = 17,
                                "Rechazado: completitud" = 15, "Rechazado: ambos" = 18),
                     name = NULL, drop = TRUE) +
  scale_x_continuous(limits = c(0, 101), breaks = seq(0, 100, 20)) +
  # Escala log(1+x): sin ella los pocos bins con >60 % de contaminacion aplastan
  # la region de interes (0-5 %), donde caen todos los genomas seleccionados.
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 1, 2, 5, 10, 25, 50, 100),
                     limits = c(0, 125)) +
  labs(x = "Completitud CheckM2 (%)", y = "Contaminacion CheckM2 (%, escala log)",
       tag = "a") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# --- Panel (b): tamano del genoma vs contenido GC -----------------------------
pb <- maestra %>%
  ggplot(aes(x = tamano_mb, y = gc_pct)) +
  geom_point(aes(colour = filo, size = n50_kb), alpha = 0.85) +
  geom_text_repel(aes(label = genoma, colour = filo), size = 2.5,
                  max.overlaps = 30, min.segment.length = 0.15,
                  segment.linewidth = 0.25, show.legend = FALSE,
                  seed = 42) +
  scale_colour_manual(values = PAL_FILO, name = "Filo") +
  scale_size_continuous(name = "N50 (kb)", range = c(1.8, 6.5),
                        breaks = c(50, 200, 500)) +
  labs(x = "Tamano del genoma (Mb)", y = "Contenido GC (%)", tag = "b") +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3)),
         size = guide_legend(order = 2))

# --- Panel (c): metricas de ensamblaje por filo -------------------------------
metricas_largo <- maestra %>%
  select(genoma, filo,
         `N50 (kb)` = n50_kb,
         `Nº de contigs` = Total_Contigs,
         `Densidad codificante (%)` = densidad_cod,
         `Completitud (%)` = Completeness) %>%
  pivot_longer(-c(genoma, filo), names_to = "metrica", values_to = "valor") %>%
  mutate(metrica = factor(metrica, levels = c("N50 (kb)", "Nº de contigs",
                                              "Densidad codificante (%)",
                                              "Completitud (%)")))

pc <- ggplot(metricas_largo, aes(x = filo, y = valor)) +
  geom_boxplot(aes(fill = filo), width = 0.55, outlier.shape = NA,
               alpha = 0.35, linewidth = 0.4, show.legend = FALSE) +
  geom_jitter(aes(colour = filo), width = 0.13, height = 0, size = 1.7,
              alpha = 0.9, show.legend = FALSE) +
  facet_wrap(~ metrica, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = PAL_FILO) +
  scale_colour_manual(values = PAL_FILO) +
  labs(x = NULL, y = NULL, tag = "c") +
  theme(axis.text.x = element_text(size = 8, angle = 22, hjust = 1))

fig2 <- (pa | pb) / pc + plot_layout(heights = c(1.25, 1))

guardar_fig(fig2, "fig2_calidad_genomica", ancho = 11, alto = 8)
