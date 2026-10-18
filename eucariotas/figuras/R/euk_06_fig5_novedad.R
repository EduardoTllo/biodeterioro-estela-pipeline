# =============================================================================
# euk_06_fig5_novedad.R — Fig. E5: ANI vs AF frente a la referencia de NCBI
# mas cercana (skani), como la Fig. 7.5
# =============================================================================

nov <- maestra %>%
  filter(!is.na(ani), !is.na(af_query)) %>%
  mutate(estado = if_else(ani >= UMBRAL_ANI & af_query >= UMBRAL_AF,
                          "Especie por ANI", "Sin especie (bajo umbral)"),
         estado = factor(estado, levels = c("Especie por ANI", "Sin especie (bajo umbral)")),
         # Como en la Fig. 7.5: los de especie confirmada solo con su codigo (el
         # taxon esta en la Tabla TE1); el resto, con su vecino mas cercano.
         etq = if_else(estado == "Especie por ANI", genome,
                       paste0(genome, "\nvecino: ", skani_organismo)))
sin_hit <- maestra %>% filter(is.na(ani))

fig5 <- ggplot(nov, aes(x = ani, y = af_query)) +
  annotate("rect", xmin = UMBRAL_ANI, xmax = 100.6, ymin = UMBRAL_AF, ymax = 102,
           fill = "#2166AC", alpha = 0.07) +
  geom_vline(xintercept = UMBRAL_ANI, linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  geom_hline(yintercept = UMBRAL_AF, linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  geom_point(aes(colour = clase, shape = clase), size = 3) +
  geom_text_repel(aes(label = etq, colour = clase), size = 2.5, lineheight = 0.9,
                  seed = 11, box.padding = 0.5, point.padding = 0.3, force = 3,
                  min.segment.length = 0, max.overlaps = 50, show.legend = FALSE) +
  annotate("text", x = 100.4, y = UMBRAL_AF - 3, hjust = 1, size = 2.7, colour = "grey30",
           label = paste0("especie: ANI >= ", UMBRAL_ANI, " % y AF >= ", UMBRAL_AF, " %")) +
  annotate("label", x = 80.3, y = 100, hjust = 0, vjust = 1, size = 2.5,
           colour = "grey30", fill = alpha("white", 0.8),
           label = paste0("Sin referencia con ANI > ~80 % (n = ", nrow(sin_hit), "):\n",
                          paste0(sin_hit$genome, " (", sin_hit$taxon_corto, ")", collapse = "\n"))) +
  scale_colour_manual(values = PAL_CLASE, name = "Clase", drop = TRUE) +
  scale_shape_manual(values = SHP_CLASE, name = "Clase", drop = TRUE) +
  scale_x_continuous(limits = c(80, 100.7), breaks = seq(80, 100, 4)) +
  scale_y_continuous(limits = c(0, 102), breaks = seq(0, 100, 20)) +
  labs(x = "ANI con la referencia de NCBI mas cercana (%)",
       y = "Fraccion del bin alineada, AF (%)")

guardar_fig(fig5, "figE5_novedad_taxonomica", ancho = 8.8, alto = 6)
