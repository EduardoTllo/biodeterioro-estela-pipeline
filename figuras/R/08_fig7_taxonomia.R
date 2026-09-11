# =============================================================================
# 08_fig7_taxonomia.R — Figura 7: composicion taxonomica del set seleccionado
#                       Figura 8 (suplementaria): calidad de la asignacion
# =============================================================================

# --- Figura 7: composicion ----------------------------------------------------
pa <- maestra %>%
  count(filo, orden, genero) %>%
  mutate(genero = fct_reorder2(genero, filo, n)) %>%
  ggplot(aes(x = n, y = genero, fill = filo)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), hjust = -0.4, size = 2.9, fontface = "bold") +
  facet_grid(filo ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = PAL_FILO, name = "Filo", guide = "none") +
  scale_x_continuous(breaks = 0:3, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Nº de genomas", y = NULL, tag = "a") +
  theme(strip.text.y = element_text(angle = 0, size = 7.5),
        axis.text.y = element_text(face = "italic", size = 8))

pb <- maestra %>%
  count(filo, clase, orden) %>%
  mutate(orden = fct_reorder(orden, n)) %>%
  ggplot(aes(x = n, y = orden, fill = filo)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), hjust = -0.4, size = 2.9, fontface = "bold") +
  scale_fill_manual(values = PAL_FILO, name = "Filo") +
  scale_x_continuous(breaks = 0:8, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Nº de genomas", y = NULL, tag = "b") +
  theme(axis.text.y = element_text(size = 8), legend.position = "bottom")

resumen_filo <- maestra %>%
  count(filo, name = "n") %>%
  mutate(pct = 100 * n / sum(n),
         etq = paste0(filo, "\n", n, " (", sprintf("%.0f", pct), " %)"))

pc <- resumen_filo %>%
  # El nombre del filo va en la leyenda comun del panel (b): dentro de los
  # segmentos estrechos no cabe.
  mutate(etq = paste0(n, " (", sprintf("%.0f", pct), " %)")) %>%
  ggplot(aes(x = n, y = "", fill = filo)) +
  geom_col(width = 0.42, colour = "white", linewidth = 0.9,
           position = position_stack(reverse = TRUE), show.legend = FALSE) +
  geom_text(aes(label = etq), position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.7, colour = "white", fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = PAL_FILO) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Nº de genomas (n = 21)", y = NULL, tag = "c") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank())

fig7 <- pa | (pb / pc + plot_layout(heights = c(3.2, 1)))
guardar_fig(fig7, "fig7_composicion_taxonomica", ancho = 10, alto = 6.2)

# --- Figura 8 (suplementaria): calidad de la asignacion taxonomica -----------
asig <- maestra %>%
  mutate(etq = fct_reorder(paste0(genoma, " | ", taxon_corto), msa_percent))

sa <- ggplot(asig, aes(x = msa_percent, y = etq, fill = filo)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.25) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "#B2182B",
             linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f", msa_percent)), hjust = -0.15, size = 2.4) +
  annotate("text", x = 50, y = nrow(asig) + 0.4, label = "minimo GTDB-Tk (50 %)",
           hjust = -0.05, size = 2.5, colour = "#B2182B") +
  scale_fill_manual(values = PAL_FILO, name = "Filo") +
  scale_x_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  labs(x = "Columnas del alineamiento bac120 retenidas (%)", y = NULL, tag = "a") +
  theme(axis.text.y = element_text(size = 7))

sb <- asig %>%
  select(etq, filo, Unicos = unicos, Faltantes = faltantes) %>%
  pivot_longer(c(Unicos, Faltantes), names_to = "tipo", values_to = "n") %>%
  mutate(tipo = factor(tipo, levels = c("Unicos", "Faltantes"))) %>%
  ggplot(aes(x = n, y = etq, fill = tipo)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Unicos" = "#4A7FB5", "Faltantes" = "grey80"),
                    name = "Marcadores bac120") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 120, 30)) +
  labs(x = "Nº de marcadores (de 120)", y = NULL, tag = "b") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

figS <- sa | sb
figS <- figS + plot_layout(widths = c(1.7, 1))
guardar_fig(figS, "figS1_calidad_asignacion", ancho = 10, alto = 5.4)
