# =============================================================================
# euk_08_fig7_taxonomia.R — Fig. E7: composicion taxonomica de los 13 genomas
#                           Fig. ES1: perfil BUSCO de los 19 bins con estimacion
# =============================================================================

comp <- maestra %>%
  mutate(grupo = case_when(!is.na(genero)  ~ genero,
                           !is.na(familia) ~ paste0(familia, " (sin genero)"),
                           TRUE            ~ paste0(orden, " (sin familia)")))

pa <- comp %>%
  count(clase, grupo) %>%
  mutate(grupo = fct_reorder(grupo, n)) %>%
  ggplot(aes(x = n, y = grupo, fill = clase)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.4, size = 2.9, fontface = "bold") +
  facet_grid(clase ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = PAL_CLASE) +
  scale_x_continuous(breaks = 0:4, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Nº de genomas", y = NULL, tag = "a") +
  theme(strip.text.y = element_text(angle = 0, size = 7.5),
        axis.text.y = element_text(face = "italic", size = 8))

pb <- comp %>%
  count(clase, orden) %>%
  mutate(orden = fct_reorder(orden, n)) %>%
  ggplot(aes(x = n, y = orden, fill = clase)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), hjust = -0.4, size = 2.9, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASE, name = "Clase", drop = TRUE) +
  scale_x_continuous(breaks = 0:6, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Nº de genomas", y = NULL, tag = "b") +
  theme(axis.text.y = element_text(size = 8), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))

pc <- comp %>%
  count(filo) %>%
  mutate(etq = paste0(n, " (", sprintf("%.0f", 100 * n / sum(n)), " %)"),
         filo = fct_reorder(filo, n)) %>%
  ggplot(aes(x = n, y = filo)) +
  geom_col(width = 0.6, fill = "grey45", colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = etq), hjust = -0.15, size = 2.9, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3)), breaks = seq(0, 12, 4)) +
  labs(x = paste0("Nº de genomas (n = ", nrow(comp), ")"), y = NULL, tag = "c") +
  theme(axis.text.y = element_text(size = 8), panel.grid.major.y = element_blank())

fig7 <- pa | (pb / pc + plot_layout(heights = c(3, 1)))
guardar_fig(fig7, "figE7_composicion_taxonomica", ancho = 10, alto = 5.6)

# --- Fig. ES1: perfil BUSCO (fungi_odb12.2) de los 19 bins con estimacion -----
busco <- sel1 %>%
  filter(!is.na(eukcc_completitud)) %>%
  transmute(bin, categoria, `Completo, copia unica` = busco_fungi_S,
            `Completo, duplicado` = busco_fungi_D, Fragmentado = busco_fungi_F,
            Ausente = busco_fungi_M) %>%
  mutate(bin = fct_reorder(bin, `Completo, copia unica`)) %>%
  pivot_longer(-c(bin, categoria), names_to = "estado", values_to = "pct") %>%
  mutate(estado = factor(estado, levels = c("Completo, copia unica", "Completo, duplicado",
                                            "Fragmentado", "Ausente")))

# Secuencial de un solo tono para "completo" + un tono aparte para duplicado
# (la senal de mezcla que importa leer) y grises para lo incompleto.
figS1 <- ggplot(busco, aes(x = pct, y = bin, fill = estado)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
  facet_grid(categoria ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Completo, copia unica" = "#2a78d6", "Completo, duplicado" = "#eb6834",
                               "Fragmentado" = "#BDBDBD", "Ausente" = "#E5E5E5"), name = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = label_percent(scale = 1)) +
  labs(x = "Genes BUSCO fungi_odb12.2", y = NULL) +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0, size = 7.5),
        axis.text.y = element_text(size = 7.5))

guardar_fig(figS1, "figES1_perfil_busco", ancho = 8.5, alto = 5.6)
