# =============================================================================
# euk_02_fig1_embudo.R — Fig. E1: embudo 97 -> 31 -> 19 -> 13 -> 10
# =============================================================================

n_crudos  <- nrow(tiara)
n_euk     <- nrow(sel1)
n_estim   <- sum(!is.na(sel1$eukcc_completitud))
n_sel     <- sum(sel1$categoria == "Seleccionado")
n_linajes <- nrow(linajes)

niveles <- c("Bins crudos", "Eucariotas\n(Tiara)", "Con estimacion\n(EukCC)",
             "Alta calidad\n(EukCC)", "Linajes\n(dRep 95% ANI)")
etapas <- tibble(
  etapa = factor(niveles, levels = niveles),
  n     = c(n_crudos, n_euk, n_estim, n_sel, n_linajes),
  fase  = c("Entrada", "Fase 1", "Fase 1", "Fase 1", "Fase 2")
) %>%
  mutate(perdida = lag(n) - n)

pa <- ggplot(etapas, aes(x = etapa, y = n, fill = fase)) +
  geom_col(width = 0.62, colour = "grey20", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.55, fontface = "bold", size = 3.6) +
  geom_text(aes(y = n / 2, label = if_else(is.na(perdida) | n < 12, "",
                                           paste0("−", perdida))),
            colour = "white", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c("Entrada" = "grey60", "Fase 1" = "#4A7FB5",
                               "Fase 2" = "#1B7837"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = pretty_breaks(5)) +
  labs(x = NULL, y = "Nº de genomas", tag = "a") +
  theme(legend.position = "top", axis.text.x = element_text(size = 8))

motivos <- sel1 %>%
  filter(categoria != "Seleccionado") %>%
  count(categoria) %>%
  mutate(categoria = fct_reorder(categoria, n))

pb <- ggplot(motivos, aes(x = n, y = categoria, fill = categoria)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.3, fontface = "bold", size = 3.4) +
  scale_fill_manual(values = PAL_DECISION_EUK) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Nº de bins eucariotas descartados", y = NULL, tag = "b")

# Panel (c): tamano de cada bin eucariota y su destino. Muestra que los 12 sin
# estimacion son fragmentos pequenos y que los 4 contaminados duplican el
# tamano tipico de los seleccionados.
pc <- sel1 %>%
  mutate(bin = fct_reorder(bin, tamano_mb)) %>%
  ggplot(aes(x = tamano_mb, y = bin, fill = categoria)) +
  geom_col(width = 0.75, colour = "grey25", linewidth = 0.2) +
  scale_fill_manual(values = PAL_DECISION_EUK, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = seq(0, 60, 10)) +
  labs(x = "Tamano del bin (Mb)", y = "Bins eucariotas (n = 31)", tag = "c") +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "top",
        legend.text = element_text(size = 7.5)) +
  guides(fill = guide_legend(nrow = 2))

fig1 <- ((pa / pb + plot_layout(heights = c(1.35, 1))) | pc) +
  plot_layout(widths = c(1, 0.8))

guardar_fig(fig1, "figE1_embudo_seleccion", ancho = 10.5, alto = 6.4)
