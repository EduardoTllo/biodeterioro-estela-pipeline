# =============================================================================
# 02_fig1_embudo.R — Figura 1: embudo de seleccion 97 -> 66 -> 21 -> 19
# =============================================================================

n_crudos  <- nrow(seleccion_f1)
n_prok    <- sum(seleccion_f1$dominio_decision == "PASA_dominio")
n_sel     <- sum(seleccion_f1$decision_final == "SELECCIONADO")
n_linajes <- nrow(linajes)

etapas <- tibble(
  etapa = factor(
    c("Bins crudos", "Procariotas\n(Tiara)", "Alta calidad\n(CheckM2)",
      "Linajes\n(dRep 95% ANI)"),
    levels = c("Bins crudos", "Procariotas\n(Tiara)", "Alta calidad\n(CheckM2)",
               "Linajes\n(dRep 95% ANI)")
  ),
  n    = c(n_crudos, n_prok, n_sel, n_linajes),
  fase = c("Entrada", "Fase 1 (OE1)", "Fase 1 (OE1)", "Fase 2 (OE2)")
) %>%
  mutate(perdida = lag(n) - n)

# --- Panel (a): embudo -------------------------------------------------------
pa <- ggplot(etapas, aes(x = etapa, y = n, fill = fase)) +
  geom_col(width = 0.62, colour = "grey20", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.55, fontface = "bold", size = 3.6) +
  geom_text(aes(y = n / 2, label = if_else(is.na(perdida), "",
                                           paste0("−", perdida))),
            colour = "white", fontface = "bold", size = 3.2, na.rm = TRUE) +
  scale_fill_manual(values = c("Entrada" = "grey60",
                               "Fase 1 (OE1)" = "#4A7FB5",
                               "Fase 2 (OE2)" = "#1B7837"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = pretty_breaks(5)) +
  labs(x = NULL, y = "Nº de genomas", tag = "a") +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 8.5))

# --- Panel (b): motivos de descarte -------------------------------------------
motivos <- seleccion_f1 %>%
  filter(categoria != "Seleccionado") %>%
  count(categoria) %>%
  mutate(categoria = fct_reorder(categoria, n))

pb <- ggplot(motivos, aes(x = n, y = categoria, fill = categoria)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.3, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.3, fontface = "bold", size = 3.4) +
  scale_fill_manual(values = PAL_DECISION) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Nº de bins descartados", y = NULL, tag = "b")

# --- Panel (c): dominio asignado por Tiara ------------------------------------
pc <- tiara %>%
  select(bin, pct_prok, pct_euk, pct_organelle, pct_unknown, majority) %>%
  pivot_longer(starts_with("pct_"), names_to = "dominio", values_to = "pct") %>%
  mutate(
    dominio = recode(dominio,
                     pct_prok = "Procariota", pct_euk = "Eucariota",
                     pct_organelle = "Organelo", pct_unknown = "Desconocido"),
    dominio = factor(dominio, levels = c("Procariota", "Eucariota",
                                         "Organelo", "Desconocido")),
    grupo   = if_else(majority == "prok", "Mayoria procariota (n = 66)",
                      "Mayoria eucariota (n = 31)"),
    bin     = fct_reorder(bin, pct * (dominio == "Procariota"), .fun = max)
  ) %>%
  ggplot(aes(x = pct, y = bin, fill = dominio)) +
  geom_col(width = 1) +
  facet_grid(grupo ~ ., scales = "free_y", space = "free_y") +
  # Paleta de alto contraste: los cuatro dominios deben distinguirse dentro de
  # una misma barra, incluidos los segmentos muy delgados de organelo.
  scale_fill_manual(values = c("Procariota" = "#08519C", "Eucariota" = "#E6550D",
                               "Organelo" = "#FFD92F", "Desconocido" = "#969696"),
                    name = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = label_percent(scale = 1)) +
  labs(x = "% de la longitud del bin clasificada por Tiara", y = "Bins (n = 97)",
       tag = "c") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top")

fig1 <- (pa / pb + plot_layout(heights = c(1.35, 1))) | pc
fig1 <- fig1 + plot_layout(widths = c(1, 0.85))

guardar_fig(fig1, "fig1_embudo_seleccion", ancho = 10.5, alto = 6.4)
