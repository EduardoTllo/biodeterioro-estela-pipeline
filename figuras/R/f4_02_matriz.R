# =============================================================================
# f4_02_matriz.R — Matriz bin x funcion y puntajes (Fase 4)
#
# Lee las tablas de f4_01_marcadores.R. Cada fila es un bin con su linaje y su
# completitud CheckM2 (la ausencia en un bin incompleto no prueba ausencia).
# =============================================================================

ETIQ_FUNCION <- c(
  acido_gluconico = "Ácido glucónico", acido_2_cetogluconico = "Ácido 2-cetoglucónico",
  eps_biopelicula = "EPS de biopelícula", curli = "Curli", nitrificacion = "Nitrificación",
  oxidacion_azufre = "Oxidación de azufre", esporulacion = "Esporulación",
  carotenoides = "Carotenoides", estres_oxidativo = "Estrés oxidativo",
  desecacion_ectoina = "Ectoína", desecacion_trehalosa = "Trehalosa",
  reparacion_uv = "Reparación UV", ciclo_tca = "Ciclo TCA", ureolisis = "Ureólisis"
)
ETIQ_CATEGORIA <- c(riesgo = "Riesgo", supervivencia = "Supervivencia", contexto = "Contexto")
PAL_PRESENCIA <- c(riesgo = "#B2182B", supervivencia = "#2166AC", contexto = "#4D4D4D",
                   ausente = "grey92")

fun <- read_tsv(file.path(DIR_TAB, "tab_f4_funciones_bin.tsv"), show_col_types = FALSE)
pun <- read_tsv(file.path(DIR_TAB, "tab_f4_puntajes.tsv"), show_col_types = FALSE)

orden_bins <- pun |>
  mutate(etiqueta = paste0(etiqueta_genoma(genoma, genero, especie), " · ", linaje_id, " · ",
                           formatC(completitud, format = "f", digits = 1, decimal.mark = ","), " %")) |>
  arrange(as.integer(sub("^L", "", linaje_id)), genoma)
orden_func <- unique(read_tsv(file.path("metadata", "marcadores_fase4.tsv"),
                              show_col_types = FALSE)$funcion)

d <- fun |>
  left_join(orden_bins |> select(genoma, etiqueta), by = "genoma") |>
  mutate(etiqueta = factor(etiqueta, levels = rev(orden_bins$etiqueta)),
         funcion = factor(funcion, levels = orden_func,
                          labels = coalesce(ETIQ_FUNCION[orden_func], orden_func)),
         categoria = factor(categoria, levels = names(ETIQ_CATEGORIA), labels = ETIQ_CATEGORIA),
         relleno = if_else(presente, names(ETIQ_CATEGORIA)[as.integer(categoria)], "ausente"))

g1 <- ggplot(d, aes(funcion, etiqueta, fill = relleno)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  facet_grid(~ categoria, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL_PRESENCIA, guide = "none") +
  labs(x = NULL, y = NULL, tag = "a") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

dp <- orden_bins |>
  select(etiqueta, starts_with("presentes_"), starts_with("evaluables_")) |>
  pivot_longer(-etiqueta, names_to = c(".value", "categoria"), names_sep = "_") |>
  mutate(etiqueta = factor(etiqueta, levels = rev(orden_bins$etiqueta)),
         categoria = factor(categoria, levels = c("riesgo", "supervivencia"),
                            labels = c("Riesgo", "Superv.")))
g2 <- ggplot(dp, aes(categoria, etiqueta, label = paste0(presentes, "/", evaluables))) +
  geom_text(size = 3) +
  labs(x = "Puntaje", y = NULL, tag = "b") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

guardar_fig(g1 + g2 + plot_layout(widths = c(6, 1.2)), "fig_f4_matriz_funciones",
            ancho = 11, alto = 0.28 * nrow(orden_bins) + 2.5)
