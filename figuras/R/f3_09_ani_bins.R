# =============================================================================
# f3_09_ani_bins.R — ANI de cada bin contra todos los genomas descargados de su
# especie (06c_ani_bins.slurm, fastANI).
#
#  (a) ANI del bin contra cada genoma, de mayor a menor; color segun si el
#      genoma esta en el pangenoma (referencia) o quedo fuera (redundante por
#      dRep o no elegido por el tope de 50). Lineas en 95 % (especie) y 99 %
#      (misma cepa o clon).
#  (b) Control de coherencia: ANI frente a la distancia patristica en el arbol
#      del core, para las referencias.
# =============================================================================

ani <- map_dfr(especies_f3$slug, function(s) {
  d <- leer_especie(s, "ani_bin_genomas.tsv")
  if (is.null(d)) return(NULL)
  refs <- colnames(leer_especie(s, "pangenoma_refs.Rtab"))[-1]
  arbol <- read.tree(file.path(RES_F3, s, "core.treefile"))
  dp <- cophenetic(arbol)
  b <- grep("^bin_", arbol$tip.label, value = TRUE)
  d |> mutate(slug = s, en_pangenoma = referencia %in% refs,
              dist_arbol = unname(dp[b, ][match(referencia, colnames(dp))]),
              fraccion_alineada = fragmentos_compartidos / fragmentos_totales) |>
    arrange(desc(ani)) |> mutate(rango = row_number())
})

if (nrow(ani)) {
  ani <- ani |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
    left_join(refs_f3 |> select(referencia = id, habitat), by = "referencia")
  top <- ani |> group_by(slug) |> slice_max(ani, n = 3, with_ties = FALSE) |> ungroup() |>
    mutate(rotulo = sprintf("%s (%.2f %%)", sub("_(\\d)$", ".\\1", referencia), ani))
  ga <- ggplot(ani, aes(rango, ani)) +
    geom_hline(yintercept = 95, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = 99, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = en_pangenoma), size = 1.6) +
    geom_text_repel(data = top, aes(label = rotulo), size = 2.3, direction = "y",
                    nudge_x = 15, hjust = 0, segment.size = 0.2, seed = 1) +
    facet_wrap(~ etiqueta, scales = "free_x", labeller = label_wrap_gen(width = 22)) +
    scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#F4A582"),
                        labels = c(`TRUE` = "referencia del pangenoma",
                                   `FALSE` = "fuera del pangenoma (dRep o tope de 50)"),
                        name = NULL) +
    labs(x = "Genomas de la especie, de mayor a menor ANI con el bin", y = "ANI con el bin (%)",
         tag = "a", caption = "Linea discontinua: 99 % (misma cepa o clon); punteada: 95 % (misma especie).") +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))
  rho <- ani |> filter(!is.na(dist_arbol)) |> group_by(etiqueta) |>
    summarise(r = cor(ani, dist_arbol, method = "spearman"), .groups = "drop") |>
    mutate(rot = sprintf("rho = %.2f", r))
  gb <- ggplot(filter(ani, !is.na(dist_arbol)), aes(dist_arbol, ani)) +
    geom_point(colour = "#2166AC", size = 1.6) +
    geom_text(data = rho, aes(x = Inf, y = Inf, label = rot), inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5, size = 2.8) +
    facet_wrap(~ etiqueta, scales = "free", labeller = label_wrap_gen(width = 22)) +
    labs(x = "Distancia en el arbol del core (sustituciones/sitio)", y = "ANI con el bin (%)",
         tag = "b", caption = "Solo referencias del pangenoma. Coherencia entre arbol del core y ANI de genoma completo.") +
    theme(strip.text = element_text(size = 8))
  guardar_fig(ga / gb, "fig_f3_ani_bins", ancho = 9, alto = 7.5)
  write_tsv(ani |> select(slug, consulta, referencia, ani, fraccion_alineada, en_pangenoma,
                          dist_arbol, habitat),
            file.path(DIR_TAB, "tab_f3_ani_bins.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_ani_bins.tsv"))
}
