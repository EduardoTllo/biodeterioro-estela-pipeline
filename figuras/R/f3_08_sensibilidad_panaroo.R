# =============================================================================
# f3_08_sensibilidad_panaroo.R — Sensibilidad del pangenoma al modo de limpieza
# de Panaroo: moderate (el usado) frente a sensitive (sin poda de extremos de
# contig). Lee results/fase3/<especie>/ (moderate) y <especie>/sensitive/.
#
#  Tabla tab_f3_sensibilidad_panaroo.tsv: familias, core, shell, cloud,
#  singletons, alpha de Heaps y, para el bin, genes fuera del pangenoma y
#  familias sin referencias (candidatos), en ambos modos.
#  Figura fig_f3_sensibilidad_panaroo: (a) particion en ambos modos;
#  (b) curvas de acumulacion del pangenoma en ambos modos.
# =============================================================================

MODOS <- c(moderate = "", sensitive = "sensitive")

leer_modo <- function(slug, sub, archivo, ...) {
  p <- file.path(RES_F3, slug, sub, archivo)
  if (!file.exists(p)) return(NULL)
  read_tsv(p, show_col_types = FALSE, ...)
}

filas <- list(); curvas <- list()
set.seed(8)
for (s in especies_f3$slug) {
  if (is.null(leer_modo(s, "sensitive", "resumen_particion.tsv"))) next
  for (m in names(MODOS)) {
    rp <- leer_modo(s, MODOS[[m]], "resumen_particion.tsv")
    rb <- leer_modo(s, MODOS[[m]], "resumen_bins.tsv")
    rt <- leer_modo(s, MODOS[[m]], "pangenoma_refs.Rtab")
    x <- as.matrix(rt[, -1]); x <- x[rowSums(x) > 0, , drop = FALSE]
    pm <- t(x)
    filas[[paste(s, m)]] <- tibble(
      slug = s, modo = m, referencias = ncol(x), familias = nrow(x),
      core = rp$core_95, shell = rp$shell_95, cloud = rp$cloud_95,
      singletons = sum(rowSums(x) == 1),
      alpha = unname(heaps_rapido(pm, n.perm = N_PERM_HEAPS)["alpha"]),
      bin_genes = rb$n_genes[1], bin_fuera_pangenoma = rb$no_en_pangenoma[1],
      bin_sin_referencias = rb$exclusivo_candidato[1])
    curvas[[paste(s, m)]] <- map_dfr(seq_len(N_PERM_CURVA), function(i) {
      cc <- curva_acumulacion(pm[sample(nrow(pm)), , drop = FALSE])
      tibble(k = seq_along(cc$pan), pan = cc$pan)
    }) |> group_by(k) |> summarise(pan = median(pan), .groups = "drop") |>
      mutate(slug = s, modo = m)
  }
}

if (length(filas)) {
  tab <- bind_rows(filas) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")
  write_tsv(tab, file.path(DIR_TAB, "tab_f3_sensibilidad_panaroo.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_sensibilidad_panaroo.tsv"))

  largo <- tab |> select(etiqueta, modo, core, shell, cloud) |>
    pivot_longer(c(core, shell, cloud), names_to = "categoria", values_to = "n") |>
    mutate(categoria = factor(categoria, levels = c("core", "shell", "cloud")))
  ga <- ggplot(largo, aes(modo, n, fill = categoria)) +
    geom_col(width = 0.65, colour = "grey30", linewidth = 0.2) +
    geom_text(aes(label = scales::comma(n)), position = position_stack(vjust = 0.5), size = 2.6) +
    facet_wrap(~ etiqueta, scales = "free_y", labeller = label_wrap_gen(width = 22)) +
    scale_fill_manual(values = PAL_CATEGORIA[c("core", "shell", "cloud")], name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Modo de limpieza de Panaroo", y = "Familias (solo referencias)", tag = "a") +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
  cv <- bind_rows(curvas) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")
  rot <- tab |> transmute(etiqueta, modo, rot = sprintf("%s: alpha = %.2f", modo, alpha)) |>
    group_by(etiqueta) |> summarise(rot = paste(rot, collapse = "\n"), .groups = "drop")
  gb <- ggplot(cv, aes(k, pan, colour = modo)) +
    geom_line(linewidth = 0.9) +
    geom_text(data = rot, aes(x = -Inf, y = Inf, label = rot), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.2, size = 2.6) +
    facet_wrap(~ etiqueta, scales = "free", labeller = label_wrap_gen(width = 22)) +
    scale_colour_manual(values = c(moderate = "#2166AC", sensitive = "#B2182B"), name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Genomas de referencia anadidos", y = "Familias (pangenoma)", tag = "b",
         caption = "Mediana de 100 ordenes aleatorios. alpha de Heaps con 500 permutaciones.") +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
  guardar_fig(ga / gb, "fig_f3_sensibilidad_panaroo", ancho = 9, alto = 8)
} else {
  message("  (sin resultados de Panaroo sensitive: se omite f3_08)")
}
