# =============================================================================
# f3_05_diagnostico_heaps.R — Diagnosticos de inflacion de la curva de Heaps
#
# Un pangenoma "inflado" parece mas abierto de lo que es (alpha mas bajo)
# porque cuenta como familias nuevas lo que son artefactos: genes partidos en
# ensamblajes fragmentados, contaminacion u ORFs espurios. Todos los
# diagnosticos usan solo las referencias de cada especie (las <= 50 del
# pangenoma; el bin no entra) y no cambian ningun criterio de seleccion.
#
#   1. Genes unicos por genoma vs. numero de contigs (Spearman). Una
#      correlacion positiva indica que la fragmentacion aporta familias falsas.
#   2. Genomas atipicos: genes unicos > mediana + 3 MAD.
#   3. alpha solo con genomas completos (Complete Genome / Chromosome), si hay
#      al menos MIN_COMPLETOS; si sube mucho respecto de todos, la
#      fragmentacion estaba inflando la curva.
#   4. Rango de alpha por jackknife (dejando fuera un genoma a la vez): muestra
#      cuanto depende alpha de genomas individuales.
#   5. Fraccion de proteinas hipoteticas por categoria (core/shell/cloud). Un
#      cloud dominado por hipoteticas sugiere ORFs espurios.
# =============================================================================

set.seed(42)
MIN_COMPLETOS  <- 5
N_PERM_JACK    <- 200
UMBRAL_MAD     <- 3
NIVELES_COMPLETOS <- c("Complete Genome", "Chromosome")

heaps_alpha <- function(pm, n_perm) {
  if (nrow(pm) < 3) return(NA_real_)
  unname(heaps_rapido(pm, n.perm = n_perm)["alpha"])
}

unicos_tab <- list()
diag_tab   <- list()
hipo_tab   <- list()

for (slug in especies_f3$slug) {
  rtab <- leer_especie(slug, "pangenoma_refs.Rtab")
  if (is.null(rtab)) next
  pm <- t(as.matrix(rtab[, -1]))          # genomas x familias
  colnames(pm) <- rtab$Gene
  n <- nrow(pm)

  meta <- refs_f3 |>
    filter(slug == !!slug) |>
    transmute(id, assembly_level,
              contigs = suppressWarnings(as.numeric(contigs)),
              n50 = suppressWarnings(as.numeric(n50)))

  # 1-2. Genes unicos por genoma y atipicos
  unicos <- tibble(id = rownames(pm),
                   unicos = as.integer(pm[, colSums(pm) == 1, drop = FALSE] |> rowSums())) |>
    left_join(meta, by = "id")
  lim <- median(unicos$unicos) + UMBRAL_MAD * mad(unicos$unicos)
  unicos <- unicos |> mutate(atipico = unicos > lim & unicos > 0, slug = slug)
  unicos_tab[[slug]] <- unicos
  ct <- suppressWarnings(cor.test(unicos$unicos, unicos$contigs,
                                  method = "spearman", exact = FALSE))

  # 3. alpha con todos vs. solo genomas completos
  alpha_todos <- heaps_alpha(pm, N_PERM_HEAPS)
  completos <- unicos$id[unicos$assembly_level %in% NIVELES_COMPLETOS]
  alpha_comp <- if (length(completos) >= MIN_COMPLETOS)
    heaps_alpha(pm[completos, , drop = FALSE], N_PERM_HEAPS) else NA_real_

  # 4. Jackknife: alpha dejando fuera un genoma a la vez
  jack <- vapply(seq_len(n), function(i) heaps_alpha(pm[-i, , drop = FALSE], N_PERM_JACK),
                 numeric(1))

  diag_tab[[slug]] <- tibble(
    slug = slug, n_referencias = n,
    alpha_todos = alpha_todos,
    alpha_jackknife_min = min(jack, na.rm = TRUE),
    alpha_jackknife_max = max(jack, na.rm = TRUE),
    n_completos = length(completos),
    alpha_solo_completos = alpha_comp,
    rho_unicos_vs_contigs = unname(ct$estimate),
    p_unicos_vs_contigs = ct$p.value,
    n_atipicos = sum(unicos$atipico),
    atipicos = paste(unicos$id[unicos$atipico], collapse = ",")
  )

  # 5. Hipoteticas por categoria (familias presentes en referencias)
  fam <- leer_especie(slug, "particion_familias.tsv", col_types = cols(.default = "c"))
  if (!is.null(fam)) {
    hipo_tab[[slug]] <- fam |>
      filter(categoria_95 %in% c("core", "shell", "cloud")) |>
      mutate(hipotetica = str_detect(str_to_lower(coalesce(anotacion, "")),
                                     "hypothetical|uncharacteri[sz]ed|duf[0-9]") |
               coalesce(anotacion, "") %in% c("", "nan")) |>
      group_by(categoria_95) |>
      summarise(familias = n(), hipoteticas = sum(hipotetica), .groups = "drop") |>
      mutate(pct_hipoteticas = 100 * hipoteticas / familias, slug = slug)
  }
}

if (length(diag_tab) > 0) {
  etiquetas <- especies_f3 |> select(slug, linaje_id, especie_gtdb, etiqueta)
  df_u <- bind_rows(unicos_tab) |> left_join(etiquetas, by = "slug")
  df_d <- bind_rows(diag_tab) |> left_join(etiquetas, by = "slug")

  # (a) genes unicos vs contigs; rho de Spearman en el titulo de cada panel
  rotulos <- df_d |>
    mutate(etiqueta_a = if_else(is.na(rho_unicos_vs_contigs),
                                paste0(etiqueta, "\nrho no calculable"),
                                sprintf("%s\nrho = %.2f, p = %.2g", etiqueta,
                                        rho_unicos_vs_contigs, p_unicos_vs_contigs))) |>
    select(slug, etiqueta_a)
  df_u <- df_u |> left_join(rotulos, by = "slug")
  ga <- ggplot(df_u, aes(contigs, unicos)) +
    geom_point(aes(colour = atipico), size = 1.8) +
    geom_text_repel(data = filter(df_u, atipico), aes(label = id), size = 2.4, seed = 1) +
    facet_wrap(~ etiqueta_a, scales = "free") +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#B2182B"),
                        labels = c("normal", "atipico (> mediana + 3 MAD)"), name = NULL) +
    labs(x = "Contigs del ensamblaje (escala log)", y = "Familias unicas del genoma",
         tag = "a") +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))

  # (b) alpha: todos (con rango jackknife) vs solo completos
  df_b <- bind_rows(
    df_d |> transmute(etiqueta, grupo = "Todas las referencias", alpha = alpha_todos,
                      lo = alpha_jackknife_min, hi = alpha_jackknife_max,
                      n = n_referencias),
    df_d |> transmute(etiqueta, grupo = "Solo genomas completos", alpha = alpha_solo_completos,
                      lo = NA_real_, hi = NA_real_, n = n_completos)
  )
  gb <- ggplot(df_b, aes(alpha, etiqueta, colour = grupo)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, na.rm = TRUE,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2.5, na.rm = TRUE, position = position_dodge(width = 0.5)) +
    geom_text(aes(x = if_else(is.na(hi), alpha, hi), label = paste0("n=", n)), size = 2.6,
              hjust = -0.3, na.rm = TRUE,
              position = position_dodge(width = 0.5), show.legend = FALSE) +
    scale_colour_manual(values = c("Todas las referencias" = "#2166AC",
                                   "Solo genomas completos" = "#E08214"), name = NULL) +
    labs(x = "alpha de Heaps (< 1 abierto; barra = rango jackknife)", y = NULL,
         tag = "b") +
    theme(legend.position = "bottom")

  graficos <- list(ga, gb)
  if (length(hipo_tab) > 0) {
    df_h <- bind_rows(hipo_tab) |> left_join(etiquetas, by = "slug") |>
      mutate(categoria_95 = factor(categoria_95, levels = c("core", "shell", "cloud")))
    graficos[[3]] <- ggplot(df_h, aes(categoria_95, pct_hipoteticas, fill = categoria_95)) +
      geom_col(width = 0.6, colour = "grey30", linewidth = 0.2) +
      geom_text(aes(label = familias), vjust = -0.4, size = 2.6) +
      facet_wrap(~ etiqueta) + theme(strip.text = element_text(size = 8)) +
      scale_fill_manual(values = PAL_CATEGORIA, guide = "none") +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
      labs(x = NULL, y = "Familias hipoteticas (%)", tag = "c",
           caption = "Numero sobre cada barra: familias de la categoria.")
  }

  guardar_fig(wrap_plots(graficos, ncol = 1, heights = c(1.2, 0.8, 0.9)[seq_along(graficos)]),
              "fig_f3_diagnostico_heaps", ancho = 8, alto = 10)

  write_tsv(df_d |> select(linaje_id, especie_gtdb, n_referencias, alpha_todos,
                           alpha_jackknife_min, alpha_jackknife_max, n_completos,
                           alpha_solo_completos, rho_unicos_vs_contigs,
                           p_unicos_vs_contigs, n_atipicos, atipicos),
            file.path(DIR_TAB, "tab_f3_diagnostico_heaps.tsv"))
  write_tsv(df_u |> select(linaje_id, especie_gtdb, id, assembly_level, contigs, n50,
                           unicos, atipico),
            file.path(DIR_TAB, "tab_f3_genes_unicos_por_genoma.tsv"))
  if (length(hipo_tab) > 0)
    write_tsv(bind_rows(hipo_tab), file.path(DIR_TAB, "tab_f3_hipoteticas_por_categoria.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_diagnostico_heaps.tsv"))
}
