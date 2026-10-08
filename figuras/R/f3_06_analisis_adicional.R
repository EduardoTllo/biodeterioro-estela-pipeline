# =============================================================================
# f3_06_analisis_adicional.R — Figuras complementarias del pangenoma (Fase 3)
#
#  fig_f3_espectro_frecuencias  Espectro de frecuencias de las familias (en
#                               cuantas referencias esta cada una): la forma en
#                               U tipica de los pangenomas bacterianos.
#  fig_f3_unicos_vecino         Familias propias de cada genoma frente a la
#                               distancia a su pariente mas cercano en el arbol
#                               del core; el bin se compara con las referencias.
#  fig_f3_pcoa_accesorio        Ordenacion (PCoA, Jaccard) del genoma accesorio,
#                               coloreada por habitat; Mantel contra la
#                               filogenia y PERMANOVA por habitat.
#  fig_f3_matriz_<especie>      Presencia/ausencia del genoma accesorio junto al
#                               arbol del core.
#  fig_f3_contexto_exclusivos   Mapa de genes alrededor de los exclusivos.
#
# El accesorio de las figuras 2-4 excluye las familias de una sola referencia
# (singletons): por definicion no aportan parecido entre genomas.
# =============================================================================

suppressPackageStartupMessages({ library(vegan); library(gggenes) })

leer_matriz <- function(slug) {
  m <- leer_especie(slug, "pangenoma_refs.Rtab")
  x <- as.matrix(m[, -1]); rownames(x) <- m$Gene
  x[rowSums(x) > 0, , drop = FALSE]
}

# Matriz genomas x familias accesorias (2..N-1 referencias) con el bin agregado.
matriz_accesoria <- function(slug) {
  x <- leer_matriz(slug); N <- ncol(x); f <- rowSums(x)
  acc <- x[f > 1 & f < N, , drop = FALSE]
  gb <- leer_especie(slug, "genes_bin.tsv")
  bin <- unique(gb$bin)
  acc <- cbind(acc, setNames(data.frame(as.integer(rownames(acc) %in% gb$familia)), bin))
  t(as.matrix(acc))
}

habitat_de <- function(ids, slug) {
  h <- refs_f3 |> filter(slug == !!slug) |> select(id, habitat) |> deframe()
  out <- unname(h[ids]); out[is.na(out)] <- "desconocido"
  out[grepl("^bin_", ids)] <- "Estela de Raimondi"
  out
}

# --- 1. Espectro de frecuencias ---------------------------------------------
espectro <- map_dfr(especies_f3$slug, function(s) {
  x <- leer_matriz(s); N <- ncol(x)
  tibble(slug = s, n = rowSums(x), N = N) |> count(slug, N, n) |>
    mutate(frac = n / N,
           categoria = case_when(frac >= UMBRAL_CORE ~ "core",
                                 frac >= UMBRAL_CLOUD ~ "shell", TRUE ~ "cloud"))
}) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")

g1 <- ggplot(espectro, aes(n, nn, fill = categoria)) +
  geom_col(width = 0.9, colour = NA) +
  facet_wrap(~ etiqueta, scales = "free", labeller = label_wrap_gen(width = 22)) +
  scale_fill_manual(values = PAL_CATEGORIA[c("core", "shell", "cloud")],
                    breaks = c("core", "shell", "cloud"), name = NULL) +
  scale_y_sqrt(labels = scales::comma) +
  labs(x = "Numero de referencias que tienen la familia",
       y = "Familias de genes (escala raiz cuadrada)",
       caption = "Cada barra cuenta las familias presentes en exactamente n referencias.") +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
guardar_fig(g1, "fig_f3_espectro_frecuencias", ancho = 8, alto = 3.8)

# --- 2. Familias propias vs distancia al pariente mas cercano ---------------
unicos <- map_dfr(especies_f3$slug, function(s) {
  x <- leer_matriz(s)
  arbol <- read.tree(file.path(RES_F3, s, "core.treefile"))
  d <- cophenetic(arbol); diag(d) <- Inf
  u <- colSums(x[rowSums(x) == 1, , drop = FALSE])
  gb <- leer_especie(s, "genes_bin.tsv"); b <- unique(gb$bin)
  bind_rows(tibble(id = names(u), unicos = as.numeric(u)),
            tibble(id = b, unicos = sum(gb$exclusivo_candidato == "si"))) |>
    mutate(slug = s, dist = apply(d[id, , drop = FALSE], 1, min), es_bin = id == b)
})
rotulos2 <- unicos |> filter(!es_bin) |> group_by(slug) |>
  summarise(rho = cor(dist, unicos, method = "spearman"),
            p = cor.test(dist, unicos, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop") |>
  left_join(unicos |> filter(es_bin) |> select(slug, ub = unicos), by = "slug") |>
  left_join(unicos |> filter(!es_bin) |> group_by(slug) |>
              summarise(med = median(unicos), .groups = "drop"), by = "slug") |>
  left_join(especies_f3 |> select(slug, etiqueta), by = "slug") |>
  mutate(rotulo = sprintf("%s\nrho = %.2f, p = %.2g", etiqueta, rho, p))
unicos <- unicos |> left_join(rotulos2 |> select(slug, rotulo), by = "slug")

g2 <- ggplot(unicos, aes(dist, unicos)) +
  geom_point(data = filter(unicos, !es_bin), colour = "grey45", size = 1.6) +
  geom_point(data = filter(unicos, es_bin), colour = "#B2182B", size = 3) +
  geom_text_repel(data = filter(unicos, es_bin), aes(label = sub("_(\\d+)_(\\d+)$", "-\\1-\\2", id)),
                  colour = "#B2182B", size = 3, seed = 1) +
  facet_wrap(~ rotulo, scales = "free") +
  scale_x_log10() +
  labs(x = "Distancia al genoma mas cercano en el arbol del core (sustituciones/sitio, log)",
       y = "Familias de genes propias del genoma",
       caption = paste("Gris: referencias (familias que no tiene ninguna otra referencia).",
                       "Rojo: bin de la Estela (exclusivos candidatos, antes de F1-F6).")) +
  theme(strip.text = element_text(size = 8))
guardar_fig(g2, "fig_f3_unicos_vecino", ancho = 9, alto = 3.8)

# --- 3. PCoA del accesorio ---------------------------------------------------
ordenes <- list(); estad <- list()
for (s in especies_f3$slug) {
  a <- matriz_accesoria(s)
  dj <- vegdist(a, method = "jaccard", binary = TRUE)
  pc <- cmdscale(dj, k = 2, eig = TRUE)
  expl <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
  refs <- !grepl("^bin_", rownames(a))
  arbol <- read.tree(file.path(RES_F3, s, "core.treefile"))
  dp <- cophenetic(arbol)[rownames(a)[refs], rownames(a)[refs]]
  djr <- vegdist(a[refs, ], method = "jaccard", binary = TRUE)
  mt <- mantel(as.dist(dp), djr, method = "spearman", permutations = 999)
  hab <- habitat_de(rownames(a)[refs], s)
  ok <- hab %in% names(which(table(hab) >= 3)) & !hab %in% c("desconocido", "sin_clasificar")
  pm <- adonis2(vegdist(a[refs, ][ok, ], "jaccard", binary = TRUE) ~ hab[ok], permutations = 999)
  estad[[s]] <- tibble(slug = s, mantel_r = mt$statistic, mantel_p = mt$signif,
                       permanova_r2 = pm$R2[1], permanova_p = pm$`Pr(>F)`[1],
                       n_permanova = sum(ok), grupos = length(unique(hab[ok])))
  ordenes[[s]] <- tibble(id = rownames(a), ej1 = pc$points[, 1], ej2 = pc$points[, 2],
                         habitat = habitat_de(rownames(a), s), slug = s,
                         ex1 = expl[1], ex2 = expl[2])
}
estad <- bind_rows(estad) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")
paneles <- map(especies_f3$slug, function(s) {
  d <- ordenes[[s]]; e <- filter(estad, slug == s)
  ggplot(d, aes(ej1, ej2)) +
    geom_point(data = filter(d, habitat != "Estela de Raimondi"), aes(fill = habitat),
               shape = 21, size = 2.2, colour = "grey30", stroke = 0.2) +
    geom_point(data = filter(d, habitat == "Estela de Raimondi"), shape = 23, size = 3.6,
               fill = "#B2182B", colour = "black") +
    scale_fill_manual(values = PAL_HABITAT, name = "Habitat", drop = FALSE,
                      limits = setdiff(names(PAL_HABITAT), "Estela de Raimondi")) +
    labs(title = e$etiqueta,
         subtitle = sprintf("Mantel vs filogenia r = %.2f (p = %.3f)\nPERMANOVA habitat R2 = %.2f (p = %.2f, n = %d)",
                            e$mantel_r, e$mantel_p, e$permanova_r2, e$permanova_p, e$n_permanova),
         x = sprintf("PCoA 1 (%.0f %%)", d$ex1[1]), y = sprintf("PCoA 2 (%.0f %%)", d$ex2[1])) +
    theme(plot.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 7))
})
paneles[-1] <- map(paneles[-1], ~ .x + guides(fill = "none"))
g3 <- wrap_plots(paneles, nrow = 1, guides = "collect") +
  plot_annotation(caption = paste("Rombo rojo: bin de la Estela. PCoA sobre distancias de Jaccard del accesorio",
                                  "(familias en 2 a N-1 referencias). PERMANOVA solo con habitats de >= 3 referencias.")) &
  theme(legend.position = "bottom")
guardar_fig(g3, "fig_f3_pcoa_accesorio", ancho = 11, alto = 4.8)
write_tsv(estad, file.path(DIR_TAB, "tab_f3_accesorio_filogenia_habitat.tsv"))
message("  -> ", file.path(DIR_TAB, "tab_f3_accesorio_filogenia_habitat.tsv"))

# --- 4. Matriz de presencia/ausencia junto al arbol ------------------------
# gheatmap() no se puede guardar con ggplot2 4 + ggtree 3.16: el mapa de calor
# se dibuja con ggplot y se alinea con el arbol por las etiquetas (aplot).
for (s in especies_f3$slug) {
  a <- matriz_accesoria(s)
  arbol <- phangorn::midpoint(read.tree(file.path(RES_F3, s, "core.treefile")))
  meta <- tibble(label = arbol$tip.label, habitat = habitat_de(arbol$tip.label, s),
                 es_bin = grepl("^bin_", arbol$tip.label))
  et <- especies_f3$etiqueta[especies_f3$slug == s]
  pa <- ggtree(arbol) %<+% meta +
    geom_tippoint(aes(colour = habitat, size = es_bin)) +
    scale_size_manual(values = c(`FALSE` = 1.2, `TRUE` = 3), guide = "none") +
    scale_colour_manual(values = PAL_HABITAT, name = "Habitat")
  orden <- order(-colSums(a))
  largo <- as_tibble(a[, orden], rownames = "label") |>
    pivot_longer(-label, names_to = "familia", values_to = "v") |>
    filter(v == 1) |> mutate(x = match(familia, colnames(a)[orden]))
  ph <- ggplot(largo, aes(x, label)) +
    geom_tile(fill = "#2166AC", width = 1, height = 0.85) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = sprintf("%d familias accesorias (2 a N-1 referencias), de mas a menos frecuentes", ncol(a)),
         y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_blank())
  comb <- aplot::as.patchwork(aplot::insert_left(ph, pa, width = 0.45)) +
    plot_annotation(title = et, caption = "Punta grande: bin de la Estela. Azul: familia presente.")
  guardar_fig(comb, paste0("fig_f3_matriz_", s), ancho = 10, alto = 6)
}

# --- 5. Contexto genomico de los exclusivos ---------------------------------
f5b <- {
  p <- file.path(RES_F3, "phase3_f5b_especie_por_ani.tsv")
  if (file.exists(p)) read_tsv(p, show_col_types = FALSE) else tibble(locus_tag = character(), resultado_f5b = character())
}
ventanas <- map_dfr(especies_f3$slug, function(s) {
  gb <- leer_especie(s, "genes_bin.tsv")
  ev <- leer_especie(s, "exclusivos_verificados.tsv", col_types = cols(.default = "c"))
  if (is.null(ev)) return(NULL)
  ani <- c(ev |> filter(grepl("F5b", motivo)) |> pull(locus_tag),
           f5b |> filter(resultado_f5b == "presente_en_la_especie") |> pull(locus_tag))
  ver <- ev |> filter(clase != "descartado" | locus_tag %in% ani) |> pull(locus_tag)
  if (!length(ver)) return(NULL)
  ctg <- gb |> filter(locus_tag %in% ver) |> group_by(contig) |>
    summarise(ini = max(0, min(inicio) - 4000), fin = max(fin) + 4000, .groups = "drop")
  gb |> inner_join(ctg, by = "contig") |> filter(fin.x >= ini, inicio <= fin.y) |>
    rename(fin = fin.x) |>
    mutate(slug = s,
           clase = case_when(locus_tag %in% ani ~ "en otra cepa de la especie (F5b, ANI)",
                             locus_tag %in% ver ~ "especifico de la cepa",
                             exclusivo_candidato == "si" ~ "candidato descartado (F1-F5)",
                             estado_panaroo != "en_familia" ~ "fuera del pangenoma",
                             TRUE ~ categoria_95),
           movil = grepl("integrase|recombinase|transposase|phage|terminase|capsid|prohead|IS[0-9]",
                         producto, ignore.case = TRUE),
           rotulo = if_else(clase %in% c("especifico de la cepa", "en otra cepa de la especie (F5b, ANI)"),
                            str_trunc(str_squish(gsub("domain-containing|family|protein|-type", "", producto)), 18), ""),
           molecula = paste0(str_extract(slug, "^L\\d+"), " ", sub("_(\\d+)_(\\d+)$", "-\\1-\\2", bin), " ", contig),
           forward = hebra == "+")
})
if (nrow(ventanas)) {
  PAL_CTX <- c("core" = "#2166AC", "shell" = "#67A9CF", "cloud" = "#D1E5F0",
               "especifico de la cepa" = "#B2182B",
               "en otra cepa de la especie (F5b, ANI)" = "#F4A582",
               "candidato descartado (F1-F5)" = "#FDDBC7",
               "fuera del pangenoma" = "#BDBDBD")
  g5 <- ggplot(ventanas, aes(xmin = inicio, xmax = fin, y = molecula, fill = clase, forward = forward)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1.2, "mm"),
                    arrow_body_height = unit(2.4, "mm"), colour = "grey30", linewidth = 0.15) +
    geom_point(data = filter(ventanas, movil), aes(x = (inicio + fin) / 2, y = molecula),
               inherit.aes = FALSE, position = position_nudge(y = -0.22), shape = 17, size = 1.1) +
    geom_text_repel(aes(x = (inicio + fin) / 2, label = rotulo), size = 1.8, nudge_y = 0.35,
                    direction = "both", max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.15, box.padding = 0.15, seed = 1) +
    facet_wrap(~ molecula, scales = "free", ncol = 1) +
    scale_fill_manual(values = PAL_CTX, name = NULL) +
    theme_genes() +
    labs(x = "Posicion en el contig (pb)", y = NULL,
         caption = "Rotulos: genes especificos y los descartados por F5b. Triangulo negro: gen movil (integrasa, recombinasa, transposasa, fago, IS).") +
    theme(legend.position = "bottom", strip.text = element_blank(),
          axis.text.y = element_text(size = 7))
  guardar_fig(g5, "fig_f3_contexto_exclusivos", ancho = 10, alto = 1.15 * n_distinct(ventanas$molecula) + 1.2)
}

# --- 6. Comparacion a igual numero de genomas --------------------------------
# El pangenoma crece con el numero de genomas: para comparar especies se usan
# 23 genomas en todas (el N de A. schindleri). (a) curvas promedio con 23
# genomas elegidos al azar; (b) alpha de Heaps en 30 submuestras de 23.
N_COMUN <- 23
set.seed(23)
curvas23 <- list(); alfas23 <- list()
for (s in especies_f3$slug) {
  pm <- t(leer_matriz(s))                       # genomas x familias
  n <- nrow(pm)
  cur <- map_dfr(seq_len(N_PERM_CURVA), function(i) {
    sub <- pm[sample(n, N_COMUN), , drop = FALSE]
    sub <- sub[, colSums(sub) > 0, drop = FALSE]
    cc <- curva_acumulacion(sub)
    tibble(k = seq_len(N_COMUN), pan = cc$pan, core = cc$core)
  })
  curvas23[[s]] <- cur |> group_by(k) |>
    summarise(across(c(pan, core), list(med = median, lo = ~ quantile(.x, 0.25),
                                        hi = ~ quantile(.x, 0.75))), .groups = "drop") |>
    mutate(slug = s)
  reps <- if (n > N_COMUN) 30 else 1
  alfas23[[s]] <- tibble(slug = s, alpha = replicate(reps, {
    sub <- pm[sample(n, N_COMUN), , drop = FALSE]
    heaps_rapido(sub[, colSums(sub) > 0, drop = FALSE], n.perm = 100)["alpha"]
  }))
}
c23 <- bind_rows(curvas23) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")
a23 <- bind_rows(alfas23) |> left_join(especies_f3 |> select(slug, etiqueta), by = "slug")
PAL_ESP <- setNames(c("#1B9E77", "#D95F02", "#7570B3"), especies_f3$etiqueta)
g6a <- ggplot(c23, aes(k, colour = etiqueta, fill = etiqueta)) +
  geom_ribbon(aes(ymin = pan_lo, ymax = pan_hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = pan_med), linewidth = 0.9) +
  geom_ribbon(aes(ymin = core_lo, ymax = core_hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = core_med), linewidth = 0.9, linetype = "dashed") +
  scale_colour_manual(values = PAL_ESP, name = NULL) +
  scale_fill_manual(values = PAL_ESP, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Genomas de referencia (submuestras de 23)", y = "Familias de genes", tag = "a",
       caption = "Linea continua: pangenoma; discontinua: core (en todos los genomas). Mediana y rango intercuartil.") +
  theme(legend.position = "bottom")
g6b <- ggplot(a23, aes(alpha, etiqueta, colour = etiqueta)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(height = 0.12, width = 0, size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = PAL_ESP, guide = "none") +
  labs(x = "alpha de Heaps con 23 genomas (< 1 = abierto)", y = NULL, tag = "b",
       caption = "30 submuestras de 23 genomas (A. schindleri tiene exactamente 23: un solo valor).")
guardar_fig(g6a / g6b + plot_layout(heights = c(2, 1)), "fig_f3_comparacion_igual_n",
            ancho = 8, alto = 7.5)
write_tsv(a23 |> group_by(slug) |> summarise(alpha_mediana = median(alpha),
                                             alpha_min = min(alpha), alpha_max = max(alpha)),
          file.path(DIR_TAB, "tab_f3_alpha_igual_n.tsv"))
