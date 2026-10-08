# =============================================================================
# f3_07_especificos.R — Analisis de los genes especificos de la cepa (Fase 3)
#
#  1. Catalogo: cada gen especifico con producto, KO/EC de Bakta, contig,
#     origen probable (mejor parecido en nr) e identidad.
#  2. Islas: genes especificos del mismo contig a <= DIST_ISLA pb entre si.
#     Para cada isla: largo, genes, marcadores de elemento movil a <= DIST_ISLA
#     pb y GC de la isla frente al GC del bin (si esta el .fna del bin).
#  3. Clase funcional por palabras clave del producto
#     (metadata/clases_funcionales.tsv; gana la primera que coincide).
#  4. Cruce con los marcadores de la Fase 4 (metadata/marcadores_fase4.tsv) por
#     numero KEGG (KO) de Bakta.
#
# Entradas por especie (results/fase3/<especie>/): exclusivos_verificados.tsv,
# genes_bin.tsv y, si se exportaron, bakta_bin/<bin>.tsv y bakta_bin/<bin>.fna.
# Acepta las clases nuevas (especifico, especifico_con_bandera) y las antiguas
# (verificado, con_bandera).
# =============================================================================

DIST_ISLA <- 5000
CLASES_ESP <- c("especifico", "especifico_con_bandera", "verificado", "con_bandera")
PATRON_MOVIL <- paste0("phage|prophage|capsid|terminase|integrase|recombinase|transposase|",
                       "IS[0-9]+|insertion element|excisionase|resolvase|relaxase|conjugal")

reglas <- read_tsv(file.path("metadata", "clases_funcionales.tsv"), show_col_types = FALSE) |>
  arrange(orden)
clase_funcional <- function(producto) {
  out <- rep("otra", length(producto))
  for (i in rev(seq_len(nrow(reglas)))) {            # la primera regla gana
    out[grepl(reglas$patron[i], producto, ignore.case = TRUE, perl = TRUE)] <- reglas$clase[i]
  }
  out
}

marcadores <- read_tsv(file.path("metadata", "marcadores_fase4.tsv"), show_col_types = FALSE) |>
  separate_rows(kos, sep = ",") |> mutate(kos = str_trim(kos)) |> filter(kos != "")

leer_fasta_simple <- function(p) {
  if (!file.exists(p)) return(NULL)
  l <- readLines(p, warn = FALSE)
  h <- grep("^>", l)
  nombres <- sub("^>(\\S+).*", "\\1", l[h])
  fin <- c(h[-1] - 1, length(l))
  setNames(vapply(seq_along(h), function(i) paste(l[(h[i] + 1):fin[i]], collapse = ""), ""), nombres)
}
gc_de <- function(x) { x <- toupper(x); n <- nchar(gsub("[^ACGT]", "", x))
  if (!n) NA_real_ else nchar(gsub("[^GC]", "", x)) / n }

catalogo <- list(); islas <- list()
for (s in especies_f3$slug) {
  ev <- leer_especie(s, "exclusivos_verificados.tsv", col_types = cols(.default = "c"))
  gb <- leer_especie(s, "genes_bin.tsv", col_types = cols(.default = "c"))
  if (is.null(ev) || is.null(gb)) next
  esp <- ev |> filter(clase %in% CLASES_ESP)
  if (!nrow(esp)) next
  b <- unique(gb$bin)
  tsv <- file.path(RES_F3, s, "bakta_bin", paste0(b, ".tsv"))
  kos <- if (file.exists(tsv)) {
    read_tsv(tsv, comment = "#", col_names = FALSE, show_col_types = FALSE,
             col_types = cols(.default = "c")) |>
      transmute(locus_tag = X6, dbx = coalesce(X9, ""),
                ko = map_chr(str_extract_all(dbx, "KEGG:K\\d{5}"), ~ paste(sub("KEGG:", "", .x), collapse = ",")),
                ec = map_chr(str_extract_all(dbx, "EC:[0-9.\\-]+"), ~ paste(sub("EC:", "", .x), collapse = ",")))
  } else tibble(locus_tag = character(), ko = character(), ec = character())
  seqs <- leer_fasta_simple(file.path(RES_F3, s, "bakta_bin", paste0(b, ".fna")))
  gc_bin <- if (is.null(seqs)) NA_real_ else gc_de(paste(seqs, collapse = ""))

  cat_s <- esp |>
    mutate(inicio = as.integer(inicio), fin = as.integer(fin)) |>
    left_join(kos |> select(locus_tag, ko, ec), by = "locus_tag") |>
    mutate(slug = s, clase_funcional = clase_funcional(producto),
           marcador_fase4 = map_chr(coalesce(ko, ""), function(k) {
             if (k == "") return("")
             m <- marcadores |> filter(kos %in% strsplit(k, ",")[[1]])
             paste(unique(paste0(m$funcion, " (", m$kos, ")")), collapse = "; ")
           })) |>
    arrange(contig, inicio)
  # islas: genes especificos consecutivos del mismo contig a <= DIST_ISLA pb
  cat_s <- cat_s |> group_by(contig) |>
    mutate(isla = cumsum(c(TRUE, diff(inicio) > 0 & (inicio - lag(fin))[-1] > DIST_ISLA))) |>
    ungroup() |> mutate(isla = paste0(sub("^L(\\d+).*", "L\\1", s), "_", contig, "_", isla))
  isl <- cat_s |> group_by(slug, isla, contig) |>
    summarise(inicio = min(inicio), fin = max(fin), genes = n(),
              con_bandera = sum(clase %in% c("especifico_con_bandera", "con_bandera")),
              clases = paste(names(sort(table(clase_funcional), decreasing = TRUE)), collapse = ","),
              .groups = "drop") |>
    rowwise() |>
    mutate(largo_kb = (fin - inicio + 1) / 1000,
           moviles_cerca = sum(gb$contig == contig &
                                 as.integer(gb$fin) >= inicio - DIST_ISLA &
                                 as.integer(gb$inicio) <= fin + DIST_ISLA &
                                 grepl(PATRON_MOVIL, gb$producto, ignore.case = TRUE)),
           gc_isla = if (is.null(seqs) || !contig %in% names(seqs)) NA_real_
                     else gc_de(substr(seqs[[contig]], inicio, fin)),
           gc_bin = gc_bin, delta_gc = gc_isla - gc_bin) |>
    ungroup()
  catalogo[[s]] <- cat_s; islas[[s]] <- isl
}
catalogo <- bind_rows(catalogo); islas <- bind_rows(islas)

if (nrow(catalogo)) {
  write_tsv(catalogo |> select(slug, isla, locus_tag, contig, inicio, fin, largo_aa, gen, producto,
                               clase_funcional, ko, ec, marcador_fase4, clase, motivo,
                               origen_candidato = any_of("origen_candidato"), f5_origen, f5_pident,
                               f5_taxon),
            file.path(DIR_TAB, "tab_f3_especificos_catalogo.tsv"))
  write_tsv(islas, file.path(DIR_TAB, "tab_f3_especificos_islas.tsv"))
  message("  -> ", file.path(DIR_TAB, "tab_f3_especificos_{catalogo,islas}.tsv"))

  et <- especies_f3 |> select(slug, etiqueta)
  PAL_CLASE <- c(movil_fago = "#7570B3", defensa = "#D95F02", resistencia_estres = "#E7298A",
                 transporte = "#1B9E77", regulacion = "#66A61E", metabolismo = "#E6AB02",
                 hipotetica = "#BDBDBD", otra = "#666666")
  ga <- catalogo |> count(slug, clase_funcional) |> left_join(et, by = "slug") |>
    mutate(clase_funcional = factor(clase_funcional, levels = names(PAL_CLASE))) |>
    ggplot(aes(n, etiqueta, fill = clase_funcional)) +
    geom_col(width = 0.6, colour = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = PAL_CLASE, name = "Clase funcional", drop = TRUE) +
    labs(x = "Genes especificos de la cepa", y = NULL, tag = "a",
         caption = "Clase asignada por palabras clave del producto de Bakta (metadata/clases_funcionales.tsv).")
  gb_ <- islas |> left_join(et, by = "slug") |>
    mutate(movil = if_else(moviles_cerca > 0, "con gen movil cerca", "sin gen movil cerca"))
  gb2 <- if (all(is.na(gb_$delta_gc))) {
    ggplot(gb_, aes(largo_kb, genes, colour = movil)) +
      geom_point(size = 2.5) + facet_wrap(~ etiqueta, scales = "free") +
      labs(x = "Largo de la isla (kb)", y = "Genes especificos en la isla", tag = "b",
           caption = "GC no disponible: falta exportar el .fna de los bins.")
  } else {
    ggplot(gb_, aes(largo_kb, 100 * delta_gc, colour = movil, size = genes)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.8) + facet_wrap(~ etiqueta, scales = "free_x") +
      scale_size_continuous(range = c(1.5, 5), name = "Genes") +
      labs(x = "Largo de la isla (kb)", y = "GC de la isla - GC del bin (puntos %)", tag = "b",
           caption = "Islas: genes especificos del mismo contig a <= 5 kb. Gen movil cerca: a <= 5 kb.")
  }
  gb2 <- gb2 + scale_colour_manual(values = c("con gen movil cerca" = "#B2182B",
                                              "sin gen movil cerca" = "#2166AC"), name = NULL) +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
  guardar_fig(ga / gb2 + plot_layout(heights = c(1, 1.4)), "fig_f3_especificos",
              ancho = 9, alto = 8)
}
