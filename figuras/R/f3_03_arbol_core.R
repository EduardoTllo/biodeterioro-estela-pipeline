# =============================================================================
# f3_03_arbol_core.R — Arbol de maxima verosimilitud del core por especie
#
# Entrada: core.treefile de IQ-TREE 3 (06b_iqtree.slurm). Enraizado en el punto
# medio solo para visualizar (phangorn::midpoint). Las puntas se colorean por
# el habitat de la referencia y el bin de la Estela se destaca. Los nodos con
# UFBoot >= 95 llevan un punto negro.
# =============================================================================

for (slug in especies_f3$slug) {
  p_arbol <- file.path(RES_F3, slug, "core.treefile")
  if (!file.exists(p_arbol)) { message("  sin core.treefile para ", slug); next }
  arbol <- phangorn::midpoint(read.tree(p_arbol))
  esp <- especies_f3 |> filter(slug == !!slug)
  bins_ids <- str_split(esp$bins_ids, ",")[[1]]
  bins_orig <- str_split(esp$bins, ",")[[1]]

  meta <- tibble(label = arbol$tip.label) |>
    left_join(refs_f3 |> filter(slug == !!slug) |>
                select(label = id, habitat, continente), by = "label") |>
    mutate(es_bin = label %in% bins_ids,
           habitat = if_else(es_bin, "Estela de Raimondi", coalesce(habitat, "desconocido")),
           etiqueta = if_else(es_bin, bins_orig[match(label, bins_ids)], NA_character_))

  soporte <- suppressWarnings(as.numeric(arbol$node.label))
  n_tips <- length(arbol$tip.label)
  nodos_ok <- tibble(node = n_tips + seq_along(soporte), ufboot = soporte) |>
    filter(!is.na(ufboot), ufboot >= 95)

  profundidad <- max(node.depth.edgelength(arbol))

  g <- ggtree(arbol) %<+% meta +
    geom_point2(aes(subset = node %in% nodos_ok$node), size = 0.9, colour = "black") +
    geom_tippoint(aes(colour = habitat, size = es_bin)) +
    geom_tiplab(aes(label = etiqueta), size = 3, fontface = "bold", colour = "#B2182B",
                offset = 0.02 * profundidad, na.rm = TRUE) +
    scale_colour_manual(values = PAL_HABITAT, name = "Habitat") +
    scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3), guide = "none") +
    geom_treescale(fontsize = 2.5, linesize = 0.3) +
    hexpand(0.15) +
    labs(title = esp$etiqueta,
         caption = "Maxima verosimilitud (IQ-TREE 3) sobre el core; punto negro: UFBoot >= 95.") +
    theme(legend.position = "right")
  guardar_fig(g, paste0("fig_f3_arbol_core_", slug), ancho = 7,
              alto = max(4, 0.12 * n_tips + 1.5))
}
