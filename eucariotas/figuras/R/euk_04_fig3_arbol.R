# =============================================================================
# euk_04_fig3_arbol.R — Fig. E3: arbol filogenomico de los 13 genomas
#
# Igual que la Fig. 7.4 (bacterias): distancias-p y Neighbor-Joining con
# bootstrap no parametrico, pero sobre la matriz concatenada de genes BUSCO de
# copia unica (fungi_odb12.2; 1 019 genes, ~511 000 posiciones) en lugar del
# user_msa bac120. Enraizado en la unica Basidiomycota (bin-2-65).
#
# Con ~511 000 columnas, ape::boot.phylo seria muy lento. Se precalcula, para
# cada par de genomas y cada columna, si la posicion es valida y si difiere; un
# replica bootstrap es entonces un vector de pesos (conteos multinomiales de
# columnas) y la distancia-p replicada es sum(w * difiere) / sum(w * valida).
# Es el mismo bootstrap por columnas, calculado con algebra matricial.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
})

set.seed(42)
N_BOOT <- 500

leer_fasta_chr <- function(path) {
  l <- readLines(path)
  cab <- grep("^>", l)
  nombres <- sub("^>(\\S+).*", "\\1", l[cab])
  fin <- c(cab[-1] - 1, length(l))
  seqs <- mapply(function(a, b) paste(l[(a + 1):b], collapse = ""), cab, fin)
  m <- do.call(rbind, strsplit(seqs, ""))
  rownames(m) <- nombres
  m
}

mat <- leer_fasta_chr(file.path(RES_E2, "phase2_euk_busco_concat.faa"))
message("Matriz BUSCO: ", nrow(mat), " genomas x ", ncol(mat), " posiciones")

valido <- !(mat %in% c("-", "X", "*", "?", "."))
dim(valido) <- dim(mat)
pares <- combn(nrow(mat), 2)
# double desde el inicio: crossprod() convertiria enteros en cada replica
dif <- matrix(0, ncol(mat), ncol(pares))
val <- matrix(0, ncol(mat), ncol(pares))
for (k in seq_len(ncol(pares))) {
  i <- pares[1, k]; j <- pares[2, k]
  ok <- valido[i, ] & valido[j, ]
  val[, k] <- ok
  dif[, k] <- ok & (mat[i, ] != mat[j, ])
}
rm(valido)

dist_desde_pesos <- function(w) {
  d <- as.vector(crossprod(w, dif)) / as.vector(crossprod(w, val))
  m <- matrix(0, nrow(mat), nrow(mat), dimnames = list(rownames(mat), rownames(mat)))
  m[t(pares)] <- d
  m[t(pares[2:1, ])] <- d
  as.dist(m)
}

outgroup <- maestra$genome[maestra$filo == "Basidiomycota"]
nj_enraizado <- function(d) {
  t <- nj(d)
  root(t, outgroup = outgroup, resolve.root = TRUE)
}

arbol <- nj_enraizado(dist_desde_pesos(rep(1, ncol(mat))))

message("Bootstrap (", N_BOOT, " replicas)...")
n <- ncol(mat)
reps <- lapply(seq_len(N_BOOT), function(b) {
  w <- as.numeric(tabulate(sample.int(n, n, replace = TRUE), nbins = n))
  nj_enraizado(dist_desde_pesos(w))
})
class(reps) <- "multiPhylo"
soporte <- prop.clades(arbol, reps, rooted = TRUE)
arbol$node.label <- ifelse(is.na(soporte), "", round(100 * soporte / N_BOOT))
arbol$node.label[1] <- ""
arbol <- ladderize(arbol)
arbol$edge.length[arbol$edge.length < 0] <- 0

anot <- maestra %>%
  select(label = genome, clase, taxon_corto, muestra, linaje_id, representante,
         Completitud = eukcc_completitud, Contaminacion = eukcc_contaminacion)

p_arbol <- ggtree(arbol, linewidth = 0.45) %<+% anot +
  geom_tippoint(aes(colour = clase, shape = clase), size = 2.3) +
  geom_tiplab(aes(label = paste0(label, "  ", taxon_corto), colour = clase),
              size = 2.8, offset = 0.004, show.legend = FALSE) +
  geom_text2(aes(subset = !isTip & label != "" & suppressWarnings(as.numeric(label)) >= 70,
                 label = label),
             hjust = 1.25, vjust = -0.45, size = 2.2, colour = "grey35") +
  geom_treescale(width = 0.05, fontsize = 2.5, linesize = 0.4, offset = 0.4) +
  scale_colour_manual(values = PAL_CLASE, name = "Clase", drop = TRUE) +
  scale_shape_manual(values = SHP_CLASE, name = "Clase", drop = TRUE) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
p_arbol <- p_arbol + xlim(0, max(p_arbol$data$x, na.rm = TRUE) * 1.6)

ylim_arbol <- ggplot_build(p_arbol)$layout$panel_params[[1]]$y.range
anot_ord <- p_arbol$data %>% filter(isTip) %>% select(label, y) %>%
  left_join(anot, by = "label")

tema_col <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
                  plot.margin = margin(5.5, 2, 5.5, 2))

p_calidad <- anot_ord %>%
  pivot_longer(c(Completitud, Contaminacion), names_to = "m", values_to = "v") %>%
  ggplot(aes(x = v, y = y, fill = clase)) +
  geom_col(width = 0.62, orientation = "y", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f", v)), hjust = -0.12, size = 2.1, colour = "grey25") +
  facet_wrap(~ m, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = PAL_CLASE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)), breaks = breaks_pretty(3)) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "EukCC (%)", y = NULL) + tema_col

p_muestra <- ggplot(anot_ord, aes(x = 1, y = y, fill = clase)) +
  geom_tile(height = 0.78, colour = "white", linewidth = 0.6, alpha = 0.45, show.legend = FALSE) +
  geom_text(aes(label = muestra), size = 2.4, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASE) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "Muestra", y = NULL) + tema_col + theme(panel.grid = element_blank())

p_linaje <- ggplot(anot_ord, aes(x = 1, y = y)) +
  geom_tile(fill = "grey93", height = 0.78, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = paste0(linaje_id, if_else(representante, "*", ""))), size = 2.4) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "Linaje", y = NULL) + tema_col + theme(panel.grid = element_blank())

fig3 <- p_arbol + p_muestra + p_linaje + p_calidad +
  plot_layout(widths = c(3.6, 0.28, 0.34, 1.5), nrow = 1)

guardar_fig(fig3, "figE3_arbol_filogenomico", ancho = 12, alto = 5.6)
write.tree(arbol, file.path(DIR_FIG, "figE3_arbol_busco_nj.newick"))
