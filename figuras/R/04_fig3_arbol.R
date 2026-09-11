# =============================================================================
# 04_fig3_arbol.R — Figura 3: arbol filogenomico anotado de los 21 genomas
#
# El arbol se construye a partir del alineamiento concatenado de los 120
# marcadores bacterianos (bac120) que produce GTDB-Tk (`user_msa`), mediante
# distancias-p y Neighbor-Joining, con soporte por bootstrap no parametrico.
# NOTA: es un arbol de distancias, no de maxima verosimilitud. El arbol de novo
# focalizado (ML) esta previsto para la Fase 4.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
})
# Biostrings se usa con espacio de nombres explicito: al adjuntarlo enmascara
# dplyr::rename / dplyr::first y rompe los scripts posteriores.
stopifnot(requireNamespace("Biostrings", quietly = TRUE))

set.seed(42)
N_BOOT <- 500

# --- Alineamiento bac120 ------------------------------------------------------
msa_path <- file.path(RES_F2, "phase2_gtdbtk_bac120_user_msa.fasta.gz")
aln <- Biostrings::readAAStringSet(msa_path)
names(aln) <- str_remove(names(aln), "\\s.*$")
mat <- as.matrix(aln)

message("Alineamiento bac120: ", nrow(mat), " genomas x ", ncol(mat), " posiciones")

# --- Distancia-p entre pares (ignorando gaps y posiciones ambiguas) -----------
dist_p <- function(m) {
  n <- nrow(m); nm <- rownames(m)
  d <- matrix(0, n, n, dimnames = list(nm, nm))
  valido <- !(m %in% c("-", "X", "*", ".")) ; dim(valido) <- dim(m)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      ok <- valido[i, ] & valido[j, ]
      d[i, j] <- d[j, i] <- if (sum(ok) == 0) NA_real_
                            else mean(m[i, ok] != m[j, ok])
    }
  }
  as.dist(d)
}

nj_desde_matriz <- function(m) ape::nj(dist_p(m))

arbol <- nj_desde_matriz(mat)

# --- Enraizado en Pseudomonadota (grupo externo respecto a Terrabacteria) ------
outgroup <- maestra$genoma[maestra$filo == "Pseudomonadota"]
outgroup <- intersect(outgroup, arbol$tip.label)
arbol <- root(arbol, outgroup = outgroup, resolve.root = TRUE)
arbol <- ladderize(arbol)

# --- Soporte por bootstrap ----------------------------------------------------
message("Calculando bootstrap (", N_BOOT, " replicas)...")
bs <- boot.phylo(arbol, mat, FUN = nj_desde_matriz, B = N_BOOT,
                 rooted = TRUE, trees = FALSE, quiet = TRUE)
arbol$node.label <- c("", round(100 * bs[-1] / N_BOOT))

# NJ puede producir ramas de longitud negativa; se fijan a 0 solo para el
# dibujo (no altera la topologia ni los valores de soporte).
arbol$edge.length[arbol$edge.length < 0] <- 0

# --- Panel del arbol ----------------------------------------------------------
anot <- maestra %>%
  select(label = genoma, filo, taxon_corto, muestra, linaje_id,
         Completeness, Contamination, ani, representante)

p_arbol <- ggtree(arbol, linewidth = 0.45) %<+% anot +
  geom_tippoint(aes(colour = filo), size = 2.1) +
  geom_tiplab(aes(label = paste0(label, "  ", taxon_corto), colour = filo),
              size = 2.7, offset = 0.004, fontface = "plain",
              show.legend = FALSE) +
  geom_text2(aes(subset = !isTip & label != "" & as.numeric(label) >= 70,
                 label = label),
             hjust = 1.25, vjust = -0.45, size = 2.1, colour = "grey35") +
  geom_treescale(width = 0.05, fontsize = 2.5, linesize = 0.4, offset = 0.4) +
  scale_colour_manual(values = PAL_FILO, name = "Filo") +
  theme(legend.position = c(0.06, 0.90), legend.justification = c(0, 1))

# Espacio a la derecha para que quepan las etiquetas
p_arbol <- p_arbol + xlim(0, max(p_arbol$data$x, na.rm = TRUE) * 1.75)

# --- Paneles de anotacion alineados a la coordenada Y del arbol ---------------
# Se usa la coordenada Y continua de cada hoja (no un factor) y los mismos
# limites del panel del arbol: asi cada fila queda exactamente a la altura de
# su rama, incluido el espacio que ocupa la barra de escala.
ylim_arbol <- ggplot_build(p_arbol)$layout$panel_params[[1]]$y.range

anot_ord <- p_arbol$data %>%
  filter(isTip) %>%
  select(label, y) %>%
  left_join(anot, by = "label")

tema_col <- theme(
  axis.title.y = element_blank(), axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
  plot.margin = margin(5.5, 2, 5.5, 2)
)

p_calidad <- anot_ord %>%
  select(y, filo, Completitud = Completeness, Contaminacion = Contamination) %>%
  pivot_longer(c(Completitud, Contaminacion), names_to = "m", values_to = "v") %>%
  mutate(m = factor(m, levels = c("Completitud", "Contaminacion"))) %>%
  ggplot(aes(x = v, y = y, fill = filo)) +
  geom_col(width = 0.62, orientation = "y", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f", v)), hjust = -0.12, size = 2.1,
            colour = "grey25") +
  facet_wrap(~ m, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = PAL_FILO) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30)),
                     breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "%", y = NULL) + tema_col

p_muestra <- anot_ord %>%
  ggplot(aes(x = 1, y = y, fill = filo)) +
  geom_tile(height = 0.78, colour = "white", linewidth = 0.6, alpha = 0.55,
            show.legend = FALSE) +
  geom_text(aes(label = muestra), size = 2.4, fontface = "bold") +
  scale_fill_manual(values = PAL_FILO) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "Muestra", y = NULL) + tema_col +
  theme(panel.grid = element_blank())

p_linaje <- anot_ord %>%
  ggplot(aes(x = 1, y = y)) +
  # Sin resaltado de linajes: la seleccion para la Fase 3 ya no se decide por
  # prevalencia sino por disponibilidad genomica publica (censo de la Fase 3).
  geom_tile(fill = "grey93", height = 0.78,
            colour = "white", linewidth = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(linaje_id, if_else(representante, "*", ""))),
            size = 2.4) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits = ylim_arbol, expand = c(0, 0)) +
  labs(x = "Linaje", y = NULL) + tema_col +
  theme(panel.grid = element_blank())

fig3 <- p_arbol + p_muestra + p_linaje + p_calidad +
  plot_layout(widths = c(3.6, 0.28, 0.34, 1.5), nrow = 1)

guardar_fig(fig3, "fig3_arbol_filogenomico", ancho = 12, alto = 6.6)

# Arbol en Newick por si se quiere reproducir en iTOL / FigTree
write.tree(arbol, file.path(DIR_FIG, "fig3_arbol_bac120_nj.newick"))
