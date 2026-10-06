# =============================================================================
# f3_00_setup.R — Configuracion comun para las figuras y tablas de la Fase 3
# Reutiliza el tema, las rutas y guardar_fig() de 00_setup.R (Fase 2).
# Entradas: results/fase3/ (traidas de Khipu con bin/export_fase3.sh + rsync).
# =============================================================================

source(file.path("figuras", "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
})
for (p in c("phangorn")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Falta el paquete '", p, "'. Instalalo con install.packages('", p, "').")
  }
}

# --- Ley de Heaps ----------------------------------------------------------------
# Reproduce exactamente micropan::heaps() (Snipen & Liland 2015; modelo de
# Tettelin et al. 2008): familias nuevas al agregar el genoma N ~ k * N^(-alpha),
# ajustado con optim (L-BFGS-B) sobre n.perm ordenes aleatorios. Misma formula,
# mismas permutaciones y mismo ajuste; da el mismo resultado que micropan con la
# misma semilla (verificado), pero es vectorizada: micropan recorre columna por
# columna y es demasiado lenta para el jackknife de f3_05.
# pm: matriz genomas x familias (0/1). alpha < 1 indica pangenoma abierto.
heaps_rapido <- function(pm, n.perm = 100) {
  pm <- (pm > 0) * 1L
  pm <- pm[, colSums(pm) > 0, drop = FALSE]
  ng <- nrow(pm)
  nmat <- matrix(0, nrow = ng - 1, ncol = n.perm)
  for (i in seq_len(n.perm)) {
    p <- pm[sample(ng), , drop = FALSE]
    primero <- max.col(t(p), ties.method = "first")   # genoma donde aparece cada familia
    nmat[, i] <- tabulate(primero, nbins = ng)[2:ng]
  }
  x <- rep(2:ng, times = n.perm)
  y <- as.numeric(nmat)
  obj <- function(p, x, y) sqrt(sum((y - p[1] * x^(-p[2]))^2)) / length(x)
  fit <- optim(c(mean(y[x == 2]), 1), obj, gr = NULL, x, y, method = "L-BFGS-B",
               lower = c(0, 0), upper = c(10000, 2))
  setNames(fit$par, c("Intercept", "alpha"))
}

# Curvas de acumulacion del pan y del core para un orden de genomas.
# pan[k]: familias vistas en los primeros k genomas; core[k]: familias presentes
# en todos los primeros k genomas.
curva_acumulacion <- function(p) {
  n <- nrow(p)
  primero <- max.col(t(p), ties.method = "first")
  primer_cero <- max.col(t(1L - p), ties.method = "first")
  primer_cero[colSums(p) == n] <- n + 1L                # presentes en todos
  list(pan = cumsum(tabulate(primero, nbins = n)),
       core = ncol(p) - cumsum(tabulate(primer_cero, nbins = n + 1))[seq_len(n)])
}

RES_F3 <- file.path("results", "fase3")
p_especies <- file.path(RES_F3, "phase3_especies_seleccionadas.tsv")
if (!file.exists(p_especies)) {
  stop("Faltan los resultados de la Fase 3 en ", RES_F3,
       ". Corre bin/export_fase3.sh en Khipu y trae export_fase3/ con rsync.")
}

especies_f3 <- read_tsv(p_especies, show_col_types = FALSE) |>
  mutate(etiqueta = paste0(especie_gtdb, " (", linaje_id, ")"))

refs_f3 <- {
  p <- file.path(RES_F3, "phase3_referencias_finales.tsv")
  if (file.exists(p)) read_tsv(p, show_col_types = FALSE, col_types = cols(.default = "c"))
  else tibble(slug = character(), id = character(), habitat = character(),
              continente = character(), nivel = character())
}

leer_especie <- function(slug, archivo, ...) {
  p <- file.path(RES_F3, slug, archivo)
  if (!file.exists(p)) return(NULL)
  read_tsv(p, show_col_types = FALSE, ...)
}

# --- Paletas ------------------------------------------------------------------
PAL_CATEGORIA <- c(
  "core"            = "#2166AC",
  "shell"           = "#67A9CF",
  "cloud"           = "#D1E5F0",
  "exclusivo"       = "#B2182B",
  "no_en_pangenoma" = "#BDBDBD"
)

PAL_HABITAT <- c(
  "petreo_arido"        = "#8C510A",
  "suelo"               = "#BF812D",
  "agua_sedimento"      = "#35978F",
  "planta"              = "#5AAE61",
  "otro_ambiental"      = "#80CDC1",
  "animal"              = "#C51B7D",
  "alimento_industrial" = "#F1B6DA",
  "clinico_humano"      = "#762A83",
  "sin_clasificar"      = "#D9D9D9",
  "desconocido"         = "#F0F0F0",
  "Estela de Raimondi"  = "#B2182B"
)

UMBRAL_CORE  <- 0.95
UMBRAL_CLOUD <- 0.15
N_PERM_CURVA <- 100
N_PERM_HEAPS <- 500
