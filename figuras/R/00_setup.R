# =============================================================================
# 00_setup.R — Configuracion comun para las figuras y tablas de la Fase 2
# Tesis: Modelo predictivo de riesgo de deterioro de la Estela de Raimondi
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# --- Rutas -------------------------------------------------------------------
# Se asume que el working directory es la raiz del repo (fase1_khipu/).
if (!dir.exists("results")) {
  stop("Ejecuta los scripts desde la raiz del repo (fase1_khipu/). WD actual: ", getwd())
}

RES_F1  <- file.path("results", "fase1")
RES_F2  <- file.path("results", "fase2")
DIR_FIG <- file.path("figuras", "figs")
DIR_TAB <- file.path("figuras", "tablas")

dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

# --- Tema grafico ------------------------------------------------------------
tema_tesis <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.25, colour = "grey90"),
      panel.border       = element_rect(linewidth = 0.4, colour = "grey40"),
      strip.background   = element_rect(fill = "grey95", colour = "grey40",
                                        linewidth = 0.4),
      strip.text         = element_text(face = "bold", size = base_size - 1),
      axis.title         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = element_text(colour = "grey30", size = base_size - 1),
      plot.tag           = element_text(face = "bold", size = base_size + 3),
      legend.key.size    = unit(0.9, "lines"),
      legend.background  = element_blank(),
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 1)
    )
}
theme_set(tema_tesis())

# --- Paletas -----------------------------------------------------------------
# Filos presentes en el set (3). Colores estables en todas las figuras.
PAL_FILO <- c(
  "Bacillota"       = "#3B7DD8",
  "Pseudomonadota"  = "#D95F02",
  "Actinomycetota"  = "#1B9E77"
)

PAL_DECISION <- c(
  "Seleccionado"            = "#1B7837",
  "Rechazado: contaminacion"= "#B2182B",
  "Rechazado: completitud"  = "#EF8A62",
  "Rechazado: ambos"        = "#762A83",
  "Descartado: dominio"     = "#999999"
)

PAL_ESPECIE <- c("Especie confirmada" = "#2166AC",
                 "Sin especie asignada" = "#B2182B")

# --- Umbrales del estudio ----------------------------------------------------
UMBRAL_COMPLETITUD   <- 70   # CheckM2, Fase 1
UMBRAL_CONTAMINACION <- 5    # CheckM2, Fase 1
UMBRAL_ANI           <- 95   # confirmacion de especie / definicion de linaje
UMBRAL_AF            <- 65   # fraccion alineada minima

# --- Guardado ----------------------------------------------------------------
# Guarda cada figura en PDF vectorial (para LaTeX/Word) y PNG 600 dpi.
guardar_fig <- function(plot, nombre, ancho, alto) {
  pdf_path <- file.path(DIR_FIG, paste0(nombre, ".pdf"))
  png_path <- file.path(DIR_FIG, paste0(nombre, ".png"))
  ggsave(pdf_path, plot, width = ancho, height = alto, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = ancho, height = alto, units = "in", dpi = 600)
  message("  -> ", pdf_path)
  message("  -> ", png_path)
  invisible(plot)
}

# --- Utilidades taxonomicas --------------------------------------------------
RANGOS <- c("dominio", "filo", "clase", "orden", "familia", "genero", "especie")

# Separa el string GTDB "d__X;p__Y;..." en 7 columnas limpias (sin prefijos).
separar_taxonomia <- function(df, col = "clasificacion", prefijo = "") {
  tax <- str_split_fixed(df[[col]], ";", 7)
  colnames(tax) <- paste0(prefijo, RANGOS)
  tax <- as_tibble(tax) %>%
    mutate(across(everything(), ~ str_remove(.x, "^[a-z]__"))) %>%
    mutate(across(everything(), ~ na_if(str_trim(.x), "")))
  bind_cols(df, tax)
}

# Etiqueta corta y legible para un genoma: "bin-6-50 (Telluria timonae)"
etiqueta_genoma <- function(genoma, genero, especie) {
  taxon <- if_else(!is.na(especie), especie, paste0(genero, " sp."))
  paste0(genoma, " (", taxon, ")")
}
