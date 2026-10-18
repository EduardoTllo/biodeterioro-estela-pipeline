# =============================================================================
# euk_00_setup.R — Configuracion de las figuras del analisis eucariota (anexo)
# Reutiliza el tema, guardar_fig() y los umbrales de figuras/R/00_setup.R.
# Uso: desde la raiz del repo (fase1_khipu/), via eucariotas/figuras/R/euk_run_all.R
# =============================================================================

source(file.path("figuras", "R", "00_setup.R"))

EUK     <- "eucariotas"
RES_E1  <- file.path(EUK, "results", "fase1", "results")
RES_E2  <- file.path(EUK, "results", "fase2", "results")
DIR_FIG <- file.path(EUK, "figuras", "figs")      # guardar_fig() usa esta global
DIR_TAB <- file.path(EUK, "figuras", "tablas")
DIR_KRAKEN <- file.path("..", "Datos", "Shotgun Analysis")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

# --- Paletas -----------------------------------------------------------------
# Color por CLASE (hay 2 filos, pero 12 de 13 genomas son Ascomycota).
# Las tres clases de Ascomycota usan los tres primeros colores de la paleta de
# referencia, validados para todos los pares (validate_palette.js, --pairs all:
# CVD dE >= 9.2, vision normal dE >= 24). La unica Basidiomycota va en gris
# neutro: un cuarto color no pasa la validacion. Todas llevan forma distinta y
# etiqueta directa (el aqua necesita etiqueta por contraste < 3:1).
PAL_CLASE <- c(
  "Eurotiomycetes"  = "#2a78d6",
  "Dothideomycetes" = "#eb6834",
  "Sordariomycetes" = "#1baf7a",
  "Tremellomycetes" = "#6b6a66"
)
SHP_CLASE <- c("Eurotiomycetes" = 16, "Dothideomycetes" = 17,
               "Sordariomycetes" = 15, "Tremellomycetes" = 18)

PAL_DECISION_EUK <- c(
  "Seleccionado"             = "#1B7837",
  "Rechazado: contaminacion" = "#B2182B",
  "Rechazado: completitud"   = "#EF8A62",
  "Sin estimacion de EukCC"  = "#999999"
)

# Etiqueta corta segun el nivel alcanzado en la cascada de la Fase 2.
taxon_corto_euk <- function(especie, genero, familia, orden, clase) {
  case_when(
    !is.na(especie) ~ especie,
    !is.na(genero)  ~ paste0(genero, " sp."),
    !is.na(familia) ~ paste0(familia, " (sin genero)"),
    !is.na(orden)   ~ paste0(orden, " (sin familia)"),
    TRUE            ~ clase
  )
}
