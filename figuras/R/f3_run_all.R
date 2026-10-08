# =============================================================================
# f3_run_all.R — Genera todas las figuras y tablas de la Fase 3
# Uso (desde la raiz del repo fase1_khipu/):
#   Rscript figuras/R/f3_run_all.R
# =============================================================================

scripts <- c(
  "f3_00_setup.R",
  "f3_01_curvas_heaps.R",
  "f3_02_posicion_bin.R",
  "f3_03_arbol_core.R",
  "f3_04_exclusivos.R",
  "f3_05_diagnostico_heaps.R",
  "f3_06_analisis_adicional.R",
  "f3_07_especificos.R",
  "f3_08_sensibilidad_panaroo.R"
)

t0 <- Sys.time()
for (s in scripts) {
  message("\n=== ", s, " ===")
  source(file.path("figuras", "R", s), local = FALSE, echo = FALSE)
}
message("\nListo en ", round(difftime(Sys.time(), t0, units = "secs")), " s.")
