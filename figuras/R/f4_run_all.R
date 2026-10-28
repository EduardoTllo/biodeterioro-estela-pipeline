# =============================================================================
# f4_run_all.R — Genera las tablas y figuras de la Fase 4
# Uso (desde la raiz del repo fase1_khipu/):
#   Rscript figuras/R/f4_run_all.R
# =============================================================================

scripts <- c(
  "f3_00_setup.R",
  "f4_01_marcadores.R",
  "f4_02_matriz.R"
)

t0 <- Sys.time()
for (s in scripts) {
  message("\n=== ", s, " ===")
  source(file.path("figuras", "R", s), local = FALSE, echo = FALSE)
}
message("\nListo en ", round(difftime(Sys.time(), t0, units = "secs")), " s.")
