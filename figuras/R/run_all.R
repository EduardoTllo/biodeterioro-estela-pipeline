# =============================================================================
# run_all.R — Genera todas las figuras y tablas de la Fase 2
# Uso (desde la raiz del repo fase1_khipu/):
#   Rscript figuras/R/run_all.R
# =============================================================================

scripts <- c(
  "00_setup.R",
  "01_load_data.R",
  "02_fig1_embudo.R",
  "03_fig2_calidad.R",
  "04_fig3_arbol.R",
  "05_fig4_ani.R",
  "06_fig5_novedad.R",
  "07_fig6_prevalencia.R",
  "08_fig7_taxonomia.R",
  "09_tablas.R"
)

t0 <- Sys.time()
for (s in scripts) {
  message("\n=== ", s, " ===")
  source(file.path("figuras", "R", s), local = FALSE, echo = FALSE)
}
message("\nListo en ", round(difftime(Sys.time(), t0, units = "secs")), " s.")
