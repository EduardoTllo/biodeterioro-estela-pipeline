# =============================================================================
# euk_run_all.R — Figuras y tablas del analisis eucariota (anexo)
# Uso (desde la raiz del repo fase1_khipu/):
#   Rscript eucariotas/figuras/R/euk_run_all.R
# Salidas: eucariotas/figuras/figs/ (PDF + PNG 600 dpi) y eucariotas/figuras/tablas/
# =============================================================================

scripts <- c("euk_00_setup.R", "euk_01_load_data.R", "euk_02_fig1_embudo.R",
             "euk_03_fig2_calidad.R", "euk_04_fig3_arbol.R", "euk_05_fig4_mash.R",
             "euk_06_fig5_novedad.R", "euk_07_fig6_prevalencia.R",
             "euk_08_fig7_taxonomia.R", "euk_09_tablas.R")

t0 <- Sys.time()
for (s in scripts) {
  message("\n=== ", s, " ===")
  source(file.path("eucariotas", "figuras", "R", s), local = FALSE, echo = FALSE)
}
message("\nListo en ", round(difftime(Sys.time(), t0, units = "secs")), " s.")
