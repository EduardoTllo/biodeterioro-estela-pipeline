# =============================================================================
# 09_tablas.R — Tablas T1-T6 de la Fase 2
# Salidas: TSV (para reutilizar), XLSX (suplementario) y HTML/PNG con `gt`
# (formateadas para la tesis).
# =============================================================================

suppressPackageStartupMessages({
  library(gt)
  library(writexl)
})

# --- T1: tabla maestra por genoma ---------------------------------------------
T1 <- maestra %>%
  transmute(
    Genoma = genoma,
    Muestra = muestra,
    Linaje = linaje_id,
    Representante = if_else(representante, "si", "no"),
    `Tamano (Mb)` = round(tamano_mb, 2),
    Contigs = Total_Contigs,
    `N50 (kb)` = round(n50_kb, 1),
    `GC (%)` = round(gc_pct, 1),
    `Completitud (%)` = Completeness,
    `Contaminacion (%)` = Contamination,
    `Calidad CheckM2` = calidad_checkm2,
    Filo = as.character(filo), Clase = clase, Orden = orden,
    Familia = familia, Genero = genero, Especie = especie,
    `Referencia GTDB` = referencia,
    `ANI (%)` = ani,
    `AF (%)` = af,
    `Especie confirmada` = especie_confirmada,
    `MSA retenido (%)` = msa_percent,
    `Marcadores bac120` = unicos
  )

# --- T2: resumen por linaje ---------------------------------------------------
T2 <- linajes %>%
  transmute(
    Linaje = linaje_id,
    `Taxon (GTDB)` = taxon,
    Filo = as.character(filo),
    `Nº genomas` = n_genomas,
    Muestras = muestras,
    Prevalencia = prevalencia,
    Representante = repres,
    `Completitud rep. (%)` = round(completitud_rep, 2),
    `Contaminacion rep. (%)` = round(contam_rep, 2),
    `Especie confirmada` = coalesce(especie_conf, "NA"),
    `Fase 3` = if_else(fase3, "SI", "no")
  )

# --- T3: estadisticas agregadas por filo --------------------------------------
resumir <- function(df, grupo) {
  df %>%
    group_by({{ grupo }}) %>%
    summarise(
      n = n(),
      `Completitud (mediana)`   = median(Completeness),
      `Completitud (rango)`     = sprintf("%.2f-%.2f", min(Completeness), max(Completeness)),
      `Contaminacion (mediana)` = median(Contamination),
      `Contaminacion (max)`     = max(Contamination),
      `Tamano Mb (mediana)`     = round(median(tamano_mb), 2),
      `Tamano Mb (rango)`       = sprintf("%.2f-%.2f", min(tamano_mb), max(tamano_mb)),
      `GC % (mediana)`          = round(median(gc_pct), 1),
      `N50 kb (mediana)`        = round(median(n50_kb), 1),
      `Contigs (mediana)`       = median(Total_Contigs),
      .groups = "drop"
    )
}

T3 <- bind_rows(
  resumir(maestra, filo) %>% rename(Grupo = filo) %>% mutate(Grupo = as.character(Grupo)),
  resumir(mutate(maestra, todos = "TOTAL"), todos) %>% rename(Grupo = todos)
)

# --- T4: candidatos a novedad (generado en 06_fig5_novedad.R) -----------------
T4 <- read_tsv(file.path(DIR_TAB, "T4_candidatos_novedad.tsv"),
               show_col_types = FALSE) %>%
  transmute(
    Genoma = genoma, Muestra = muestra, Linaje = linaje_id,
    Filo = filo, Familia = familia, Genero = genero,
    `Referencia mas cercana` = referencia,
    `ANI (%)` = ani, `AF (%)` = af,
    `Metodo de clasificacion` = classification_method,
    Nota = note
  )

# --- T5: ANI intracluster (generada en 05_fig4_ani.R) y T6: versiones --------
T5 <- read_tsv(file.path(DIR_TAB, "T5_ani_intracluster.tsv"), show_col_types = FALSE)

T6 <- tibble(linea = read_lines(file.path(RES_F2, "phase2_versions.txt"))) %>%
  filter(str_trim(linea) != "")

# --- Escritura ----------------------------------------------------------------
write_tsv(T1, file.path(DIR_TAB, "T1_genomas.tsv"))
write_tsv(T2, file.path(DIR_TAB, "T2_linajes.tsv"))
write_tsv(T3, file.path(DIR_TAB, "T3_estadisticas.tsv"))

write_xlsx(list(T1_genomas = T1, T2_linajes = T2, T3_estadisticas = T3,
                T4_novedad = T4, T5_ani_intracluster = T5, T6_versiones = T6),
           file.path(DIR_TAB, "tablas_suplementarias_fase2.xlsx"))

# --- Versiones formateadas con gt ---------------------------------------------
gt_T1 <- T1 %>%
  gt() %>%
  tab_header(
    title = md("**Tabla 1.** Genomas seleccionados en la Fase 2 (n = 21)"),
    subtitle = "Metricas de ensamblaje, calidad CheckM2 y asignacion taxonomica GTDB-Tk (R220)"
  ) %>%
  fmt_number(columns = c("ANI (%)", "AF (%)", "MSA retenido (%)",
                         "Completitud (%)", "Contaminacion (%)"), decimals = 2) %>%
  sub_missing(missing_text = "—") %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(columns = c(Genero, Especie))) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Genoma)) %>%
  tab_options(table.font.size = px(11), data_row.padding = px(3))

gt_T2 <- T2 %>%
  gt() %>%
  tab_header(
    title = md("**Tabla 2.** Linajes definidos por desreplicacion a 95 % de ANI (n = 19)"),
    subtitle = "Ordenados por prevalencia espacial; los tres primeros pasan a la Fase 3"
  ) %>%
  tab_style(style = cell_fill(color = "#FCEFC7"),
            locations = cells_body(rows = `Fase 3` == "SI")) %>%
  tab_style(style = cell_text(style = "italic"),
            locations = cells_body(columns = `Taxon (GTDB)`)) %>%
  tab_options(table.font.size = px(11), data_row.padding = px(3))

gt_T3 <- T3 %>%
  gt() %>%
  tab_header(title = md("**Tabla 3.** Estadisticas agregadas por filo")) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(rows = Grupo == "TOTAL")) %>%
  tab_options(table.font.size = px(11), data_row.padding = px(3))

gtsave(gt_T1, file.path(DIR_TAB, "T1_genomas.html"))
gtsave(gt_T2, file.path(DIR_TAB, "T2_linajes.html"))
gtsave(gt_T3, file.path(DIR_TAB, "T3_estadisticas.html"))

message("Tablas escritas en ", DIR_TAB)
