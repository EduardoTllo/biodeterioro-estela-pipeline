# =============================================================================
# f4_01_marcadores.R — Cruce de DRAM con los marcadores y el pangenoma (Fase 4)
#
#  1. Funciones por bin: aplica las reglas de metadata/marcadores_fase4.tsv
#     (D35) a los KO de KOfam que DRAM asigno a cada bin. Una funcion esta
#     presente si alguna alternativa cumple todos sus grupos (cada grupo exige
#     un minimo de KO presentes). Los modulos (M00528, M00595, M00097) ya estan
#     escritos en la tabla como "n de m KO".
#  2. Puntajes por bin: funciones presentes / funciones de la categoria, para
#     riesgo y supervivencia, junto a la completitud CheckM2. Provisionales
#     hasta acordarlos con el asesor (D33).
#  3. Fisher por especie del pangenoma (D34, A.1). Unidad: gen del bin que esta
#     en el pangenoma (sin exclusivos ni no_en_pangenoma). Tabla 2x2: marcador
#     si/no x accesorio/core. Marcador = gen con un KO de una alternativa
#     cumplida (la funcion esta presente en el bin). BH dentro de cada
#     definicion de accesorio:
#       principal   accesorio = shell + cloud (D27)
#       solo_cloud  accesorio = cloud; el shell queda fuera de la tabla
#       core_90     particion con core >= 90 %
#  4. Genes especificos de la cepa (D34, A.2) con la anotacion de DRAM.
#  5. Contexto: CAZymes de dbCAN por familia y bin.
#
# Entradas: results/fase4/dram/annotations.tsv[.gz] (DRAM annotate_genes sobre
# las proteinas de Bakta, ver README_fase4.md), results/fase1, results/fase2 y
# results/fase3/<especie>/{genes_bin,exclusivos_verificados}.tsv.
# =============================================================================

RES_F4 <- file.path("results", "fase4")
CLASES_ESP <- c("especifico", "especifico_con_bandera", "verificado", "con_bandera")
DEFINICIONES <- list(
  principal  = list(col = "categoria_95", acc = c("shell", "cloud")),
  solo_cloud = list(col = "categoria_95", acc = "cloud"),
  core_90    = list(col = "categoria_90", acc = c("shell", "cloud"))
)

id_de <- function(x) gsub("[^A-Za-z0-9]", "_", x)    # bin-4-52 -> bin_4_52, como en la Fase 3
txt <- function(l) map_chr(l, ~ paste(.x, collapse = ","))

# --- Reglas de la tabla de marcadores -----------------------------------------
reglas_f4 <- read_tsv(file.path("metadata", "marcadores_fase4.tsv"), show_col_types = FALSE) |>
  mutate(ko_lista = strsplit(kos, ","))

# Una fila por grupo: KO hallados, si el grupo llega a su minimo y si su
# alternativa (todos sus grupos) se cumple.
evaluar_grupos <- function(kos_genoma) {
  reglas_f4 |>
    mutate(hallados = map(ko_lista, ~ intersect(.x, kos_genoma)),
           grupo_ok = lengths(hallados) >= minimo) |>
    group_by(funcion, alternativa) |> mutate(alt_ok = all(grupo_ok)) |> ungroup()
}
evaluar_funciones <- function(kos_genoma) {
  evaluar_grupos(kos_genoma) |>
    group_by(categoria, funcion) |>
    summarise(presente = any(alt_ok), kos_hallados = paste(unique(unlist(hallados)), collapse = ","),
              .groups = "drop")
}

local({
  f <- function(k) with(evaluar_funciones(k), setNames(presente, funcion))
  stopifnot(
    !f(c("K00117", "K06135", "K06136"))[["acido_gluconico"]],            # gcd con 2 de 5 pqq
    f(c("K00117", "K06135", "K06136", "K06137"))[["acido_gluconico"]],
    !f(c("K00117", "K06135", "K06136", "K06137"))[["acido_2_cetogluconico"]],
    f(c("K06044", "K01236"))[["desecacion_trehalosa"]],                  # solo la ruta TreYZ
    !f("K06044")[["desecacion_trehalosa"]],
    f(sprintf("K%05d", 19419:19427))[["eps_biopelicula"]],               # 9 de 13 eps
    !f(sprintf("K%05d", 19419:19426))[["eps_biopelicula"]],
    f(c("K06336", "K19433"))[["eps_biopelicula"]],
    f(c("K28504", "K10945", "K10535"))[["nitrificacion"]],               # 3 de 4
    !f(c("K28504", "K10945"))[["nitrificacion"]]
  )
})

# --- DRAM -------------------------------------------------------------------------
p_dram <- file.path(RES_F4, "dram", c("annotations.tsv", "annotations.tsv.gz"))
p_dram <- p_dram[file.exists(p_dram)][1]
if (is.na(p_dram)) {
  stop("Falta la anotacion de DRAM en ", file.path(RES_F4, "dram"), ". Ver README_fase4.md.")
}
dram <- read_tsv(p_dram, show_col_types = FALSE, col_types = cols(.default = "c"),
                 name_repair = "minimal") |>
  rename(gen_dram = 1)
for (k in c("fasta", "ko_id")) {
  if (!k %in% names(dram)) stop("La anotacion de DRAM no tiene la columna '", k, "': ", p_dram)
}
dram <- dram |>
  mutate(id = id_de(fasta),
         # annotate_genes antepone el nombre del archivo al ID de Bakta
         locus_tag = str_remove(gen_dram, paste0("^\\Q", fasta, "_\\E")),
         ko = map(str_extract_all(coalesce(ko_id, ""), "K\\d{5}"), unique),
         cazy = map(str_extract_all(coalesce(if ("cazy_ids" %in% names(dram)) cazy_ids else NA, ""),
                                    "(GH|GT|PL|CE|AA|CBM)\\d+"), unique))

for (i in seq_len(nrow(especies_f3))) {
  gb <- leer_especie(especies_f3$slug[i], "genes_bin.tsv", col_types = cols(.default = "c"))
  if (is.null(gb)) stop("Falta genes_bin.tsv de ", especies_f3$slug[i])
  frac <- mean(gb$locus_tag %in% dram$locus_tag[dram$id == especies_f3$bins_ids[i]])
  if (frac < 0.95) {
    stop(sprintf("Solo %.1f %% de los genes de %s calzan con DRAM. Revisa los IDs de annotations.tsv.",
                 100 * frac, especies_f3$bins_ids[i]))
  }
}

bins <- read_tsv(file.path(RES_F2, "phase2_selection.tsv"), show_col_types = FALSE) |>
  separar_taxonomia() |>
  left_join(read_tsv(file.path(RES_F1, "phase1_selection.tsv"), show_col_types = FALSE) |>
              select(genome = bin, completitud = completeness, contaminacion = contamination),
            by = "genome") |>
  transmute(id = id_de(genome), genoma = genome, linaje_id, genero, especie, completitud, contaminacion)
faltan <- setdiff(unique(dram$id), bins$id)
if (length(faltan)) stop("Bins de DRAM sin datos de la Fase 2: ", paste(faltan, collapse = ", "))
message("  DRAM: ", n_distinct(dram$id), " bins, ", nrow(dram), " genes, ",
        sum(lengths(dram$ko) > 0), " con KO")

# --- 1-2. Funciones y puntajes por bin -------------------------------------------------
kos_bin <- dram |> group_by(id) |> summarise(kos = list(unique(unlist(ko))), .groups = "drop")
funciones <- kos_bin |>
  mutate(f = map(kos, evaluar_funciones)) |>
  select(id, f) |> unnest(f) |>
  left_join(bins, by = "id")

n_func <- reglas_f4 |> distinct(categoria, funcion) |> count(categoria, name = "evaluables")
puntajes <- funciones |>
  filter(categoria != "contexto") |>
  group_by(id, categoria) |> summarise(presentes = sum(presente), .groups = "drop") |>
  left_join(n_func, by = "categoria") |>
  mutate(puntaje = presentes / evaluables) |>
  pivot_wider(id_cols = id, names_from = categoria, values_from = c(presentes, evaluables, puntaje)) |>
  left_join(bins, by = "id") |>
  mutate(puntaje_estado = "provisional (D33)") |>
  arrange(as.integer(sub("^L", "", linaje_id)), genoma)

# --- 3. Fisher por especie del pangenoma ----------------------------------------------
fisher <- list(); genes_marcador <- list()
for (i in seq_len(nrow(especies_f3))) {
  s <- especies_f3$slug[i]; b <- especies_f3$bins_ids[i]
  gb <- leer_especie(s, "genes_bin.tsv", col_types = cols(.default = "c")) |>
    left_join(dram |> filter(id == b) |> select(locus_tag, ko), by = "locus_tag")
  activos <- evaluar_grupos(kos_bin$kos[[match(b, kos_bin$id)]]) |> filter(alt_ok)
  for (cat in c("riesgo", "supervivencia")) {
    a <- activos |> filter(categoria == cat)
    kos_cat <- unique(unlist(a$ko_lista))
    marcador <- map_lgl(gb$ko, ~ any(.x %in% kos_cat))
    genes_marcador[[paste(s, cat)]] <- gb[marcador, ] |>
      transmute(slug = s, bin = b, categoria = cat, locus_tag, gen, producto, ko = txt(ko),
                funcion = map_chr(ko, function(k) paste(unique(a$funcion[map_lgl(a$ko_lista, ~ any(k %in% .x))]),
                                                        collapse = ",")),
                categoria_95, categoria_90, f_ref)
    for (d in names(DEFINICIONES)) {
      part <- gb[[DEFINICIONES[[d]]$col]]
      acc <- part %in% DEFINICIONES[[d]]$acc
      core <- part == "core"
      tab <- matrix(c(sum(marcador & acc), sum(marcador & core),
                      sum(!marcador & acc), sum(!marcador & core)), 2, byrow = TRUE)
      ft <- if (sum(tab[1, ]) > 0) fisher.test(tab) else NULL      # sin genes marcador no hay prueba
      fisher[[length(fisher) + 1]] <- tibble(
        slug = s, bin = b, categoria = cat, definicion = d,
        marcador_accesorio = tab[1, 1], marcador_core = tab[1, 2],
        otro_accesorio = tab[2, 1], otro_core = tab[2, 2],
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        ic95_inf = if (is.null(ft)) NA_real_ else ft$conf.int[1],
        ic95_sup = if (is.null(ft)) NA_real_ else ft$conf.int[2],
        p = if (is.null(ft)) NA_real_ else ft$p.value)
    }
  }
}
fisher <- bind_rows(fisher) |>
  group_by(definicion) |> mutate(p_bh = p.adjust(p, method = "BH")) |> ungroup()
genes_marcador <- bind_rows(genes_marcador)

# --- 4. Genes especificos de la cepa con DRAM ------------------------------------------
funcion_de_ko <- function(k) {
  paste(unique(reglas_f4$funcion[map_lgl(reglas_f4$ko_lista, ~ any(k %in% .x))]), collapse = ",")
}
especificos <- map_dfr(especies_f3$slug, function(s) {
  ev <- leer_especie(s, "exclusivos_verificados.tsv", col_types = cols(.default = "c"))
  if (is.null(ev)) return(NULL)
  ev |> filter(clase %in% CLASES_ESP) |> mutate(slug = s, id = bin)
})
if (nrow(especificos)) {
  especificos <- especificos |>
    left_join(dram |> select(id, locus_tag, ko, any_of(c("kegg_hit", "cazy_ids", "peptidase_family"))),
              by = c("id", "locus_tag")) |>
    mutate(marcador_fase4 = map_chr(ko, funcion_de_ko), ko = txt(ko)) |>
    select(slug, bin, locus_tag, clase, contig, gen, producto, ko, any_of(c("kegg_hit", "cazy_ids",
           "peptidase_family")), marcador_fase4, f5_origen, f5_pident, f5_taxon)
}

# --- 5. Contexto: CAZymes por familia y bin ----------------------------------------------
cazy <- dram |>
  select(id, locus_tag, cazy) |> unnest_longer(cazy, values_to = "familia") |>
  distinct() |> count(id, familia, name = "genes") |>
  mutate(clase_cazy = str_extract(familia, "^[A-Z]+")) |>
  left_join(bins |> select(id, genoma, linaje_id), by = "id") |>
  select(genoma, linaje_id, clase_cazy, familia, genes)

# --- Salidas -------------------------------------------------------------------------
write_tsv(funciones |> select(genoma, linaje_id, genero, especie, completitud, categoria, funcion,
                              presente, kos_hallados),
          file.path(DIR_TAB, "tab_f4_funciones_bin.tsv"))
write_tsv(puntajes |> select(genoma, linaje_id, genero, especie, completitud, contaminacion,
                             starts_with("presentes_"), starts_with("evaluables_"),
                             starts_with("puntaje_")),
          file.path(DIR_TAB, "tab_f4_puntajes.tsv"))
write_tsv(fisher, file.path(DIR_TAB, "tab_f4_fisher.tsv"))
write_tsv(genes_marcador, file.path(DIR_TAB, "tab_f4_genes_marcador.tsv"))
write_tsv(especificos, file.path(DIR_TAB, "tab_f4_especificos_dram.tsv"))
write_tsv(cazy, file.path(DIR_TAB, "tab_f4_cazy_bin.tsv"))
message("  -> ", file.path(DIR_TAB, paste0("tab_f4_{funciones_bin,puntajes,fisher,genes_marcador,",
                                         "especificos_dram,cazy_bin}.tsv")))
print(fisher |> filter(definicion == "principal") |>
        select(slug, categoria, marcador_accesorio, marcador_core, odds_ratio, p, p_bh))
