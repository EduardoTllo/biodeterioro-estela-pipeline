# =============================================================================
# f4_01_marcadores.R — Cruce de DRAM con los marcadores y el pangenoma (Fase 4)
#
#  1. Funciones por bin (D42): aplica las reglas de metadata/marcadores_fase4.tsv
#     a los KO de KOfam que DRAM asigno a cada bin. Cada alternativa es una
#     lista de pasos. Los pasos de un modulo de KEGG copian su definicion
#     (",": alternativa; "+": subunidad obligatoria; "-": subunidad opcional;
#     espacio: pasos encadenados). Los de bibliografia usan KO sueltos o
#     "n/K1,K2,..." (al menos n de la lista). Una alternativa se cumple si:
#       - estan todos sus pasos esenciales,
#       - le faltan como maximo T pasos, y
#       - le falta menos de la mitad de sus pasos (asi una via de 2 pasos
#         nunca se cumple con 1 solo).
#     T = 1 en el criterio principal, 0 en el estricto y 2 en el flexible.
#     Una funcion esta presente si alguna alternativa se cumple.
#  2. Puntajes por bin: funciones presentes / funciones de la categoria, para
#     riesgo y supervivencia, junto a la completitud CheckM2, con los tres
#     criterios. Provisionales hasta acordarlos con el asesor (D33).
#  3. Fisher por especie del pangenoma (D34, A.1). Unidad: gen del bin que esta
#     en el pangenoma (sin exclusivos ni no_en_pangenoma). Tabla 2x2: marcador
#     si/no x accesorio/core. Marcador = gen con un KO de un paso marcado como
#     marcador = "si" en una alternativa cumplida (los pasos universales de
#     M00033 y M00565 cuentan para la presencia pero no como marcador). Se
#     repite con los tres criterios. BH dentro de cada criterio y definicion
#     de accesorio:
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
CRITERIOS <- c(estricto = 0L, principal = 1L, flexible = 2L)   # pasos que pueden faltar

# Evalua un paso. "n/K1,K2,...": al menos n de la lista. Si no, una expresion
# de KEGG: espacio = pasos encadenados (todos), "," = alternativas,
# "+" = subunidad obligatoria, "-" = subunidad opcional, parentesis = grupo.
evaluar_def <- function(def, kos) {
  if (grepl("^\\d+/", def)) {
    n <- as.integer(sub("/.*", "", def))
    return(sum(strsplit(sub("^\\d+/", "", def), ",")[[1]] %in% kos) >= n)
  }
  s <- gsub("\\s+", " ", trimws(def)); n <- nchar(s); pos <- 1L
  ver <- function() if (pos <= n) substr(s, pos, pos) else ""
  secuencia <- function() {
    v <- alternativas()
    while (ver() == " ") { pos <<- pos + 1L; v <- alternativas() & v }
    v
  }
  alternativas <- function() {
    v <- complejo()
    while (ver() == ",") { pos <<- pos + 1L; v <- complejo() | v }
    v
  }
  complejo <- function() {
    v <- factor_()
    while (ver() %in% c("+", "-")) {
      op <- ver(); pos <<- pos + 1L
      f <- factor_()
      if (op == "+") v <- v & f                     # "-": subunidad opcional, no cuenta
    }
    v
  }
  factor_ <- function() {
    if (ver() == "(") {
      pos <<- pos + 1L; v <- secuencia()
      if (ver() != ")") stop("Falta ')' en: ", def)
      pos <<- pos + 1L
      return(v)
    }
    ko <- regmatches(substring(s, pos), regexpr("^K\\d{5}", substring(s, pos)))
    if (!length(ko)) stop("No se entiende la definicion en la posicion ", pos, ": ", def)
    pos <<- pos + 6L
    ko %in% kos
  }
  v <- secuencia()
  if (pos <= n) stop("Sobra texto en la definicion: ", def)
  v
}

reglas_f4 <- read_tsv(file.path("metadata", "marcadores_fase4.tsv"), show_col_types = FALSE) |>
  mutate(ko_lista = str_extract_all(def, "K\\d{5}"),
         esencial = esencial == "si", marcador = marcador == "si")
stopifnot(all(!is.na(reglas_f4$esencial)), all(!is.na(reglas_f4$marcador)))

# Una fila por paso: si esta, KO hallados y si su alternativa se cumple con
# cada criterio.
evaluar_pasos <- function(kos_genoma) {
  r <- reglas_f4 |>
    mutate(paso_ok = map_lgl(def, evaluar_def, kos = kos_genoma),
           hallados = map(ko_lista, ~ intersect(.x, kos_genoma))) |>
    group_by(funcion, alternativa) |>
    mutate(n_pasos = n(), faltan = sum(!paso_ok), esenciales_ok = all(paso_ok | !esencial)) |>
    ungroup()
  for (k in names(CRITERIOS)) {
    r[[paste0("alt_ok_", k)]] <- r$esenciales_ok & r$faltan <= CRITERIOS[[k]] & r$faltan < r$n_pasos / 2
  }
  r
}
evaluar_funciones <- function(kos_genoma) {
  evaluar_pasos(kos_genoma) |>
    group_by(categoria, funcion) |>
    summarise(presente_estricto = any(alt_ok_estricto),
              presente = any(alt_ok_principal),
              presente_flexible = any(alt_ok_flexible),
              # la alternativa con mas pasos cumplidos, como contexto
              mejor = {
                a <- tibble(alternativa, fuente, n_pasos, faltan) |> distinct() |>
                  arrange(faltan / n_pasos, faltan)
                sprintf("%s %s: %d/%d", a$alternativa[1], a$fuente[1],
                        a$n_pasos[1] - a$faltan[1], a$n_pasos[1])
              },
              kos_hallados = paste(unique(unlist(hallados)), collapse = ","),
              .groups = "drop")
}

local({
  e <- function(d, k) evaluar_def(d, k)
  f <- function(k, crit = "presente") with(evaluar_funciones(k), setNames(get(crit), funcion))
  sox <- c("K17222", "K17223", "K17224", "K17226", "K17227")
  stopifnot(
    e("K17222+K17223-K17225", c("K17222", "K17223")),                     # "-" es opcional
    !e("K17222+K17223", "K17222"),
    e("((K02293 K15744 K00514 K09835),K10027,K15745)", "K10027"),
    !e("((K02293 K15744 K00514 K09835),K10027,K15745)", c("K02293", "K00514")),
    e("3/K06135,K06136,K06137,K06138,K06139", c("K06135", "K06136", "K06139")),
    !f(c("K00117", "K06135", "K06136"))[["acido_gluconico"]],            # gcd con 2 de 5 pqq
    f(c("K00117", "K06135", "K06136", "K06137"))[["acido_gluconico"]],
    f("K19813")[["acido_gluconico"]],                                     # GDH GMC
    !f(c("K00117", "K06135", "K06136", "K06137"))[["acido_2_cetogluconico"]],  # gad es esencial
    !f("K10535")[["nitrificacion"]],                                      # 1 de 2 pasos no basta
    f(c("K28504", "K10945", "K10946", "K10535"))[["nitrificacion"]],
    !f(c("K00370", "K00371"), "presente_flexible")[["nitrificacion"]],    # narGH solo: 1 de 3
    !f(c("K19713"))[["oxidacion_azufre"]],                                # 1 de 2 pasos (S4I)
    f(sox)[["oxidacion_azufre"]],
    f(c("K03897", "K03896", "K03894"))[["sideroforos"]],                  # aerobactina 3 de 4
    !f(c("K03897", "K03896", "K03894"), "presente_estricto")[["sideroforos"]],
    !f(c("K03897", "K03896"), "presente_flexible")[["sideroforos"]],      # 2 de 4: la mitad
    f(c("K02361", "K01252", "K00216", "K02363"))[["sideroforos"]],        # catecolato sin NRPS
    !f("K10208")[["carotenoides"]],                                       # crtN es esencial
    f(c("K10208", "K10209"))[["carotenoides"]],
    f(c("K10208", "K10209"))[["pigmentos"]],
    f(c("K00928", "K00133", "K00836", "K06718"))[["desecacion_ectoina"]], # 4 de 5 (sin ectC)
    !f(c("K00928", "K00133", "K00836", "K06718"), "presente_estricto")[["desecacion_ectoina"]],
    !f(c("K06044", "K01236"))[["desecacion_trehalosa"]],                  # treYZ: 2 de 6 pasos
    f(c("K00697", "K01087"))[["desecacion_trehalosa"]],
    !f("K00697")[["desecacion_trehalosa"]],
    f("K03716")[["reparacion_uv"]],
    f(c("K07699", "K06378"))[["esporulacion"]],                           # 2 de 3
    !f("K07699", "presente_flexible")[["esporulacion"]]
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
col_presente <- c(estricto = "presente_estricto", principal = "presente", flexible = "presente_flexible")
conteo <- map_dfr(names(CRITERIOS), function(k) {
  funciones |>
    filter(categoria != "contexto") |>
    group_by(id, categoria) |> summarise(presentes = sum(.data[[col_presente[[k]]]]), .groups = "drop") |>
    mutate(criterio = k)
}) |>
  left_join(n_func, by = "categoria") |>
  mutate(puntaje = presentes / evaluables)
# Puntajes con el criterio principal (los usa f4_02_matriz.R)
puntajes <- conteo |>
  filter(criterio == "principal") |>
  pivot_wider(id_cols = id, names_from = categoria, values_from = c(presentes, evaluables, puntaje)) |>
  left_join(bins, by = "id") |>
  mutate(puntaje_estado = "provisional (D33)") |>
  arrange(as.integer(sub("^L", "", linaje_id)), genoma)
# Sensibilidad: puntajes con los tres criterios y funciones que cambian (anexo)
sensibilidad <- conteo |>
  pivot_wider(id_cols = c(id, categoria), names_from = criterio, values_from = presentes,
              names_prefix = "presentes_") |>
  left_join(n_func, by = "categoria") |>
  left_join(bins |> select(id, genoma, linaje_id, genero, especie, completitud), by = "id") |>
  arrange(as.integer(sub("^L", "", linaje_id)), genoma, categoria)
cambios <- funciones |>
  filter(presente_estricto != presente_flexible) |>
  arrange(as.integer(sub("^L", "", linaje_id)), genoma, categoria, funcion)

# --- 3. Fisher por especie del pangenoma ----------------------------------------------
fisher <- list(); genes_marcador <- list()
for (i in seq_len(nrow(especies_f3))) {
  s <- especies_f3$slug[i]; b <- especies_f3$bins_ids[i]
  gb <- leer_especie(s, "genes_bin.tsv", col_types = cols(.default = "c")) |>
    left_join(dram |> filter(id == b) |> select(locus_tag, ko), by = "locus_tag")
  pasos <- evaluar_pasos(kos_bin$kos[[match(b, kos_bin$id)]])
  for (k in names(CRITERIOS)) {
    activos <- pasos[pasos[[paste0("alt_ok_", k)]] & pasos$marcador, ]
    for (cat in c("riesgo", "supervivencia")) {
      a <- activos |> filter(categoria == cat)
      kos_cat <- unique(unlist(a$ko_lista))
      marcador <- map_lgl(gb$ko, ~ any(.x %in% kos_cat))
      if (k == "principal") {
        genes_marcador[[paste(s, cat)]] <- gb[marcador, ] |>
          transmute(slug = s, bin = b, categoria = cat, locus_tag, gen, producto, ko = txt(ko),
                    funcion = map_chr(ko, function(x) {
                      paste(unique(a$funcion[map_lgl(a$ko_lista, ~ any(x %in% .x))]), collapse = ",")
                    }),
                    categoria_95, categoria_90, f_ref)
      }
      for (d in names(DEFINICIONES)) {
        part <- gb[[DEFINICIONES[[d]]$col]]
        acc <- part %in% DEFINICIONES[[d]]$acc
        core <- part == "core"
        tab <- matrix(c(sum(marcador & acc), sum(marcador & core),
                        sum(!marcador & acc), sum(!marcador & core)), 2, byrow = TRUE)
        ft <- if (sum(tab[1, ]) > 0) fisher.test(tab) else NULL      # sin genes marcador no hay prueba
        fisher[[length(fisher) + 1]] <- tibble(
          slug = s, bin = b, categoria = cat, criterio = k, definicion = d,
          marcador_accesorio = tab[1, 1], marcador_core = tab[1, 2],
          otro_accesorio = tab[2, 1], otro_core = tab[2, 2],
          odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
          ic95_inf = if (is.null(ft)) NA_real_ else ft$conf.int[1],
          ic95_sup = if (is.null(ft)) NA_real_ else ft$conf.int[2],
          p = if (is.null(ft)) NA_real_ else ft$p.value)
      }
    }
  }
}
fisher <- bind_rows(fisher) |>
  group_by(criterio, definicion) |> mutate(p_bh = p.adjust(p, method = "BH")) |> ungroup()
genes_marcador <- bind_rows(genes_marcador)

# --- 4. Genes especificos de la cepa con DRAM ------------------------------------------
reglas_marcador <- reglas_f4 |> filter(marcador)
funcion_de_ko <- function(k) {
  paste(unique(reglas_marcador$funcion[map_lgl(reglas_marcador$ko_lista, ~ any(k %in% .x))]), collapse = ",")
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
                              presente, presente_estricto, presente_flexible, mejor, kos_hallados),
          file.path(DIR_TAB, "tab_f4_funciones_bin.tsv"))
write_tsv(puntajes |> select(genoma, linaje_id, genero, especie, completitud, contaminacion,
                             starts_with("presentes_"), starts_with("evaluables_"),
                             starts_with("puntaje_")),
          file.path(DIR_TAB, "tab_f4_puntajes.tsv"))
write_tsv(sensibilidad |> select(genoma, linaje_id, genero, especie, completitud, categoria,
                                 presentes_estricto, presentes_principal, presentes_flexible, evaluables),
          file.path(DIR_TAB, "tab_f4_sensibilidad.tsv"))
write_tsv(cambios |> select(genoma, linaje_id, genero, especie, completitud, categoria, funcion,
                            presente_estricto, presente, presente_flexible, mejor),
          file.path(DIR_TAB, "tab_f4_sensibilidad_funciones.tsv"))
write_tsv(fisher, file.path(DIR_TAB, "tab_f4_fisher.tsv"))
write_tsv(genes_marcador, file.path(DIR_TAB, "tab_f4_genes_marcador.tsv"))
write_tsv(especificos, file.path(DIR_TAB, "tab_f4_especificos_dram.tsv"))
write_tsv(cazy, file.path(DIR_TAB, "tab_f4_cazy_bin.tsv"))
message("  -> ", file.path(DIR_TAB, paste0("tab_f4_{funciones_bin,puntajes,sensibilidad,",
                                         "sensibilidad_funciones,fisher,genes_marcador,",
                                         "especificos_dram,cazy_bin}.tsv")))
print(fisher |> filter(criterio == "principal", definicion == "principal") |>
        select(slug, categoria, marcador_accesorio, marcador_core, odds_ratio, p, p_bh))
