# =============================================================================
# 06_fig5_novedad.R — Figura 5: ANI vs fraccion alineada frente a la referencia
# GTDB mas cercana. Identifica candidatos a novedad taxonomica.
#
# El reparto por rango taxonomico (16 genomas con especie, 5 solo hasta genero)
# se comenta en el texto: son dos cifras y no justifican un panel propio.
# =============================================================================

nov <- maestra %>%
  filter(!is.na(ani), !is.na(af)) %>%
  mutate(
    estado = if_else(ani >= UMBRAL_ANI & af >= UMBRAL_AF,
                     "Especie confirmada", "Candidato a novedad"),
    estado = factor(estado, levels = c("Especie confirmada", "Candidato a novedad"))
  )

sin_ani <- maestra %>% filter(is.na(ani) | is.na(af))

pal_estado <- c("Especie confirmada" = "#2166AC", "Candidato a novedad" = "#B2182B")

fig5 <- ggplot(nov, aes(x = ani, y = af)) +
  annotate("rect", xmin = UMBRAL_ANI, xmax = 100.6, ymin = UMBRAL_AF, ymax = 102,
           fill = "#2166AC", alpha = 0.07) +
  geom_vline(xintercept = UMBRAL_ANI, linetype = "dashed", colour = "grey35",
             linewidth = 0.45) +
  geom_hline(yintercept = UMBRAL_AF, linetype = "dashed", colour = "grey35",
             linewidth = 0.45) +
  geom_point(aes(colour = estado, shape = filo), size = 2.8) +
  # Los genomas dentro del cuadrante de especie confirmada se etiquetan solo con
  # su codigo (el taxon esta en la Tabla 1); los candidatos a novedad, completos.
  geom_text_repel(data = filter(nov, estado == "Especie confirmada"),
                  aes(label = genoma), colour = "#2166AC",
                  size = 2.4, seed = 7, box.padding = 0.28,
                  min.segment.length = 0.1, max.overlaps = 40,
                  show.legend = FALSE) +
  geom_text_repel(data = filter(nov, estado == "Candidato a novedad"),
                  aes(label = paste0(genoma, "\n", taxon_corto)),
                  colour = "#B2182B", size = 2.6, lineheight = 0.95,
                  seed = 7, nudge_y = -3, box.padding = 0.4,
                  min.segment.length = 0.1, show.legend = FALSE) +
  annotate("text", x = 100.4, y = UMBRAL_AF + 1.5, hjust = 1, size = 2.7,
           colour = "grey30",
           label = paste0("especie confirmada: ANI >= ", UMBRAL_ANI,
                          " % y AF >= ", UMBRAL_AF, " %")) +
  annotate("label", x = 89.1, y = 100, hjust = 0, size = 2.5, label.size = 0,
           colour = "grey30", fill = alpha("white", 0.75),
           label = paste0("Sin ANI calculado (n = ", nrow(sin_ani), "): ",
                          paste(sin_ani$genoma, collapse = ", "),
                          "\nGTDB-Tk los clasifico solo por topologia del arbol")) +
  scale_colour_manual(values = pal_estado, name = NULL) +
  scale_shape_manual(values = c("Bacillota" = 16, "Pseudomonadota" = 17,
                                "Actinomycetota" = 15), name = "Filo") +
  scale_x_continuous(limits = c(89, 100.7), breaks = seq(90, 100, 2)) +
  scale_y_continuous(limits = c(38, 102), breaks = seq(40, 100, 10)) +
  labs(x = "ANI con la referencia GTDB mas cercana (%)",
       y = "Fraccion alineada, AF (%)") +
  theme(legend.position = "right")

guardar_fig(fig5, "fig5_novedad_taxonomica", ancho = 8.6, alto = 5.6)

# --- Cifras para el texto (no van en figura) ---------------------------------
prof <- maestra %>%
  mutate(
    rango = case_when(
      !is.na(especie) ~ "Especie",
      !is.na(genero)  ~ "Genero",
      !is.na(familia) ~ "Familia",
      TRUE            ~ "Superior"
    )
  ) %>%
  count(filo, rango) %>%
  arrange(rango, filo)

message("Rango taxonomico mas profundo asignado (para el texto):")
print(as.data.frame(prof))
print(as.data.frame(count(mutate(maestra, r = if_else(is.na(especie), "Genero", "Especie")), r)))

# Tabla de apoyo: candidatos a novedad
candidatos <- maestra %>%
  filter(is.na(especie)) %>%
  select(genoma, muestra, linaje_id, filo, familia, genero,
         referencia, ani, af, classification_method, note) %>%
  arrange(desc(is.na(ani)), ani)

write_tsv(candidatos, file.path(DIR_TAB, "T4_candidatos_novedad.tsv"))
message("  -> ", file.path(DIR_TAB, "T4_candidatos_novedad.tsv"))
