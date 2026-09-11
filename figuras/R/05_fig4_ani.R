# =============================================================================
# 05_fig4_ani.R — Figura 4: matriz de distancias genomicas por pares
#
# IMPORTANTE: dRep solo ejecuta fastANI DENTRO de cada cluster primario (MASH),
# por lo que NO existe una matriz ANI 21x21 completa. La matriz completa
# disponible es la de similitud MASH (todos contra todos). Los valores exactos
# de fastANI de las dos parejas comparadas se exportan a una tabla y se
# comentan en el texto, no en una figura propia (son solo dos numeros).
# =============================================================================

# --- Matriz MASH todos contra todos ------------------------------------------
mash <- drep_mdb %>%
  transmute(g1 = genome1, g2 = genome2, similitud = 100 * (1 - dist))

orden_genomas <- maestra %>%
  arrange(filo, genero, genoma) %>%
  pull(genoma)

mash_mat <- mash %>%
  filter(g1 %in% orden_genomas, g2 %in% orden_genomas) %>%
  pivot_wider(names_from = g2, values_from = similitud) %>%
  column_to_rownames("g1") %>%
  as.matrix()
mash_mat <- mash_mat[orden_genomas, orden_genomas]

etiquetas <- maestra %>%
  transmute(genoma, etq = paste0(genoma, " | ", taxon_corto))

mash_largo <- as.data.frame(mash_mat) %>%
  rownames_to_column("g1") %>%
  pivot_longer(-g1, names_to = "g2", values_to = "similitud") %>%
  left_join(rename(etiquetas, g1 = genoma, etq1 = etq), by = "g1") %>%
  left_join(rename(etiquetas, g2 = genoma, etq2 = etq), by = "g2") %>%
  mutate(etq1 = factor(etq1, levels = etiquetas$etq[match(orden_genomas, etiquetas$genoma)]),
         etq2 = factor(etq2, levels = rev(levels(etq1))))

fig4 <- ggplot(mash_largo, aes(x = etq1, y = etq2, fill = similitud)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(data = filter(mash_largo, similitud >= 80),
            aes(label = sprintf("%.0f", similitud),
                colour = similitud >= 92),
            size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    name = "Similitud\nMASH (%)", limits = c(76, 100), oob = squish,
    breaks = c(80, 85, 90, 95, 100)
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6.8),
        axis.text.y = element_text(size = 6.8),
        panel.grid = element_blank())

guardar_fig(fig4, "fig4_matriz_mash", ancho = 8.6, alto = 7.4)

# --- Valores de fastANI para el texto (no van en figura) ---------------------
# fastANI reporta cada pareja en los dos sentidos -> promediamos por pareja.
ani_intracluster <- drep_ndb %>%
  filter(reference != querry) %>%
  mutate(a = pmin(reference, querry), b = pmax(reference, querry)) %>%
  group_by(a, b) %>%
  summarise(ani = 100 * mean(ani),
            cobertura = 100 * mean(alignment_coverage), .groups = "drop") %>%
  left_join(select(maestra, a = genoma, linaje_id, taxon_corto, muestra),
            by = "a") %>%
  transmute(Linaje = linaje_id, `Taxon (GTDB)` = taxon_corto,
            `Genoma A` = a, `Genoma B` = b,
            `ANI (%)` = round(ani, 2),
            `Cobertura del alineamiento (%)` = round(cobertura, 1),
            `Supera 95 % ANI` = if_else(ani >= UMBRAL_ANI, "si", "no"))

write_tsv(ani_intracluster, file.path(DIR_TAB, "T5_ani_intracluster.tsv"))
message("  -> ", file.path(DIR_TAB, "T5_ani_intracluster.tsv"))
message("ANI intracluster (para el texto):")
print(as.data.frame(ani_intracluster))
