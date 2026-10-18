# =============================================================================
# euk_05_fig4_mash.R — Fig. E4: similitud Mash 13 x 13 (dRep) y tabla de
# fastANI dentro de los clusteres (como la Fig. 7.6)
# =============================================================================

orden <- maestra %>% arrange(clase, orden, familia, genero, genome) %>% pull(genome)
etq <- maestra %>% transmute(genome, etq = paste0(genome, " | ", taxon_corto))

mash_largo <- drep_mdb %>%
  transmute(g1 = genome1, g2 = genome2, similitud = 100 * (1 - dist)) %>%
  filter(g1 %in% orden, g2 %in% orden) %>%
  left_join(rename(etq, g1 = genome, etq1 = etq), by = "g1") %>%
  left_join(rename(etq, g2 = genome, etq2 = etq), by = "g2") %>%
  mutate(etq1 = factor(etq1, levels = etq$etq[match(orden, etq$genome)]),
         etq2 = factor(etq2, levels = rev(levels(etq1))))

fig4 <- ggplot(mash_largo, aes(x = etq1, y = etq2, fill = similitud)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(data = filter(mash_largo, similitud >= 80),
            aes(label = sprintf("%.0f", similitud), colour = similitud >= 92),
            size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       name = "Similitud\nMash (%)", limits = c(76, 100), oob = squish,
                       breaks = c(80, 85, 90, 95, 100)) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7), panel.grid = element_blank())

guardar_fig(fig4, "figE4_matriz_mash", ancho = 8.2, alto = 7)

ani_intra <- drep_ndb %>%
  filter(reference != querry) %>%
  mutate(a = pmin(reference, querry), b = pmax(reference, querry)) %>%
  group_by(a, b) %>%
  summarise(ani = 100 * mean(ani), cobertura = 100 * mean(alignment_coverage), .groups = "drop") %>%
  left_join(select(maestra, a = genome, linaje_id, taxon_corto), by = "a") %>%
  transmute(`Linaje de A` = linaje_id, `Taxon de A` = taxon_corto, `Genoma A` = a, `Genoma B` = b,
            `ANI (%)` = round(ani, 2), `Cobertura (%)` = round(cobertura, 1),
            `Mismo linaje (>= 95 %)` = if_else(ani >= UMBRAL_ANI, "si", "no"))
write_tsv(ani_intra, file.path(DIR_TAB, "TE4_ani_intracluster.tsv"))
print(as.data.frame(ani_intra))
