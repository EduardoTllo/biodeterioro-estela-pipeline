# Figuras y tablas de la Fase 2 (OE2)

Genera las figuras y tablas de la asignacion taxonomica y la definicion de
linajes a partir de los resultados versionados en `results/fase1/` y
`results/fase2/`. Todo corre en el laptop (segundos a minutos); no requiere
Khipu.

## Uso

Desde la raiz del repo (`fase1_khipu/`):

```bash
Rscript figuras/R/run_all.R
```

Para regenerar una sola figura, carga primero el setup y los datos:

```bash
Rscript -e 'source("figuras/R/00_setup.R"); source("figuras/R/01_load_data.R"); source("figuras/R/04_fig3_arbol.R")'
```

## Requisitos

R >= 4.5 con: `tidyverse` (readr, dplyr, tidyr, stringr, tibble, purrr,
forcats, ggplot2), `ggrepel`, `patchwork`, `scales`, `gt`, `writexl`, y de
Bioconductor `ggtree`, `treeio`, `Biostrings`, ademas de `ape`.

```r
install.packages(c("tidyverse", "ggrepel", "patchwork", "scales", "gt", "writexl", "ape", "BiocManager"))
BiocManager::install(c("ggtree", "treeio", "Biostrings"))
```

## Estructura

| Script | Produce |
|---|---|
| `00_setup.R` | Tema grafico, paletas por filo, umbrales, funcion `guardar_fig()` |
| `01_load_data.R` | Objetos `maestra` (21 genomas) y `linajes` (19 linajes) |
| `02_fig1_embudo.R` | Fig. 1 — embudo 97 → 66 → 21 → 19 |
| `03_fig2_calidad.R` | Fig. 2 — completitud/contaminacion, tamano/GC, metricas por filo |
| `04_fig3_arbol.R` | Fig. 3 — arbol bac120 anotado + `.newick` |
| `05_fig4_ani.R` | Fig. 4 — matriz MASH; Tabla S2 (ANI intracluster) |
| `06_fig5_novedad.R` | Fig. 5 — ANI vs AF, novedad taxonomica; Tabla S1 |
| `07_fig6_prevalencia.R` | Fig. 6 — matriz linaje × muestra y priorizacion |
| `08_fig7_taxonomia.R` | Fig. 7 — composicion taxonomica; Fig. S1 — calidad de la asignacion |
| `09_tablas.R` | Tablas 1-3 (TSV + HTML con `gt`) y XLSX con las 6 tablas |

Salidas en `figuras/figs/` (PDF vectorial + PNG 600 dpi) y `figuras/tablas/`.
Los pies de figura en espanol estan en `figuras/pies_de_figura.md`.

## Datos de entrada

Todos los archivos que consumen los scripts viven en el repo:

- `results/fase1/phase1_tiara_bin_summary.tsv`, `phase1_selection.tsv`,
  `phase1_checkm2_quality_report.tsv`
- `results/fase2/phase2_gtdbtk_bac120_summary.tsv`, `phase2_selection.tsv`,
  `phase2_drep_{Cdb,Ndb,Mdb,Wdb}.csv`, `phase2_gtdbtk_markers_summary.tsv`,
  `phase2_gtdbtk_bac120_user_msa.fasta.gz`, `phase2_versions.txt`

## Advertencias metodologicas

1. El arbol de la Fig. 3 es de distancias (NJ), no de maxima verosimilitud, y
   solo incluye los genomas de consulta. El arbol de novo focalizado por linaje
   (query, referencias GTDB del taxon y outgroup) esta previsto para la Fase 4.
2. No existe una matriz de ANI 21x21: dRep ejecuta fastANI solo dentro de cada
   cluster primario de MASH. La Fig. 4a usa similitud MASH.
3. La calidad se reporta con los criterios de CheckM2, no con MIMAG completo.
   MIMAG exige ademas la presencia de rRNA 5S/16S/23S y al menos 18 tRNAs; esos
   conteos son extraibles de los `.gbk` de Bakta y quedan pendientes.
4. La prevalencia es por muestra y no por zona fisica de la estela; falta el
   mapeo de muestra a zona.

## Fase 3 (OE3)

Las figuras de la Fase 3 leen `results/fase3/`, que se trae de Khipu con
`bin/export_fase3.sh` y `rsync` (ver [README_fase3.md](../README_fase3.md)).

```bash
Rscript figuras/R/f3_run_all.R
```

Requieren ademas `micropan` (ley de Heaps) y `phangorn` (enraizado en el punto
medio):

```r
install.packages(c("micropan", "phangorn"))
```

| Script | Produce |
|---|---|
| `f3_00_setup.R` | Carga `00_setup.R`, las especies seleccionadas, las referencias finales y las paletas por categoria y habitat |
| `f3_01_curvas_heaps.R` | `fig_f3_curvas_acumulacion` (pan y core, solo referencias) y `tab_f3_heaps.tsv` (alpha de Heaps; < 1 = abierto) |
| `f3_02_posicion_bin.R` | `fig_f3_posicion_bin` (genes del bin por categoria; recuperacion del core vs completitud) y `tab_f3_resumen_pangenoma.tsv` |
| `f3_03_arbol_core.R` | `fig_f3_arbol_core_<especie>` (arbol ML del core con el habitat de cada referencia y el bin destacado) |
| `f3_04_exclusivos.R` | `fig_f3_exclusivos` (embudo F1-F6 y origen segun nr) y `tab_f3_embudo_exclusivos.tsv` |
