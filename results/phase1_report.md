# Fase 1 — Control de calidad y seleccion de genomas (OE1)

**Proyecto:** Modelo predictivo de riesgo de deterioro de la Estela de Raimondi

**Cluster:** Khipu (UTEC) — job SLURM `49286`

**Fecha de ejecucion:** 2026-08-11 11:51:36

---

## 1. Resumen del proceso

La Fase 1 aplica dos filtros secuenciales sobre los bins genomicos crudos entregados por el centro de secuenciacion, para quedarnos con genomas procariotas de calidad suficiente para el analisis pangenomico posterior.

| Etapa | Herramienta | Criterio | Entran | Pasan |
|---|---|---|---:|---:|
| 1. Filtrado por dominio | Tiara | contig >= 3000 bp; se descarta el bin si su longitud mayoritaria es eucariota u organelar | 97 | 66 |
| 2. Calidad genomica | CheckM2 | completitud > 70% y contaminacion < 5% | 66 | 21 |

**Resultado final: 21 de 97 genomas seleccionados** para la Fase 2.

## 2. Genomas seleccionados (OE1)

| Bin | Contigs | Long (bp) | Dom. mayoritario | Completitud (%) | Contaminacion (%) | Decision | Motivo |
|---|---:|---:|---|---:|---:|---|---|
| bin-1-51 | 51 | 4727435 | prok | 100.00 | 0.00 | SELECCIONADO | - |
| bin-1-68 | 50 | 2465224 | prok | 100.00 | 1.79 | SELECCIONADO | - |
| bin-2-58 | 17 | 4771531 | prok | 100.00 | 0.06 | SELECCIONADO | - |
| bin-6-50 | 591 | 5543214 | prok | 100.00 | 0.65 | SELECCIONADO | - |
| bin-6-61 | 30 | 4792746 | prok | 100.00 | 0.74 | SELECCIONADO | - |
| bin-1-57 | 22 | 4032412 | prok | 99.99 | 3.23 | SELECCIONADO | - |
| bin-2-53 | 39 | 5164804 | prok | 99.99 | 0.57 | SELECCIONADO | - |
| bin-2-69 | 61 | 2481071 | prok | 99.99 | 0.18 | SELECCIONADO | - |
| bin-1-59 | 98 | 5907869 | prok | 99.97 | 3.07 | SELECCIONADO | - |
| bin-2-57 | 28 | 5172263 | prok | 99.96 | 1.01 | SELECCIONADO | - |
| bin-4-52 | 96 | 3497259 | prok | 99.86 | 0.25 | SELECCIONADO | - |
| bin-1-53 | 30 | 2476631 | prok | 99.84 | 0.12 | SELECCIONADO | - |
| bin-1-64 | 163 | 3944265 | prok | 99.84 | 0.59 | SELECCIONADO | - |
| bin-3-52 | 204 | 5184467 | prok | 98.97 | 0.58 | SELECCIONADO | - |
| bin-2-63 | 150 | 4282638 | prok | 92.39 | 0.34 | SELECCIONADO | - |
| bin-1-54 | 49 | 4312038 | prok | 92.25 | 0.01 | SELECCIONADO | - |
| bin-9-71 | 346 | 1832072 | prok | 85.09 | 2.13 | SELECCIONADO | - |
| bin-5-60 | 97 | 2849856 | prok | 77.53 | 0.10 | SELECCIONADO | - |
| bin-5-50 | 1041 | 3472002 | prok | 76.56 | 0.44 | SELECCIONADO | - |
| bin-5-63 | 90 | 2758939 | prok | 72.46 | 0.05 | SELECCIONADO | - |
| bin-3-64 | 459 | 2352198 | prok | 70.14 | 0.64 | SELECCIONADO | - |

## 3. Detalle completo (todos los bins)

| Bin | Contigs | Long (bp) | Dom. mayoritario | Completitud (%) | Contaminacion (%) | Decision | Motivo |
|---|---:|---:|---|---:|---:|---|---|
| bin-1-49 | 87 | 545053 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-1-50 | 224 | 2609893 | prok | 62.29 | 2.66 | RECHAZADO | completitud 62.29<=70 |
| bin-1-51 | 51 | 4727435 | prok | 100.00 | 0.00 | SELECCIONADO | - |
| bin-1-52 | 26 | 415732 | prok | 6.05 | 0.20 | RECHAZADO | completitud 6.05<=70 |
| bin-1-53 | 30 | 2476631 | prok | 99.84 | 0.12 | SELECCIONADO | - |
| bin-1-54 | 49 | 4312038 | prok | 92.25 | 0.01 | SELECCIONADO | - |
| bin-1-55 | 247 | 34242081 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-1-56 | 12 | 230366 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-1-57 | 22 | 4032412 | prok | 99.99 | 3.23 | SELECCIONADO | - |
| bin-1-58 | 4725 | 12125594 | unknown | 0.00 | 17.71 | RECHAZADO | completitud 0.00<=70; contaminacion 17.71>=5 |
| bin-1-59 | 98 | 5907869 | prok | 99.97 | 3.07 | SELECCIONADO | - |
| bin-1-60 | 302 | 8273570 | prok | 100.00 | 24.87 | RECHAZADO | contaminacion 24.87>=5 |
| bin-1-61 | 448 | 40131343 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-1-62 | 320 | 7390811 | prok | 99.72 | 72.67 | RECHAZADO | contaminacion 72.67>=5 |
| bin-1-63 | 1626 | 5185746 | prok | 84.85 | 39.63 | RECHAZADO | contaminacion 39.63>=5 |
| bin-1-64 | 163 | 3944265 | prok | 99.84 | 0.59 | SELECCIONADO | - |
| bin-1-65 | 322 | 29636962 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-1-68 | 50 | 2465224 | prok | 100.00 | 1.79 | SELECCIONADO | - |
| bin-1-69 | 67 | 251086 | prok | 6.60 | 0.14 | RECHAZADO | completitud 6.60<=70 |
| bin-1-70 | 38054 | 94251631 | unknown | 31.70 | 28.68 | RECHAZADO | completitud 31.70<=70; contaminacion 28.68>=5 |
| bin-1-71 | 371 | 1498934 | prok | 67.71 | 15.05 | RECHAZADO | completitud 67.71<=70; contaminacion 15.05>=5 |
| bin-1-72 | 167 | 408652 | unknown | 5.87 | 0.19 | RECHAZADO | completitud 5.87<=70 |
| bin-10-71 | 506 | 2845557 | prok | 97.57 | 30.13 | RECHAZADO | contaminacion 30.13>=5 |
| bin-11-71 | 41610 | 101946556 | unknown | 65.68 | 17.72 | RECHAZADO | completitud 65.68<=70; contaminacion 17.72>=5 |
| bin-12-71 | 383 | 1015904 | unknown | 43.42 | 3.52 | RECHAZADO | completitud 43.42<=70 |
| bin-2-49 | 134 | 27557437 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-50 | 464 | 53406620 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-52 | 270 | 5347688 | prok | 99.93 | 8.10 | RECHAZADO | contaminacion 8.10>=5 |
| bin-2-53 | 39 | 5164804 | prok | 99.99 | 0.57 | SELECCIONADO | - |
| bin-2-54 | 106 | 7273841 | prok | 100.00 | 15.52 | RECHAZADO | contaminacion 15.52>=5 |
| bin-2-56 | 816 | 49647189 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-57 | 28 | 5172263 | prok | 99.96 | 1.01 | SELECCIONADO | - |
| bin-2-58 | 17 | 4771531 | prok | 100.00 | 0.06 | SELECCIONADO | - |
| bin-2-59 | 424 | 23212147 | prok | 91.94 | 115.89 | RECHAZADO | contaminacion 115.89>=5 |
| bin-2-60 | 290 | 10390751 | prok | 100.00 | 100.33 | RECHAZADO | contaminacion 100.33>=5 |
| bin-2-61 | 9 | 614253 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-62 | 5 | 227679 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-63 | 150 | 4282638 | prok | 92.39 | 0.34 | SELECCIONADO | - |
| bin-2-64 | 9242 | 26083286 | unknown | 0.02 | 16.03 | RECHAZADO | completitud 0.02<=70; contaminacion 16.03>=5 |
| bin-2-65 | 141 | 17056380 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-2-68 | 70 | 246182 | prok | 4.99 | 0.03 | RECHAZADO | completitud 4.99<=70 |
| bin-2-69 | 61 | 2481071 | prok | 99.99 | 0.18 | SELECCIONADO | - |
| bin-2-70 | 87 | 278635 | unknown | 5.38 | 0.09 | RECHAZADO | completitud 5.38<=70 |
| bin-2-71 | 618 | 1368774 | unknown | 56.83 | 26.01 | RECHAZADO | completitud 56.83<=70; contaminacion 26.01>=5 |
| bin-2-72 | 667 | 2485725 | prok | 79.16 | 54.30 | RECHAZADO | contaminacion 54.30>=5 |
| bin-3-50 | 302 | 607902 | unknown | 13.03 | 0.08 | RECHAZADO | completitud 13.03<=70 |
| bin-3-52 | 204 | 5184467 | prok | 98.97 | 0.58 | SELECCIONADO | - |
| bin-3-53 | 301 | 29244092 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-56 | 160 | 26628065 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-57 | 818 | 25745728 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-58 | 265 | 29647785 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-60 | 959 | 2321293 | unknown | 39.27 | 0.00 | RECHAZADO | completitud 39.27<=70 |
| bin-3-61 | 50 | 10715439 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-62 | 117 | 449659 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-63 | 4589 | 19740558 | prok | 99.91 | 102.18 | RECHAZADO | contaminacion 102.18>=5 |
| bin-3-64 | 459 | 2352198 | prok | 70.14 | 0.64 | SELECCIONADO | - |
| bin-3-65 | 5256 | 27711612 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-3-70 | 1696 | 4261799 | unknown | 77.48 | 75.48 | RECHAZADO | contaminacion 75.48>=5 |
| bin-3-71 | 3657 | 12324162 | prok | 97.39 | 70.70 | RECHAZADO | contaminacion 70.70>=5 |
| bin-3-72 | 1741 | 6597747 | prok | 98.53 | 98.35 | RECHAZADO | contaminacion 98.35>=5 |
| bin-4-50 | 97 | 1654013 | prok | 43.00 | 0.23 | RECHAZADO | completitud 43.00<=70 |
| bin-4-52 | 96 | 3497259 | prok | 99.86 | 0.25 | SELECCIONADO | - |
| bin-4-56 | 122 | 1739035 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-4-57 | 977 | 39075173 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-4-60 | 4482 | 25724799 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-4-61 | 1610 | 11012084 | prok | 100.00 | 58.51 | RECHAZADO | contaminacion 58.51>=5 |
| bin-4-62 | 2707 | 9759046 | prok | 90.25 | 65.70 | RECHAZADO | contaminacion 65.70>=5 |
| bin-4-63 | 856 | 56508797 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-4-64 | 428 | 29693972 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-4-65 | 2506 | 20150702 | prok | 99.99 | 106.28 | RECHAZADO | contaminacion 106.28>=5 |
| bin-4-70 | 484 | 1565460 | prok | 64.00 | 2.55 | RECHAZADO | completitud 64.00<=70 |
| bin-4-71 | 216 | 488562 | unknown | 32.24 | 2.32 | RECHAZADO | completitud 32.24<=70 |
| bin-4-72 | 72 | 248801 | prok | 4.69 | 0.01 | RECHAZADO | completitud 4.69<=70 |
| bin-5-50 | 1041 | 3472002 | prok | 76.56 | 0.44 | SELECCIONADO | - |
| bin-5-52 | 1099 | 3740736 | prok | 95.51 | 12.95 | RECHAZADO | contaminacion 12.95>=5 |
| bin-5-56 | 146 | 552204 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-5-57 | 12 | 224138 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-5-60 | 97 | 2849856 | prok | 77.53 | 0.10 | SELECCIONADO | - |
| bin-5-61 | 273 | 5811268 | prok | 97.82 | 10.64 | RECHAZADO | contaminacion 10.64>=5 |
| bin-5-62 | 222 | 27886446 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-5-63 | 90 | 2758939 | prok | 72.46 | 0.05 | SELECCIONADO | - |
| bin-5-64 | 2173 | 9450268 | prok | 100.00 | 67.99 | RECHAZADO | contaminacion 67.99>=5 |
| bin-5-65 | 24 | 222213 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-5-70 | 282 | 628283 | unknown | 48.48 | 2.40 | RECHAZADO | completitud 48.48<=70 |
| bin-5-71 | 118 | 438388 | prok | 12.74 | 0.54 | RECHAZADO | completitud 12.74<=70 |
| bin-6-50 | 591 | 5543214 | prok | 100.00 | 0.65 | SELECCIONADO | - |
| bin-6-56 | 116 | 9578458 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-6-57 | 126 | 448524 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-6-61 | 30 | 4792746 | prok | 100.00 | 0.74 | SELECCIONADO | - |
| bin-6-63 | 21 | 367167 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-6-64 | 71 | 239057 | euk | - | - | RECHAZADO | dominio mayoritario: euk |
| bin-6-70 | 3194 | 8190971 | unknown | 68.24 | 17.99 | RECHAZADO | completitud 68.24<=70; contaminacion 17.99>=5 |
| bin-6-71 | 34 | 453921 | prok | 11.79 | 0.16 | RECHAZADO | completitud 11.79<=70 |
| bin-7-64 | 1246 | 3769741 | prok | 60.17 | 0.22 | RECHAZADO | completitud 60.17<=70 |
| bin-7-71 | 95 | 258968 | unknown | 7.21 | 0.05 | RECHAZADO | completitud 7.21<=70 |
| bin-8-71 | 908 | 2450575 | unknown | 81.92 | 35.87 | RECHAZADO | contaminacion 35.87>=5 |
| bin-9-71 | 346 | 1832072 | prok | 85.09 | 2.13 | SELECCIONADO | - |

## 4. Parametros y trazabilidad

- Longitud minima de contig (Tiara): **3000 bp**
- Umbral de completitud: **> 70 %**
- Umbral de contaminacion: **< 5 %**
- Tabla completa en formato TSV: `phase1_selection.tsv`
- Clasificacion por contig (Tiara): `01_tiara/`
- Reporte de calidad (CheckM2): `02_checkm2/quality_report.tsv`
- Genomas seleccionados (FASTA): `phase1_selected_genomes/`

> Nota: el umbral (>70%%, <5%%) no corresponde a una categoria MIMAG estandar; es un compromiso propio del proyecto que flexibiliza la completitud pero mantiene la contaminacion estricta para no sesgar la reconstruccion del genoma accesorio (tesis, seccion 4.1.3.2).
