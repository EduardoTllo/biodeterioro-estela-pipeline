# Fase 1 eucariota — Control de calidad de los bins eucariotas

**Job SLURM:** `57148` — **Fecha:** 2026-10-17 15:51

Entrada: bins con mayoria de longitud eucariota segun Tiara (Fase 1). Calidad con EukCC (completitud > 70 % y contaminacion < 5 %). BUSCO como medida complementaria (no decide).

| Metrica | Valor |
|---|---:|
| Bins eucariotas (Tiara) | 31 |
| Sin estimacion de EukCC | 12 |
| Rechazo solo por completitud | 2 |
| Rechazo solo por contaminacion | 4 |
| Rechazo por ambos | 0 |
| **Seleccionados** | **13** |

## Detalle por bin

| Bin | Mb | Contigs | N50 | EukCC compl. | EukCC cont. | BUSCO fungi C/D | Linaje EukCC | Decision |
|---|---:|---:|---:|---:|---:|---|---|---|
| bin-1-55 | 34.24 | 247 | 403592 | 100.0 | 0.0 | 98.9/0.0 | Fungi;Ascomycota;Dothideomycetes;Pleosporales;Pleosporaceae;Alternaria | SELECCIONADO |
| bin-1-65 | 29.64 | 322 | 155267 | 99.25 | 2.26 | 99.4/0.8 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | SELECCIONADO |
| bin-2-49 | 27.56 | 134 | 337091 | 99.62 | 0.0 | 98.8/0.7 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | SELECCIONADO |
| bin-2-65 | 17.06 | 141 | 257073 | 91.8 | 0.2 | 94.9/0.2 | Fungi;Basidiomycota;Tremellomycetes;Tremellales;Cuniculitremaceae | SELECCIONADO |
| bin-3-53 | 29.24 | 301 | 172947 | 99.62 | 0.0 | 99.4/0.0 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | SELECCIONADO |
| bin-3-56 | 26.63 | 160 | 275458 | 98.5 | 0.0 | 97.0/0.5 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | SELECCIONADO |
| bin-3-57 | 25.75 | 818 | 76707 | 98.5 | 1.5 | 97.4/0.5 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | SELECCIONADO |
| bin-3-58 | 29.65 | 265 | 493691 | 98.12 | 0.0 | 97.7/0.7 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | SELECCIONADO |
| bin-3-65 | 27.71 | 5256 | 6529 | 83.46 | 0.38 | 81.5/0.5 | Fungi;Ascomycota;Sordariomycetes;Microascales;Microascaceae | SELECCIONADO |
| bin-4-57 | 39.08 | 977 | 84533 | 98.5 | 1.5 | 98.4/0.9 | Fungi;Ascomycota;Dothideomycetes;Pleosporales;Didymosphaeriaceae | SELECCIONADO |
| bin-4-60 | 25.72 | 4482 | 7360 | 89.47 | 1.13 | 84.8/0.1 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | SELECCIONADO |
| bin-4-64 | 29.69 | 428 | 164221 | 100.0 | 0.38 | 99.1/0.6 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | SELECCIONADO |
| bin-5-62 | 27.89 | 222 | 456475 | 99.25 | 0.0 | 97.3/0.8 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | SELECCIONADO |
| bin-1-49 | 0.55 | 87 | 9806 |  |  | 0.0/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-1-56 | 0.23 | 12 | 22274 |  |  | 0.1/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-1-61 | 40.13 | 448 | 148710 | 99.62 | 51.13 | 99.1/32.2 | Fungi;Ascomycota | contaminacion 51.13>=5 |
| bin-2-50 | 53.41 | 464 | 190157 | 100.0 | 96.99 | 100.0/61.0 | Fungi;Ascomycota | contaminacion 96.99>=5 |
| bin-2-56 | 49.65 | 816 | 108462 | 98.5 | 56.39 | 98.9/43.5 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | contaminacion 56.39>=5 |
| bin-2-61 | 0.61 | 9 | 124060 |  |  | 1.4/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-2-62 | 0.23 | 5 | 136727 |  |  | 1.2/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-3-61 | 10.72 | 50 | 315453 | 46.99 | 0.0 | 43.3/0.0 | Fungi;Ascomycota;Sordariomycetes;Hypocreales | completitud 46.99<=70 |
| bin-3-62 | 0.45 | 117 | 4403 |  |  | 0.0/0.0 | - | sin salida de EukCC |
| bin-4-56 | 1.74 | 122 | 29480 |  |  | 1.6/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-4-63 | 56.51 | 856 | 112022 | 100.0 | 99.62 | 99.6/70.9 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | contaminacion 99.62>=5 |
| bin-5-56 | 0.55 | 146 | 4220 |  |  | 0.0/0.0 | - | sin salida de EukCC |
| bin-5-57 | 0.22 | 12 | 29714 |  |  | 0.1/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-5-65 | 0.22 | 24 | 13496 |  |  | 0.0/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-6-56 | 9.58 | 116 | 132098 | 42.86 | 0.75 | 38.7/0.2 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Teratosphaeriaceae | completitud 42.86<=70 |
| bin-6-57 | 0.45 | 126 | 3993 |  |  | 0.0/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-6-63 | 0.37 | 21 | 26065 |  |  | 0.5/0.0 | - | EukCC sin estimacion (sin set de marcadores) |
| bin-6-64 | 0.24 | 71 | 3755 |  |  | 0.0/0.0 | - | EukCC sin estimacion (sin set de marcadores) |

Tabla completa: `phase1_euk_selection.tsv`
