# Fase 2 eucariota — Taxonomia y linajes de los bins eucariotas

**Job SLURM:** `57311` — **Fecha:** 2026-10-18 01:44

| Metrica | Valor |
|---|---:|
| Genomas (Fase 1 eucariota) | 13 |
| Linajes (dRep, 0.95 ANI) | 10 |
| Asignados a nivel: especie | 9 |
| Asignados a nivel: familia | 3 |
| Asignados a nivel: orden (consenso skani/EukCC) | 1 |
| Con ITS recuperado | 1 |
| Conflictos de genero skani vs ITS | 1 |
| Especie por ANI con otras especies sobre el umbral | 2 |

## Linajes

| Linaje | Prev. | Genomas | Muestras | Representante | Taxonomia | Nivel |
|---|---:|---:|---|---|---|---|
| E1 | 3 | 3 | 49,56,62 | bin-2-49 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium olsonii | especie |
| E2 | 2 | 2 | 57,60 | bin-3-57 | Fungi;Ascomycota;Dothideomycetes;Cladosporiales;Cladosporiaceae;Cladosporium;Cladosporium dominicanum | especie |
| E3 | 1 | 1 | 65 | bin-1-65 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | familia |
| E4 | 1 | 1 | 53 | bin-3-53 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | familia |
| E5 | 1 | 1 | 57 | bin-4-57 | Fungi;Ascomycota;Dothideomycetes;Pleosporales | orden (consenso skani/EukCC) |
| E6 | 1 | 1 | 55 | bin-1-55 | Fungi;Ascomycota;Dothideomycetes;Pleosporales;Pleosporaceae;Alternaria;Alternaria alternata | especie |
| E7 | 1 | 1 | 64 | bin-4-64 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | familia |
| E8 | 1 | 1 | 58 | bin-3-58 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium coffeae | especie |
| E9 | 1 | 1 | 65 | bin-3-65 | Fungi;Ascomycota;Sordariomycetes;Microascales;Microascaceae;Scopulariopsis;Scopulariopsis brevicaulis | especie |
| E10 | 1 | 1 | 65 | bin-2-65 | Fungi;Basidiomycota;Tremellomycetes;Tremellales;Rhynchogastremaceae;Papiliotrema;Papiliotrema fonsecae | especie |

## Detalle por genoma

| Genoma | Linaje | Asignacion | Nivel | Metodo | Vecino skani | ANI | AF | ITS UNITE (%) | ITS RefSeq (%) | Conflicto | Otras especies >= umbral |
|---|---|---|---|---|---|---:|---:|---|---|---|---|
| bin-2-49 | E1 | Penicillium olsonii | especie | ANI (skani) | Penicillium olsonii | 97.65 | 94.8 | - | - | - | - |
| bin-3-56 | E1 | Penicillium olsonii | especie | ANI (skani) | Penicillium olsonii | 99.39 | 98.19 | - | - | - | - |
| bin-5-62 | E1 | Penicillium olsonii | especie | ANI (skani) | Penicillium olsonii | 99.32 | 96.2 | - | - | - | - |
| bin-2-65 | E10 | Papiliotrema fonsecae | especie | ANI (skani) | Papiliotrema fonsecae | 96.82 | 95.7 | - | - | - | - |
| bin-3-57 | E2 | Cladosporium dominicanum | especie | ANI (skani) | Cladosporium dominicanum | 96.32 | 89.54 | - | - | - | - |
| bin-4-60 | E2 | Cladosporium dominicanum | especie | ANI (skani) | Cladosporium dominicanum | 96.12 | 82.8 | - | - | - | - |
| bin-1-65 | E3 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - | - |
| bin-3-53 | E4 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - | - |
| bin-4-57 | E5 | Pleosporales | orden (consenso skani/EukCC) | vecino skani Periconia (ANI 83.30 %, AF 15.9 % < 50) y EukCC coinciden hasta orden | Periconia digitata | 83.3 | 15.94 | - | - | - | - |
| bin-1-55 | E6 | Alternaria alternata | especie | ANI (skani) | Alternaria alternata | 98.26 | 89.45 | - | - | - | Alternaria alstroemeriae; Alternaria alternarina; Alternaria angustiovoidea; Alternaria arborescens; Alternaria brassicae; Alternaria burnsii; Alternaria gaisen; Alternaria longipes; Alternaria porri; Alternaria postmessia; Alternaria tenuissima |
| bin-4-64 | E7 | Aspergillaceae | familia | colocacion EukCC | - | - | - | - | - | - | - |
| bin-3-58 | E8 | Penicillium coffeae | especie | ANI (skani) | Penicillium coffeae | 96.46 | 90.7 | - | - | - | - |
| bin-3-65 | E9 | Scopulariopsis brevicaulis | especie | ANI (skani) | Scopulariopsis brevicaulis | 99.59 | 97.22 | - | Cladosporium sphaerospermum (100.0) | genero skani=Scopulariopsis vs ITS=Cladosporium | Scopulariopsis alboflavescens |

## Reglas

1. skani: ANI >= 95 % y AF del bin >= 65 % -> especie.
2. ITS (cobertura >= 90 %), umbrales de Vu et al. 2019: especie >= 98.41, genero >= 94.3, familia >= 88.5, orden >= 81.2 % de identidad.
3. Sin lo anterior: genero del vecino mas cercano de skani (ANI < 95 %) si alinea >= 50 % del bin; si alinea menos, el rango mas profundo que comparten ese vecino y EukCC.
4. Ultimo recurso: linaje de la colocacion de EukCC (orientativo).

El 95 % de ANI es un umbral operativo: no hay un umbral de especie validado para hongos.

> Resumen recalculado en la Mac el 18-10-2026 con `summarize_euk.py` (regla: genero por vecino skani solo con AF >= 50 %). Mismos datos del job 57311; solo cambia bin-4-57.
