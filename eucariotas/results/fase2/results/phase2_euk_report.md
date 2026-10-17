# Fase 2 eucariota — Taxonomia y linajes de los bins eucariotas

**Job SLURM:** `57295` — **Fecha:** 2026-10-17 17:59

| Metrica | Valor |
|---|---:|
| Genomas (Fase 1 eucariota) | 13 |
| Linajes (dRep, 0.95 ANI) | 10 |
| Asignados a nivel: especie | 4 |
| Asignados a nivel: familia | 9 |
| Con ITS recuperado | 1 |
| Conflictos de genero skani vs ITS | 1 |

## Linajes

| Linaje | Prev. | Genomas | Muestras | Representante | Taxonomia | Nivel |
|---|---:|---:|---|---|---|---|
| E1 | 3 | 3 | 49,56,62 | bin-2-49 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | familia |
| E2 | 2 | 2 | 57,60 | bin-3-57 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | familia |
| E3 | 1 | 1 | 65 | bin-1-65 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | familia |
| E4 | 1 | 1 | 53 | bin-3-53 | Fungi;Ascomycota;Dothideomycetes;Mycosphaerellales;Mycosphaerellaceae | familia |
| E5 | 1 | 1 | 57 | bin-4-57 | Fungi;Ascomycota;Dothideomycetes;Pleosporales;Didymosphaeriaceae | familia |
| E6 | 1 | 1 | 55 | bin-1-55 | Fungi;Ascomycota;Dothideomycetes;Pleosporales;Pleosporaceae;Alternaria;Alternaria alternata | especie |
| E7 | 1 | 1 | 58 | bin-3-58 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | familia |
| E8 | 1 | 1 | 64 | bin-4-64 | Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae | familia |
| E9 | 1 | 1 | 65 | bin-3-65 | Fungi;Ascomycota;Sordariomycetes;Microascales;Microascaceae;Scopulariopsis;Scopulariopsis brevicaulis | especie |
| E10 | 1 | 1 | 65 | bin-2-65 | Fungi;Basidiomycota;Tremellomycetes;Tremellales;Cuniculitremaceae | familia |

## Detalle por genoma

| Genoma | Linaje | Asignacion | Nivel | Metodo | Vecino skani | ANI | AF | ITS UNITE (%) | ITS RefSeq (%) | Conflicto |
|---|---|---|---|---|---|---:|---:|---|---|---|
| bin-2-49 | E1 | Aspergillaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-3-56 | E1 | Penicillium olsonii | especie | ANI (skani) | Penicillium olsonii | 99.39 | 98.19 | - | - | - |
| bin-5-62 | E1 | Penicillium olsonii | especie | ANI (skani) | Penicillium olsonii | 99.32 | 96.2 | - | - | - |
| bin-2-65 | E10 | Cuniculitremaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-3-57 | E2 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-4-60 | E2 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-1-65 | E3 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-3-53 | E4 | Mycosphaerellaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-4-57 | E5 | Didymosphaeriaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-1-55 | E6 | Alternaria alternata | especie | ANI (skani) | Alternaria alternata | 98.26 | 89.45 | - | - | - |
| bin-3-58 | E7 | Aspergillaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-4-64 | E8 | Aspergillaceae | familia | colocacion EukCC | - | - | - | - | - | - |
| bin-3-65 | E9 | Scopulariopsis brevicaulis | especie | ANI (skani) | Scopulariopsis brevicaulis | 99.59 | 97.22 | - | sphaerospermum CBS (100.0) | genero skani=Scopulariopsis vs ITS=sphaerospermum |

## Reglas

1. skani: ANI >= 95 % y AF del bin >= 65 % -> especie.
2. ITS (cobertura >= 90 %), umbrales de Vu et al. 2019: especie >= 98.41, genero >= 94.3, familia >= 88.5, orden >= 81.2 % de identidad.
3. Sin lo anterior: genero del vecino mas cercano de skani (ANI < 95 %).
4. Ultimo recurso: linaje de la colocacion de EukCC (orientativo).

El 95 % de ANI es un umbral operativo: no hay un umbral de especie validado para hongos.
