# Fase 2 — Asignacion taxonomica y definicion de linajes (OE2)

**Proyecto:** Modelo predictivo de riesgo de deterioro de la Estela de Raimondi

**Cluster:** Khipu (UTEC) — job SLURM `49316`

**Fecha de ejecucion:** 2026-08-12 00:15:08

---

## 1. Resumen

Partiendo de los 21 genomas seleccionados en la Fase 1, se asigno taxonomia con GTDB-Tk (datos R220), se confirmo la especie por ANI y se agruparon los genomas redundantes en **linajes** (clusteres de especie a 95 % de ANI con dRep). Cada linaje se prioriza por su **prevalencia espacial** (nº de muestras distintas que lo aportan).

| Metrica | Valor |
|---|---:|
| Genomas de entrada (Fase 1) | 21 |
| Linajes definidos (clusteres 95% ANI) | 19 |
| Genomas con especie confirmada (ANI>=95%, AF>=65%) | 16 |
| Linajes priorizados para la Fase 3 | 3 |

## 2. Linajes priorizados (top-3 por prevalencia)

| Linaje | Prevalencia (nº muestras) | Genomas | Muestras | Representante | Taxonomia (GTDB) | Fase 3 |
|---|---:|---:|---|---|---|---|
| L1 | 2 | 2 | 50,52 | bin-6-50 | g__Telluria s__Telluria timonae | SI |
| L2 | 2 | 2 | 68,69 | bin-2-69 | f__Mycobacteriaceae g__Corynebacterium | SI |
| L3 | 1 | 1 | 61 | bin-6-61 | g__Bacillus_AB s__Bacillus_AB infantis | SI |
| L4 | 1 | 1 | 64 | bin-3-64 | f__Planococcaceae g__Sporosarcina | no |
| L5 | 1 | 1 | 51 | bin-1-51 | f__Enterobacteriaceae g__Pantoea | no |
| L6 | 1 | 1 | 57 | bin-1-57 | g__Brevibacterium s__Brevibacterium sp021023135 | no |
| L7 | 1 | 1 | 59 | bin-1-59 | g__Paenibacillus s__Paenibacillus sp001955535 | no |
| L8 | 1 | 1 | 52 | bin-4-52 | g__Acinetobacter s__Acinetobacter schindleri | no |
| L9 | 1 | 1 | 60 | bin-5-60 | g__Priestia s__Priestia flexa | no |
| L10 | 1 | 1 | 71 | bin-9-71 | f__Acetobacteraceae g__Aristophania | no |
| L11 | 1 | 1 | 53 | bin-1-53 | g__Staphylococcus s__Staphylococcus warneri_A | no |
| L12 | 1 | 1 | 54 | bin-1-54 | g__Peribacillus s__Peribacillus frigoritolerans | no |
| L13 | 1 | 1 | 63 | bin-2-63 | g__Rossellomorea s__Rossellomorea marisflavi | no |
| L14 | 1 | 1 | 64 | bin-1-64 | g__Bacillus s__Bacillus licheniformis | no |
| L15 | 1 | 1 | 53 | bin-2-53 | g__Metabacillus s__Metabacillus halosaccharovorans | no |
| L16 | 1 | 1 | 63 | bin-5-63 | g__Bacillus s__Bacillus altitudinis | no |
| L17 | 1 | 1 | 57 | bin-2-57 | g__Cytobacillus s__Cytobacillus oceanisediminis | no |
| L18 | 1 | 1 | 50 | bin-5-50 | g__Cytobacillus s__Cytobacillus firmus_B | no |
| L19 | 1 | 1 | 58 | bin-2-58 | g__Mesobacillus s__Mesobacillus sp014856545 | no |

## 3. Detalle por genoma

| Genoma | Muestra | Linaje | Clasificacion GTDB | ANI (%) | AF (%) | Especie confirmada |
|---|---|---|---|---:|---:|---|
| bin-3-52 | 52 | L1 | g__Telluria s__Telluria timonae | 99.96 | 99.2 | si |
| bin-6-50 | 50 | L1 | g__Telluria s__Telluria timonae | 99.92 | 95.3 | si |
| bin-9-71 | 71 | L10 | f__Acetobacteraceae g__Aristophania | - | - | NA |
| bin-1-53 | 53 | L11 | g__Staphylococcus s__Staphylococcus warneri_A | 99.02 | 90.1 | si |
| bin-1-54 | 54 | L12 | g__Peribacillus s__Peribacillus frigoritolerans | 97.2 | 87.9 | si |
| bin-2-63 | 63 | L13 | g__Rossellomorea s__Rossellomorea marisflavi | 98.53 | 92.6 | si |
| bin-1-64 | 64 | L14 | g__Bacillus s__Bacillus licheniformis | 99.75 | 98.5 | si |
| bin-2-53 | 53 | L15 | g__Metabacillus s__Metabacillus halosaccharovorans | 98.9 | 92.0 | si |
| bin-5-63 | 63 | L16 | g__Bacillus s__Bacillus altitudinis | 98.19 | 93.7 | si |
| bin-2-57 | 57 | L17 | g__Cytobacillus s__Cytobacillus oceanisediminis | 95.61 | 79.7 | si |
| bin-5-50 | 50 | L18 | g__Cytobacillus s__Cytobacillus firmus_B | 98.08 | 84.7 | si |
| bin-2-58 | 58 | L19 | g__Mesobacillus s__Mesobacillus sp014856545 | 98.01 | 85.8 | si |
| bin-1-68 | 68 | L2 | f__Mycobacteriaceae g__Corynebacterium | 90.55 | 44.3 | no |
| bin-2-69 | 69 | L2 | f__Mycobacteriaceae g__Corynebacterium | 90.5 | 44.0 | no |
| bin-6-61 | 61 | L3 | g__Bacillus_AB s__Bacillus_AB infantis | 98.22 | 89.6 | si |
| bin-3-64 | 64 | L4 | f__Planococcaceae g__Sporosarcina | 91.25 | 60.4 | no |
| bin-1-51 | 51 | L5 | f__Enterobacteriaceae g__Pantoea | - | - | NA |
| bin-1-57 | 57 | L6 | g__Brevibacterium s__Brevibacterium sp021023135 | 98.31 | 91.5 | si |
| bin-1-59 | 59 | L7 | g__Paenibacillus s__Paenibacillus sp001955535 | 99.99 | 94.6 | si |
| bin-4-52 | 52 | L8 | g__Acinetobacter s__Acinetobacter schindleri | 98.03 | 80.7 | si |
| bin-5-60 | 60 | L9 | g__Priestia s__Priestia flexa | 99.44 | 92.9 | si |

## 4. Parametros y trazabilidad

- Clasificacion: GTDB-Tk `classify_wf` (marcadores bac120/ar53, datos R220).
- Confirmacion de especie: ANI >= 95 % y AF >= 65 % (FastANI interno de GTDB-Tk).
- Definicion de linaje: cluster de especie a 95 % de ANI (dRep, `-pa 0.90 -sa 0.95`), reutilizando la calidad de CheckM2.
- Prevalencia espacial: nº de muestras distintas por cluster.
- Arbol filogenomico: **arbol de colocacion de `classify_wf`** (`classify/*.classify.tree`). El **arbol de novo focalizado** (query + referencias del taxon asignado + outgroup) se realizara en la **Fase 4**, solo con los linajes seleccionados.
- Tabla por genoma (TSV): `phase2_selection.tsv`
- Genomas de los linajes top-3 (para Fase 3): `phase2_selected_lineages/`

> Nota: solo se analizan los genomas que superaron la Fase 1. La prevalencia no depende de las lecturas crudas (no disponibles) y por eso es el criterio primario de importancia de un linaje.
