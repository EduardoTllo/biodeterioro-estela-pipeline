# Fase 3 (OE3) - Censo de disponibilidad genomica y seleccion de especies diana

Generado: 2026-10-04 19:34

## Criterio de seleccion

Disponibilidad genomica publica. Umbral de viabilidad: **>= 15 genomas** de
calidad adecuada, aplicado sobre el conteo **post-QC** y no sobre el conteo
bruto de NCBI.

Admision de cada genoma publico al set de referencia:

- completitud >= 95.0 %, contaminacion <= 5.0 %
- <= 300 contigs (la fragmentacion parte genes e infla el genoma accesorio)
- aislados unicamente; se excluyen MAGs y genomas de muestra ambiental
- tope de 50 genomas por especie (computo y sesgo por clonalidad)

## Tabla de viabilidad

| # | Linaje | Especie GTDB | Especie NCBI / TaxID | GenBank | RefSeq | Cluster GTDB | Post-QC | Prev. esp. | Viable |
|---|---|---|---|---|---|---|---|---|---|
| 1 | L14 | *Bacillus licheniformis* | Bacillus licheniformis / 1402 | 1 | 1 | 307 | 279 | 1 | **Si** |
| 2 | L16 | *Bacillus altitudinis* | Bacillus altitudinis / 293387 | 1 | 1 | 217 | 209 | 1 | **Si** |
| 3 | L12 | *Peribacillus frigoritolerans* | Peribacillus frigoritolerans / 450367 | 1 | 1 | 82 | 72 | 1 | **Si** |
| 4 | L8 | *Acinetobacter schindleri* | Acinetobacter schindleri / 108981 | 1 | 1 | 37 | 29 | 1 | **Si** |
| 5 | L9 | *Priestia flexa* | Priestia flexa / 86664 | 1 | 1 | 31 | 29 | 1 | **Si** |
| 6 | L11 | *Staphylococcus warneri_A* | Staphylococcus warneri / 1292 | 1 | 1 | 27 | 24 | 1 | **Si** |
| 7 | L18 | *Cytobacillus firmus_B* | Cytobacillus firmus / 1399 | 1 | 1 | 18 | 18 | 1 | **Si** |
| NA | L1 | *Telluria timonae* | Massilia timonae / 47229 | 1 | 1 | 8 | 3 | 2 | **No** |
| NA | L3 | *Bacillus_AB infantis* | Bacillus infantis / 324767 | 1 | 1 | 12 | 9 | 1 | **No** |
| NA | L13 | *Rossellomorea marisflavi* | Rossellomorea marisflavi / 189381 | 1 | 1 | 4 | 4 | 1 | **No** |
| NA | L15 | *Metabacillus halosaccharovorans* | Metabacillus halosaccharovorans / 930124 | 1 | 1 | 4 | 4 | 1 | **No** |
| NA | L17 | *Cytobacillus oceanisediminis* | Cytobacillus oceanisediminis / 665099 | 1 | 1 | 1 | 1 | 1 | **No** |

## Especies que pasan a la descarga (ranking preliminar)

El top-3 definitivo se fija tras la desreplicacion al 99 % ANI
(04_clasificar_refs.py seleccionar).

1. **Bacillus licheniformis** (L14) - 279 genomas post-QC, 50 para el pangenoma; 83 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
2. **Bacillus altitudinis** (L16) - 209 genomas post-QC, 50 para el pangenoma; 90 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
3. **Peribacillus frigoritolerans** (L12) - 72 genomas post-QC, 50 para el pangenoma; 27 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
4. **Acinetobacter schindleri** (L8) - 29 genomas post-QC, 29 para el pangenoma; 19 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
5. **Priestia flexa** (L9) - 29 genomas post-QC, 29 para el pangenoma; 14 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
6. **Staphylococcus warneri_A** (L11) - 24 genomas post-QC, 24 para el pangenoma; 15 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).
7. **Cytobacillus firmus_B** (L18) - 18 genomas post-QC, 18 para el pangenoma; 5 habitat(s) distinto(s); prevalencia espacial 1; 1 bin(s) propio(s).

## Advertencias para la redaccion

- Linajes con prevalencia espacial > 1 que NO alcanzan el umbral: L1 (Telluria timonae).
  Se excluyen del pangenoma por factibilidad analitica, no por falta de
  relevancia ecologica; deben reportarse igualmente como resultado de OE2.
- Conteos NCBI NO comparables por sufijo GTDB (NCBI no separa la especie
  escindida, el total esta inflado): L11 (Staphylococcus warneri_A), L18 (Cytobacillus firmus_B), L3 (Bacillus_AB infantis). Usar el tamano del cluster GTDB.
- El ranking es PRELIMINAR: el orden definitivo se fija tras la desreplicacion
  a 99 % ANI (04_drep_refs.slurm y 04_clasificar_refs.py seleccionar), que
  mide diversidad de cepas no redundante.
- Los bins propios son MAGs (70-100 % de completitud) frente a aislados al
  >= 95 %. Las conclusiones sobre PRESENCIA de genes son validas; las de
  AUSENCIA no lo son sin controlar por la completitud del MAG.
