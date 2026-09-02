# Fase 3 en Khipu: analisis pangenomico comparativo (OE3)

Toma los linajes definidos en la Fase 2, mide cuantos genomas publicos de calidad
existe para cada especie, selecciona las especies viables y construye un
pangenoma por especie que ubica a los bins propios dentro de la diversidad
conocida del taxon.

## Criterio de seleccion

**Disponibilidad genomica publica**, con un umbral de **>= 10 genomas** de calidad
adecuada por especie.

Este criterio reemplaza al de prevalencia espacial que la Fase 2 usaba para
priorizar linajes. El motivo es de factibilidad: un pangenoma necesita un
conjunto de referencia suficiente, y varios de los linajes mas prevalentes no lo
tienen. La prevalencia espacial se sigue calculando y reportando como resultado
de la Fase 2, y aparece como columna en la tabla de viabilidad para dejar
constancia de los casos en que ambos criterios discrepan.

El umbral se aplica sobre el conteo **post-QC**, no sobre el conteo bruto de
NCBI. Una especie con 12 MAGs fragmentados no es viable; una con 11 aislados
completos si lo es.

## Etapas

| Etapa | Herramienta | Criterio |
|---|---|---|
| 3.0 Censo de disponibilidad | `03_censo_genomas.py censo` + metadata de GTDB R220 | Mapeo GTDB<->NCBI, tamano del cluster de especie, conteo RefSeq/GenBank |
| 3.1 Filtro de calidad | metadata de GTDB (CheckM2) | Completitud >= 95 %, contaminacion <= 5 %, <= 300 contigs, solo aislados |
| 3.2 Seleccion de especies diana | `03_censo_genomas.py tabla` | Umbral >= 10 genomas post-QC; ranking por disponibilidad; tope de 50 genomas |
| 3.3 Anotacion homogenea | Bakta | Misma version y misma BD para referencias y bins propios |
| 3.4 Pangenoma | Panaroo (`--clean-mode moderate`) + PPanGGOLiN | Core >= 95 %, shell 15-95 %, cloud < 15 % |
| 3.5 Analisis | scripts | Curvas de acumulacion, posicion del MAG propio, genes exclusivos |

Los genes accesorios y exclusivos que salgan de la etapa 3.5 son la entrada de la
Fase 4 (DRAM, marcadores KO de biodeterioro y supervivencia).

## Solo se censan los linajes con binomio real

De los 19 linajes de la Fase 2 se censan **12**. Quedan fuera 7:

- **Sin especie asignada** (`s__` vacio): L2 *Corynebacterium*, L4 *Sporosarcina*,
  L5 *Pantoea*, L10 *Aristophania*.
- **Nombre placeholder de GTDB** (`sp<digitos>`): L6 *Brevibacterium* sp021023135,
  L7 *Paenibacillus* sp001955535, L19 *Mesobacillus* sp014856545. Son clusteres
  validos a nivel de especie, pero sin binomio publicado ni equivalente estable en
  NCBI, asi que no admiten un censo por nombre.

L2 merece mencion aparte: con ANI 90,55 % y AF 44,3 % es candidato a especie
nueva, y ademas es el segundo linaje mas prevalente. Se excluye del pangenoma por
falta de conjunto de referencia, no por falta de interes, y debe reportarse como
hallazgo de novedad taxonomica.

## Dos trampas taxonomicas que el script resuelve

**Nombres reasignados.** GTDB llama *Telluria timonae* a lo que NCBI Taxonomy
sigue llamando *Massilia timonae*. Buscar el nombre GTDB en NCBI devuelve cero
genomas y llevaria a declarar L1 inviable por error. El script deriva el nombre y
el TaxID de NCBI por voto mayoritario de `ncbi_species_taxid` dentro del cluster
de GTDB, asi que el caso se resuelve solo y queda avisado en pantalla.

**Sufijos GTDB.** Los sufijos `_A`, `_B`, `_AB` marcan especies que GTDB escinde y
NCBI no separa (L3 *Bacillus_AB infantis*, L11 *Staphylococcus warneri_A*,
L18 *Cytobacillus firmus_B*). En esos linajes el conteo de NCBI esta inflado
porque agrupa varias especies GTDB bajo un solo nombre. El script los marca con
`conteo_ncbi_comparable = no`; ahi el numero valido es el tamano del cluster de
GTDB.

El script tambien reporta `acuerdo_mapeo`, la fraccion de genomas del cluster que
comparten el TaxID mayoritario. Por debajo de 0,90 avisa: significa que NCBI
reparte ese cluster entre varias especies y conviene revisar el mapeo a mano.

## Requisitos

La Fase 2 debe estar ejecutada, con `results/fase2/phase2_selection.tsv`
disponible.

Hace falta el metadata del mismo release de GTDB usado en la Fase 2, que **no
viene** en el tarball de GTDB-Tk y se descarga aparte (~120 MB comprimido) en el
nodo de login:

```bash
wget -P ~/gtdbtk_data/release220/ https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz
```

```bash
gunzip ~/gtdbtk_data/release220/bac120_metadata_r220.tsv.gz
```

Para la etapa 3.0 basta con eso y con Python 3; no se necesita ningun entorno
conda. El CLI de NCBI Datasets solo hace falta para los conteos vigentes de
GenBank y RefSeq.

## Paso a paso

### 1) Subir el script a Khipu

```bash
scp 03_censo_genomas.py USUARIO@khipu.utec.edu.pe:~/estela/fase3/
```

### 2) Censo de disponibilidad (no requiere internet)

Desde `~/estela/fase3`, con el `phase2_selection.tsv` de la Fase 2 a mano:

```bash
python 03_censo_genomas.py censo --selection ~/estela/fase2/results/phase2_selection.tsv --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv --outdir results/fase3
```

Recorre el metadata en streaming (son cientos de MB) y tarda un par de minutos.
Fijate en la salida: imprime la particion 12 / 7 de los linajes, avisa de cada
nombre NCBI distinto al esperado y marca los mapeos con acuerdo bajo.

Los umbrales de calidad se pueden mover con `--min-completitud`,
`--max-contaminacion` y `--max-contigs`.

### 3) Conteos vigentes de NCBI (nodo de LOGIN, con internet)

Este paso es opcional: sin el, la tabla se arma solo con los conteos de GTDB, que
son mas conservadores porque el metadata corresponde a un release congelado.

Instalar el CLI una sola vez, en el nodo de login, respetando la restriccion de
canales de Khipu:

```bash
conda create -y -n ncbi-datasets --override-channels -c conda-forge -c bioconda python=3.11 ncbi-datasets-cli
```

```bash
conda activate ncbi-datasets && python 03_censo_genomas.py ncbi --outdir results/fase3
```

Son tres consultas por linaje (GenBank, RefSeq, RefSeq completos). Si una falla,
el script avisa y sigue con las demas en vez de abortar.

### 4) Tabla de viabilidad y seleccion

```bash
python 03_censo_genomas.py tabla --outdir results/fase3 --umbral 10 --top 3
```

`--umbral` fija el minimo de genomas post-QC, `--tope` el maximo de genomas por
pangenoma (50 por defecto) y `--top` cuantas especies pasan a la etapa 3.3.

El ranking que sale aqui es **preliminar**: ordena por numero de genomas post-QC,
que premia a las especies clinicas o industriales sobresecuenciadas. El orden
definitivo se fija tras la desreplicacion a 99 % ANI de la etapa siguiente, que
mide diversidad de cepas no redundante.

## Salidas

Todo va a `results/fase3/`:

| Archivo | Contenido |
|---|---|
| `phase3_censo_genomas.tsv` | Censo por linaje, 30 columnas |
| `phase3_genomas_candidatos.tsv` | Accesiones que pasan el QC, una por fila. Es la entrada de la etapa 3.3 |
| `phase3_mapeo_ncbi.tsv` | Mapeo GTDB<->NCBI con la fraccion de acuerdo |
| `phase3_ncbi_counts.tsv` | Conteos GenBank/RefSeq vigentes y fecha de consulta |
| `phase3_viabilidad.tsv` | Tabla comparativa con la decision de viabilidad |
| `phase3_censo_report.md` | Reporte legible, con la tabla ya en markdown |

La columna `descartes_qc` desglosa por que se cayo cada genoma (`no_aislado`,
`fragmentado`, `completitud_baja`, `contaminacion_alta`, `sin_metricas`), que es
lo que permite justificar el filtro en la tesis.

## Advertencia de interpretacion

Los bins propios son MAGs, con 70-100 % de completitud, mientras que los genomas
de referencia entran al set con >= 95 %. Por construccion, un MAG va a *parecer*
que le faltan genes del core.

En consecuencia, **las conclusiones sobre presencia de genes son validas; las
conclusiones sobre ausencia no lo son** sin controlar por la completitud del MAG.
Hay que reportar la completitud CheckM2 de cada bin junto a su perfil y limitar
las afirmaciones fuertes a los genes presentes. Esto condiciona tambien la Fase 4:
un marcador KO ausente en el MAG no prueba ausencia funcional.

## Estado

Etapas 3.0-3.2 implementadas en `03_censo_genomas.py`, probadas end-to-end contra
un metadata sintetico. Pendientes: `04_fetch_refs.sh` (descarga, QC y dRep al
99 %), `05_bakta_all.slurm` y `06_panaroo.slurm`.
