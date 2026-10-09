# Pipeline de biodeterioro - Estela de Raimondi

Flujo bioinformatico para caracterizar el potencial de biodeterioro y las senales
de adaptacion genomica de los linajes bacterianos de la fraccion viable de la
Estela de Raimondi. Se ejecuta en el cluster Khipu (UTEC).

Proyecto final de carrera, Bioingenieria, UTEC. Datos de metagenomica shotgun de
la fraccion cultivable, provistos por el Laboratorio de Genomica Microbiana
(UPCH).

## Estado del proyecto

| Fase | Objetivo | Estado | Resultados |
|---|---|---|---|
| 1 | Control de calidad y seleccion de genomas (OE1) | Ejecutada | [`results/fase1/`](results/fase1/) |
| 2 | Asignacion taxonomica y definicion de linajes (OE2) | Ejecutada | [`results/fase2/`](results/fase2/) |
| 3 | Analisis pangenomico comparativo (OE3) | Ejecutada; falta la parte remota de la v2 de genes especificos | [`results/fase3/phase3_resultados.md`](results/fase3/phase3_resultados.md) |
| 4 | Potencial metabolico y de riesgo (OE4) | Pendiente | - |

En la Fase 1, de 97 bins crudos 66 pasan el filtro por dominio (Tiara) y quedan
21 genomas seleccionados por CheckM2 (completitud > 70 %, contaminacion < 5 %).
En la Fase 2, esos 21 genomas se agrupan en 19 linajes (clusteres de especie a
95 % ANI), 16 de ellos con especie confirmada. Los tres linajes mas prevalentes
son L1 *Telluria timonae* (muestras 50 y 52), L2 *Corynebacterium* sp. (68 y 69)
y L3 *Bacillus_AB infantis* (61).

La Fase 3 selecciona las especies para el analisis pangenomico por
**disponibilidad genomica publica**, con un umbral de 15 genomas de calidad, y no
por prevalencia espacial. El cambio es de factibilidad: un pangenoma necesita un
conjunto de referencia suficiente, y varios de los linajes mas prevalentes no lo
tienen. La prevalencia espacial se mantiene como resultado de la Fase 2 y viaja
como columna en la tabla de viabilidad, para dejar constancia de los casos en que
los dos criterios discrepan. Cada pangenoma se arma con hasta 50 referencias
publicas (GTDB R220 + NCBI) mas el bin de la Estela, todo anotado con Bakta y
agrupado con Panaroo; la particion core/shell/cloud se calcula solo con las
referencias, y los genes candidatos a exclusivos de la Estela pasan por seis
filtros de verificacion. El detalle esta en [README_fase3.md](README_fase3.md).

Pasaron el umbral 3 especies: *Bacillus altitudinis* (L16, 50 referencias),
*Peribacillus frigoritolerans* (L12, 42) y *Acinetobacter schindleri* (L8, 23).
Los tres pangenomas son abiertos (7 059, 16 125 y 7 599 familias). El flujo
completo con sus parametros, la ejecucion en Khipu y la interpretacion de cada
figura estan en
[results/fase3/phase3_resultados.md](results/fase3/phase3_resultados.md).

## Estructura del repositorio

```
.
├── README.md                    # este archivo
│
├── README_fase1.md              # guia de la Fase 1
├── 01_prep_bins.sh              #   prepara bins_clean/ + tabla de procedencia
├── 00_setup_khipu.sh            #   entornos conda (tiara, checkm2) + BD
├── run_fase1.slurm              #   job SLURM de la Fase 1
│
├── README_fase2.md              # guia de la Fase 2
├── 02_setup_fase2_khipu.sh      #   entornos conda (gtdbtk, drep) + BD GTDB R220
├── run_fase2.slurm              #   job SLURM de la Fase 2
│
├── README_fase3.md              # guia de la Fase 3
├── 03_setup_fase3_khipu.sh      #   entornos (bakta, panaroo, iqtree, blast...) + BD
├── 03_censo_genomas.py          #   censo de disponibilidad genomica (GTDB R220)
├── 04_fetch_refs.sh             #   descarga de referencias (NCBI Datasets)
├── 04_clasificar_refs.py        #   verificacion, prioridad/habitat, ranking, seleccion
├── 04_drep_refs.slurm           #   dRep 99 % sobre las referencias
├── 05_bakta_all.slurm           #   anotacion homogenea con Bakta (array)
├── 06_panaroo.slurm             #   pangenoma con Panaroo + particion
├── 06b_iqtree.slurm             #   arbol del core por especie
├── 07_particion.py              #   control de Bakta y particion solo con referencias
├── 08_exclusivos.py             #   filtros F1-F6 de genes exclusivos
├── 08_exclusivos_local.slurm    #   F1-F4 (tblastn/blastn locales)
├── 08_exclusivos_remoto.sh      #   F5-F6 (BLAST remoto contra nr/core_nt)
│
├── bin/
│   ├── summarize_phase1.py      # tablas y reporte de la Fase 1
│   ├── summarize_phase2.py      # tablas y reporte de la Fase 2
│   └── export_fase3.sh          # reune los resultados versionables de la Fase 3
│
├── figuras/                     # scripts en R para las figuras y tablas
│
├── metadata/
│   ├── bin_provenance.tsv       # los 97 bins con su muestra de origen (49-72)
│   ├── habitat_keywords.tsv     # patrones para clasificar la fuente de aislamiento
│   └── paises_continentes.tsv   # pais -> continente y marca de Latinoamerica
│
├── envs/                        # lock files de conda (versiones exactas)
│
└── results/
    ├── fase1/
    │   ├── phase1_report.md
    │   ├── phase1_selection.tsv              # decision por bin
    │   ├── phase1_checkm2_quality_report.tsv # completitud/contaminacion
    │   ├── phase1_tiara_bin_summary.tsv      # composicion de dominio
    │   └── phase1_versions.txt
    ├── fase2/
    │   ├── phase2_report.md
    │   ├── phase2_selection.tsv              # taxonomia + linaje por genoma
    │   ├── phase2_gtdbtk_bac120_summary.tsv  # salida completa de GTDB-Tk
    │   ├── phase2_drep_{Cdb,Wdb,Ndb,Sdb}.csv # clusteres/linajes de dRep
    │   ├── phase2_gtdbtk_backbone_bac120.classify.tree
    │   ├── phase2_gtdbtk_tree_mapping.tsv
    │   ├── phase2_genomeInfo.csv
    │   ├── phase2_versions.txt
    │   └── figuras_drep/                     # dendrogramas y graficos de dRep
    └── fase3/
        ├── phase3_resultados.md              # flujo, parametros y resultados de la fase
        ├── phase3_viabilidad.tsv             # decision de viabilidad por linaje
        ├── phase3_ranking_especies.tsv       # ranking por genomas no redundantes
        ├── phase3_referencias_finales.tsv    # referencias de cada pangenoma
        └── <especie>/                        # particion, arbol del core y exclusivos
```

En `results/` se versionan solo tablas, reportes y arboles, que son la evidencia
reproducible del analisis. Los genomas (FASTA), las bases de datos y los archivos
pesados, como el MSA de 230 MB de GTDB-Tk, quedan excluidos por `.gitignore`:
viven en el cluster y en el almacenamiento local.

## Por donde empezar

Para entender que se hizo, lee los reportes
[`results/fase1/phase1_report.md`](results/fase1/phase1_report.md),
[`results/fase2/phase2_report.md`](results/fase2/phase2_report.md) y
[`results/fase3/phase3_resultados.md`](results/fase3/phase3_resultados.md).

Para reproducir una fase, sigue su README ([Fase 1](README_fase1.md),
[Fase 2](README_fase2.md), [Fase 3](README_fase3.md)), con el paso a paso de
subida a Khipu, setup y envio del job.

## Convenciones

Los genomas se nombran `bin-<N>-<muestra>` (por ejemplo `bin-6-50` es el bin 6 de
la libreria 50). El nombre codifica la muestra de origen, que es lo que permite
calcular la prevalencia espacial dado que no se dispone de las lecturas crudas.
La prevalencia mide la distribucion de un linaje en la estela; la seleccion de
especies para el pangenoma se decide aparte, por disponibilidad genomica
(ver [Fase 3](README_fase3.md)).

Los jobs usan la particion `standard` dentro de los limites de la cuenta
(32 cores, 98 GB de RAM, 24 h). Los entornos conda llevan versiones fijadas,
incluido `python=3.11`, y estan exportados a `envs/*.lock.yml`.
