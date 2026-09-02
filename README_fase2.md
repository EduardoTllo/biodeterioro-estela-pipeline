# Fase 2 en Khipu: asignacion taxonomica y definicion de linajes (OE2)

Toma los 21 genomas seleccionados en la Fase 1, los clasifica taxonomicamente,
confirma su especie por ANI, los agrupa en linajes y los ordena por prevalencia
espacial. Produce un reporte del proceso.

Solo se analizan los genomas que superaron la Fase 1 (completitud > 70 %,
contaminacion < 5 %). Los nombres `bin-<N>-<muestra>` conservan la muestra de
origen, que es lo que permite calcular la prevalencia espacial.

> **Nota.** El top 3 por prevalencia que produce esta fase ya **no** determina
> que especies pasan al pangenoma. Desde la Fase 3 la seleccion se hace por
> disponibilidad genomica publica (>= 10 genomas de calidad); ver
> [README_fase3.md](README_fase3.md). La prevalencia sigue siendo un resultado de
> la Fase 2 y se reporta como tal.

## Que hace (tesis, seccion 4.1.4)

| Paso | Herramienta | Criterio |
|------|-------------|----------|
| 4.1.4.1 Clasificacion taxonomica | GTDB-Tk `classify_wf` (bac120/ar53, datos R220) | Ubica cada genoma en la taxonomia GTDB |
| 4.1.4.2 Confirmacion de especie | FastANI (interno de GTDB-Tk) | ANI >= 95 % y AF >= 65 % |
| 4.1.4.3 Desreplicacion y linajes | dRep (`-pa 0.90 -sa 0.95`) | Linaje = cluster de especie a 95 % de ANI; reutiliza la calidad de CheckM2 de la Fase 1 via `--genomeInfo`, sin re-correr CheckM |
| 4.1.4.4 Confirmacion filogenomica | GTDB-Tk (arbol de colocacion) | Arbol de referencia con los genomas injertados (`classify/*.classify.tree`) |
| 4.1.4.5 Seleccion de linajes | script | Ranking por prevalencia espacial y top 3 para la Fase 3 |

Aqui se usa el arbol de colocacion de `classify_wf`, que ya viene incluido y es
ligero, para verificar la ubicacion de cada linaje. El arbol de novo focalizado
(genomas propios, referencias del taxon asignado y outgroup) se hara en la Fase 4
y solo con los linajes seleccionados, para evitar el costo de un arbol de novo de
dominio completo, que no cabe en los 98 GB de RAM de la cuenta.

## Requisitos

La Fase 1 debe estar ejecutada, con los 21 genomas en
`~/estela/fase1/results/phase1_selected_genomes/` y el reporte de CheckM2 en
`~/estela/fase1/results/02_checkm2/quality_report.tsv`.

Hace falta ademas la base GTDB-Tk R220 (`~/gtdbtk_data/gtdbtk_r220_data.tar.gz`,
~102 GB, unos 107 GB descomprimida en `~/gtdbtk_data/release220/`). En `/home`
hay 18 TB libres.

## Versiones fijadas

El setup instala versiones pinneadas:

| Componente | Version | Motivo |
|---|---|---|
| GTDB-Tk | `2.6.1` | Debe ser compatible con el release de los datos. Segun la tabla oficial, R220 admite 2.4.0-2.6.1 y R232 requiere 2.7.0 o superior. |
| Python | `3.11` | Sin pin, conda resuelve Python 3.14, con el que GTDB-Tk falla ([issue #669](https://github.com/Ecogenomics/GTDBTk/issues/669)). |
| Canales conda | solo `conda-forge` y `bioconda` | En Khipu el canal `defaults` (`repo.anaconda.com`) no resuelve. El setup reescribe `~/.condarc` y guarda un backup. |

Se pueden sobreescribir por variable de entorno al migrar a otro release, por
ejemplo para datos R232:

```bash
GTDBTK_VERSION=2.7.2 bash 02_setup_fase2_khipu.sh
```

## Paso a paso

### 1) Subir los scripts a Khipu

```bash
ssh USUARIO@khipu.utec.edu.pe "mkdir -p ~/estela/fase2/scripts"
```

```bash
scp 02_setup_fase2_khipu.sh run_fase2.slurm USUARIO@khipu.utec.edu.pe:~/estela/fase2/
```

```bash
scp bin/summarize_phase2.py USUARIO@khipu.utec.edu.pe:~/estela/fase2/scripts/
```

### 2) Preparar el entorno (nodo de login, con internet)

```bash
ssh USUARIO@khipu.utec.edu.pe
```

```bash
cd ~/estela/fase2 && bash 02_setup_fase2_khipu.sh
```

Crea los entornos conda `gtdbtk` y `drep` con las versiones fijadas, verifica el
espacio y descomprime la base R220 en `~/gtdbtk_data/release220` (una sola vez,
unos 107 GB). Tambien fija `GTDBTK_DATA_PATH` dentro del entorno `gtdbtk` y
verifica la integridad de la base.

Para comprobar que quedo bien antes de seguir:

```bash
conda activate gtdbtk && gtdbtk --version && python --version && echo "DATA=$GTDBTK_DATA_PATH"
```

Debe imprimir GTDB-Tk 2.6.1, Python 3.11.x y la ruta a `release220`.

La descompresion es larga y conviene lanzarla con `sbatch` o dentro de `tmux`,
para que no dependa de la sesion SSH.

### 3) Enviar el job y monitorear

```bash
cd ~/estela/fase2 && sbatch run_fase2.slurm
```

```bash
squeue --me
```

```bash
tail -f fase2_<JOBID>.out
```

### 4) Resultados (`~/estela/fase2/results/`)

| Archivo | Contenido |
|---------|-----------|
| `phase2_report.md` | Reporte con taxonomia, linajes, prevalencia y top 3 |
| `phase2_selection.tsv` | Tabla por genoma (taxonomia, ANI/AF, linaje) |
| `03_gtdbtk/classify/*summary.tsv` | Clasificacion GTDB-Tk con ANI y AF |
| `03_gtdbtk/classify/*.classify.tree` | Arbol de colocacion |
| `04_drep/data_tables/Cdb.csv` | Asignacion de cada genoma a su cluster o linaje |
| `phase2_selected_lineages/L1,L2,L3/` | Genomas de los 3 linajes priorizados (entrada de la Fase 3) |
| `versions.txt` | Versiones y umbrales usados |

```bash
scp USUARIO@khipu.utec.edu.pe:~/estela/fase2/results/phase2_report.md .
```

En el repositorio, los resultados de texto de esta fase estan en
[`results/fase2/`](results/fase2/). Las rutas de la tabla de arriba son las del
cluster; en el repo llevan el prefijo `phase2_`. El MSA de GTDB-Tk (~230 MB) y
los FASTA no se versionan.

## Recursos SLURM

`run_fase2.slurm` pide particion `standard`, 32 CPU, 90 GB de RAM y 12 h, con
`--pplacer_cpus 1` para acotar la memoria de GTDB-Tk (pplacer bac120 usa unos
55 GB).

## Reproducibilidad

Las versiones exactas de GTDB-Tk y dRep quedan en `envs/gtdbtk.lock.yml` y
`envs/drep.lock.yml`, generadas por el setup, y en el `results/versions.txt` de
cada corrida. Los umbrales se pasan por parametro: ANI de especie 95 %, AF 65 %,
dRep `-pa 0.90 -sa 0.95` y top 3 linajes.

dRep reutiliza la calidad de CheckM2 de la Fase 1 via `genomeInfo.csv`; su filtro
interno se desactiva (`-comp 0 -con 100`) porque el QC ya se aplico antes.

`classify_wf` corre con `--skip_ani_screen`, que omite el prefiltro Mash: no
aporta a la confirmacion de especie y evita construir una base Mash extra. La
confirmacion de especie se hace igual con el FastANI interno de GTDB-Tk.
