# Fase 2 en Khipu — Asignación taxonómica y definición de linajes (OE2)

Toma los **21 genomas seleccionados en la Fase 1** y los clasifica
taxonómicamente, confirma su especie por ANI, los agrupa en **linajes** y
prioriza los **3 linajes más prevalentes** para la Fase 3. Produce un
**reporte legible** de todo el proceso.

> Solo se analizan los genomas que **superaron la Fase 1** (completitud > 70 %,
> contaminación < 5 %). Los nombres `bin-<N>-<muestra>` conservan la muestra de
> origen, que es lo que permite calcular la prevalencia espacial.

## Qué hace (según la tesis, sección 4.1.4)

| Paso | Herramienta | Criterio |
|------|-------------|----------|
| 4.1.4.1 Clasificación taxonómica | **GTDB-Tk `classify_wf`** (bac120/ar53, datos **R220**) | Ubica cada genoma en la taxonomía GTDB |
| 4.1.4.2 Confirmación de especie | **FastANI** (interno de GTDB-Tk) | ANI ≥ 95 % y AF ≥ 65 % |
| 4.1.4.3 Desreplicación y linajes | **dRep** (`-pa 0.90 -sa 0.95`) | Linaje = clúster de especie a 95 % de ANI; reutiliza la calidad de CheckM2 (Fase 1) vía `--genomeInfo`, sin re-correr CheckM |
| 4.1.4.4 Confirmación filogenómica | **GTDB-Tk** (árbol de colocación) | Árbol de referencia con tus genomas injertados (`classify/*.classify.tree`) |
| 4.1.4.5 Selección de linajes | script | Ranking por **prevalencia espacial** → top 3 para la Fase 3 |

### 🌳 Nota sobre el árbol filogenómico

En la **Fase 2** se usa el **árbol de colocación** de `classify_wf` (ligero, ya
incluido) para verificar visualmente la ubicación de cada linaje. El **árbol de
novo focalizado** (tus genomas + referencias del taxón asignado + outgroup,
inferido de novo) se realizará en la **Fase 4**, y **solo con los linajes
seleccionados** — evitando el costo de un árbol de novo de dominio completo, que
no cabe en los 98 GB de RAM de la cuenta.

## Requisitos

- Fase 1 ejecutada: los 21 genomas en `~/estela/fase1/results/phase1_selected_genomes/`
  y el reporte de CheckM2 en `~/estela/fase1/results/02_checkm2/quality_report.tsv`.
- Base **GTDB-Tk R220** (`~/gtdbtk_data/gtdbtk_r220_data.tar.gz`, ~102 GB; ~110 GB
  descomprimida). En `/home` hay 18 TB libres, así que entra sin problema.

## Paso a paso

### 1) Subir los scripts de Fase 2 a Khipu

```bash
ssh USUARIO@khipu.utec.edu.pe "mkdir -p ~/estela/fase2/scripts"
scp 02_setup_fase2_khipu.sh run_fase2.slurm USUARIO@khipu.utec.edu.pe:~/estela/fase2/
scp bin/summarize_phase2.py                 USUARIO@khipu.utec.edu.pe:~/estela/fase2/scripts/
```

### 2) Preparar el entorno (nodo de login, con internet)

```bash
ssh USUARIO@khipu.utec.edu.pe
cd ~/estela/fase2
bash 02_setup_fase2_khipu.sh
```

Crea los entornos conda `gtdbtk` y `drep`, verifica el espacio y **descomprime la
base R220** en `~/gtdbtk_data/release220` (una sola vez; ~110 GB). Fija
`GTDBTK_DATA_PATH` automáticamente en el entorno `gtdbtk`.

### 3) Enviar el job

```bash
cd ~/estela/fase2
sbatch run_fase2.slurm
```

### 4) Monitorear

```bash
squeue --me
tail -f fase2_<JOBID>.out
```

### 5) Resultados (`~/estela/fase2/results/`)

| Archivo | Contenido |
|---------|-----------|
| `phase2_report.md` | **Reporte legible**: taxonomía, linajes, prevalencia, top-3 |
| `phase2_selection.tsv` | Tabla por genoma (taxonomía, ANI/AF, linaje) |
| `03_gtdbtk/classify/*summary.tsv` | Clasificación GTDB-Tk + ANI/AF |
| `03_gtdbtk/classify/*.classify.tree` | Árbol de colocación |
| `04_drep/data_tables/Cdb.csv` | Asignación de cada genoma a su clúster/linaje |
| `phase2_selected_lineages/L1,L2,L3/` | Genomas de los 3 linajes priorizados (entrada de la Fase 3) |
| `versions.txt` | Versiones y umbrales usados |

```bash
scp USUARIO@khipu.utec.edu.pe:~/estela/fase2/results/phase2_report.md .
```

## Recursos SLURM (dentro de tus límites: 32 cores, 98 GB, 24 h)

`run_fase2.slurm` pide **partición `standard`, 32 CPU, 90 GB RAM, 12 h**, con
`--pplacer_cpus 1` para acotar la memoria de GTDB-Tk (pplacer bac120 ≈ 55 GB).

## Reproducibilidad

- Versiones exactas de GTDB-Tk/dRep en `envs/gtdbtk.lock.yml`, `envs/drep.lock.yml`
  (generadas por el setup) y en `results/versions.txt` de cada corrida.
- Umbrales por parámetros: ANI especie 95 %, AF 65 %, dRep `-pa 0.90 -sa 0.95`,
  top 3 linajes.
- dRep reutiliza la calidad de CheckM2 (Fase 1) vía `genomeInfo.csv`; se desactiva
  su filtro interno (`-comp 0 -con 100`) porque el QC ya se aplicó en la Fase 1.

> Nota: `classify_wf` usa `--skip_ani_screen` (omite el pre-filtro Mash, que no
> aporta a la confirmación de especie y evita construir una base Mash extra); la
> confirmación de especie igual se hace con el FastANI interno de GTDB-Tk.
