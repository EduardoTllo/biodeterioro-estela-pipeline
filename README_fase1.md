# Fase 1 en Khipu: control de calidad y seleccion de genomas (OE1)

Toma los bins genomicos crudos y devuelve los genomas procariotas con calidad
suficiente para la Fase 2, junto con un reporte del proceso.

## Que hace (tesis, seccion 4.1.3)

| Paso | Herramienta | Criterio |
|------|-------------|----------|
| 1. Filtrado por dominio | Tiara | Longitud minima de contig de 3000 bp. Se descarta un bin si su fraccion mayoritaria de longitud es eucariota u organelar. |
| 2. Calidad genomica | CheckM2 | Se conservan genomas con completitud > 70 % y contaminacion < 5 %. |

La salida es la lista de genomas seleccionados, las tablas intermedias y
`phase1_report.md`.

## Datos de entrada

La fuente es `C:\Tesis-EstelaRaimondi\Datos\Shotgun Analysis\bins completos\`,
con los 97 bins individuales nombrados `bin-<N>-<muestra>.fasta` (por ejemplo
`bin-1-49.fasta`, bin 1 de la muestra 49). El nombre codifica la muestra de
origen, asi que la trazabilidad para la Fase 2 queda preservada. El conteo por
muestra coincide con la Tabla 4.1 de la tesis: 97 bins, con las muestras 66 y 67
sin bins.

No se usan `MAGs_v2/`, `MAGs_v2.zip` ni los `<muestra>_mag.fasta`, que son 24
ensamblajes por muestra (contigs `k141_*`) y no bins; un QC sobre ellos daria
contaminacion altisima. Tampoco `Datos\MAGs_all\`, una copia obsoleta de
`bins completos` con nombres corruptos (`bin-1.fasta 10`).

`01_prep_bins.sh` copia los 97 `bin-<N>-<muestra>.fasta` a `bins_clean/` y genera
`bin_rename_map.tsv` con la muestra de origen, el numero de contigs y la longitud
de cada bin.

## Archivos

```
fase1_khipu/
├── 00_setup_khipu.sh          # entornos conda + DB de CheckM2 (nodo login)
├── 01_prep_bins.sh            # normaliza nombres de bins + tabla de mapeo
├── run_fase1.slurm            # job SLURM que corre toda la Fase 1
├── bin/summarize_phase1.py    # genera tablas y el reporte
├── envs/                      # lock files de conda (se crean en el setup)
└── README_fase1.md            # este archivo
```

## Paso a paso

Los `scp` y la preparacion de bins se corren en el laptop (WSL/Ubuntu), en una
terminal que no este conectada a Khipu; en WSL el disco `C:` esta en `/mnt/c/`.
El `ssh`, el `tar`, el setup y el `sbatch` se corren dentro de Khipu. Reemplaza
`USUARIO` por tu usuario del cluster (por ejemplo `eduardo.tello`).

### 1) Preparar los bins (una vez, en el laptop)

```bash
bash 01_prep_bins.sh "/mnt/c/Tesis-EstelaRaimondi/Datos/Shotgun Analysis/bins completos" "/mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean"
```

```bash
tar czf /mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz -C "/mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean" .
```

El tarball pesa unos 335 MB (los FASTA comprimen bien: ~1.1 GB sin comprimir).

### 2) Subir a Khipu (en el laptop)

```bash
ssh USUARIO@khipu.utec.edu.pe "mkdir -p ~/estela/fase1/bins ~/estela/fase1/scripts"
```

```bash
scp /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/00_setup_khipu.sh /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/run_fase1.slurm USUARIO@khipu.utec.edu.pe:~/estela/fase1/
```

```bash
scp /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/bin/summarize_phase1.py USUARIO@khipu.utec.edu.pe:~/estela/fase1/scripts/
```

```bash
scp /mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz USUARIO@khipu.utec.edu.pe:~/estela/fase1/
```

### 3) Preparar el entorno (nodo de login de Khipu)

El nodo de login es el unico con salida a internet, que hace falta para crear los
entornos y descargar la base de CheckM2 (~3 GB).

```bash
ssh USUARIO@khipu.utec.edu.pe
```

```bash
cd ~/estela/fase1 && tar xzf bins_clean.tar.gz -C bins && bash 00_setup_khipu.sh
```

### 4) Enviar el job y monitorear

```bash
cd ~/estela/fase1 && sbatch run_fase1.slurm
```

```bash
squeue --me
```

```bash
tail -f fase1_<JOBID>.out
```

### 5) Resultados (`~/estela/fase1/results/`)

| Archivo | Contenido |
|---------|-----------|
| `phase1_report.md` | Reporte del proceso completo |
| `phase1_selection.tsv` | Decision final por bin |
| `01_tiara/tiara_bin_summary.tsv` | Composicion de dominio por bin |
| `02_checkm2/quality_report.tsv` | Completitud y contaminacion (CheckM2) |
| `phase1_selected_genomes/` | FASTA de los genomas seleccionados |
| `versions.txt` | Versiones de herramientas y umbrales usados |

Para traer el reporte al laptop:

```bash
scp USUARIO@khipu.utec.edu.pe:~/estela/fase1/results/phase1_report.md /mnt/c/Tesis-EstelaRaimondi/
```

En el repositorio, los resultados de texto de esta fase estan en
[`results/fase1/`](results/fase1/). Las rutas de la tabla de arriba son las del
cluster, donde el job escribe; en el repo llevan el prefijo `phase1_`.

## Recursos SLURM

`run_fase1.slurm` pide particion `standard`, 32 CPU, 90 GB de RAM y 24 h, sin
GPU, dentro de los limites de la cuenta (32 cores, 98 GB, 24 h). Bajar `--time`
reduce la espera en cola.

## Reproducibilidad

Las versiones exactas de Tiara y CheckM2 quedan en `envs/*.lock.yml`, generadas
por el setup, y en el `results/versions.txt` de cada corrida. Los umbrales se
pasan por parametro y no estan hardcodeados: contig minimo de 3000 bp,
completitud > 70 %, contaminacion < 5 %. `bin_rename_map.tsv` documenta cada bin
con su muestra de origen (49-72), numero de contigs y longitud, que es el insumo
directo para la prevalencia espacial de la Fase 2.

El umbral de >70 % y <5 % es criterio propio del proyecto y no una categoria
MIMAG estandar (tesis 4.1.3.2). Tiara no separa bacteria de arquea; eso se
resuelve en la Fase 2 con GTDB-Tk.
