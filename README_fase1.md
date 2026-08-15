# Fase 1 en Khipu — Control de calidad y selección de genomas (OE1)

Pipeline reproducible para ejecutar la **Fase 1** de la tesis de la Estela de
Raimondi en el clúster **Khipu (UTEC)**. Toma los *bins* genómicos crudos y
devuelve los genomas procariotas de calidad suficiente para la Fase 2, más un
**reporte legible** de todo el proceso.

## Qué hace (según la tesis, sección 4.1.3)

| Paso | Herramienta | Criterio |
|------|-------------|----------|
| 1. Filtrado por dominio | **Tiara** | Longitud mínima de contig = 3000 bp. Se descarta un bin si su fracción **mayoritaria de longitud** es eucariota u organelar. |
| 2. Calidad genómica | **CheckM2** | Se conservan genomas con **completitud > 70 %** y **contaminación < 5 %**. |

Salida: lista de genomas seleccionados + tablas + `phase1_report.md`.

---

## Datos de entrada

- **Fuente correcta:** `C:\Tesis-EstelaRaimondi\Datos\Shotgun Analysis\bins completos\`
  → los **97 bins individuales** con nombres `bin-<N>-<muestra>.fasta`
  (p.ej. `bin-1-49.fasta` = bin 1 de la muestra 49). El nombre **codifica la
  muestra de origen**, así que la trazabilidad para la Fase 2 queda preservada.
  El conteo por muestra coincide exactamente con la Tabla 4.1 de la tesis
  (97 bins; muestras 66 y 67 sin bins).
- **No usar** `MAGs_v2/`, `MAGs_v2.zip` ni los `<muestra>_mag.fasta`: son 24
  ensamblajes **por muestra** (contigs `k141_*`), no bins; un QC sobre ellos
  daría contaminación altísima.
- **No usar** `Datos\MAGs_all\`: es una copia con nombres corruptos
  (`bin-1.fasta 10`) de `bins completos`; quedó obsoleta.

El script `01_prep_bins.sh` copia los 97 `bin-<N>-<muestra>.fasta` a `bins_clean/`
y genera `bin_rename_map.tsv` con la muestra de origen, nº de contigs y longitud
de cada bin.

---

## Archivos de este repositorio

```
fase1_khipu/
├── 00_setup_khipu.sh          # entornos conda + DB de CheckM2 (nodo login)
├── 01_prep_bins.sh            # normaliza nombres de bins + tabla de mapeo
├── run_fase1.slurm            # job SLURM que corre toda la Fase 1
├── bin/summarize_phase1.py    # genera tablas y el reporte legible
├── envs/                      # lock files de conda (se crean en el setup)
├── .gitignore
└── README_fase1.md            # este archivo
```

---

## Paso a paso

> **Dónde se corre cada cosa (¡importante!):**
> - Los `scp` y la preparación de bins → **en tu laptop (WSL/Ubuntu)**, en una
>   terminal que **NO** esté conectada a Khipu. En WSL, tu disco `C:` está en
>   `/mnt/c/`.
> - `ssh`, `tar`, `bash 00_setup...`, `sbatch` → **dentro de Khipu** (por SSH).
> - Reemplaza `USUARIO` por tu usuario de Khipu (p.ej. `eduardo.tello`).

### 0) Preparar los bins (una vez, en tu laptop / WSL)

```bash
cd /mnt/c/Tesis-EstelaRaimondi/fase1_khipu
bash 01_prep_bins.sh "/mnt/c/Tesis-EstelaRaimondi/Datos/Shotgun Analysis/bins completos" "/mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean"
# Empaquetar para subir (FASTA comprime bien: ~1.1 GB -> ~335 MB)
tar czf /mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz -C "/mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean" .
```

### 1) Subir a Khipu (en tu laptop / WSL, NO dentro de Khipu)

```bash
ssh USUARIO@khipu.utec.edu.pe "mkdir -p ~/estela/fase1/bins ~/estela/fase1/scripts"
```
```bash
# scripts
scp /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/00_setup_khipu.sh /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/run_fase1.slurm USUARIO@khipu.utec.edu.pe:~/estela/fase1/
scp /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/bin/summarize_phase1.py USUARIO@khipu.utec.edu.pe:~/estela/fase1/scripts/
# datos (tarball)
scp /mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz USUARIO@khipu.utec.edu.pe:~/estela/fase1/
```

### 2) Preparar el entorno (en el NODO DE LOGIN de Khipu, requiere internet)

```bash
ssh USUARIO@khipu.utec.edu.pe
cd ~/estela/fase1
tar xzf bins_clean.tar.gz -C bins      # deja bin-<N>-<muestra>.fasta en bins/
bash 00_setup_khipu.sh                  # crea entornos + descarga DB (~3 GB)
```

> ⚠️ El nodo de acceso necesita salida a internet para este paso. Si `conda`
> no resuelve `conda.anaconda.org` (p.ej. mantenimiento de red), espera a que
> el internet del clúster esté disponible y reintenta.

### 3) Enviar el job

```bash
cd ~/estela/fase1
sbatch run_fase1.slurm
```

### 4) Monitorear

```bash
squeue --me
tail -f fase1_<JOBID>.out
```

### 5) Resultados (`~/estela/fase1/results/`)

| Archivo | Contenido |
|---------|-----------|
| `phase1_report.md` | **Reporte legible** del proceso completo |
| `phase1_selection.tsv` | Decisión final por bin |
| `01_tiara/tiara_bin_summary.tsv` | Composición de dominio por bin |
| `02_checkm2/quality_report.tsv` | Completitud/contaminación (CheckM2) |
| `phase1_selected_genomes/` | FASTA de los genomas seleccionados |
| `versions.txt` | Versiones de herramientas y umbrales usados |

Para traerlo a tu laptop (en tu laptop / WSL):

```bash
scp USUARIO@khipu.utec.edu.pe:~/estela/fase1/results/phase1_report.md /mnt/c/Tesis-EstelaRaimondi/
```

> **En este repositorio**, los resultados de texto de esta fase están versionados
> en [`results/fase1/`](results/fase1/). Las rutas de la tabla de arriba son las
> del clúster (donde el job escribe); en el repo llevan el prefijo `phase1_`.

---

## Recursos SLURM (dentro de tus límites: 32 cores, 98 GB, 24 h)

`run_fase1.slurm` pide **partición `standard`, 32 CPU, 90 GB RAM, 24 h** (sin GPU).
Puedes bajar `--time` para reducir la espera en cola.

## Reproducibilidad

- Versiones exactas de Tiara/CheckM2 en `envs/*.lock.yml` (generadas por el setup)
  y en `results/versions.txt` de cada corrida.
- Umbrales fijados por parámetros (no hardcodeados): min contig 3000 bp,
  completitud > 70 %, contaminación < 5 %.
- `bin_rename_map.tsv` documenta cada bin con su **muestra de origen** (49–72),
  nº de contigs y longitud — insumo directo para la prevalencia espacial de la Fase 2.

> Nota: el umbral (>70 %, <5 %) es criterio propio del proyecto, no una categoría
> MIMAG estándar (tesis 4.1.3.2). Tiara no separa bacteria/arquea; eso se resuelve
> en la Fase 2 con CheckM2 y GTDB-Tk.
