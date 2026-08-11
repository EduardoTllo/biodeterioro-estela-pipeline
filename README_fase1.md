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

- **Input correcto:** los **97 bins individuales** (`C:\Tesis-EstelaRaimondi\Datos\MAGs_all`),
  con contigs `contig_*`. El QC de CheckM2 se hace **por bin**.
- ⚠️ **No usar** `MAGs_v2/` ni `MAGs_v2.zip`: esos son 24 ensamblajes **por muestra**
  (contigs `k141_*`), no bins; un QC sobre ellos daría contaminación altísima.

### Dos avisos importantes sobre `MAGs_all`

1. **Nombres rotos.** Al consolidar los bins, Windows renombró duplicados como
   `bin-1.fasta 10`, `bin-1.fasta 22`, etc. El script `01_prep_bins.sh` los
   normaliza a `bin_0001.fasta … bin_0097.fasta` y guarda `bin_rename_map.tsv`.
2. **Trazabilidad a la muestra perdida.** Al aplanar la carpeta se perdió de qué
   librería (49–72) vino cada bin. **No afecta la Fase 1** (QC independiente por
   bin), pero **sí la Fase 2** (prevalencia espacial). Antes de la Fase 2 hay que
   recuperar el mapeo bin→muestra del binning original (MetaBAT2, UPCH).

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

### 0) Preparar los bins (una vez, en tu laptop con Git Bash)

```bash
cd /c/Tesis-EstelaRaimondi/fase1_khipu
bash 01_prep_bins.sh "/c/Tesis-EstelaRaimondi/Datos/MAGs_all" "/c/Tesis-EstelaRaimondi/Datos/bins_clean"
# Empaquetar para subir (FASTA comprime bien: ~1.1 GB -> ~300 MB)
tar czf /c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz -C "/c/Tesis-EstelaRaimondi/Datos/bins_clean" .
```

### 1) Subir a Khipu (PowerShell/Git Bash; reemplaza `USUARIO`)

```bash
ssh USUARIO@khipu.utec.edu.pe "mkdir -p ~/estela/fase1/bins ~/estela/fase1/scripts"
# scripts
scp 00_setup_khipu.sh run_fase1.slurm USUARIO@khipu.utec.edu.pe:~/estela/fase1/
scp bin/summarize_phase1.py           USUARIO@khipu.utec.edu.pe:~/estela/fase1/scripts/
# datos (tarball)
scp /c/Tesis-EstelaRaimondi/Datos/bins_clean.tar.gz USUARIO@khipu.utec.edu.pe:~/estela/fase1/
```

### 2) Preparar el entorno (en el NODO DE LOGIN, con internet)

```bash
ssh USUARIO@khipu.utec.edu.pe
cd ~/estela/fase1
tar xzf bins_clean.tar.gz -C bins      # deja bin_0001.fasta ... en bins/
bash 00_setup_khipu.sh                  # crea entornos + descarga DB (~3 GB)
```

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

```bash
scp USUARIO@khipu.utec.edu.pe:~/estela/fase1/results/phase1_report.md .
```

---

## Recursos SLURM (dentro de tus límites: 32 cores, 98 GB, 24 h)

`run_fase1.slurm` pide **partición `standard`, 32 CPU, 90 GB RAM, 24 h** (sin GPU).
Puedes bajar `--time` para reducir la espera en cola.

## Reproducibilidad

- Versiones exactas de Tiara/CheckM2 en `envs/*.lock.yml` (generadas por el setup)
  y en `results/versions.txt` de cada corrida.
- Umbrales fijados por parámetros (no hardcodeados): min contig 3000 bp,
  completitud > 70 %, contaminación < 5 %.
- `bin_rename_map.tsv` documenta el renombrado de cada bin.

> Nota: el umbral (>70 %, <5 %) es criterio propio del proyecto, no una categoría
> MIMAG estándar (tesis 4.1.3.2). Tiara no separa bacteria/arquea; eso se resuelve
> en la Fase 2 con CheckM2 y GTDB-Tk.
