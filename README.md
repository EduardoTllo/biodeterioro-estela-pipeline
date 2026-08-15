# Pipeline de biodeterioro — Estela de Raimondi

Flujo bioinformático reproducible para caracterizar el **potencial de biodeterioro
y las señales de adaptación genómica** de los linajes bacterianos de la fracción
viable de la Estela de Raimondi, ejecutado en el clúster **Khipu (UTEC)**.

Proyecto final de carrera — Bioingeniería, UTEC. Datos de metagenómica *shotgun*
de la fracción cultivable, provistos por el Laboratorio de Genómica Microbiana (UPCH).

---

## Estado del proyecto

| Fase | Objetivo | Estado | Resultados |
|---|---|---|---|
| **1** | Control de calidad y selección de genomas (OE1) | ✅ Ejecutada | [`results/fase1/`](results/fase1/) |
| **2** | Asignación taxonómica y definición de linajes (OE2) | ✅ Ejecutada | [`results/fase2/`](results/fase2/) |
| **3** | Análisis pangenómico comparativo (OE3) | ⏳ Pendiente | — |
| **4** | Potencial metabólico y de riesgo (OE4) | ⏳ Pendiente | — |

### Resumen de resultados

- **Fase 1:** de **97 bins** crudos → 66 pasan el filtro por dominio (Tiara) →
  **21 genomas seleccionados** (CheckM2: completitud > 70 %, contaminación < 5 %).
- **Fase 2:** los 21 genomas → **19 linajes** (clústeres de especie a 95 % ANI),
  16 con especie confirmada. **Top-3 por prevalencia espacial:**
  **L1** *Telluria timonae* (muestras 50, 52) · **L2** *Corynebacterium* sp. (68, 69) ·
  **L3** *Bacillus_AB infantis* (61).

---

## Estructura del repositorio

```
.
├── README.md                    # este archivo
│
├── README_fase1.md              # guía completa de la Fase 1
├── 01_prep_bins.sh              #   prepara bins_clean/ + tabla de procedencia
├── 00_setup_khipu.sh            #   entornos conda (tiara, checkm2) + BD
├── run_fase1.slurm              #   job SLURM de la Fase 1
│
├── README_fase2.md              # guía completa de la Fase 2 (+ troubleshooting)
├── 02_setup_fase2_khipu.sh      #   entornos conda (gtdbtk, drep) + BD GTDB R220
├── run_fase2.slurm              #   job SLURM de la Fase 2
│
├── bin/
│   ├── summarize_phase1.py      # tablas y reporte legible de la Fase 1
│   └── summarize_phase2.py      # tablas y reporte legible de la Fase 2
│
├── metadata/
│   └── bin_provenance.tsv       # los 97 bins con su muestra de origen (49–72)
│
├── envs/                        # lock files de conda (versiones exactas)
│
└── results/                     # SOLO evidencia de texto (ver nota abajo)
    ├── fase1/
    │   ├── phase1_report.md                  # ← reporte legible
    │   ├── phase1_selection.tsv              # decisión por bin
    │   ├── phase1_checkm2_quality_report.tsv # completitud/contaminación
    │   ├── phase1_tiara_bin_summary.tsv      # composición de dominio
    │   └── phase1_versions.txt
    └── fase2/
        ├── phase2_report.md                  # ← reporte legible
        ├── phase2_selection.tsv              # taxonomía + linaje por genoma
        ├── phase2_gtdbtk_bac120_summary.tsv  # salida completa de GTDB-Tk
        ├── phase2_drep_{Cdb,Wdb,Ndb,Sdb}.csv # clústeres/linajes de dRep
        ├── phase2_gtdbtk_backbone_bac120.classify.tree
        ├── phase2_gtdbtk_tree_mapping.tsv
        ├── phase2_genomeInfo.csv
        └── phase2_versions.txt
```

> **Qué se versiona y qué no.** En `results/` se guardan **solo tablas, reportes y
> árboles** (texto), que son la evidencia reproducible del análisis. Los genomas
> (FASTA), las bases de datos y los archivos pesados —como el MSA de 230 MB de
> GTDB-Tk— quedan excluidos por `.gitignore`: viven en el clúster y en el
> almacenamiento local, no en git.

---

## Por dónde empezar

1. **Para entender qué se hizo:** lee los reportes
   [`results/fase1/phase1_report.md`](results/fase1/phase1_report.md) y
   [`results/fase2/phase2_report.md`](results/fase2/phase2_report.md).
2. **Para reproducir una fase:** sigue su README
   ([Fase 1](README_fase1.md) · [Fase 2](README_fase2.md)), que incluye el
   paso a paso de subida a Khipu, setup y envío del job.
3. **Si algo falla:** el [README de la Fase 2](README_fase2.md) tiene una sección
   de **resolución de problemas** con los errores reales encontrados en Khipu
   (versiones de GTDB-Tk, Python/pydantic, canales de conda, `HASH MISMATCH`).

## Convenciones

- **Nombres de genoma:** `bin-<N>-<muestra>` (p. ej. `bin-6-50` = bin 6 de la
  librería 50). El nombre **codifica la muestra de origen**, que es lo que permite
  calcular la *prevalencia espacial* — el criterio primario de importancia de un
  linaje, dado que no se dispone de las lecturas crudas.
- **Recursos SLURM:** partición `standard`, dentro de los límites de la cuenta
  (32 cores, 98 GB RAM, 24 h).
- **Entornos conda:** versiones **fijadas** (incluido `python=3.11`) y exportadas
  a `envs/*.lock.yml`. Ver la tabla de versiones en el README de cada fase.
