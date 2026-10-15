# Fases 1 y 2 con los bins eucariotas (anexo)

Repite el control de calidad (OE1) y la taxonomía y los linajes (OE2) con los
bins que Tiara clasificó como eucariotas. Todo lo del análisis eucariota vive en esta
carpeta (`eucariotas/`): código, lock files (`envs/`), resultados (`results/`) y
figuras (`figuras/`). Las tablas de las fases bacterianas no se tocan. El fin es comparar con las bacterias
(pedido del asesor) y dejar el análisis listo para un anexo o un paper.

## Qué cambia respecto a las bacterias

| Paso | Bacterias (Fases 1–2) | Eucariotas (este anexo) |
|---|---|---|
| Dominio | Tiara: se quedan los procariotas | Tiara: se quedan los 31 bins con mayoría eucariota (misma corrida, no se repite) |
| Calidad | CheckM2 > 70 % / < 5 % | **EukCC 2.2.0** > 70 % / < 5 %; **BUSCO 6.1.0** como complemento |
| Especie | GTDB-Tk + ANI (skani interno) ≥ 95 % / AF ≥ 65 % | **skani 0.3.2** contra ~5 600 genomas de referencia de hongos de NCBI, mismos umbrales |
| Segunda línea | — | **ITS** (barrnap + ITSx) contra UNITE y RefSeq ITS, umbrales de Vu et al. 2019 |
| Linaje | dRep 95 % ANI | dRep 95 % ANI (umbral operativo en hongos) |
| Árbol | NJ sobre `user_msa` bac120 | NJ sobre la matriz concatenada de genes BUSCO de copia única |
| Prevalencia | n.º de muestras por linaje | Igual |

Los 5 bins de mayoría `unknown` con fracción clasificada eucariota (bin-2-64,
1-58, 1-70, 11-71, 6-70) quedan fuera. Tienen entre el 52 y el 68 % de su
longitud en contigs de menos de 3 kb y contaminación bacteriana según CheckM2.

## Archivos

```
eucariotas/
├── 00_setup_euk_khipu.sh  # entornos y bases (etapas 'base' y 'refs'), nodo de login
├── run_fase1_euk.slurm    # Fase 1 eucariota: EukCC + BUSCO + selección
├── run_fase2_euk.slurm    # Fase 2 eucariota: skani + ITS + árbol BUSCO + dRep
├── bin/summarize_euk.py   # tablas y reportes (solo biblioteca estándar)
```

## Paso a paso

Los comandos marcados **Mac** van en la terminal de la Mac. Los marcados
**Khipu** van dentro de Khipu, en una sesión de `tmux`.

### 1) Subir el código (Mac)

```bash
cd ~/Documents/Tesis-EstelaRaimondi/fase1_khipu
```

```bash
git add eucariotas && git commit -m "Analisis eucariota: Fases 1 y 2 (EukCC, BUSCO, skani, ITS)" && git push
```

### 2) Bajar UNITE (Mac, a mano)

UNITE no tiene un enlace de descarga fijo.

1. Entra a <https://unite.ut.ee/repository.php>.
2. En **General FASTA release**, baja la versión más reciente para **Fungi** (la 10.0 o posterior).
3. Descomprime y quédate con el archivo `sh_general_release_dynamic_*.fasta`.
4. Súbelo a Khipu con este nombre exacto:

```bash
ssh eduardo.tello@khipu.utec.edu.pe "mkdir -p ~/dbs/its"
```

```bash
scp ~/Downloads/sh_general_release_dynamic_*.fasta eduardo.tello@khipu.utec.edu.pe:~/dbs/its/unite_fungi.fasta
```

Anota la versión y la fecha de UNITE, porque van en el anexo. Si te saltas este
paso, la Fase 2 usa solo RefSeq ITS.

### 3) Preparación base (Khipu, nodo de login, ~1 h)

```bash
tmux new -s euk
```

```bash
cd ~/estela/fase3/scripts && git pull && bash eucariotas/00_setup_euk_khipu.sh base
```

Crea los entornos `eukcc`, `busco` y `euktax`, y baja la BD de EukCC (v1.2),
los datasets `fungi_odb12.2` y `eukaryota_odb12.2` de BUSCO y RefSeq ITS. Al
final imprime las versiones y deja los lock files en `eucariotas/envs/`.

### 4) Enviar la Fase 1 (Khipu)

```bash
mkdir -p ~/estela/euk && cd ~/estela/euk && sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" ~/estela/fase3/scripts/eucariotas/run_fase1_euk.slurm
```

Tarda entre 6 y 15 h para los 31 bins; los grandes (40–56 Mb) son los más
lentos. Si el job se corta por tiempo, se reenvía con el mismo comando y
retoma donde quedó.

### 5) Bajar las referencias mientras corre la Fase 1 (Khipu, nodo de login)

En otra ventana de tmux (`Ctrl-b c`):

```bash
cd ~/estela/fase3/scripts && bash eucariotas/00_setup_euk_khipu.sh refs
```

Baja un genoma de referencia por especie de hongo de NCBI (~5 600 genomas,
~200 GB comprimidos; hay espacio de sobra en `/home`). Puede tardar varias
horas. Si se corta, vuelve a correr el mismo comando y continúa. Al terminar
revisa el conteo:

```bash
wc -l ~/dbs/euk_refs/refs_fasta.txt ~/dbs/euk_refs/refs_fungi.tsv
```

### 6) Revisar la Fase 1 (checkpoint)

```bash
cat ~/estela/euk/fase1/results/phase1_euk_report.md
```

Qué mirar:

- **Cuántos pasan.** Hacen falta al menos 3 para el árbol y 2 para dRep.
- **`01_eukcc/codigos_salida.tsv`.** Los códigos 201 (sin set de marcadores) y
  204 (sin proteínas) son esperables en los 12 bins pequeños (0,2–1,7 Mb).
- **El linaje de EukCC de los seleccionados.** Si alguno no es un hongo
  (`reino` distinto de Fungi), edita `BUSCO_ARBOL=eukaryota_odb12.2` en
  `run_fase2_euk.slurm` antes del paso 7.
- **BUSCO D (duplicados) alto con EukCC bajo en contaminación.** Anótalo:
  puede ser una mezcla de dos cepas o un genoma con duplicaciones reales.

### 7) Enviar la Fase 2 (Khipu)

Requiere la Fase 1 terminada y la etapa `refs` completa.

```bash
cd ~/estela/euk && sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" ~/estela/fase3/scripts/eucariotas/run_fase2_euk.slurm
```

La primera vez arma el sketch de skani de las referencias (~1–2 h). El resto
tarda menos de 2 h.

### 8) Revisar la Fase 2 (checkpoint)

```bash
cat ~/estela/euk/fase2/results/phase2_euk_report.md
```

Qué mirar:

- **Nivel de asignación por bin:** `especie` (ANI), `especie (hipotesis ITS)`,
  `genero`, etc.
- **Columna `Conflicto`.** Marca los casos en que el género de skani y el del
  ITS no coinciden; hay que revisarlos a mano.
- **Cuántos bins tienen ITS.** Es normal que muchos no tengan, porque el ADNr
  multicopia suele quedar fuera de los bins.
- **En *Fusarium*, *Penicillium* y *Aspergillus*,** una «especie» por ITS es en
  realidad un complejo de especies. Para el anexo, vale más la del ANI.

### 9) Traer los resultados (Mac)

Copia solo tablas y reportes; los FASTA y las salidas intermedias se quedan
en Khipu.

```bash
rsync -avh --prune-empty-dirs --exclude='bins_euk/' --exclude='phase1_euk_selected_genomes/' --exclude='01_eukcc/*/' --exclude='02_busco/' --exclude='gff/' --exclude='genes/' --exclude='aln/' --exclude='trim/' --exclude='04_drep/data/' --include='*/' --include='*.tsv' --include='*.md' --include='*.txt' --include='*.csv' --include='*.faa' --include='itsx.full.fasta' --exclude='*' eduardo.tello@khipu.utec.edu.pe:~/estela/euk/ ~/Documents/Tesis-EstelaRaimondi/fase1_khipu/eucariotas/results/
```

### 10) Figuras y tablas (Mac, con Claude)

Con los resultados en `eucariotas/results/`, se adaptan los scripts de
`figuras/R/` y la versión eucariota de cada figura queda en `eucariotas/figuras/`:

| Bacterias | Eucariotas |
|---|---|
| Fig. 7.1 embudo | 97 → 31 eucariotas → seleccionados por EukCC → linajes |
| Fig. 7.2 calidad | Completitud vs contaminación (EukCC); tamaño vs GC; BUSCO por bin |
| Fig. 7.3 composición | Filo, orden y género (taxonomía NCBI/UNITE) |
| Fig. 7.4 árbol | NJ sobre la matriz BUSCO (`phase2_euk_busco_concat.faa`) |
| Fig. 7.5 novedad | ANI vs AF de skani |
| Fig. 7.6 matriz | Similitud Mash de dRep (`Mdb.csv`) |
| Fig. 7.7 prevalencia | Linaje × muestra |
| S1 | Completitud y duplicación BUSCO (en lugar de `msa_percent`) |

Además, una tabla de contraste con el Kraken2 del centro de secuenciación
(local, `Datos/Shotgun Analysis/<muestra>/<bin>/*_taxonomy.tsv`).

## Archivos de salida principales

| Archivo | Contenido |
|---|---|
| `fase1/results/phase1_euk_selection.tsv` | Una fila por bin: tamaño, N50, GC, EukCC, BUSCO, linaje de EukCC, decisión y motivo |
| `fase1/results/phase1_euk_report.md` | Embudo y detalle legible |
| `fase2/results/phase2_euk_selection.tsv` | Una fila por genoma: taxonomía en 8 rangos, nivel y método de asignación, skani, ITS, linaje |
| `fase2/results/phase2_euk_linajes.tsv` | Una fila por linaje (E1, E2…): prevalencia, muestras, representante |
| `fase2/results/phase2_euk_busco_concat.faa` | Matriz del árbol |
| `fase2/results/04_drep/data_tables/` | `Cdb`, `Wdb`, `Mdb`, `Ndb` de dRep |

## Problemas previsibles

- **EukCC tarda mucho en un bin.** El job es reanudable: reenvíalo. Si un bin se
  cuelga siempre, revisa su `.log` en `01_eukcc/`.
- **`QOSMaxSubmitJobPerUserLimit`.** Hay otro job de la tesis en cola; espera o
  usa la cuenta `pregrado` (máximo 8 h, no alcanza para la Fase 1 completa).
- **skani no reporta un bin.** No hay ninguna referencia con ANI mayor de ~80 %.
  Esa especie probablemente no tiene genoma en NCBI. Queda el ITS o el linaje de
  EukCC.
- **ITSx no encuentra nada.** Es normal en MAGs. El reporte lo indica.
