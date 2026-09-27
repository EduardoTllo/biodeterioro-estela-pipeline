# Fase 3 en Khipu: analisis pangenomico comparativo (OE3)

Toma los linajes definidos en la Fase 2, mide cuantos genomas publicos de calidad
existen para cada especie, selecciona las especies viables y construye un
pangenoma por especie que ubica a los bins propios dentro de la diversidad
conocida del taxon. Termina con la lista de genes candidatos a exclusivos de la
Estela, verificados, que es la entrada de la Fase 4.

## Criterio de seleccion

**Disponibilidad genomica publica**, con un umbral de **>= 15 genomas no
redundantes** de calidad adecuada por especie, y las **3 primeras** especies
viables en el ranking (si pasan menos, se trabaja con las que haya).

Este criterio reemplaza al de prevalencia espacial que la Fase 2 usaba para
priorizar linajes. El motivo es de factibilidad: con 17 de 19 linajes
representados por un solo genoma no hay pangenoma "interno" posible, asi que cada
pangenoma se arma con referencias publicas mas el bin de la Estela. La
prevalencia espacial se sigue reportando como resultado de la Fase 2 y aparece
como columna en la tabla de viabilidad.

### De donde sale el 15

Gautreau et al. 2020 (*PLoS Comput Biol* 16(3):e1007732) recomiendan al menos 15
genomas con variacion genomica real para particionar un pangenoma de forma
confiable. Se usa como referencia bibliografica del umbral: la herramienta de ese
trabajo (PPanGGOLiN) no forma parte del flujo, que usa solo Panaroo. Por debajo
de esa cifra el paradigma core/accesorio sigue siendo calculable desde 5 genomas
(Tettelin et al. 2005, *PNAS* 102(39):13950-13955), pero pierde soporte. Una
version previa fijaba 10 genomas sin respaldo y fue corregida.

El umbral se aplica al conteo **post-QC y desreplicado**, no al bruto de NCBI:
los genomas casi identicos inflan el core y pueden hacer parecer cerrado un
pangenoma abierto (Guerra 2026, *Bioinform Adv* 6(1):vbag069).

## Decisiones de diseno

| Tema | Decision |
|---|---|
| Fuente de referencias | Lista de genomas: metadata de GTDB **R220** (mismo release que la Fase 2). FASTA: NCBI GenBank/RefSeq via NCBI Datasets, por accesion con version exacta, sin sustituciones. No se agregan genomas posteriores a R220 |
| QC de referencias | Completitud >= 95 %, contaminacion <= 5 % (CheckM2 del metadata), <= 300 contigs, solo aislados (sin MAGs) |
| Linajes censados | Solo los 12 con binomio real; quedan fuera 4 sin especie (L2, L4, L5, L10) y 3 con nombre placeholder de GTDB (L6, L7, L19) |
| Desreplicacion | dRep al 99 % ANI solo sobre las referencias (los bins no entran) |
| Tope de 50 | Por niveles: T0 representante de GTDB y cepa tipo; T1 Latinoamerica y el Caribe; T2 sustrato petreo o ambiente arido; T3 el resto, por rotacion entre habitats y continentes. Reglas fijas, sin azar |
| Anotacion | Bakta 1.12.1 + BD v6.0 full para TODO (referencias y bins), mismas opciones. Nunca GFF de PGAP |
| Pangenoma | Panaroo 1.8.0, `--clean-mode moderate` (strict elimina los genes presentes en un solo genoma, que es justamente lo que son los exclusivos) |
| Particion | Core >= 95 %, shell 15-95 %, cloud < 15 %, calculada **solo con las referencias**; sensibilidad con core >= 90 % |
| Accesorio | Shell + cloud; "solo cloud" queda como analisis de sensibilidad para la Fase 4 |
| Arbol | IQ-TREE 3 sobre el alineamiento del core de cada especie |
| Exclusivos | Seis filtros (F1-F6), ver mas abajo |
| Bins de baja completitud | Entran con advertencia; se reporta su completitud junto a cada perfil |

## Etapas y scripts

| Etapa | Script | Donde | Salida principal |
|---|---|---|---|
| 0 Setup | `03_setup_fase3_khipu.sh` | login | entornos, BD Bakta, conversor de Panaroo, taxdump, metadata GTDB, lock files |
| 3.0-3.2 Censo | `03_censo_genomas.py` (`censo`, `ncbi`, `tabla`) | login | `results/00_censo/phase3_viabilidad.tsv` |
| 3.2b Descarga | `04_fetch_refs.sh` (usa `04_clasificar_refs.py listas` y `verificar`) | login | `data/refs/<especie>/*.fna`, `results/01_descarga/` |
| 3.2b Clasificacion | `04_clasificar_refs.py clasificar` | login | `results/02_drep/refs_clasificadas.tsv` + entradas de dRep |
| 3.2b dRep 99 % | `04_drep_refs.slurm` | SLURM | `results/02_drep/<especie>/data_tables/` |
| 3.2b Seleccion | `04_clasificar_refs.py seleccionar` | login | `results/03_seleccion/` (ranking, referencias, manifiestos) |
| 3.3 Anotacion | `05_bakta_all.slurm` + `07_particion.py bakta` | SLURM | `results/04_bakta/<genoma>/` |
| 3.4 Pangenoma | `06_panaroo.slurm` (llama a `07_particion.py particion`) | SLURM | `results/05_panaroo/<especie>/` |
| 3.4 Arbol | `06b_iqtree.slurm` | SLURM | `results/06_iqtree/<especie>/core.treefile` |
| 3.5 Exclusivos F1-F4 | `08_exclusivos_local.slurm` (usa `08_exclusivos.py`) | SLURM | `results/07_exclusivos/<especie>/` |
| 3.5 Exclusivos F5-F6 | `08_exclusivos_remoto.sh` | login | `exclusivos_verificados.tsv` |
| Exportar | `bin/export_fase3.sh` | login | `export_fase3/` -> `results/fase3/` del repo |
| Figuras | `figuras/R/f3_run_all.R` | laptop | `figuras/figs/fig_f3_*`, `figuras/tablas/tab_f3_*` |

Todos los scripts de Python usan solo la biblioteca estandar. Los parametros
estan como variables al inicio de cada script.

## Paso a paso

Rutas en Khipu:

```
~/estela/fase3/
|-- scripts/            # copia de este repo (incluye bin/, metadata/)
|-- data/
|   |-- bins/           # bins de los linajes censados (de la Fase 1)
|   |-- bins_libreria/  # los 97 bins crudos con nombre original (filtro F3)
|   `-- refs/<especie>/ # FASTA descargados
|-- results/            # 00_censo ... 07_exclusivos
`-- logs/               # salidas de SLURM
~/dbs/bakta/db/         # BD de Bakta v6.0 full
~/dbs/taxdump/          # NCBI taxdump
```

Las reglas de Khipu de siempre aplican: solo el nodo de login tiene internet;
operaciones largas del login dentro de `tmux`; los `sbatch` se envian desde
`~/estela/fase3` y la carpeta `logs/` debe existir antes.

### 0) Subir el codigo y los datos propios (desde el laptop, WSL)

```bash
rsync -avh --exclude .git --exclude results /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/ eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/scripts/
```

```bash
rsync -avh /mnt/c/Tesis-EstelaRaimondi/Datos/bins_clean/ eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/data/bins_libreria/
```

En Khipu, los bins seleccionados de la Fase 1:

```bash
mkdir -p ~/estela/fase3/data/bins ~/estela/fase3/logs && cp ~/estela/fase1/results/phase1_selected_genomes/*.fasta ~/estela/fase3/data/bins/
```

### 1) Setup (login, en tmux)

```bash
bash ~/estela/fase3/scripts/03_setup_fase3_khipu.sh
```

Checkpoint C0: los cinco entornos responden, `~/dbs/bakta/db/version.json` dice
`"major": 6` y hay lock files en `scripts/envs/` (copiarlos al repo, carpeta
`envs/`).

### 2) Censo (login)

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

Los pasos de Python de login (2, 4 y 6) usan el python de este entorno.

```bash
python scripts/03_censo_genomas.py censo --selection ~/estela/fase2/results/phase2_selection.tsv --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py ncbi --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py tabla --outdir results/00_censo --umbral 15
```

`tabla` deja pasar a la descarga **todas** las especies viables (`--top 0`, por
defecto): el top-3 se decide despues del dRep, sobre genomas no redundantes.

Checkpoint C1: 12 linajes censados; L1 mapea a NCBI *Massilia timonae*; L3, L11 y
L18 con `conteo_ncbi_comparable = no`; revisar a mano todo `acuerdo_mapeo < 0.90`.
Llevar `phase3_censo_report.md` al asesor.

### 3) Descarga y verificacion (login, en tmux)

```bash
bash scripts/04_fetch_refs.sh
```

Descarga todas las candidatas post-QC de cada especie viable (el conjunto
completo se usa luego en el filtro F4), las renombra a `<ID>.fna` y controla
que la longitud de cada FASTA coincida con `genome_size` de GTDB (+/- 1 %).

Checkpoint C2: en `results/01_descarga/descarga_estado.tsv`, descargados +
fallidos = solicitados; revisar `tamano_discrepante` y `version_distinta`.

### 4) Clasificacion por prioridad, habitat y continente (login)

```bash
python scripts/04_clasificar_refs.py clasificar --censo-dir results/00_censo --descarga-dir results/01_descarga --refs-dir data/refs --keywords scripts/metadata/habitat_keywords.tsv --paises scripts/metadata/paises_continentes.tsv --outdir results/02_drep
```

Checkpoint C3 (revision manual): leer `results/02_drep/fuentes_unicas.tsv`
(cada fuente de aislamiento con su categoria) y `paises_sin_mapear.tsv`. Si hace
falta, ampliar `metadata/habitat_keywords.tsv` (patrones regex; gana la primera
categoria que coincide, en el orden del archivo) o
`metadata/paises_continentes.tsv`, y volver a correr el comando.

### 5) dRep al 99 % (SLURM)

```bash
N=$(tail -n +2 results/01_descarga/especies_para_drep.tsv | wc -l); sbatch --array=1-$N scripts/04_drep_refs.slurm
```

### 6) Ranking definitivo y seleccion de referencias (login)

```bash
python scripts/04_clasificar_refs.py seleccionar --descarga-dir results/01_descarga --drep-dir results/02_drep --refs-dir data/refs --bins-dir data/bins --checkm2 ~/estela/fase1/results/02_checkm2/quality_report.tsv --outdir results/03_seleccion
```

El script se detiene si dRep no aplico la tabla de pesos extra (dRep la ignora
sin avisar si no puede leerla) o si ninguna especie llega a 15.

Checkpoint C4: `phase3_seleccion_report.md` con el top-3 y la composicion de
las referencias por nivel, habitat y continente. **Informar al asesor antes de
seguir.**

### 7) Piloto P0 (antes de anotar todo)

Tres referencias y los bins de la especie n.o 1 pasan por Bakta, la conversion
de GFF y Panaroo:

```bash
N=$(tail -n +2 results/03_seleccion/piloto_manifest.tsv | wc -l); sbatch --array=1-$N --export=ALL,MANIFEST=results/03_seleccion/piloto_manifest.tsv scripts/05_bakta_all.slurm
```

```bash
sbatch --export=ALL,PILOTO=1 scripts/06_panaroo.slurm
```

Y una prueba del BLAST remoto desde el login (1 proteina contra nr y 1 contig
contra core_nt). El piloto aprueba si Panaroo termina sin `Error reading prokka
input!` y el `gene_presence_absence.csv` tiene una columna por genoma.

### 8) Anotacion completa (SLURM)

```bash
N=$(tail -n +2 results/03_seleccion/bakta_manifest.tsv | wc -l); sbatch --array=1-$N%4 scripts/05_bakta_all.slurm
```

```bash
python scripts/07_particion.py bakta --manifest results/03_seleccion/bakta_manifest.tsv --bakta-dir results/04_bakta --out results/04_bakta_resumen.tsv
```

Checkpoint C5: ningun genoma sin anotacion; revisar las referencias marcadas
`cds_fuera_de_rango` (+/- 20 % de la mediana de CDS de su especie).

### 9) Pangenoma y arbol (SLURM)

```bash
N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l); sbatch --array=1-$N%1 scripts/06_panaroo.slurm
```

```bash
sbatch --array=1-$N%1 scripts/06b_iqtree.slurm
```

Si IQ-TREE no termina en 24 h para una especie, repetir esa tarea con
`--export=ALL,MODO=snps` (alineamiento de SNPs + `+ASC`).

Checkpoint C6: `results/05_panaroo/<especie>/particion/particion_report.md`.
Cada referencia debe recuperar >= 97 % del core de referencias y cada bin debe
quedar a <= 10 puntos de su completitud CheckM2.

### 10) Genes exclusivos

```bash
sbatch --array=1-$N scripts/08_exclusivos_local.slurm
```

Luego, en el login (dentro de tmux):

```bash
bash scripts/08_exclusivos_remoto.sh
```

Checkpoint C7: `results/07_exclusivos/<especie>/exclusivos_report.md`; revisar a
mano una muestra de 10 exclusivos verificados.

### 11) Traer los resultados y generar las figuras

```bash
bash ~/estela/fase3/scripts/bin/export_fase3.sh
```

Desde el laptop (WSL):

```bash
rsync -avh eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/export_fase3/ /mnt/c/Tesis-EstelaRaimondi/fase1_khipu/results/fase3/
```

```bash
Rscript figuras/R/f3_run_all.R
```

## Particion y posicion del bin

`07_particion.py particion` lee `gene_presence_absence.csv` de Panaroo y calcula
la frecuencia de cada familia **solo entre las referencias** (`f_ref`):

| Categoria | Regla |
|---|---|
| core | f_ref >= 0.95 (sensibilidad: 0.90) |
| shell | 0.15 <= f_ref < 0.95 |
| cloud | 0 < f_ref < 0.15 |
| exclusivo (candidato) | f_ref = 0 y presente en >= 1 bin de la Estela |

Asi, un gen que le falta al bin por incompletitud no saca a su familia del core.
El `--core_threshold 0.95` de Panaroo solo define el alineamiento del core para
el arbol. `genes_bin.tsv` lista cada gen del bin (locus tag de Bakta) con su
familia y su categoria: es la interfaz con la Fase 4, donde DRAM se corre sobre
las proteinas de Bakta del bin para que los identificadores coincidan.

## Verificacion de genes exclusivos

| Filtro | Regla | Resultado |
|---|---|---|
| F1 estructura | pseudogen de Bakta, < 100 aa o CDS a < 100 pb del borde del contig | descartado |
| F4 especie ampliada | `tblastn` contra todos los genomas post-QC de la especie; identidad >= 80 % y cobertura >= 80 % | descartado (presente en la especie) |
| F2 anclaje | el contig lleva al menos un gen de una familia presente en referencias; si no, es huerfano | huerfano -> F6 |
| F3 misma libreria | `blastn` contra los demas bins de la misma muestra; identidad >= 95 % y cobertura >= 80 % | bandera |
| F6 contig huerfano | `blastn` remoto contra core_nt; mejor hit de otro genero | descartado (contaminacion); sin hit = bandera |
| F5 origen | `blastp` remoto contra nr; clasifica el mejor hit (mismo genero, misma familia, mismo filo, otro filo, sin hit) | descartado si hay hit de la misma especie con >= 80 % / >= 80 % |

Clases finales en `exclusivos_verificados.tsv`: `verificado`, `con_bandera` o
`descartado` con su motivo. Bakta se corre sin `--partial`, asi que no predice
genes que se salen del borde de un contig; un fragmento con inicio alternativo si
puede aparecer, y eso es lo que caza F1.

## Salidas versionadas (`results/fase3/`)

| Archivo | Contenido |
|---|---|
| `phase3_viabilidad.tsv`, `phase3_censo_report.md` | tabla de viabilidad del censo |
| `phase3_genomas_candidatos.tsv` | accesiones post-QC con pais, fuente de aislamiento y designacion de tipo |
| `phase3_descarga_estado.tsv`, `phase3_descarga_resumen.tsv` | control de la descarga |
| `phase3_refs_clasificadas.tsv`, `phase3_fuentes_unicas.tsv` | nivel, habitat y continente de cada referencia |
| `phase3_ranking_especies.tsv`, `phase3_seleccion_report.md` | ranking definitivo por genomas no redundantes |
| `phase3_referencias_finales.tsv`, `phase3_composicion_referencias.tsv` | referencias del pangenoma y su composicion |
| `phase3_bakta_resumen.tsv` | control de la anotacion |
| `<especie>/drep_*.csv` | tablas de dRep |
| `<especie>/particion_*`, `genes_bin.tsv`, `recuperacion_core.tsv`, `resumen_*.tsv`, `pangenoma_refs.Rtab` | particion y posicion del bin |
| `<especie>/core.treefile`, `iqtree_*.txt` | arbol del core |
| `<especie>/exclusivos_verificados.tsv`, `embudo_exclusivos.tsv`, `exclusivos_report.md` | genes exclusivos |

## Advertencia de interpretacion

Los bins propios son MAGs, con 70-100 % de completitud, mientras que las
referencias entran con >= 95 %. Por construccion, a un MAG le van a *faltar*
genes del core. Por eso **las conclusiones sobre presencia de genes son validas y
las conclusiones sobre ausencia no lo son** sin controlar por la completitud del
MAG. Hay que reportar la completitud CheckM2 de cada bin junto a su perfil (la
tabla `recuperacion_core.tsv` lo hace). Esto condiciona tambien la Fase 4: un
marcador KO ausente en el MAG no prueba ausencia funcional.

## Estado

Codigo completo de las etapas 3.0 a 3.5, probado de punta a punta con datos
sinteticos (censo, descarga, clasificacion, seleccion, particion, filtros de
exclusivos y figuras). Pendiente de ejecucion con datos reales en Khipu.
