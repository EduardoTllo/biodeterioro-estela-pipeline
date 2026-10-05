# Fase 3: analisis pangenomico comparativo (OE3)

Esta fase toma los linajes bacterianos definidos en la Fase 2 (bins de
metagenomica de la Estela de Raimondi, clasificados con GTDB-Tk y agrupados con
dRep) y, para cada especie con suficientes genomas publicos, construye un
pangenoma con esos genomas de referencia mas el bin propio. El resultado es la
posicion del bin dentro de la diversidad conocida de su especie y una lista
verificada de genes candidatos a exclusivos, que es la entrada de la Fase 4.

## Flujo

```
Fase 2: phase2_selection.tsv          metadata GTDB R220
                \                          /
   [1] Censo de genomas publicos por especie (QC)       03_censo_genomas.py
                          |
   [2] Descarga de referencias desde NCBI + verificacion 04_fetch_refs.sh
                          |
   [3] Clasificacion por prioridad, habitat y continente 04_clasificar_refs.py clasificar
                          |
   [4] Desreplicacion al 99 % ANI                        04_drep_refs.slurm
                          |
   [5] Ranking de especies y seleccion de referencias    04_clasificar_refs.py seleccionar
                          |
   [6] Anotacion homogenea (referencias + bins)          05_bakta_all.slurm
                          |
   [7] Pangenoma por especie + particion core/accesorio  06_panaroo.slurm
                          |
   [8] Arbol del core por especie                        06b_iqtree.slurm
                          |
   [9] Verificacion de genes exclusivos (F1-F6)          08_exclusivos_local.slurm
                          |                              08_exclusivos_remoto.sh
   [10] Exportacion de resultados y figuras              bin/export_fase3.sh, figuras/R/f3_*.R
```

Los pasos 1, 2, 3, 5 y la parte remota del 9 necesitan internet y se corren en
el nodo de login. Los pasos 4, 6, 7, 8 y la parte local del 9 se envian a SLURM.

## Requisitos

### Entorno de computo

- Cluster con **SLURM** y **conda** (o mamba). Los scripts cargan conda con
  `module load miniconda/3.0` (Lmod).
- Un nodo con salida a internet (normalmente el de login) para instalar
  programas, descargar bases de datos y consultar NCBI. Los nodos de computo
  pueden no tener internet.
- Recursos por trabajo: hasta 32 nucleos, 90 GB de RAM y 24 h (Panaroo e
  IQ-TREE). Bakta usa 8 nucleos y 24 GB por genoma.
- Disco: ~130 GB libres para la base de datos de Bakta, mas ~30 GB para
  genomas y resultados.

### Software

Lo instala `03_setup_fase3_khipu.sh` desde conda-forge y bioconda, con versiones
fijas y `python=3.11`:

| Entorno | Paquetes |
|---|---|
| `ncbi-datasets` | ncbi-datasets-cli |
| `bakta` | bakta 1.12.1 (base de datos v6.0, tipo full) |
| `panaroo` | panaroo 1.8.0 (incluye mafft) |
| `iqtree` | iqtree 3.1.3, snp-sites 2.5.1 |
| `blast` | blast 2.17.0, taxonkit 0.20.0 |
| `drep` | dRep, fastANI y Mash; **se reutiliza el entorno de la Fase 2** |

Ademas descarga:
- el conversor de GFF3 de Bakta a formato Prokka de Panaroo (`scripts/convert_bakta_to_prokka_gff.py`, tag v1.8.0);
- el NCBI taxdump;
- el metadata de GTDB R220 (`bac120_metadata_r220.tsv`).

Las figuras se generan en R (>= 4.5) con los paquetes de
[figuras/README.md](figuras/README.md), mas `micropan` y `phangorn`.

### Entradas de las fases anteriores

| Archivo | Lo produce |
|---|---|
| `phase2_selection.tsv` (taxonomia y linaje de cada genoma) | Fase 2 (`run_fase2.slurm`) |
| Bins seleccionados, `bin-<N>-<muestra>.fasta` | Fase 1 (`results/phase1_selected_genomes/`) |
| Bins crudos de todas las muestras, `bin-<N>-<muestra>.fasta` | Fase 1 (carpeta `bins/`) |
| `quality_report.tsv` de CheckM2 | Fase 1 (`results/02_checkm2/`) |

El nombre `bin-<N>-<muestra>` es obligatorio: se usa para saber de que muestra
viene cada bin.

## Configuracion

Todos los scripts trabajan sobre una carpeta de trabajo, por defecto
`$HOME/estela/fase3`. Para usar otra, exporta la variable antes de correr los
scripts y pasala a SLURM con `--export=ALL`:

```bash
export WORKDIR=/ruta/a/mi/fase3
```

Otras variables que se pueden cambiar sin editar los scripts:

| Variable | Valor por defecto | Usada por |
|---|---|---|
| `DBS_DIR` | `$HOME/dbs` | `03_setup_fase3_khipu.sh` (Bakta y taxdump) |
| `GTDB_RELEASE_DIR` | `$HOME/gtdbtk_data/release220` | `03_setup_fase3_khipu.sh` |
| `CHECKM2_REPORT` | `$HOME/estela/fase1/results/02_checkm2/quality_report.tsv` | `06_panaroo.slurm` |
| `TAXDUMP_DIR` | `$HOME/dbs/taxdump` | `08_exclusivos_remoto.sh` |
| `NR_DB`, `NT_DB` | `nr`, `core_nt` | `08_exclusivos_remoto.sh` |

Lo que es propio del cluster donde se desarrollo (Khipu, UTEC) y hay que editar
a mano en otro cluster:

- la linea `#SBATCH --partition=standard` de cada `.slurm`;
- las lineas `module purge` / `module load miniconda/3.0` de cada script, si tu
  cluster carga conda de otra forma.

**Aviso:** `03_setup_fase3_khipu.sh` reescribe `~/.condarc` para usar solo
conda-forge y bioconda si encuentra el canal `defaults`. Antes guarda una copia
como `~/.condarc.bak.<fecha>`.

## Instalacion

Estructura de trabajo, codigo y datos de entrada:

```bash
mkdir -p ~/estela/fase3/data/bins ~/estela/fase3/data/bins_libreria ~/estela/fase3/logs
```

```bash
git clone https://github.com/EduardoTllo/biodeterioro-estela-pipeline.git ~/estela/fase3/scripts
```

```bash
cp <fase1>/results/phase1_selected_genomes/*.fasta ~/estela/fase3/data/bins/
```

```bash
cp <fase1>/bins/*.fasta ~/estela/fase3/data/bins_libreria/
```

`<fase1>` es la carpeta de trabajo de la Fase 1 (por defecto `~/estela/fase1`).

Instalacion de entornos y bases de datos, en el nodo con internet. Tarda varias
horas por la base de datos de Bakta y conviene correrlo en `tmux` o `screen`:

```bash
bash ~/estela/fase3/scripts/03_setup_fase3_khipu.sh 2>&1 | tee ~/estela/fase3/logs/setup.log
```

El script es idempotente: si se interrumpe, se vuelve a correr y omite lo que ya
esta hecho. Al terminar imprime las versiones instaladas y deja los lock files
de todos los entornos en `scripts/envs/`.

## Ejecucion

Los comandos se corren desde `~/estela/fase3` (los `.slurm` escriben su registro
en `logs/`, relativo a esa carpeta). Los pasos de login necesitan conda activado:

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

### 1. Censo de genomas publicos (login)

```bash
python scripts/03_censo_genomas.py censo --selection <fase2>/results/phase2_selection.tsv --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py ncbi --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py tabla --outdir results/00_censo --umbral 15
```

`censo` no necesita internet; `ncbi` consulta los conteos vigentes de NCBI y es
opcional. El reporte queda en `results/00_censo/phase3_censo_report.md`.

### 2. Descarga de referencias (login)

```bash
bash scripts/04_fetch_refs.sh 2>&1 | tee logs/descarga.log
```

Descarga todas las referencias de las especies viables y verifica que la
longitud de cada FASTA coincida con la de GTDB. Resumen:
`results/01_descarga/descarga_resumen.tsv`.

Tambien aplica la lista de exclusiones manuales `metadata/exclusiones.tsv`
(columnas `accession`, `linaje_id`, `motivo`): esos genomas quedan con estado
`excluido` y no entran a la desreplicacion, al pangenoma ni a la busqueda del
filtro F4. Si se edita la lista despues de descargar, basta con repetir la
verificacion:

```bash
python scripts/04_clasificar_refs.py verificar --censo-dir results/00_censo --descarga-dir results/01_descarga --refs-dir data/refs --exclusiones scripts/metadata/exclusiones.tsv
```

### 3. Clasificacion de referencias (login)

```bash
python scripts/04_clasificar_refs.py clasificar --censo-dir results/00_censo --descarga-dir results/01_descarga --refs-dir data/refs --keywords scripts/metadata/habitat_keywords.tsv --paises scripts/metadata/paises_continentes.tsv --outdir results/02_drep
```

El habitat se asigna con las expresiones regulares de
`metadata/habitat_keywords.tsv` (gana la primera categoria que coincide) y el
continente con `metadata/paises_continentes.tsv`. Las fuentes de aislamiento y
paises que no se pudieron clasificar quedan en `results/02_drep/fuentes_unicas.tsv`
y `results/02_drep/paises_sin_mapear.tsv`.

### 4. Desreplicacion (SLURM)

```bash
N=$(tail -n +2 results/01_descarga/especies_para_drep.tsv | wc -l); sbatch --array=1-$N scripts/04_drep_refs.slurm
```

### 5. Ranking y seleccion de referencias (login)

```bash
python scripts/04_clasificar_refs.py seleccionar --descarga-dir results/01_descarga --drep-dir results/02_drep --refs-dir data/refs --bins-dir data/bins --checkm2 <fase1>/results/02_checkm2/quality_report.tsv --outdir results/03_seleccion
```

Elige las especies y sus referencias y escribe el manifiesto de Bakta. Reporte:
`results/03_seleccion/phase3_seleccion_report.md`. El script se detiene si dRep
no aplico la tabla de pesos (`extraW.tsv`) o si ninguna especie es viable.

### 6. Anotacion (SLURM)

Se recomienda una prueba piloto antes de la corrida completa: Bakta y Panaroo
sobre 3 referencias y el bin de la primera especie
(`results/03_seleccion/piloto_manifest.tsv`).

```bash
N=$(tail -n +2 results/03_seleccion/piloto_manifest.tsv | wc -l); sbatch --array=1-$N --export=ALL,MANIFEST=results/03_seleccion/piloto_manifest.tsv scripts/05_bakta_all.slurm
```

```bash
sbatch --export=ALL,PILOTO=1 scripts/06_panaroo.slurm
```

Corrida completa (como maximo 4 genomas a la vez; los ya anotados se omiten):

```bash
N=$(tail -n +2 results/03_seleccion/bakta_manifest.tsv | wc -l); sbatch --array=1-$N%4 scripts/05_bakta_all.slurm
```

Resumen de las anotaciones:

```bash
python scripts/07_particion.py bakta --manifest results/03_seleccion/bakta_manifest.tsv --bakta-dir results/04_bakta --out results/04_bakta_resumen.tsv
```

### 7. Pangenoma (SLURM)

```bash
N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l); sbatch --array=1-$N%1 scripts/06_panaroo.slurm
```

Convierte los GFF3 de Bakta, corre Panaroo y calcula la particion. Reporte por
especie: `results/05_panaroo/<especie>/particion/particion_report.md`.

### 8. Arbol del core (SLURM)

```bash
N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l); sbatch --array=1-$N%1 scripts/06b_iqtree.slurm
```

Si una especie supera las 24 h, se repite solo esa tarea sobre el alineamiento
de SNPs:

```bash
sbatch --array=<n> --export=ALL,MODO=snps scripts/06b_iqtree.slurm
```

### 9. Verificacion de genes exclusivos

Parte local (SLURM):

```bash
N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l); sbatch --array=1-$N scripts/08_exclusivos_local.slurm
```

Parte remota (login; BLAST contra NCBI, puede tardar horas y se puede retomar si
se interrumpe):

```bash
bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/exclusivos_remoto.log
```

Si `core_nt` no esta disponible en el BLAST remoto, se usa `nt`:

```bash
NT_DB=nt bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/exclusivos_remoto.log
```

### 10. Exportacion y figuras

En el cluster:

```bash
bash scripts/bin/export_fase3.sh
```

Luego se copia `~/estela/fase3/export_fase3/` a `results/fase3/` del repositorio
y, desde la raiz del repositorio:

```bash
Rscript figuras/R/f3_run_all.R
```

## Salidas

En `~/estela/fase3/results/`:

| Carpeta | Contenido |
|---|---|
| `00_censo/` | Censo por especie, accesiones candidatas y tabla de viabilidad |
| `01_descarga/` | Estado de la descarga por genoma y resumen por especie |
| `02_drep/` | Clasificacion de referencias y resultados de dRep por especie |
| `03_seleccion/` | Ranking de especies, referencias elegidas y manifiestos de Bakta |
| `04_bakta/<genoma>/` | Anotacion de Bakta |
| `05_panaroo/<especie>/` | Pangenoma de Panaroo y `particion/` (clasificacion de cada gen) |
| `06_iqtree/<especie>/` | Arbol del core |
| `07_exclusivos/<especie>/` | Verificacion de genes exclusivos |

`bin/export_fase3.sh` reune las tablas, reportes y arboles en
`export_fase3/` (sin FASTA, GFF ni alineamientos); esa carpeta es la que se
versiona en `results/fase3/`. Archivos principales:

| Archivo | Contenido |
|---|---|
| `phase3_viabilidad.tsv` | Especies con suficientes genomas publicos |
| `phase3_ranking_especies.tsv` | Ranking por genomas no redundantes |
| `phase3_referencias_finales.tsv` | Referencias de cada pangenoma, con habitat y continente |
| `<especie>/genes_bin.tsv` | Cada gen del bin con su familia y categoria (entrada de la Fase 4) |
| `<especie>/recuperacion_core.tsv` | Fraccion del core de referencias presente en cada genoma |
| `<especie>/core.treefile` | Arbol de maxima verosimilitud del core |
| `<especie>/exclusivos_verificados.tsv` | Genes exclusivos con su clase final y el motivo |

Figuras (`figuras/figs/`): curvas de acumulacion, posicion del bin en el
pangenoma, arbol del core por especie y embudo de genes exclusivos.

## Criterios y parametros

| Etapa | Criterio | Referencia |
|---|---|---|
| Especies censadas | Solo linajes con nombre de especie en GTDB; se excluyen los sin especie y los de nombre provisional (`sp<digitos>`) | |
| QC de referencias | Completitud >= 95 %, contaminacion <= 5 % (CheckM2), <= 300 contigs, solo aislados (sin MAGs) | |
| Fuente | Genomas listados en GTDB R220 y descargados de NCBI por accesion con version exacta; sin sustituciones | |
| Exclusiones manuales | Genomas del cluster GTDB cuyo nombre en NCBI pertenece a otro filo (riesgo de ensamblaje quimerico o contaminado); listados con su motivo en `metadata/exclusiones.tsv` | |
| Desreplicacion | dRep, ANI >= 99 % (fastANI), solo sobre referencias | |
| Viabilidad | >= 15 genomas no redundantes por especie; se toman las 3 especies con mas genomas | Gautreau et al. 2020; Guerra 2026 |
| Tope por especie | 50 referencias. Entran primero el representante de GTDB y la cepa tipo (peso +1000 en dRep), luego los aislados de sustratos petreos o aridos (+500); el resto se reparte rotando entre habitats y, dentro de cada habitat, entre continentes en orden alfabetico. Sin preferencia por region | |
| Anotacion | Bakta 1.12.1, BD v6.0, mismas opciones para todos los genomas; no se usan anotaciones de NCBI | Schwengers et al. 2021 |
| Pangenoma | Panaroo 1.8.0, `--clean-mode moderate`, demas parametros por defecto | Tonkin-Hill et al. 2020 |
| Particion | Frecuencias calculadas solo con referencias: core >= 95 % (sensibilidad 90 %), shell 15-95 %, cloud < 15 %, exclusivo = ausente en todas las referencias | Tettelin et al. 2005 |
| Apertura del pangenoma | Ley de Heaps (`micropan::heaps`); alpha < 1 indica pangenoma abierto | Tettelin et al. 2008 |
| Arbol | IQ-TREE 3, ModelFinder, 1000 replicas UFBoot, sobre el alineamiento del core | |
| Genes exclusivos | F1: descarta pseudogenes, < 100 aa y CDS a < 100 pb del borde de contig. F4: descarta si aparece (>= 80 % identidad y cobertura) en cualquier genoma de la especie. F2/F6: contigs sin genes de la especie se contrastan con core_nt y se descartan si su mejor hit es de otro genero. F3: marca genes presentes en otro bin de la misma muestra. F5: BLASTp contra nr para el origen probable | |

Referencias:

- Gautreau G, et al. (2020). PPanGGOLiN: depicting microbial diversity via a partitioned pangenome graph. *PLoS Comput Biol* 16(3):e1007732.
- Guerra A. (2026). The pangenome: a statistical model, not a fixed biological property. *Bioinform Adv* 6(1):vbag069.
- Schwengers O, et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 7(11):000685.
- Tettelin H, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*. *PNAS* 102(39):13950-13955.
- Tettelin H, et al. (2008). Comparative genomics: the bacterial pan-genome. *Curr Opin Microbiol* 11(5):472-477.
- Tonkin-Hill G, et al. (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol* 21:180.

## Limitaciones

- Los bins son MAGs (70-100 % de completitud) y las referencias tienen >= 95 %.
  La **presencia** de un gen en el bin es confiable; su **ausencia** no, salvo
  que se controle por la completitud del bin (`recuperacion_core.tsv`).
- Las referencias se limitan a GTDB R220; no incluyen genomas depositados
  despues.
- El arbol del core no corrige por recombinacion.

## Solucion de problemas

| Sintoma | Causa probable |
|---|---|
| Un trabajo de SLURM termina sin dejar registro | La carpeta `logs/` no existia al enviarlo, o no se envio desde `~/estela/fase3` |
| `module: command not found` | El cluster no usa Lmod: ajustar las lineas `module` de los scripts |
| `Error reading prokka input!` en Panaroo | Algun GFF no paso por el conversor; revisar `conversion_descartes.tsv` y el registro de esa especie |
| `seleccionar` se detiene por `extraW.tsv` | dRep no leyo la tabla de pesos; borrar `results/02_drep/<especie>/data_tables` y reenviar esa especie |
| IQ-TREE con estado `TIMEOUT` | Repetir esa especie con `MODO=snps` |

Todos los scripts omiten lo que ya esta hecho, asi que un paso fallido se puede
reenviar sin repetir el trabajo anterior. Para reenviar una sola especie o un
solo genoma, se usa `--array=<n>` con el numero de fila del archivo
correspondiente.
