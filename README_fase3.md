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
                          |                              06c_ani_bins.slurm (ANI bin vs especie)
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
  IQ-TREE). Bakta usa 10 nucleos y 32 GB por lote.
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
[figuras/README.md](figuras/README.md), mas `phangorn`.

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

### Carpetas

Todo se hace dentro de una carpeta de trabajo, `~/estela/fase3`. Los comandos
de esta guia suponen ademas que la Fase 1 y la Fase 2 se corrieron en
`~/estela/fase1` y `~/estela/fase2`. Si usas otras rutas, cambialas en los
comandos donde aparecen.

Variables que se pueden cambiar sin editar los scripts:

| Variable | Valor por defecto | Usada por |
|---|---|---|
| `WORKDIR` | `$HOME/estela/fase3` | todos los scripts (pasarla a SLURM con `--export=ALL`) |
| `DBS_DIR` | `$HOME/dbs` | `03_setup_fase3_khipu.sh` (Bakta y taxdump) |
| `GTDB_RELEASE_DIR` | `$HOME/gtdbtk_data/release220` | `03_setup_fase3_khipu.sh` |
| `CHECKM2_REPORT` | `$HOME/estela/fase1/results/02_checkm2/quality_report.tsv` | `06_panaroo.slurm` |
| `TAXDUMP_DIR` | `$HOME/dbs/taxdump` | `08_exclusivos_remoto.sh` |
| `NR_DB`, `NT_DB` | `nr`, `core_nt` | `08_exclusivos_remoto.sh` |
| `ESPERA_MAX_H` | `24` | `08_exclusivos_remoto.sh` (horas maximas de espera por busqueda en NCBI) |

### Cuenta de SLURM, correo y limites de cola

Los comandos de esta guia son los que se usaron en Khipu (UTEC): todos los
envios llevan `--account=tesis` y piden un correo al terminar o fallar. En otro
cluster cambia `tesis` por tu cuenta y la linea `#SBATCH --partition=standard`
de cada `.slurm`; si tu cluster no usa Lmod, ajusta tambien las lineas
`module purge` / `module load miniconda/3.0` de los scripts.

El correo se guarda una sola vez en `~/.bashrc` como la variable `CORREO`
(reemplaza la direccion por la tuya):

```bash
echo 'export CORREO=nombre.apellido@universidad.edu' >> ~/.bashrc && source ~/.bashrc
```

La cuenta `tesis` de Khipu admite **5 trabajos enviados y 3 corriendo** a la
vez (32 CPU, 98 GB, 24 h por trabajo). En un *job array* cada tarea cuenta como
un trabajo; por eso algunos pasos se envian por tramos o con `%1` (una tarea a
la vez). Si al enviar aparece `QOSMaxSubmitJobPerUserLimit`, espera a que
terminen trabajos anteriores. Para ver cuantos tienes en cola:

```bash
squeue --me
```

Para ver los limites de tu cuenta:

```bash
sacctmgr -n show qos format=Name%20,MaxSubmitPU,MaxJobsPU,MaxTRESPU%40,MaxWall
```

**Aviso:** `03_setup_fase3_khipu.sh` reescribe `~/.condarc` para usar solo
conda-forge y bioconda si encuentra el canal `defaults`. Antes guarda una copia
como `~/.condarc.bak.<fecha>`.

## Instalacion

Todo en el nodo de login.

**1. Carpetas, codigo y bins.**

```bash
mkdir -p ~/estela/fase3/data/bins ~/estela/fase3/data/bins_libreria ~/estela/fase3/logs
```

```bash
git clone https://github.com/EduardoTllo/biodeterioro-estela-pipeline.git ~/estela/fase3/scripts
```

```bash
cp ~/estela/fase1/results/phase1_selected_genomes/*.fasta ~/estela/fase3/data/bins/
```

```bash
cp ~/estela/fase1/bins/*.fasta ~/estela/fase3/data/bins_libreria/
```

`data/bins/` recibe los bins seleccionados en la Fase 1 (21 en este proyecto) y
`data/bins_libreria/` todos los bins crudos (97), que se usan en el filtro F3.
Comprobacion:

```bash
ls ~/estela/fase3/data/bins | wc -l; ls ~/estela/fase3/data/bins_libreria | wc -l
```

Para actualizar el codigo mas adelante:

```bash
cd ~/estela/fase3/scripts && git pull
```

**2. Comprobar la tabla de la Fase 2.** El censo necesita que
`phase2_selection.tsv` tenga la clasificacion de GTDB de cada bin. Una version
antigua del parser de la Fase 2 la dejaba vacia y el censo devolvia 0 linajes.
Revisa que la columna de clasificacion tenga texto (`d__Bacteria;p__...`):

```bash
head -n 3 ~/estela/fase2/results/phase2_selection.tsv
```

Si esta vacia, usa la version del repositorio (guardando la anterior):

```bash
mv ~/estela/fase2/results/phase2_selection.tsv ~/estela/fase2/results/phase2_selection.sin_taxonomia.bak.tsv && cp ~/estela/fase3/scripts/results/fase2/phase2_selection.tsv ~/estela/fase2/results/
```

**3. Entornos y bases de datos.** Tarda varias horas por la base de datos de
Bakta; correrlo dentro de `tmux` para que no se corte si se cierra la conexion
(`tmux new -s fase3`; salir sin cortarlo con `Ctrl+b` y luego `d`; volver con
`tmux attach -t fase3`):

```bash
bash ~/estela/fase3/scripts/03_setup_fase3_khipu.sh 2>&1 | tee ~/estela/fase3/logs/setup.log
```

El script se puede volver a correr si se interrumpe: omite lo que ya esta
hecho. Deja los lock files de los entornos en `scripts/envs/`.

**4. Comprobar los entornos.** Cada linea debe imprimir una version. En Khipu
el entorno `blast` llego a quedar solo con python, y el error recien aparecio
en el paso 9 (`makeblastdb: command not found`):

```bash
module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda run -n ncbi-datasets datasets --version && conda run -n bakta bakta --version && conda run -n panaroo panaroo --version && conda run -n iqtree iqtree3 --version | head -n 1 && conda run -n blast makeblastdb -version | head -n 1 && conda run -n blast taxonkit version && conda run -n drep dRep -h | head -n 2
```

Si falta alguno, se recrea ese entorno. Por ejemplo, `blast`:

```bash
conda env remove -n blast -y; conda create -y -n blast --override-channels -c conda-forge -c bioconda blast=2.17.0 taxonkit=0.20.0 python=3.11
```

## Ejecucion

Todos los comandos se corren desde `~/estela/fase3`: los `.slurm` escriben su
registro en `logs/` relativo a esa carpeta. Los pasos 1, 2, 3 y 5 corren en el
login con el entorno `ncbi-datasets` activado:

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

Para revisar un trabajo de SLURM terminado (estado, tiempo y memoria), cambia
el nombre segun el paso (`f3_drep`, `f3_bakta`, `f3_panaroo`, `f3_iqtree`,
`f3_exclus`):

```bash
sacct -u $USER --name=f3_drep -X --format=JobID%14,State,ExitCode,Elapsed,MaxRSS
```

### 1. Censo de genomas publicos (login)

```bash
python scripts/03_censo_genomas.py censo --selection ~/estela/fase2/results/phase2_selection.tsv --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py ncbi --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py tabla --outdir results/00_censo --umbral 15
```

`censo` no necesita internet; `ncbi` consulta los conteos vigentes de NCBI y es
opcional. Reporte: `results/00_censo/phase3_censo_report.md`. Antes de
descargar, conviene revisar los nombres de NCBI dentro de cada especie: los
genomas rotulados con otro filo se agregan a `scripts/metadata/exclusiones.tsv`
(columnas `accession`, `linaje_id`, `motivo`).

### 2. Descarga de referencias (login)

```bash
bash scripts/04_fetch_refs.sh 2>&1 | tee logs/descarga.log
```

Descarga las referencias de las especies viables, verifica que la longitud de
cada FASTA coincida con la de GTDB y aplica `metadata/exclusiones.tsv` (esos
genomas quedan como `excluido` y no entran a la desreplicacion, al pangenoma
ni al filtro F4). Resumen: `results/01_descarga/descarga_resumen.tsv`. Si se
edita la lista de exclusiones despues de descargar, basta con repetir la
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
continente con `metadata/paises_continentes.tsv`. Lo que no se pudo clasificar
queda en `results/02_drep/fuentes_unicas.tsv` y
`results/02_drep/paises_sin_mapear.tsv`.

### 4. Desreplicacion (SLURM)

Una tarea por especie viable (16 CPU, 32 GB, hasta 6 h). Numero de especies:

```bash
tail -n +2 results/01_descarga/especies_para_drep.tsv | wc -l
```

Con 7 especies y el limite de 5 envios, va en dos tramos. Primero las 5
primeras, de a 2 a la vez:

```bash
sbatch --account=tesis --array=1-5%2 --mail-type=END,FAIL --mail-user="$CORREO" scripts/04_drep_refs.slurm
```

Cuando terminen, las restantes:

```bash
sbatch --account=tesis --array=6-7 --mail-type=END,FAIL --mail-user="$CORREO" scripts/04_drep_refs.slurm
```

Con 5 especies o menos basta un envio (`--array=1-N`).

### 5. Ranking y seleccion de referencias (login)

```bash
python scripts/04_clasificar_refs.py seleccionar --descarga-dir results/01_descarga --drep-dir results/02_drep --refs-dir data/refs --bins-dir data/bins --checkm2 ~/estela/fase1/results/02_checkm2/quality_report.tsv --outdir results/03_seleccion
```

Elige las 3 especies con mas genomas no redundantes y hasta 50 referencias por
especie, y escribe los manifiestos de Bakta. Reporte:
`results/03_seleccion/phase3_seleccion_report.md`. El script se detiene si dRep
no aplico la tabla de pesos (`extraW.tsv`) o si ninguna especie es viable.
Desde aqui los pasos usan 3 especies (`--array=1-3`).

### 6. Anotacion con Bakta (SLURM)

`05_bakta_all.slurm` trabaja por lotes: cada tarea anota varios genomas
seguidos (10 CPU y 32 GB por lote). El numero de tareas del array debe ser
igual a `LOTES`.

**Piloto** (3 referencias y el bin de la primera especie, un lote):

```bash
sbatch --account=tesis --array=1 --export=ALL,LOTES=1,MANIFEST=results/03_seleccion/piloto_manifest.tsv --mail-type=END,FAIL --mail-user="$CORREO" scripts/05_bakta_all.slurm
```

Cuando termine, Panaroo sobre el piloto (unos minutos):

```bash
sbatch --account=tesis --export=ALL,PILOTO=1 --mail-type=END,FAIL --mail-user="$CORREO" scripts/06_panaroo.slurm
```

El piloto esta bien si el registro de Panaroo no muestra `Error reading prokka
input`, `conversion_descartes.tsv` no tiene CDS descartados y la recuperacion
del core en el bin (`particion_report.md`) es parecida a su completitud.

**Corrida completa** (118 genomas en 3 lotes en paralelo, unas 4-5 h; los
genomas del piloto se saltan):

```bash
sbatch --account=tesis --array=1-3 --export=ALL,LOTES=3 --mail-type=END,FAIL --mail-user="$CORREO" scripts/05_bakta_all.slurm
```

Resumen y control de las anotaciones (alerta si un genoma se aleja mas de 20 %
de la mediana de CDS de su especie):

```bash
python scripts/07_particion.py bakta --manifest results/03_seleccion/bakta_manifest.tsv --bakta-dir results/04_bakta --out results/04_bakta_resumen.tsv
```

### 7. Pangenoma con Panaroo (SLURM)

Una especie a la vez (32 CPU y 90 GB cada una; en este proyecto, entre 12 y
50 min por especie):

```bash
sbatch --account=tesis --array=1-3%1 --mail-type=END,FAIL --mail-user="$CORREO" scripts/06_panaroo.slurm
```

Convierte los GFF3 de Bakta, corre Panaroo y calcula la particion. Reporte por
especie: `results/05_panaroo/<especie>/particion/particion_report.md`. Revisar
las alertas de `recuperacion_core.tsv`.

### 8. Arbol del core con IQ-TREE (SLURM)

Despues de que termine Panaroo (Panaroo e IQ-TREE juntos serian 6 trabajos,
mas que el limite de 5):

```bash
sbatch --account=tesis --array=1-3%1 --mail-type=END,FAIL --mail-user="$CORREO" scripts/06b_iqtree.slurm
```

Si una especie termina en `TIMEOUT`, se repite solo esa tarea (su numero de
fila, por ejemplo 2) sobre el alineamiento de SNPs:

```bash
sbatch --account=tesis --array=2 --export=ALL,MODO=snps --mail-type=END,FAIL --mail-user="$CORREO" scripts/06b_iqtree.slurm
```

### 9. Verificacion de genes exclusivos

**Parte local, F1-F4 (SLURM).** Puede correr al mismo tiempo que IQ-TREE si en
la cola hay 2 trabajos o menos (`squeue --me`). Tarda menos de un minuto por
especie:

```bash
sbatch --account=tesis --array=1-3 --mail-type=END,FAIL --mail-user="$CORREO" scripts/08_exclusivos_local.slurm
```

El registro de cada especie (`logs/fase3_exclusivos_*.out`) indica cuantos
candidatos descarto cada filtro y cuantas proteinas pasan a F5.

**Parte remota, F5-F6 (login, dentro de `tmux`).** Necesita internet, por eso
no va a SLURM:

```bash
bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/fase3_exclusivos_remoto.log
```

Las busquedas se envian a NCBI por su URL API con `curl`, el script consulta
el estado cada minuto y recoge el resultado con `blast_formatter`. No se usa
`blastp -remote`, que desde Khipu se quedaba esperando y en otras redes fallaba
con `Connection stream is in bad state`. Segun la cola de NCBI, cada especie
puede tardar de minutos a horas; mientras espera, el registro escribe
`... sigue en cola` cada 30 min. Si se corta, se vuelve a correr el mismo
comando: las proteinas ya buscadas no se repiten (aunque cambie la lista de
candidatos) y las busquedas en curso se retoman por su numero (RID).

Despues de F5, el paso F5b consulta en NCBI que genomas tienen la proteina de
cada parecido (registro IPG) y la especie de cada genoma por su ANI contra las
cepas tipo. Si alguno es de la especie del bin, el gen se descarta aunque el
genoma este depositado como "sp." o con un sinonimo. Tarda minutos; con la
variable `F5B=0` se omite. Al final imprime `[ok] embudo: ...` por especie: el
paso `pasa_F5b` son los genes especificos de la cepa y `con_bandera` cuantos
de ellos tienen un parecido casi identico (>= 99 %) en otra especie.

Si la cola de NCBI no avanza (puede pasar horas en `WAITING`), F5 se puede
correr en el servicio BLAST del EBI contra UniProtKB, que suele responder en
minutos. Usa la misma tabla de resultados y la misma clasificacion (los taxid
de UniProt son los de NCBI) y escribe en `results/07_exclusivos_ebi/`, sin
tocar la corrida de NCBI, asi que ambas pueden correr a la vez en ventanas
distintas de `tmux`. El EBI exige un correo; se toma de `$CORREO`. No cubre
F6 (contigs huerfanos contra `core_nt`):

```bash
MOTOR=ebi bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/fase3_exclusivos_remoto_ebi.log
```

Si `core_nt` no estuviera disponible, se usa `nt`:

```bash
NT_DB=nt bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/fase3_exclusivos_remoto.log
```

### 9b. Analisis complementarios

ANI de cada bin contra todos los genomas descargados de su especie (un solo
trabajo, minutos):

```bash
sbatch --account=tesis --mail-type=END,FAIL --mail-user="$CORREO" scripts/06c_ani_bins.slurm
```

Categoria funcional (COG) de cada familia del pangenoma, a partir de las
anotaciones de Bakta (login, segundos por especie):

```bash
for d in results/05_panaroo/*/; do python scripts/07_particion.py cog --panaroo-dir "$d" --particion-dir "$d/particion" --bakta-dir results/04_bakta --out "$d/particion/familias_cog.tsv"; done
```

### 10. Exportacion y figuras

En el cluster, reunir las tablas, reportes y arboles en `export_fase3/`:

```bash
bash scripts/bin/export_fase3.sh
```

En la computadora local (Linux, macOS o WSL), desde la raiz del repositorio,
traer los resultados y los lock files de los entornos (cambia `usuario` y el
nombre del cluster por los tuyos):

```bash
rsync -avh usuario@khipu.utec.edu.pe:~/estela/fase3/export_fase3/ results/fase3/
```

```bash
rsync -avh usuario@khipu.utec.edu.pe:~/estela/fase3/scripts/envs/ envs/
```

Y generar las figuras:

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
| `<especie>/exclusivos_verificados.tsv` | Cada candidato con su clase final (`especifico`, `especifico_con_bandera` o `descartado`) y el motivo |
| `<especie>/f5b_ani.tsv` | Especie por ANI de los genomas que tienen la proteina de cada hit de F5 |
| `<especie>/ani_bin_genomas.tsv` | ANI del bin contra todos los genomas descargados de su especie |
| `<especie>/familias_cog.tsv` | Categoria COG de cada familia del pangenoma |

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
| Apertura del pangenoma | Ley de Heaps (misma implementacion que `micropan::heaps`, vectorizada); alpha < 1 indica pangenoma abierto. Diagnosticos de inflacion: genes unicos vs. contigs, genomas atipicos, alpha solo con genomas completos, rango jackknife y fraccion de hipoteticas por categoria | Tettelin et al. 2008; Snipen y Liland 2015 |
| Arbol | IQ-TREE 3, ModelFinder, 1000 replicas UFBoot, sobre el alineamiento del core | |
| Genes especificos de la cepa | Candidatos: genes del bin cuya familia no esta en ninguna referencia y genes del bin que Panaroo elimino (en modo moderate poda de forma recursiva los genes de extremo de contig presentes en < 2 genomas). F1: descarta pseudogenes, < 100 aa y CDS a < 100 pb del borde de contig. F4: descarta si aparece (>= 80 % identidad y cobertura) en cualquier genoma de la especie. F2/F6: contigs sin genes de la especie se contrastan con core_nt y se descartan si su mejor hit es de otro genero. F3: marca genes presentes en otro bin de la misma muestra. F5: BLASTp contra nr; descarta si hay un hit >= 80/80 rotulado con la especie y da el origen probable. F5b: descarta si alguno de los genomas con la proteina de un hit >= 80/80 es de la especie segun su ANI (NCBI). Advertencia si el mejor hit tiene >= 99 % de identidad y >= 90 % de cobertura en otra especie (posible transferencia reciente) | Tonkin-Hill et al. 2020 |

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
| `QOSMaxSubmitJobPerUserLimit` al enviar | El array tiene mas tareas que los envios permitidos; esperar a que terminen trabajos o enviar por tramos |
| El censo devuelve 0 linajes | `phase2_selection.tsv` sin clasificacion de GTDB (ver Instalacion, paso 2) |
| `makeblastdb: command not found` (exit 127) en el paso 9 | El entorno `blast` quedo sin BLAST; recrearlo (Instalacion, paso 4) |
| Un trabajo de Bakta termina por memoria (`OUT_OF_MEMORY`) | Un lote usa hasta ~31 GB; no bajar `--mem` de 32G |
| `-bash: ...: No such file or directory` al enviar con correo | Se copio un marcador con `<...>`; usar `--mail-user="$CORREO"` |

Todos los scripts omiten lo que ya esta hecho, asi que un paso fallido se puede
reenviar sin repetir el trabajo anterior. Para reenviar una sola especie se usa
`--array=<n>` con su numero de fila en el archivo correspondiente. En Bakta se
reenvia el lote completo (`--array=<k> --export=ALL,LOTES=3`): los genomas que
ya estan anotados se saltan.
