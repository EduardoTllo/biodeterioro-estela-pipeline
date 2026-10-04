# Fase 3: pangenoma comparativo (OE3)

## Que hace esta fase

Para cada especie que tiene suficientes genomas publicos, arma un pangenoma con
esos genomas (las "referencias") mas el bin de la Estela, y responde: **que
genes tiene el bin de la Estela que no tienen las demas cepas conocidas de su
especie**. Esos genes, ya verificados, pasan a la Fase 4.

En resumen:

1. Cuenta cuantos genomas publicos de calidad hay para cada especie de la Fase 2.
2. Descarga esos genomas de NCBI.
3. Quita los genomas casi identicos y elige hasta 50 por especie.
4. Anota todos los genomas (referencias y bins) con Bakta.
5. Construye un pangenoma por especie con Panaroo y un arbol con IQ-TREE.
6. Revisa los genes que solo aparecen en el bin de la Estela para descartar errores.

Las decisiones de diseno (umbrales, criterios) estan al final, en
[Como funciona](#como-funciona).

---

## Antes de empezar

### Que necesitas

- Acceso SSH a Khipu (usuario `eduardo.tello` en `khipu.utec.edu.pe`).
- Acceso al repositorio de GitHub `EduardoTllo/biodeterioro-estela-pipeline`. Es
  **privado**: para clonarlo en Khipu necesitas un *personal access token* de
  GitHub (GitHub > Settings > Developer settings > Personal access tokens). Git
  lo pide como "Password". No es la contrasena de GitHub.
- Las Fases 1 y 2 ya corridas **en Khipu**. Esta fase usa sus resultados.

Esta fase no necesita subir nada desde una laptop: todo lo que usa ya esta en
Khipu o se descarga desde ahi.

### Comprobar que las Fases 1 y 2 estan en Khipu

Conectate a Khipu:

```bash
ssh eduardo.tello@khipu.utec.edu.pe
```

Corre estos cuatro comandos. Cada uno debe mostrar archivos, no un error:

```bash
ls ~/estela/fase1/bins/*.fasta | wc -l
```

Debe dar **97** (los bins crudos de la Fase 1).

```bash
ls ~/estela/fase1/results/phase1_selected_genomes/*.fasta | wc -l
```

Debe dar **21** (los bins que pasaron el control de calidad).

```bash
ls ~/estela/fase1/results/02_checkm2/quality_report.tsv
```

```bash
ls ~/estela/fase2/results/phase2_selection.tsv
```

Si alguno falla, esas carpetas estan en otra ruta: encuentrala con
`find ~ -name phase2_selection.tsv` (o el archivo que falte) y usa esa ruta en
los comandos de abajo.

### Cuatro cosas que conviene saber

**Login y nodos de computo.** Al entrar por SSH quedas en el *nodo de login*. Es
el unico con internet, asi que ahi se instala todo y se descargan datos. Los
calculos pesados se mandan a los *nodos de computo* con `sbatch`; esos no tienen
internet.

**`sbatch`, `squeue`, `sacct`.** `sbatch archivo.slurm` envia un trabajo a la
cola y te devuelve un numero (JobID). `squeue --me` muestra tus trabajos en cola
o corriendo; cuando un trabajo desaparece de ahi, termino. `sacct` dice si
termino bien (`COMPLETED`) o mal (`FAILED`, `TIMEOUT`). Lo que imprime cada
trabajo queda en la carpeta `logs/`.

**`tmux`.** Si ejecutas algo largo directamente en el login y se corta la
conexion, el proceso muere. `tmux` crea una sesion que sigue viva aunque te
desconectes:

| Que quieres | Comando |
|---|---|
| Crear una sesion llamada `fase3` | `tmux new -s fase3` |
| Salir dejandola corriendo | `Ctrl+b`, soltar, y luego `d` |
| Volver a entrar | `tmux attach -t fase3` |
| Ver las sesiones que existen | `tmux ls` |

Los trabajos de `sbatch` no necesitan `tmux`: corren solos en los nodos de
computo.

**`2>&1 | tee archivo.log`.** Agregado al final de un comando, muestra la salida
en pantalla y ademas la guarda en `archivo.log`, incluidos los errores. Sirve
para revisar despues que paso.

### Activar conda

Varios pasos usan Python o programas instalados con conda. En cada sesion nueva
de Khipu (o cada vez que abras `tmux`), antes de esos pasos:

```bash
module load miniconda/3.0 && eval "$(conda shell.bash hook)"
```

---

## Paso a paso

Todos los comandos se corren **en Khipu**, salvo el paso 14. Los `sbatch` se
envian **siempre desde `~/estela/fase3`**, porque los trabajos guardan su
registro en `logs/` relativo a esa carpeta.

### Paso 1. Crear las carpetas

```bash
mkdir -p ~/estela/fase3/data/bins ~/estela/fase3/data/bins_libreria ~/estela/fase3/logs
```

### Paso 2. Descargar el codigo

Clona el repositorio en `~/estela/fase3/scripts`. La carpeta `scripts` no debe
existir todavia (git la crea):

```bash
git clone https://github.com/EduardoTllo/biodeterioro-estela-pipeline.git ~/estela/fase3/scripts
```

Git pide `Username` (tu usuario de GitHub) y `Password` (el token, no la
contrasena). Si dice que la carpeta ya existe y no esta vacia, revisa que hay
con `ls -A ~/estela/fase3/scripts` antes de borrarla.

Comprueba que tienes la ultima version:

```bash
cd ~/estela/fase3/scripts && git log --oneline -1
```

Si mas adelante se corrige algo en GitHub, actualiza con:

```bash
cd ~/estela/fase3/scripts && git pull
```

### Paso 3. Copiar los bins de la Fase 1

Los 21 bins seleccionados (los que entran al pangenoma):

```bash
cp ~/estela/fase1/results/phase1_selected_genomes/*.fasta ~/estela/fase3/data/bins/
```

Los 97 bins crudos (se usan en un control de contaminacion, el filtro F3):

```bash
cp ~/estela/fase1/bins/*.fasta ~/estela/fase3/data/bins_libreria/
```

Comprueba los conteos y los nombres:

```bash
ls ~/estela/fase3/data/bins/*.fasta | wc -l; ls ~/estela/fase3/data/bins_libreria/*.fasta | wc -l; ls ~/estela/fase3/data/bins_libreria | head -3
```

Debe dar 21 y 97, y los nombres deben ser como `bin-1-49.fasta`. Si los nombres
son distintos (por ejemplo `bin_0001.fasta`), no sigas: los scripts necesitan el
nombre `bin-<N>-<muestra>`.

### Paso 4. Instalar el entorno (en `tmux`, varias horas)

Instala los programas (Bakta, Panaroo, IQ-TREE, BLAST, NCBI Datasets) y descarga
la base de datos de Bakta (30 GB comprimida, 84 GB descomprimida), la taxonomia
de NCBI y el metadata de GTDB.

```bash
tmux new -s fase3
```

```bash
bash ~/estela/fase3/scripts/03_setup_fase3_khipu.sh 2>&1 | tee ~/estela/fase3/logs/setup.log
```

Sal con `Ctrl+b`, `d` y vuelve mas tarde con `tmux attach -t fase3`. Termina
cuando aparece `SETUP DE LA FASE 3 COMPLETO`.

**Comprobar:**

```bash
grep -n -i "error\|aviso" ~/estela/fase3/logs/setup.log
```

No deberia haber lineas `ERROR`. Los `AVISO` sobre entornos `tiara`, `checkm2`
o `gtdbtk` solo significan que no se pudo exportar su lock file.

```bash
cat $(cat ~/estela/fase3/bakta_db_path.txt)/version.json
```

Debe contener `"major": 6`.

Si el script se corto o fallo, vuelve a correrlo: salta lo que ya esta hecho.

### Paso 5. Censo de genomas publicos (login)

Cuenta, para cada especie de la Fase 2, cuantos genomas publicos de calidad hay.

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

```bash
python scripts/03_censo_genomas.py censo --selection ~/estela/fase2/results/phase2_selection.tsv --gtdb-metadata ~/gtdbtk_data/release220/bac120_metadata_r220.tsv --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py ncbi --outdir results/00_censo
```

```bash
python scripts/03_censo_genomas.py tabla --outdir results/00_censo --umbral 15
```

**Revisar:**

```bash
less results/00_censo/phase3_censo_report.md
```

(`q` para salir de `less`.) Debe haber 12 linajes en la tabla. Las especies con
**Viable = Si** pasan al paso siguiente. Si **ninguna** es viable, no sigas y
consultalo con el asesor.

**Pausa: mostrar este reporte al asesor antes de seguir.**

### Paso 6. Descargar las referencias (en `tmux`, 1 a 3 horas)

Descarga de NCBI todos los genomas que pasaron el censo y comprueba que cada uno
llego completo.

```bash
tmux attach -t fase3
```

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && bash scripts/04_fetch_refs.sh 2>&1 | tee logs/descarga.log
```

**Revisar:**

```bash
column -t -s$'\t' results/01_descarga/descarga_resumen.tsv
```

La columna `ok` dice cuantos genomas llegaron bien por especie. Si hay valores
en `tamano_discrepante` o `version_distinta`, el detalle esta en
`results/01_descarga/descarga_estado.tsv`. Esos genomas quedan fuera; no es un
error del script.

### Paso 7. Clasificar las referencias (login)

Asigna a cada genoma su prioridad, su habitat (segun la fuente de aislamiento) y
su continente.

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

```bash
python scripts/04_clasificar_refs.py clasificar --censo-dir results/00_censo --descarga-dir results/01_descarga --refs-dir data/refs --keywords scripts/metadata/habitat_keywords.tsv --paises scripts/metadata/paises_continentes.tsv --outdir results/02_drep
```

**Revisar a mano:**

```bash
column -t -s$'\t' results/02_drep/fuentes_unicas.tsv | less -S
```

Cada fila es una fuente de aislamiento con el habitat asignado. Si una esta mal
clasificada (por ejemplo, una roca marcada como "suelo"), corrige
`metadata/habitat_keywords.tsv` en el repo, actualiza con `git pull` y repite
este paso.

```bash
cat results/02_drep/paises_sin_mapear.tsv
```

Los paises que aparezcan aqui se agregan a `metadata/paises_continentes.tsv`.

### Paso 8. Quitar genomas casi identicos (SLURM, menos de 6 h)

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/01_descarga/especies_para_drep.tsv | wc -l) && echo "$N especies" && sbatch --array=1-$N scripts/04_drep_refs.slurm
```

Espera a que `squeue --me` ya no muestre trabajos `f3_drep`. Luego:

```bash
sacct -u $USER --name=f3_drep -S today --format=JobID,State,Elapsed
```

Todas las filas deben decir `COMPLETED`.

### Paso 9. Elegir las especies y las referencias (login)

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets
```

```bash
python scripts/04_clasificar_refs.py seleccionar --descarga-dir results/01_descarga --drep-dir results/02_drep --refs-dir data/refs --bins-dir data/bins --checkm2 ~/estela/fase1/results/02_checkm2/quality_report.tsv --outdir results/03_seleccion
```

**Revisar:**

```bash
less results/03_seleccion/phase3_seleccion_report.md
```

Muestra las especies elegidas (hasta 3) y de que habitats y continentes vienen
sus referencias.

Si el script se detiene con "dRep no leyo extraW.tsv": borra
`results/02_drep/<especie>/data_tables`, vuelve a enviar solo esa especie con
`sbatch --array=<numero de fila> scripts/04_drep_refs.slurm` y repite este paso.

**Pausa: mostrar este reporte al asesor antes de seguir.**

### Paso 10. Prueba piloto (obligatoria)

Antes de anotar unos 150 genomas, se prueba todo con 3 referencias y el bin de
la primera especie.

**10a. Anotar con Bakta** (SLURM, menos de 1 h):

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/03_seleccion/piloto_manifest.tsv | wc -l) && sbatch --array=1-$N --export=ALL,MANIFEST=results/03_seleccion/piloto_manifest.tsv scripts/05_bakta_all.slurm
```

**10b. Panaroo**, cuando `squeue --me` ya no muestre `f3_bakta`:

```bash
cd ~/estela/fase3 && sbatch --export=ALL,PILOTO=1 scripts/06_panaroo.slurm
```

Cuando termine, comprueba estas tres cosas:

```bash
grep -l "Error reading prokka input" logs/fase3_panaroo_*
```

No debe mostrar nada.

```bash
head -n1 results/05_panaroo/_piloto_*/gene_presence_absence.csv
```

Debe listar los nombres de los 4 o 5 genomas del piloto.

```bash
cat results/05_panaroo/_piloto_*/particion/particion_report.md
```

Debe existir y mostrar una tabla.

**10c. Probar el BLAST por internet** (login):

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate blast && ID=$(tail -n +2 results/03_seleccion/piloto_manifest.tsv | head -n1 | cut -f1) && head -n2 results/04_bakta/$ID/$ID.faa > /tmp/prueba.faa && blastp -remote -db nr -query /tmp/prueba.faa -max_target_seqs 3 -outfmt "6 qseqid sseqid pident"
```

```bash
head -n2 results/04_bakta/$ID/$ID.fna | cut -c1-2000 > /tmp/prueba.fna && blastn -remote -db core_nt -task megablast -query /tmp/prueba.fna -max_target_seqs 3 -outfmt "6 qseqid sseqid pident"
```

Cada uno puede tardar unos minutos y debe devolver lineas con resultados. Si el
segundo falla, prueba cambiando `core_nt` por `nt`. Si `nt` funciona, en el
paso 13 usaras `NT_DB=nt`.

Si algo del piloto falla, **no sigas**: guarda el archivo `.err` de `logs/` y
revisalo antes de continuar.

### Paso 11. Anotar todos los genomas (SLURM, 15 a 25 h)

Los genomas del piloto no se repiten.

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/03_seleccion/bakta_manifest.tsv | wc -l) && echo "$N genomas" && sbatch --array=1-$N%4 scripts/05_bakta_all.slurm
```

(`%4` limita a 4 trabajos a la vez, para no pasar los 32 nucleos de la cuenta.)

Cuando termine, el resumen:

```bash
cd ~/estela/fase3 && module load miniconda/3.0 && eval "$(conda shell.bash hook)" && conda activate ncbi-datasets && python scripts/07_particion.py bakta --manifest results/03_seleccion/bakta_manifest.tsv --bakta-dir results/04_bakta --out results/04_bakta_resumen.tsv
```

Si dice que algun genoma esta "sin anotacion", busca su error en
`logs/fase3_bakta_<JobID>_<n>.err` y reenvia solo ese con
`sbatch --array=<n> scripts/05_bakta_all.slurm`.

### Paso 12. Pangenoma y arbol (SLURM, hasta 24 h cada uno)

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l) && sbatch --array=1-$N%1 scripts/06_panaroo.slurm
```

Cuando terminen todos los `f3_panaroo`:

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l) && sbatch --array=1-$N%1 scripts/06b_iqtree.slurm
```

**Revisar:**

```bash
cat ~/estela/fase3/results/05_panaroo/*/particion/particion_report.md
```

La seccion "Checkpoint C6" debe decir "Sin alertas". Si hay alertas, revisalas
antes de interpretar resultados.

Si `sacct` muestra `TIMEOUT` para IQ-TREE en una especie, reenvia solo esa
especie con la version rapida:

```bash
cd ~/estela/fase3 && sbatch --array=<numero de fila> --export=ALL,MODO=snps scripts/06b_iqtree.slurm
```

### Paso 13. Revisar los genes exclusivos

**13a. Parte local** (SLURM):

```bash
cd ~/estela/fase3 && N=$(tail -n +2 results/03_seleccion/especies_seleccionadas.tsv | wc -l) && sbatch --array=1-$N scripts/08_exclusivos_local.slurm
```

**13b. Parte por internet** (login, en `tmux`, puede tardar horas):

```bash
tmux attach -t fase3
```

```bash
cd ~/estela/fase3 && bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/exclusivos_remoto.log
```

Si en el paso 10c solo funciono `nt`, usa en su lugar
`NT_DB=nt bash scripts/08_exclusivos_remoto.sh 2>&1 | tee logs/exclusivos_remoto.log`.
Si se corta, vuelve a correrlo: continua donde quedo.

**Revisar:**

```bash
cat ~/estela/fase3/results/07_exclusivos/*/exclusivos_report.md
```

### Paso 14. Llevar los resultados al repositorio

**En Khipu**, junta en una carpeta los archivos que se guardan en el repo:

```bash
bash ~/estela/fase3/scripts/bin/export_fase3.sh
```

**En la computadora donde tienes el repositorio** (en Windows, desde WSL). Primero
entra a la carpeta del repo. La ruta depende de donde lo clonaste; por ejemplo:

```bash
cd /mnt/c/Tesis-EstelaRaimondi/fase1_khipu
```

Trae los resultados y los lock files de los entornos:

```bash
rsync -avh eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/export_fase3/ results/fase3/
```

```bash
rsync -avh eduardo.tello@khipu.utec.edu.pe:~/estela/fase3/scripts/envs/ envs/
```

Genera las figuras (necesita R con `micropan` y `phangorn`; ver
[figuras/README.md](figuras/README.md)):

```bash
Rscript figuras/R/f3_run_all.R
```

Luego commitea `results/fase3/` y `envs/`.

---

## Si algo falla

| Sintoma | Que hacer |
|---|---|
| `No such file or directory` al correr un script | Revisa que estas en la carpeta correcta y que el paso 2 se completo (`ls ~/estela/fase3/scripts`) |
| `conda: command not found` o `python: command not found` | Falta activar conda (seccion [Activar conda](#activar-conda)) |
| Un trabajo de `sbatch` no aparece en `squeue --me` y no hizo nada | Mira su error en `logs/`. Si `logs/` no existia al enviarlo, SLURM no pudo escribir: crea la carpeta (paso 1) y reenvia |
| `sacct` dice `FAILED` | Lee el `.err` de ese trabajo en `logs/`. Corrige y reenvia solo esa tarea con `--array=<n>` |
| `sacct` dice `TIMEOUT` | El trabajo paso de su tiempo maximo. Para IQ-TREE usa `MODO=snps` (paso 12) |
| Se corto la conexion durante un paso largo | Vuelve a entrar y usa `tmux attach -t fase3`. Si no usaste `tmux`, vuelve a correr el paso: todos los scripts saltan lo que ya esta hecho |

Para ver los trabajos de hoy y como terminaron:

```bash
sacct -u $USER -S today --format=JobID,JobName,State,Elapsed,MaxRSS
```

---

## Como funciona

### Criterio para elegir especies

Una especie entra al pangenoma si tiene **al menos 15 genomas publicos de buena
calidad y no redundantes**. Se eligen las 3 con mas genomas; si hay menos de 3,
se usan las que haya.

- **Por que 15**: Gautreau et al. 2020 (*PLoS Comput Biol* 16(3):e1007732)
  recomiendan al menos 15 genomas para dividir un pangenoma de forma confiable.
- **Por que no redundantes**: los genomas casi identicos inflan el core y pueden
  hacer parecer cerrado un pangenoma abierto (Guerra 2026, *Bioinform Adv*
  6(1):vbag069). Por eso se cuenta despues de agrupar al 99 % de identidad (ANI).
- Este criterio reemplaza al de prevalencia espacial de la Fase 2, porque 17 de
  los 19 linajes tienen un solo genoma propio y el pangenoma depende de genomas
  publicos. La prevalencia se sigue reportando.

### Que genomas se usan como referencia

| Tema | Regla |
|---|---|
| Fuente | La lista sale del metadata de GTDB R220 (el mismo release de la Fase 2); los archivos se descargan de NCBI. No se agregan genomas posteriores a R220 |
| Calidad | Completitud >= 95 %, contaminacion <= 5 %, <= 300 contigs y solo aislados (sin MAGs) |
| Especies censadas | Solo las 12 con nombre de especie real. Quedan fuera L2, L4, L5 y L10 (sin especie) y L6, L7 y L19 (nombre provisional de GTDB) |
| Descarga | Version exacta de cada genoma; si no baja o su tamano no coincide con GTDB (+/- 1 %), queda fuera |
| Maximo 50 por especie | Entran primero el genoma representativo de GTDB y la cepa tipo; luego los de Latinoamerica y los de roca o ambientes aridos; el resto se reparte entre habitats y continentes |

### Anotacion y pangenoma

- **Bakta 1.12.1** con la base de datos v6.0 para todos los genomas, con las
  mismas opciones. No se usan las anotaciones de NCBI, para que las diferencias
  entre genomas no vengan del programa de anotacion.
- **Panaroo 1.8.0** en modo `moderate`. El modo `strict` borra los genes que
  aparecen en un solo genoma, que es justamente lo que se busca.
- Panaroo no lee bien los archivos de Bakta, asi que antes se convierten con un
  script oficial de Panaroo (lo descarga el paso 4).

### Como se clasifican los genes

La frecuencia de cada familia de genes se calcula **solo entre las
referencias**, para que lo que le falta al bin por estar incompleto no cambie el
resultado.

| Categoria | Presente en |
|---|---|
| core | 95 % o mas de las referencias |
| shell | entre 15 % y 95 % |
| cloud | menos de 15 % |
| exclusivo (candidato) | ninguna referencia, pero si el bin |

"Accesorio" = shell + cloud. Tambien se calcula con el core al 90 %, como
comprobacion.

### Como se revisan los genes exclusivos

Un gen "exclusivo" puede ser un error. Cada candidato pasa por seis filtros:

| Filtro | Descarta o marca si... |
|---|---|
| F1 | es un pseudogen, mide menos de 100 aminoacidos o esta a menos de 100 pb del borde de un contig (puede ser un gen cortado) |
| F4 | aparece en algun otro genoma de la especie, aunque no este entre los 50 elegidos |
| F2 | su contig no tiene ningun gen de la especie (contig "huerfano"): pasa a F6 |
| F3 | aparece casi igual en otro bin de la misma muestra (se marca: puede ser contaminacion) |
| F6 | su contig huerfano se parece mas a otro genero en NCBI (se descarta: contaminacion) |
| F5 | se busca en la base `nr` de NCBI para saber de donde podria venir |

Resultado final de cada gen: `verificado`, `con_bandera` o `descartado` (con el
motivo).

### Advertencia para interpretar

Los bins tienen entre 70 y 100 % de completitud y las referencias 95 % o mas. Por
eso **que un gen este en el bin es confiable, pero que falte no lo es**. Cada
resultado se reporta junto a la completitud del bin. Lo mismo vale para la Fase 4.

---

## Archivos

### Scripts

| Archivo | Que hace | Donde corre |
|---|---|---|
| `03_setup_fase3_khipu.sh` | Instala programas y bases de datos | login |
| `03_censo_genomas.py` | Censo de genomas publicos | login |
| `04_fetch_refs.sh` | Descarga y verifica las referencias | login |
| `04_clasificar_refs.py` | Clasifica referencias y elige especies y genomas | login |
| `04_drep_refs.slurm` | Agrupa genomas casi identicos (dRep) | SLURM |
| `05_bakta_all.slurm` | Anota con Bakta | SLURM |
| `06_panaroo.slurm` | Pangenoma con Panaroo y clasificacion de genes | SLURM |
| `06b_iqtree.slurm` | Arbol del core con IQ-TREE | SLURM |
| `07_particion.py` | Resumen de Bakta y clasificacion de genes | lo llaman otros pasos |
| `08_exclusivos.py` | Filtros de genes exclusivos | lo llaman otros pasos |
| `08_exclusivos_local.slurm` | Filtros F1 a F4 | SLURM |
| `08_exclusivos_remoto.sh` | Filtros F5 y F6 (por internet) | login |
| `bin/export_fase3.sh` | Junta los resultados para el repo | login |
| `figuras/R/f3_*.R` | Figuras y tablas | computadora local |
| `metadata/habitat_keywords.tsv` | Palabras para clasificar el habitat | se edita a mano si hace falta |
| `metadata/paises_continentes.tsv` | Pais -> continente | se edita a mano si hace falta |

### Carpetas en Khipu

```
~/estela/fase3/
|-- scripts/                 el repositorio (paso 2)
|-- data/bins/               21 bins seleccionados (paso 3)
|-- data/bins_libreria/      97 bins crudos (paso 3)
|-- data/refs/<especie>/     genomas descargados (paso 6)
|-- results/00_censo/        paso 5
|-- results/01_descarga/     paso 6
|-- results/02_drep/         pasos 7 y 8
|-- results/03_seleccion/    paso 9
|-- results/04_bakta/        pasos 10 y 11
|-- results/05_panaroo/      pasos 10 y 12
|-- results/06_iqtree/       paso 12
|-- results/07_exclusivos/   paso 13
|-- export_fase3/            paso 14
`-- logs/                    registros de cada trabajo
~/dbs/bakta/db/              base de datos de Bakta (paso 4)
~/dbs/taxdump/               taxonomia de NCBI (paso 4)
~/gtdbtk_data/release220/    metadata de GTDB (paso 4)
```

### Resultados que se guardan en el repo (`results/fase3/`)

| Archivo | Contenido |
|---|---|
| `phase3_censo_report.md`, `phase3_viabilidad.tsv` | Censo: que especies tienen suficientes genomas |
| `phase3_descarga_resumen.tsv` | Cuantos genomas se descargaron bien |
| `phase3_seleccion_report.md`, `phase3_referencias_finales.tsv` | Especies elegidas y sus referencias |
| `<especie>/particion_report.md`, `genes_bin.tsv` | Clasificacion de cada gen del bin |
| `<especie>/core.treefile` | Arbol del core |
| `<especie>/exclusivos_verificados.tsv`, `exclusivos_report.md` | Genes exclusivos revisados |

## Estado

Codigo completo, probado de punta a punta con datos sinteticos. Pendiente de
ejecucion con datos reales en Khipu.
