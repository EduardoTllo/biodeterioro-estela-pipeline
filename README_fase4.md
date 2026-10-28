# Fase 4: potencial metabolico y riesgo (OE4)

Caracteriza el potencial metabolico de los bins de la Estela con DRAM y lo
cruza con una tabla de marcadores de biodeterioro y supervivencia fijada antes
de ver los datos, y con la particion del pangenoma de la Fase 3.

Estado: codigo completo (instalacion, anotacion y cruce), probado solo con una
anotacion de DRAM sintetica. Falta correr Bakta en los 18 bins restantes y DRAM.

## Pasos disponibles

Los comandos se corren en el cluster, desde `~/estela/fase3` (donde estan los
scripts y las anotaciones de Bakta).

### 1. Instalar DRAM y sus bases de datos (nodo con internet, en tmux)

```bash
bash scripts/09_setup_dram_khipu.sh 2>&1 | tee logs/setup_dram.log
```

Instala DRAM 1.5.0 con el `environment.yaml` oficial y prepara solo las bases
que usa esta fase: KOfam (identificadores KO de KEGG), dbCAN (CAZymes), MEROPS
(peptidasas) y los formularios de resumen de DRAM. Se omiten UniRef90, Pfam,
VOGDB y RefSeq viral, que no hacen falta para los marcadores y son las que mas
memoria requieren al procesarse.

### 2. Anotar con Bakta los bins que no entraron al pangenoma

DRAM se corre sobre las proteinas de Bakta para que cada gen conserve el mismo
identificador que en la Fase 3. Los 3 bins del pangenoma ya estan anotados;
los demas se anotan con el mismo script de la Fase 3:

```bash
{ printf "id\tfasta\ttipo\tslug\tgenoma_original\n"; for f in data/bins/*.fasta; do b=$(basename "$f" .fasta); printf "%s\t%s\tbin\tfase4\t%s\n" "$(echo "$b" | sed 's/[^A-Za-z0-9]/_/g')" "$PWD/$f" "$b"; done; } > results/03_seleccion/bins21_manifest.tsv
```

```bash
sbatch --array=1 --export=ALL,LOTES=1,MANIFEST=results/03_seleccion/bins21_manifest.tsv scripts/05_bakta_all.slurm
```

### 3. Anotacion funcional con DRAM (SLURM)

```bash
sbatch scripts/10_dram.slurm
```

Corre `DRAM.py annotate_genes` sobre las proteinas de Bakta de todos los bins
de `data/bins/` y luego `DRAM.py distill`. Salida en
`~/estela/fase4/results/dram/`: `anotacion/annotations.tsv` (una fila por gen)
y `resumen/` (modulos KEGG, ciclos C/N/S y CAZymes por bin).

Si el cluster exige cuenta, se agrega `--account=<cuenta>` a cada `sbatch`
(ver [README_fase3.md](README_fase3.md)).

### 4. Traer DRAM a la Mac

```bash
mkdir -p results/fase4/dram && rsync -avh eduardo.tello@khipu.utec.edu.pe:~/estela/fase4/results/dram/{anotacion/annotations.tsv,versions.txt,resumen} results/fase4/dram/
```

Si `annotations.tsv` pasa de 50 MB, se comprime con `gzip` antes del commit; el
script lee `annotations.tsv` o `annotations.tsv.gz`.

### 5. Cruce con los marcadores y el pangenoma (en la Mac, R)

```bash
MAMBA_ROOT_PREFIX=~/micromamba ~/micromamba/bin/micromamba run -n r-fase3 Rscript figuras/R/f4_run_all.R
```

`f4_01_marcadores.R` aplica las reglas de `metadata/marcadores_fase4.tsv` (D35)
a los KO de KOfam de cada bin. Una funcion esta presente si alguna alternativa
cumple todos sus grupos. Se detiene si menos del 95 % de los genes de un bin del
pangenoma calza con DRAM (los IDs no coinciden). Salidas en `figuras/tablas/`:

| Tabla | Contenido |
|---|---|
| `tab_f4_funciones_bin.tsv` | bin x funcion: presente y KO hallados |
| `tab_f4_puntajes.tsv` | funciones presentes / evaluables por bin, con completitud (provisional, D33) |
| `tab_f4_fisher.tsv` | Fisher por especie (D34, A.1): 2x2, OR, IC 95 %, p y BH por definicion de accesorio |
| `tab_f4_genes_marcador.tsv` | genes marcador de los bins del pangenoma y su particion |
| `tab_f4_especificos_dram.tsv` | genes especificos de la cepa con su anotacion DRAM (D34, A.2) |
| `tab_f4_cazy_bin.tsv` | genes por familia CAZy y bin (contexto) |

En Fisher, marcador = gen con un KO de una alternativa cumplida de una funcion
presente en el bin. Solo entran los genes del bin que estan en el pangenoma
(sin exclusivos ni `no_en_pangenoma`). Definiciones de accesorio: `principal`
(shell + cloud, D27), `solo_cloud` (cloud frente a core; el shell queda fuera) y
`core_90` (particion con core >= 90 %). Si un bin no tiene genes marcador de una
categoria, la prueba queda en `NA`.

`f4_02_matriz.R` dibuja `figuras/figs/fig_f4_matriz_funciones` (bin x funcion,
con linaje, completitud y puntaje).
