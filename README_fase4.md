# Fase 4: potencial metabolico y riesgo (OE4)

Caracteriza el potencial metabolico de los bins de la Estela con DRAM y lo
cruza con una tabla de marcadores de biodeterioro y supervivencia fijada antes
de ver los datos, y con la particion del pangenoma de la Fase 3.

Estado: en construccion. Por ahora estan listos la instalacion de DRAM y la
anotacion de los bins; el cruce con los marcadores y el pangenoma se agregara
cuando se cierre el diseno.

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
