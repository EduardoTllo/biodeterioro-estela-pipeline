# Fase 3 - Analisis critico de los resultados

Documento de trabajo (08-10-2026). Revisa cada resultado de la Fase 3, lo pone a
prueba con analisis adicionales y senala que se sostiene, que hay que corregir
y que queda abierto. Complementa la guia de lectura
(`phase3_guia_resultados.md`), que explica las figuras 1-7.

**Estado (08-10, tarde).** Las correcciones de la seccion 6 se decidieron y se
implementaron (commit `2dfbe1f`): A, B, D, E y F si; C no. Ya corrieron la
parte local de la verificacion (v2, con los genes que Panaroo habia eliminado)
y el ANI de cada bin. Falta la parte remota (F5 en `nr` de las proteinas
nuevas, F5b y F6): **los numeros finales de genes especificos de la cepa estan
pendientes** y las cifras de la v1 (0 / 0 / 27) se mantienen solo como
referencia. Desde esta version, "exclusivo verificado" pasa a llamarse **gen
especifico de la cepa**.

Especies: *Bacillus altitudinis* (L16, bin-5-63, 72,5 % completo),
*Peribacillus frigoritolerans* (L12, bin-1-54, 92,3 %) y
*Acinetobacter schindleri* (L8, bin-4-52, 99,9 %).

---

## 0. Resumen

**Lo que se sostiene**

1. Los tres pangenomas son abiertos y la apertura no es un artefacto: se
   mantiene con el mismo numero de genomas (23) y con solo genomas completos
   (salvo una atenuacion en *P. frigoritolerans*).
2. Calcular la particion solo con referencias funciona: el core recuperado en
   cada bin coincide con su completitud.
3. Cada bin tiene un pariente cercano identificable en el arbol del core, con
   soporte maximo.
4. El contenido accesorio sigue a la filogenia (Mantel r = 0,25-0,73,
   p = 0,001) y no al habitat de aislamiento (PERMANOVA no significativa).

**Lo que hay que corregir**

5. **Los "exclusivos verificados" eran 2, 0 y 30; con una comprobacion por ANI
   pasan a 0, 0 y 27.** Cinco genes que F5 dio por exclusivos estan, identicos,
   en genomas de la misma especie publicados en 2025 con otro nombre en NCBI
   (*Bacillus* sp. 22483, ANI 98,3 % con *B. altitudinis*; *Acinetobacter*
   sp. 22323, ANI 97,75 % con *A. schindleri*). El filtro F5 decidia "misma
   especie" por el nombre de NCBI, no por ANI.
6. **F6 no vio contigs huerfanos porque Panaroo ya los habia eliminado.**
   Panaroo quito 140 genes de los bins, entre ellos 4 contigs completos
   (40 kb) del bin de *A. schindleri*, 3 con GC de 0,35 frente a 0,42 del
   genoma. Esos genes nunca llegaron a ser candidatos ni pasaron por F6. La
   causa esta confirmada en el codigo de Panaroo (seccion 3.7). En la v2 entran
   como candidatos: 46 / 41 / 171 candidatos (antes 14 / 34 / 70), 14 / 1 / 97
   proteinas a F5 (antes 2 / 1 / 36) y 6 / 0 / 4 contigs a F6 (antes 0).
7. En la guia se dijo que la rama de bin-5-63 era larga por su incompletitud.
   No es asi: mide 0,0083, casi lo mismo que la de su pariente (0,0077).

**Lo que cambia la interpretacion**

8. El numero de exclusivos no mide "adaptacion": cualquier genoma de estas
   especies tiene genes propios (mediana 18, 120 y 104 por referencia), y los
   bins estan dentro de lo tipico (percentiles 34, 2 y 26).
9. bin-1-54 pertenece a un clon cosmopolita de *P. frigoritolerans* (8 genomas
   a >= 99 % ANI de 4 continentes y 6 habitats). Confirmado con fastANI:
   99,98 % con GCA_900000145.1, es practicamente la misma cepa. Por eso no
   tiene ningun gen propio. bin-5-63 tambien pertenece a un clon ya
   secuenciado (99,06 %); bin-4-52 no (98,39 % con su pariente mas cercano).
10. Los exclusivos de *A. schindleri* son casi todos ADN movil (profagos,
    islas con integrasas y transposasas); 19 de 27 son identicos (>= 99 %) a
    proteinas de otras especies de *Acinetobacter*.

---

## 1. Conceptos

| Termino | Significado en este trabajo |
|---|---|
| Familia de genes | Grupo de genes que Panaroo considera el mismo gen en distintos genomas (por similitud de secuencia y por los genes vecinos). Es la unidad que se cuenta en el pangenoma |
| Pangenoma | Todas las familias presentes en al menos un genoma de la especie |
| Core / shell / cloud | Familias en >= 95 %, en 15-95 % y en < 15 % de las referencias. Shell + cloud = genoma accesorio |
| Singleton | Familia presente en una sola referencia: lo mas propio de cada genoma |
| Pangenoma abierto / cerrado | Abierto: cada genoma nuevo sigue aportando familias nuevas. Cerrado: tras cierto numero de genomas ya no aparecen familias nuevas |
| Ley de Heaps, alpha | Ajuste de cuantas familias nuevas aporta el genoma n: nuevas ~ k * n^(-alpha). alpha < 1 = abierto. Es un resumen estadistico de la muestra, no una propiedad fija de la especie (Guerra 2026) |
| Referencia | Genoma publico de la misma especie (aislado, completitud >= 95 %, contaminacion <= 5 %) usado para construir el pangenoma |
| ANI | Identidad nucleotidica promedio entre dos genomas. >= 95 % = misma especie; >= 99 % = practicamente la misma cepa o clon |
| Desreplicacion (dRep) | Agrupar los genomas con >= 99 % de ANI y quedarse con uno por grupo, para que una cepa muy secuenciada no pese de mas |
| MAG / bin | Genoma reconstruido a partir de un metagenoma. Puede estar incompleto: que un gen falte en el bin no prueba que la bacteria no lo tenga |
| Distancia patristica | Suma de las ramas del arbol entre dos genomas (sustituciones por sitio del core). 0,0017 ~ 99,8 % de identidad en el core |
| UFBoot | Soporte de cada agrupamiento del arbol (0-100). >= 95 se considera confiable |
| Exclusivo candidato | Gen del bin cuya familia no esta en ninguna referencia |
| Gen especifico de la cepa (antes "exclusivo verificado") | Candidato que sobrevive a los filtros F1-F6 y F5b: bien formado y sin una copia parecida (>= 80 % de identidad en >= 80 % del largo) en ningun otro genoma conocido de la especie |
| Candidato eliminado por Panaroo | Gen del bin que Panaroo saco de la grafica por estar en el extremo de un contig y presente en un solo genoma. Desde la v2 tambien es candidato |
| F5b | Comprobacion por ANI de la especie de todos los genomas que tienen la proteina de cada parecido de F5 (IPG + ANI de NCBI) |
| Advertencia (bandera) | Gen especifico cuyo mejor parecido tiene >= 99 % de identidad y >= 90 % de cobertura en otra especie: posible transferencia reciente. No se descarta |
| Transferencia horizontal (HGT) | Adquisicion de ADN de otro organismo, no heredado del ancestro. Suele viajar en elementos moviles |
| Elemento genetico movil | ADN que se mueve entre genomas: profagos (virus integrados), transposones, secuencias de insercion (IS), islas con integrasa |
| IPG (Identical Protein Groups) | Registro de NCBI que lista todos los genomas que contienen una proteina identica |
| PCoA (Jaccard) | Mapa en 2D donde genomas con contenido accesorio parecido quedan cerca. Jaccard = proporcion de familias no compartidas |
| Mantel | Prueba si dos matrices de distancia se correlacionan (aqui: accesorio frente a filogenia) |
| PERMANOVA | Prueba si un factor (aqui: habitat) explica parte de las diferencias de contenido accesorio. R2 = fraccion explicada |

---

## 2. Que significa "especifico de la cepa" (antes "verificado") y por que 0, 2 o 30

Un gen especifico de la cepa responde a una sola pregunta: **este gen del bin de la
Estela, existe en algun otro genoma conocido de su especie?** Si la respuesta
es no (con los umbrales de F4 y F5), el gen es especifico de la cepa de la
Estela frente a todo lo secuenciado hasta hoy.

No significa que el gen sea nuevo para la ciencia, ni que sea una adaptacion a
la piedra, ni que sea raro en la naturaleza: puede estar identico en otra
especie del mismo genero.

| | v1: candidatos | v1: tras F1-F6 | v1 corregido por ANI | v2: candidatos | v2: a F5 | Lectura |
|---|---:|---:|---:|---:|---:|---|
| *B. altitudinis* | 14 | 2 | **0** | 46 | 14 | Los 2 de la v1 estan en una cepa de rizosfera china publicada en 2025 |
| *P. frigoritolerans* | 34 | 0 | **0** | 41 | 1 | El bin es practicamente la misma cepa que otro genoma (ANI 99,98 %) |
| *A. schindleri* | 70 | 30 | **27** | 171 | 97 | 3 de la v1 estan en otra cepa china de 2025; el resto es ADN movil compartido con otras *Acinetobacter* |

La v2 agrega los genes que Panaroo habia eliminado (seccion 3.7). El resultado
final de la v2 esta pendiente de F5, F5b y F6.

Por que tan distintos:

- **Cuantos genomas hay de la especie.** F4 busco en 208, 69 y 29 genomas.
  Con 29 genomas de *A. schindleri* es mucho mas probable que un gen no
  aparezca en ninguno.
- **Que tan cerca esta el pariente secuenciado mas proximo.** bin-1-54 tiene
  un pariente a 0,0017 (casi identico); los otros dos, a 0,013-0,016.
- **Completitud.** A bin-5-63 le falta ~25 % del genoma: menos genes que
  comparar.
- **Biologia del genero.** *Acinetobacter* intercambia mucho ADN movil
  (profagos, plasmidos, IS); en la figura 12 se ve.

---

## 3. Cuestionamientos

Cada punto: la pregunta, la evidencia y el veredicto.

### 3.1 Se pueden comparar tres pangenomas con 50, 42 y 23 referencias?

No directamente: el pangenoma crece con el numero de genomas. Se submuestrearon
23 genomas de cada especie (200 repeticiones):

| | Pangenoma con 23 | Core (100 %) con 23 | alpha con 23 (rango) |
|---|---:|---:|---|
| *B. altitudinis* | 5 994 (5 628-6 310) | 3 126 | 0,67 (0,59-0,76) |
| *P. frigoritolerans* | 12 905 (12 117-13 573) | 3 580 | 0,58 (0,51-0,60) |
| *A. schindleri* | 7 599 | 2 237 | 0,56 |

**Veredicto:** el orden se mantiene a igual N. *P. frigoritolerans* tiene de
verdad un pangenoma mas grande. alpha cambia poco con N, pero depende de que
genomas entren: en *B. altitudinis* va de 0,59 a 0,76 segun los 23 elegidos.
En la tesis conviene reportar el rango, no solo el valor puntual.

### 3.2 El pangenoma de *P. frigoritolerans* esta inflado?

- 38 % de sus familias estan en una sola referencia (22 % en *B. altitudinis*).
- La fragmentacion no lo explica (rho = -0,08 entre genes unicos y contigs).
- alpha sube de 0,57 a 0,74 si se usan solo los 7 genomas completos.
- Dos referencias aportan muchos genes propios: GCA_001636475.1 (~650) y la
  cepa divergente GCF_018613135.1 (~350).
- El accesorio esta muy estructurado por la filogenia (Mantel r = 0,73): la
  especie tiene linajes profundos (ANI minimo 95,8 %, cerca del limite de
  especie) con repertorios distintos (figura 11, matriz).

**Veredicto:** el tamano es en parte real (especie heterogenea, con linajes
profundos) y en parte atribuible a pocos genomas atipicos y borradores. No
cambia ninguna conclusion sobre el bin.

### 3.3 Tiene el bin "muchos" genes propios?

Se conto, para cada referencia, cuantas familias tiene que ninguna otra
referencia tiene, y se comparo con los candidatos del bin (figura 9).

| | Referencias: mediana (rango) | Bin | Percentil del bin |
|---|---|---:|---:|
| *B. altitudinis* | 18 (0-263) | 14 | 34 |
| *P. frigoritolerans* | 120 (16-655) | 34 | 2 |
| *A. schindleri* | 104 (20-308) | 70 | 26 |

Los genomas con un pariente mas lejano tienen algo mas de genes propios
(rho = 0,21-0,35).

**Veredicto:** ningun bin tiene un contenido propio inusual. Tener genes
exclusivos es lo normal en estas especies; por si mismo no indica adaptacion
al ambiente de la Estela.

### 3.4 De donde venian los candidatos que F4 descarto?

| | Candidatos descartados por F4 | Genoma donde se encontraron |
|---|---:|---|
| *B. altitudinis* | 9 | 6 en representantes no elegidos por el tope de 50; 3 en genomas eliminados por dRep |
| *P. frigoritolerans* | 26 | Todos en genomas eliminados por dRep del grupo del pariente del bin; 24 en uno solo (GCA_900000145.1) |
| *A. schindleri* | 20 | 17 en referencias **del propio pangenoma**; 3 en genomas eliminados por dRep |

- La desreplicacion y el tope de 50 dejan fuera genomas casi identicos que si
  comparten genes accesorios con el bin. F4 los recupera.
- En *A. schindleri*, 17 "candidatos" estaban en referencias del pangenoma,
  con 80-97 % de identidad: integrasas, transposasas y proteinas de fago que
  Panaroo separo en familias distintas por ser variantes divergentes o por su
  contexto.

**Veredicto:** el numero de candidatos esta inflado por la desreplicacion y
por como Panaroo agrupa los genes moviles. F4 era imprescindible; reportar los
candidatos sin F4 habria sido un error.

### 3.5 bin-1-54 es una cepa nueva?

Su pariente (GCA_024160055.1, suelo, Corea) representa un grupo de 8 genomas a
>= 99 % de ANI: manzana podrida y suelo agricola (Alemania), semilla (Francia),
suelo (Corea), nieve (Antartida), sala limpia y piel humana (EE. UU.) y uno sin
datos. bin-1-54 esta a 0,0017 sustituciones/sitio de ese representante.

**ANI (fastANI, bin contra los 69 genomas descargados de la especie):**
99,98 % con GCA_900000145.1 (sin datos de origen), 99,75 % con GCF_025142885.1
(piel humana, EE. UU.) y 99,74 % con GCF_024159205.1 (sala limpia, EE. UU.).

**Veredicto:** confirmado. La cepa de la Estela es practicamente la misma que
GCA_900000145.1 y pertenece a un clon cosmopolita y generalista. Eso explica
que no tenga genes propios.

Para comparar, los otros bins: bin-5-63 tiene 99,06 % con GCF_900119345.1 y
GCF_900188195.1, dos genomas que no estan en el pangenoma (los dejo fuera la
desreplicacion o el tope de 50; son los mismos donde F4 encontro candidatos),
asi que tambien es miembro de un clon ya secuenciado. bin-4-52 tiene como
maximo 98,39 % (GCF_025514435.1, su hermana en el arbol): es una cepa distinta
de todas las secuenciadas.

### 3.6 Los exclusivos verificados son de verdad exclusivos?

F5 descarta un gen si su proteina aparece en `nr` en un organismo con el
**nombre** de la misma especie. Los genomas nuevos suelen depositarse como
"sp." o con sinonimos. Para cada verificado se consulto el IPG de su mejor
parecido (todos los genomas con la proteina identica, 341 ensamblajes) y el ANI
que NCBI calcula de cada uno contra las cepas tipo
(`phase3_f5b_especie_por_ani.tsv`).

- Los 2 de *B. altitudinis* (serina proteasa y DUF4342) estan en *Bacillus*
  sp. 22483 (GCF_052044605.1), cuyo mejor ANI es 98,3 % con *B. aerius*,
  sinonimo de *B. altitudinis* (su cepa de referencia es nuestra T0
  GCF_029894105.1).
- 3 de *A. schindleri* (transportador de potasio, regulador AraC, polifosfato
  quinasa) estan en *Acinetobacter* sp. 22323 (GCF_052044895.1): ANI 97,75 %
  con la cepa tipo de *A. schindleri*.
- Ambas cepas: rizosfera, China (2023), BioProject PRJNA1299421, publicadas el
  13-08-2025, despues de GTDB R220.
- La sospecha sobre "*Acinetobacter* sp. UBA2581" (un MAG, 5 genes) se
  descarto: su ANI es 97,6 % con *A. variabilis*, otra especie.

**Veredicto:** 5 de 32 verificados no son exclusivos. Corregido: 0 / 0 / 27.
Limitacion que sigue: solo se reviso el mejor parecido de cada gen y F5 pidio
10 parecidos por proteina; una proteina compartida por cientos de genomas
(por ejemplo MgtC, 210 ensamblajes) puede tener una copia en la especie fuera
de esos 10.

### 3.7 Que quedo sin examinar?

Panaroo, en su limpieza, elimina genes poco respaldados, que trata como
posibles errores o contaminacion (Tonkin-Hill et al. 2020). Segun su codigo
(v1.8.0, `set_default_args.py`), en modo `moderate` elimina todo gen que este
en el extremo de un contig y aparezca en menos de 2 genomas
(`min_trailing_support = max(2, 1 % de los genomas)`), y repite la poda hacia
adentro hasta encontrar un gen compartido (`trailing_recursive` ilimitado). En
los bins, 136 de los 140 genes eliminados forman justamente tramos continuos
desde el borde de un contig; los otros 4 son proteinas muy cortas (31-45 aa) y
una transposasa. Esos genes no entran a ninguna familia y por eso nunca fueron
candidatos.

| | Genes del bin eliminados | Contigs eliminados enteros |
|---|---:|---|
| *B. altitudinis* | 32 | 13 contigs cortos (1,5-3,4 kb; 26 kb), GC normal: fragmentos |
| *P. frigoritolerans* | 7 | 1 contig de 3,3 kb |
| *A. schindleri* | 101 | **4 contigs (40 kb, 34 genes)**: contig_19 (14 kb), contig_73 (19 kb), contig_1 y contig_68; 3 con GC 0,35 frente a 0,42 del bin |

Los 101 genes de *A. schindleri* tienen largo normal (mediana 267 aa) y muchos
son de sistemas de defensa o elementos moviles (metiltransferasas de
restriccion-modificacion, helicasas DEAD/DEAH, dinamina, dominios WYL).

**Veredicto:** no es una falla de Panaroo sino un filtro pensado para
pangenomas de aislados, que va en contra de lo que aqui se busca: un gen
especifico esta, por definicion, en un solo genoma, y un MAG tiene muchos
extremos de contig. El "0 contigs huerfanos" de F6 era un artefacto del orden
de los pasos. Un GC 7 puntos mas bajo que el del genoma es tipico de ADN
adquirido (plasmidos, fagos) o de contaminacion de otro organismo.

**v2 (corrida local del 08-10):** estos genes ya son candidatos.

| | Candidatos (de ellos, eliminados por Panaroo) | Pasan F1 | Pasan F4 | Contigs a F6 |
|---|---:|---:|---:|---:|
| *B. altitudinis* | 46 (32) | 25 | 14 | 6 |
| *P. frigoritolerans* | 41 (7) | 29 | 1 | 0 |
| *A. schindleri* | 171 (101) | 137 | 97 | 4 |

Pendiente: F5, F5b y F6 de la parte remota.

### 3.8 El habitat explica el contenido accesorio?

| | Mantel accesorio ~ filogenia | PERMANOVA por habitat |
|---|---|---|
| *B. altitudinis* | r = 0,25, p = 0,001 | R2 = 0,17, p ~ 0,35 (38 genomas, 7 habitats) |
| *P. frigoritolerans* | r = 0,73, p = 0,001 | R2 = 0,04, p ~ 0,3 (28, 2) |
| *A. schindleri* | r = 0,54, p = 0,001 | R2 = 0,15, p ~ 0,06 (17, 3) |

**Veredicto:** los genomas emparentados comparten accesorio; el habitat no
agrega una senal detectable, igual que en el grupo *B. pumilus* (Fu et al.
2021, "el linaje tiene prioridad sobre el nicho"). Advertencia: pocos genomas
por habitat y rotulos gruesos; "sin efecto detectable" no es "sin efecto". El
pangenoma por si solo no permite afirmar adaptacion al sustrato petreo.

### 3.9 Es confiable el arbol?

- Los tres bins se ubican con UFBoot 100.
- La rama propia de bin-5-63 (0,0083) es comparable a la de su pariente
  (0,0077): la incompletitud no la alarga.
- En *B. altitudinis* 18 de 48 nodos tienen UFBoot < 95: el esqueleto
  profundo de la especie esta mal resuelto (ramas cortas, posible
  recombinacion). No afecta la posicion del bin.
- No se corrigio la recombinacion; la raiz es de punto medio (solo para
  dibujar).

**Veredicto:** la posicion de cada bin es solida; las relaciones profundas de
*B. altitudinis*, no.

### 3.10 Las referencias estan sesgadas?

- *A. schindleri*: 11 de 23 clinicas; el bin cae entre cepas clinicas. Es lo
  que se ha secuenciado, no una afinidad clinica del bin.
- *P. frigoritolerans*: suelo y plantas; *B. altitudinis*: variado.
- Se excluyeron MAGs de las referencias y de F4 (8 en *A. schindleri*, 3 en
  cada *Bacillus*). Ninguno resulto ser el origen de un gen especifico de la v1
  (ver 3.6), pero F4 no los reviso. Se decidio no ampliar F4 (propuesta C).
- GTDB R220 no incluye genomas posteriores; justamente ahi aparecieron las dos
  cepas chinas de 2025.

### 3.11 Que son los exclusivos que quedan?

- *A. schindleri*: 27 genes en 8 contigs, la mayoria junto a integrasas,
  recombinasas, transposasas IS y proteinas de fago (figura 12). Identidad del
  mejor parecido: >= 99 % en 16, 90-99 % en 7, < 90 % en 4. Hipoteticas:
  24 % de los candidatos frente a 8 % del genoma.
- Origen: 24 en otras especies de *Acinetobacter* (*A. lwoffii*,
  *A. radioresistens*, *A. variabilis*, *A. johnsonii*...), 1 en un fago, 1 en
  *Cloacibacterium* (Bacteroidota, 99 %), 1 en *Thauera* (55 %).

**Veredicto:** son genes adquiridos por transferencia horizontal, en buena
parte recientes (identidad casi total con otras especies). Su relevancia
funcional para el biodeterioro es una pregunta para la Fase 4. Los tres genes
con funcion de estres mas sugerente (transporte de potasio, polifosfato
quinasa, AraC) son justamente los que tambien tiene la cepa china de
rizosfera, asi que no son propios de la Estela.

---

## 4. Figuras

### 4.1 Existentes (1-7)

Curvas de acumulacion, diagnostico de Heaps, posicion del bin, tres arboles del
core y embudo de exclusivos. Ver `phase3_guia_resultados.md`.

### 4.2 Nuevas (8-12)

Generadas con `figuras/R/f3_06_analisis_adicional.R`.

**Figura 8. Espectro de frecuencias.** Cuantas familias estan en exactamente
1, 2, ..., N referencias. Los pangenomas bacterianos tienen forma de U: muchas
familias raras (izquierda), muchas universales (derecha) y pocas intermedias
(Horesh et al. 2021). Se cumple en las tres especies; *P. frigoritolerans*
tiene la rama izquierda mucho mas alta (6 103 familias en un solo genoma).

![Espectro](../../figuras/figs/fig_f3_espectro_frecuencias.png)

**Figura 9. Genes propios frente a distancia al pariente mas cercano.** Cada
punto gris es una referencia; el rojo, el bin. Responde si el bin tiene mas
genes propios de lo esperable (seccion 3.3).

![Unicos vs vecino](../../figuras/figs/fig_f3_unicos_vecino.png)

**Figura 10. Ordenacion del genoma accesorio.** Genomas con accesorio parecido
quedan cerca; color = habitat. No se ven grupos por habitat; Mantel y
PERMANOVA en el subtitulo (seccion 3.8).

![PCoA](../../figuras/figs/fig_f3_pcoa_accesorio.png)

**Figura 11. Matriz de presencia/ausencia junto al arbol.** Cada fila es un
genoma (orden del arbol) y cada columna una familia accesoria. Bloques de
familias compartidas por clados vecinos muestran que el accesorio sigue a la
filogenia. Es la figura estandar de Roary/Panaroo.

![Matriz L16](../../figuras/figs/fig_f3_matriz_L16_Bacillus_altitudinis.png)
![Matriz L12](../../figuras/figs/fig_f3_matriz_L12_Peribacillus_frigoritolerans.png)
![Matriz L8](../../figuras/figs/fig_f3_matriz_L8_Acinetobacter_schindleri.png)

**Figura 12. Contexto genomico de los exclusivos.** Genes alrededor de cada
exclusivo, coloreados por categoria; triangulo = gen movil. Muestra que los
exclusivos estan en islas junto a integrasas, transposasas y fagos.

![Contexto](../../figuras/figs/fig_f3_contexto_exclusivos.png)

### 4.3 Posibles, pero requieren datos que estan en Khipu

| Figura | Para que | Que falta |
|---|---|---|
| ANI del bin contra todos los genomas de su especie | Confirmar el clon cosmopolita de *P. frigoritolerans* y medir que tan nueva es cada cepa | **Hecho (08-10, `06c_ani_bins.slurm`):** resultados en 3.5; la figura se hace con la exportacion final |
| Categorias funcionales (COG) por core/shell/cloud/exclusivo | Figura habitual en los estudios de pangenoma (por ejemplo, Fu et al. 2021) | **Descartada en la Fase 3 (08-10):** Bakta anota la categoria COG en 4-31 % de los CDS segun el genoma y KEGG/EC/GO en 17-26 %; una figura con esa cobertura mostraria cuanto se sabe de cada gen, no que hace. Una anotacion uniforme (DRAM o eggNOG-mapper sobre las familias) se dejo para la Fase 4 |
| Ganancia y perdida de genes sobre el arbol | Ubicar cuando entro cada isla | Software adicional (Count, GLOOME); prescindible |
| Correccion por recombinacion (Gubbins, ClonalFrameML) | Mejorar el esqueleto del arbol de *B. altitudinis* | Prescindible para la posicion del bin |
| Asociacion gen-habitat (pan-GWAS, Scoary) | Buscar genes ligados a ambientes petreos | **No recomendable:** 1-2 genomas petreos por especie, sin potencia |

---

## 5. Lo que se puede afirmar en la tesis

1. Los pangenomas de las tres especies son abiertos, con un core de 2 400 a
   4 000 familias.
2. Cada cepa de la Estela se ubica con soporte maximo junto a una cepa
   conocida de otro continente y otro ambiente (cereales fermentados, suelo,
   cepa clinica). Ninguna se agrupa con las cepas de ambientes petreos.
3. El contenido accesorio sigue al parentesco y no al habitat de aislamiento.
4. Las cepas de la Estela no tienen mas genes propios que cualquier otra cepa
   de su especie. Los genes especificos de la cepa (v1: 0, 0 y 27; cifra final
   pendiente de la v2) son en su mayoria ADN movil adquirido de otras especies
   del mismo genero.
4b. La cepa de *P. frigoritolerans* de la Estela es practicamente identica
   (ANI 99,98 %) a una cepa ya secuenciada y pertenece a un clon cosmopolita;
   la de *B. altitudinis* tambien pertenece a un clon conocido (99,06 %); la de
   *A. schindleri* es una cepa distinta de las secuenciadas (98,39 %).
5. La exclusividad depende de la base de datos y de la definicion de especie:
   con UniProtKB en lugar de `nr` habrian pasado 7 genes mas, y con nombres de
   NCBI en lugar de ANI, 5 mas.

---

## 6. Correcciones propuestas y decisiones

Decididas el 08-10 e implementadas en el commit `2dfbe1f`.

| # | Propuesta | Costo | Recomendacion | Decision |
|---|---|---|---|---|
| A | Incorporar a F5 la comprobacion por ANI de los genomas portadores (IPG + ANI de NCBI) para todos los parecidos >= 80/80, no solo el mejor | Bajo (minutos, desde el login o la laptop) | Si: alinea F5 con la definicion de especie del resto del flujo (GTDB/ANI) | **Si**: paso F5b (`08_exclusivos.py ani`) |
| B | Pasar por F1-F6 los 140 genes que Panaroo elimino de los bins (incluye F6 para los 4 contigs de *A. schindleri*) | Medio (una corrida local + BLAST remoto) | Si, al menos los contigs eliminados enteros | **Si**, todos (tras confirmar la causa en el codigo de Panaroo) |
| C | Ampliar F4 a todos los genomas del cluster GTDB (incluidos MAGs y los que no pasaron el QC: 8, 10 y 8) | Bajo | Si: para buscar presencia la calidad importa menos | **No** |
| D | `fastANI` de cada bin contra sus referencias | Muy bajo | Si: confirma 3.5 y da la figura de ANI | **Si**: hecho (`06c_ani_bins.slurm`) |
| E | Renombrar "verificado" a "especifico de la cepa" en tablas y figuras | Nulo | Opcional | **Si**: clases `especifico`, `especifico_con_bandera`, `descartado` |
| F | Marca de advertencia para exclusivos con >= 99 % de identidad con otra especie ("posible transferencia reciente") | Nulo | Opcional; afectaria a 16 de 27 | **Si**: >= 99 % de identidad y >= 90 % de cobertura |

---

## 7. Referencias

- Brockhurst MA, et al. (2019). The ecology and evolution of pangenomes. *Curr Biol* 29(20):R1094-R1103. doi:10.1016/j.cub.2019.08.012
- Fu X, et al. (2021). *Bacillus pumilus* group comparative genomics: toward pangenome features, diversity, and marine environmental adaptation. *Front Microbiol* 12:571212. doi:10.3389/fmicb.2021.571212
- Guerra A. (2026). The pangenome: a statistical model, not a fixed biological property. *Bioinform Adv* 6(1):vbag069.
- Horesh G, et al. (2021). Different evolutionary trends form the twilight zone of the bacterial pan-genome. *Microb Genom* 7(9):000670. doi:10.1099/mgen.0.000670
- Li T, Yin Y. (2022). Critical assessment of pan-genomic analysis of metagenome-assembled genomes. *Brief Bioinform* 23(6):bbac413. doi:10.1093/bib/bbac413
- McInerney JO, McNally A, O'Connell MJ. (2017). Why prokaryotes have pangenomes. *Nat Microbiol* 2:17040. doi:10.1038/nmicrobiol.2017.40
- Tettelin H, et al. (2008). Comparative genomics: the bacterial pan-genome. *Curr Opin Microbiol* 11(5):472-477.
- Tonkin-Hill G, et al. (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol* 21:180. doi:10.1186/s13059-020-02090-4
