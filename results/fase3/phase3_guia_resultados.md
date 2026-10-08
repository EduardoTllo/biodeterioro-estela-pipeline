# Fase 3 - Guia de lectura de los resultados

Esta guia explica el reporte de seleccion (`phase3_seleccion_report.md`) y las
figuras de la Fase 3: que muestra cada una, como se lee y que dicen nuestros
datos. Corresponde a la exportacion del 08-10-2026 (Fase 3 completa).

**Importante:** el analisis critico (`phase3_analisis_critico.md`) corrige
dos resultados de esta guia: los exclusivos verificados pasan de 2/0/30 a
**0/0/27** al comprobar por ANI la especie de los genomas que tienen cada
proteina, y F6 no vio contigs huerfanos porque Panaroo ya los habia eliminado.
Tambien agrega las figuras 8-12.

Especies analizadas:

| Linaje | Especie (GTDB) | Bin de la Estela | Completitud del bin |
|---|---|---|---|
| L16 | *Bacillus altitudinis* | bin-5-63 | 72,5 % |
| L12 | *Peribacillus frigoritolerans* | bin-1-54 | 92,3 % |
| L8 | *Acinetobacter schindleri* | bin-4-52 | 99,9 % |

---

## 1. El reporte de seleccion (`phase3_seleccion_report.md`)

Responde dos preguntas: **que especies** entran al analisis pangenomico y **con
que genomas de referencia** se construye el pangenoma de cada una.

### 1.1 Tabla de ranking

| Columna | Significado |
|---|---|
| Rango | Posicion de la especie entre las viables, ordenadas por genomas no redundantes. `NA` = no viable |
| Linaje | Linaje de la Fase 2 al que pertenece el bin |
| Especie GTDB | Especie asignada por GTDB-Tk (R220) |
| Post-QC | Genomas publicos de la especie en GTDB que pasan el control de calidad: completitud >= 95 %, contaminacion <= 5 %, <= 300 contigs, solo aislados |
| Descargados | Genomas descargados de NCBI y verificados (longitud igual a la de GTDB), ya sin los excluidos a mano (D30) |
| No redundantes | Genomas que quedan despues de agrupar con dRep los que comparten >= 99 % de identidad (ANI) y quedarse con uno por grupo |
| Viable | `si` si tiene >= 15 genomas no redundantes |
| Seleccionada | `si` para las 3 viables con mas genomas no redundantes |

La diferencia entre **Post-QC** y **Descargados** son las 4 exclusiones de D30:
genomas que GTDB pone en la especie pero que NCBI rotula como bacterias de
otro filo (3 en L12, 1 en L16). Se sacaron por riesgo de ensamblaje
contaminado.

La diferencia entre **Descargados** y **No redundantes** muestra cuanta
diversidad real tiene cada especie en las bases publicas:

| Especie | Descargados | No redundantes | Lectura |
|---|---:|---:|---|
| *B. altitudinis* | 208 | 137 | Mucha diversidad; se usan 50 (tope) |
| *P. frigoritolerans* | 69 | 42 | Se usan las 42 |
| *A. schindleri* | 29 | 23 | Se usan las 23 |
| *C. firmus_B* | 18 | 12 | No llega a 15 |
| *S. warneri_A* | 24 | 4 | Casi todos los genomas son variantes de pocas cepas |
| *P. flexa* | 29 | 4 | Igual |
| *B. licheniformis* | 279 | **1** | Los 279 genomas son practicamente una sola cepa al 99 % (especie casi clonal en las bases publicas) |

Por eso *B. licheniformis* queda fuera aunque es la especie con mas genomas:
279 copias casi identicas no permiten estimar un pangenoma. El analisis de
sensibilidad (99 / 99,5 / 99,9 %) confirmo que el resultado no depende de un
umbral demasiado estricto (decision D32).

### 1.2 Composicion de las referencias

Para cada especie seleccionada el reporte resume las referencias elegidas en
tres ejes.

**Nivel de prioridad** (como se eligieron, decision D31):

| Nivel | Que es | Efecto |
|---|---|---|
| T0 | Genoma representante de la especie en GTDB o cepa tipo | Entra siempre (peso +1000 en dRep: si esta en un grupo de genomas casi identicos, es el elegido) |
| T1 | Aislado de sustrato petreo o arido (roca, mineria, ceramica, arqueologico, desierto) | Entra despues de T0 (peso +500); es lo mas parecido al ambiente de la Estela |
| T2 | El resto | Completa hasta 50, rotando entre habitats y, dentro de cada habitat, entre continentes, para maximizar la diversidad |

Solo en *B. altitudinis* hubo que elegir (137 > 50). En las otras dos entran
todas las no redundantes y el nivel es solo informativo.

Las T0 y T1 seleccionadas:

| Especie | Nivel | Genoma | Origen |
|---|---|---|---|
| *B. altitudinis* | T0 | GCF_029894105.1 | EE. UU., sin fuente declarada |
| *B. altitudinis* | T0 | GCF_000691145.1 | India, tubos para muestras de aire de gran altitud |
| *B. altitudinis* | T1 | GCF_010747395.1 | Portugal, residuo minero |
| *P. frigoritolerans* | T0 | GCF_001636405.1 | Alemania, suelo |
| *P. frigoritolerans* | T0 | GCA_021012855.1 | Espana, heces humanas |
| *P. frigoritolerans* | T0 | GCF_024169475.1 | Marruecos, sin fuente declarada |
| *P. frigoritolerans* | T1 | GCF_022603155.1 | Italia, anfora romana |
| *A. schindleri* | T0 | GCF_000368625.1 | Republica Checa, orina |

**Habitat**: categoria asignada a partir de la fuente de aislamiento
(`isolation_source`) con las reglas de `metadata/habitat_keywords.tsv`.
`desconocido` = la fuente esta vacia o solo nombra una coleccion de cultivos;
`sin_clasificar` = hay texto, pero no coincide con ninguna categoria.

**Continente**: a partir del pais de aislamiento. `Oceano` = muestras de mar
abierto; `desconocido` = sin pais.

Lectura por especie:

- ***B. altitudinis*** tiene las referencias mas variadas: 10 habitats y
  8 continentes, sin ninguno dominante. Es la mejor muestra de la diversidad
  global de la especie.
- ***P. frigoritolerans*** es sobre todo de suelo y plantas (28 de 42) y de
  Asia y Europa (30 de 42). Es una especie de suelo; la muestra refleja eso.
- ***A. schindleri*** tiene **11 de 23 referencias clinicas** (orina, sangre,
  etc.). Es un sesgo de las bases publicas (se secuencian mas cepas de
  hospital), no de la seleccion: se usaron todas las disponibles. Hay que
  tenerlo en cuenta al interpretar que el bin se parezca a cepas clinicas.

---

## 2. Figuras

Todas estan en `figuras/figs/` (PNG y PDF) y sus tablas en `figuras/tablas/`.
Las frecuencias del pangenoma se calculan **solo con las referencias** (el bin
no cuenta), porque el bin es incompleto y haria pasar genes del core por
accesorios (decision D22).

### Figura 1. Curvas de acumulacion del pangenoma

![Curvas de acumulacion](../../figuras/figs/fig_f3_curvas_acumulacion.png)

**Que muestra.** Como crecen el pangenoma (rojo: todas las familias de genes
vistas en al menos un genoma) y el genoma core (azul: familias presentes en
todos) a medida que se agregan genomas de referencia. Se repite con 100
ordenes aleatorios; la linea es la mediana y la banda el rango intercuartil.

**Como se lee.**
- Si la curva roja se aplana, nuevos genomas ya no aportan genes nuevos:
  pangenoma **cerrado**. Si sigue subiendo, cada cepa trae genes propios:
  pangenoma **abierto**.
- La curva azul baja y se estabiliza: el core es lo que comparten todas las
  cepas.

**Que dicen los datos.** En las tres especies la curva roja sigue subiendo con
el ultimo genoma: los tres pangenomas son abiertos. *P. frigoritolerans*
destaca: con 42 genomas acumula 16 125 familias, el doble que las otras dos con
un numero parecido de genomas.

| | Core (>= 95 %) | Shell (15-95 %) | Cloud (< 15 %) | Pangenoma |
|---|---:|---:|---:|---:|
| *B. altitudinis* (50 refs) | 3 303 | 847 | 2 909 | 7 059 |
| *P. frigoritolerans* (42) | 3 961 | 2 142 | 10 022 | 16 125 |
| *A. schindleri* (23) | 2 389 | 1 114 | 4 096 | 7 599 |

*Core:* genes de casi todas las cepas (funciones basicas). *Shell:* genes
frecuentes pero no universales. *Cloud:* genes raros, de pocas cepas
(adaptaciones especificas, elementos moviles, o errores de anotacion).

### Figura 2. Diagnostico de la ley de Heaps

![Diagnostico de Heaps](../../figuras/figs/fig_f3_diagnostico_heaps.png)

La ley de Heaps (Tettelin et al. 2008) resume la curva roja en un numero,
**alpha**: con alpha < 1 el pangenoma es abierto, con alpha > 1 es cerrado;
cuanto mas bajo, mas abierto. Esta figura comprueba si esa apertura es real o
esta inflada por genomas de mala calidad.

**Panel a. Genes unicos frente a fragmentacion.** Cada punto es una
referencia: cuantas familias tiene que ninguna otra referencia tiene (eje y)
frente a en cuantos contigs esta partido su ensamblaje (eje x). Si los
ensamblajes mas fragmentados tuvieran mas genes unicos, la apertura vendria de
genes partidos en pedazos (artefacto). En rojo, genomas atipicos (muchos mas
genes unicos que el resto).

- *B. altitudinis* (rho = -0,05) y *P. frigoritolerans* (rho = -0,08): sin
  relacion; la fragmentacion no infla.
- *A. schindleri* (rho = -0,52, p = 0,01): relacion **negativa**, lo contrario
  de la inflacion.
- Atipicos: 4 en *B. altitudinis*, 2 en *P. frigoritolerans* (uno es
  GCF_018613135.1, la cepa divergente revisada en el control C6; el otro,
  GCA_001636475.1, aporta ~650 familias propias) y 1 en *A. schindleri*.

**Panel b. alpha con y sin genomas en borrador.** Azul: alpha con todas las
referencias; la barra es el rango al quitar un genoma cada vez (jackknife):
si es estrecha, ningun genoma solo cambia el resultado. Naranja: alpha solo
con genomas completos (cerrados), que no tienen problemas de fragmentacion.

| Especie | alpha (todas) | alpha (solo completos) |
|---|---|---|
| *B. altitudinis* | 0,66 (n = 50) | 0,64 (n = 22) |
| *P. frigoritolerans* | 0,57 (n = 42) | **0,74 (n = 7)** |
| *A. schindleri* | 0,55 (n = 23) | 0,55 (n = 7) |

En *B. altitudinis* y *A. schindleri* los dos valores coinciden: la apertura es
robusta. En *P. frigoritolerans* alpha sube a 0,74 con solo completos: sigue
abierto, pero menos. Parte de su apertura puede venir de genomas en borrador y
de los dos atipicos; con 7 genomas completos la estimacion es imprecisa. Se
reportan los dos valores.

**Panel c. Proteinas hipoteticas por categoria.** Porcentaje de familias sin
funcion conocida en core, shell y cloud. Es normal que el cloud tenga mas
(genes raros, poco estudiados); si fuera casi todo hipotetico, sugeriria genes
espurios. Aqui el core tiene 13-17 % y el cloud 29-45 %: el cloud no esta
dominado por hipoteticas.

**Conclusion de la figura.** Los tres pangenomas son abiertos y la apertura no
es un artefacto de la calidad de los ensamblajes; en *P. frigoritolerans* es
algo menor de lo que indica el alpha con todas las referencias.

### Figura 3. Posicion del bin en el pangenoma

![Posicion del bin](../../figuras/figs/fig_f3_posicion_bin.png)

**Panel a. A que categoria pertenecen los genes del bin.** Cada barra es el
100 % de los genes del bin de la Estela, repartidos segun la categoria que
tiene su familia en el pangenoma de referencias:

- **Core / shell / cloud:** genes que el bin comparte con otras cepas.
- **Exclusivo (candidato):** genes de familias que no tiene ninguna
  referencia. Son candidatos a adaptaciones propias de la Estela y pasan por
  los filtros F1-F6 antes de aceptarse (14, 34 y 70 candidatos).
- **Fuera del pangenoma:** genes que Panaroo descarto en su limpieza
  (fragmentos o genes dudosos); 32, 7 y 101.

La mayor parte de cada bin es core (como se espera de una cepa de la especie).
*A. schindleri* tiene mas shell y cloud: el bin esta casi completo y sus
referencias son pocas (23).

**Panel b. Recuperacion del core frente a completitud.** Eje x: completitud
del bin segun CheckM2 (Fase 1). Eje y: porcentaje del core de las referencias
que se encontro en el bin. La diagonal es la igualdad y la banda gris +/- 10
puntos. Si el bin cae dentro de la banda, sus ausencias se explican por
incompletitud y no por un problema del analisis.

| Bin | Completitud | Core recuperado |
|---|---:|---:|
| bin-5-63 (*B. altitudinis*) | 72,5 % | 75,2 % |
| bin-1-54 (*P. frigoritolerans*) | 92,3 % | 85,7 % |
| bin-4-52 (*A. schindleri*) | 99,9 % | 99,0 % |

Los tres estan dentro de la banda. Consecuencia practica: que un gen **este**
en el bin es confiable; que **falte** no lo es, sobre todo en bin-5-63, al que
le falta ~25 % del genoma.

### Figuras 4-6. Arbol del core de cada especie

![Arbol B. altitudinis](../../figuras/figs/fig_f3_arbol_core_L16_Bacillus_altitudinis.png)

![Arbol P. frigoritolerans](../../figuras/figs/fig_f3_arbol_core_L12_Peribacillus_frigoritolerans.png)

![Arbol A. schindleri](../../figuras/figs/fig_f3_arbol_core_L8_Acinetobacter_schindleri.png)

**Que muestran.** El parentesco entre el bin (rojo, grande) y las referencias,
calculado con IQ-TREE sobre el alineamiento de los genes core (2-3 millones de
posiciones). El color de cada punta es el habitat de origen de la referencia.

**Como se leen.**
- Dos puntas unidas por un nodo cercano son cepas muy parecidas; la longitud
  horizontal de las ramas es la cantidad de cambios (la barra de escala indica
  0,001 o 0,002 sustituciones por posicion).
- Punto negro en un nodo: soporte UFBoot >= 95 (el agrupamiento es confiable).
- La raiz se puso en el punto medio solo para dibujar; no indica el ancestro
  real.

**Que dicen los datos.** En las tres especies el bin se agrupa con una
referencia concreta con soporte 100:

| Bin | Pariente mas cercano | Distancia | Lectura |
|---|---|---:|---|
| bin-5-63 (*B. altitudinis*) | GCF_000828455.1, cereales fermentados, Paises Bajos | 0,016 | Forma con ella un grupo separado del resto de la especie. La rama del bin (0,0083) es comparable a la de su pariente (0,0077) |
| bin-1-54 (*P. frigoritolerans*) | GCA_024160055.1, suelo, Corea del Sur | **0,0017** | Practicamente la misma cepa en el core |
| bin-4-52 (*A. schindleri*) | GCF_025514435.1, clinico (endovascular), EE. UU. | 0,013 | Dentro de un grupo de cepas clinicas, coherente con el sesgo de las referencias |

Ninguno de los bins se agrupa con las referencias de habitat petreo o arido
(T1): la cepa de la Estela no esta emparentada especificamente con otras
cepas de piedra. El habitat de la cepa mas cercana no indica el habitat de
origen del bin: refleja que cepas se han secuenciado.

---

### Figura 7. Embudo de verificacion de genes exclusivos

![Exclusivos](../../figuras/figs/fig_f3_exclusivos.png)

Un **candidato a exclusivo** es un gen del bin cuya familia no aparece en
ninguna de las referencias del pangenoma. Puede ser una adaptacion propia de la
cepa de la Estela, pero tambien un gen mal predicho, un gen que otras cepas si
tienen pero que no estaban entre las referencias, o contaminacion. Los filtros
F1-F6 separan esos casos (decision D26).

**Panel a. Cuantos candidatos sobreviven a cada filtro.**

| Filtro | Que descarta |
|---|---|
| F1 estructura | Pseudogenes, genes de < 100 aminoacidos y genes a < 100 pb del borde del contig (suelen ser genes partidos) |
| F4 especie ampliada | Genes presentes (>= 80 % identidad y cobertura) en cualquiera de los genomas descargados de la especie (hasta 208), no solo en los <= 50 del pangenoma |
| F6 contigs huerfanos | Genes en contigs sin genes de la especie cuyo mejor parecido en `core_nt` es de otro genero (contaminacion). No hubo contigs huerfanos: ningun candidato se descarto aqui |
| F5 nr (misma especie) | Genes cuya proteina aparece en otra cepa de la misma especie en `nr` de NCBI (>= 80 % identidad y cobertura): no son exclusivos, solo faltaban en GTDB |
| Verificados | Lo que queda. F3 (presencia en otro bin de la misma muestra) solo marca con advertencia y no marco ninguno |

| | Candidatos | F1 | F4 | F5 | Verificados |
|---|---:|---:|---:|---:|---:|
| *B. altitudinis* | 14 | -3 | -9 | 0 | **2** |
| *P. frigoritolerans* | 34 | -7 | -26 | -1 | **0** |
| *A. schindleri* | 70 | -14 | -20 | -6 | **30** |

F4 es el filtro que mas descarta en las *Bacillus*: casi todos sus "exclusivos"
eran genes de cepas de la especie que no entraron entre las referencias.
*P. frigoritolerans* termina sin ningun exclusivo verificado.
*A. schindleri* conserva 30, coherente con un bin casi completo (99,9 %) frente
a solo 23 referencias.

**Panel b. Origen probable de los verificados**, segun la especie de su mejor
parecido en `nr`:

| Origen | *B. altitudinis* | *A. schindleri* |
|---|---:|---:|
| Mismo genero | 2 | 27 |
| Mismo filo | 0 | 1 |
| Otro filo | 0 | 2 |

Casi todos los exclusivos tienen parientes en otras especies del mismo
genero, muchos con 99-100 % de identidad: no son genes "nuevos", sino genes
que circulan entre especies cercanas y que esta cepa adquirio y sus
conespecificas no. Los dos de "otro filo" en *A. schindleri* son una proteina
de fago (mejor parecido: un virus, *Caudoviricetes*) y un regulador TetR
(99 % con *Cloacibacterium*, Bacteroidota).

**Contexto genomico.** Los exclusivos no estan dispersos: se agrupan en pocos
contigs junto a integrasas, proteinas de fago y transportadores. Son islas
genomicas y profagos, es decir, ADN movil:
- *B. altitudinis*: una serina proteasa y una proteina con dominio DUF4342,
  consecutivas al inicio de un contig, justo despues de una integrasa, e
  identicas (100 %) a proteinas de *B. subtilis*. Indica una transferencia
  horizontal reciente desde el grupo *B. subtilis*. El contig continua con
  genes core de *B. altitudinis*, por lo que no es contaminacion.
- *A. schindleri*: grupos en 6 contigs (contig_40, 63, 69, 71, 91 y 93) con
  integrasas, proteinas de fago (HK97, P1), transportadores de potasio y
  magnesio (MgtC) y una polifosfato quinasa; muchos identicos a proteinas de
  *A. lwoffii*, *A. radioresistens* y otras especies del genero.

**Sensibilidad a la base de datos (NCBI frente a EBI).** F5 se corrio tambien
contra UniProtKB en el EBI, porque la cola de NCBI tardo hasta 7 h por
especie. UniProtKB no sirve como reemplazo: elimina proteomas redundantes y
le faltan muchos genomas de generos muy secuenciados.
- No encontro en la especie 7 genes que `nr` si encontro (1 en
  *P. frigoritolerans*, 6 en *A. schindleri*): con UniProtKB habrian pasado
  como exclusivos.
- En *B. altitudinis* asigno los 2 genes a otra familia y a otro genero del
  filo (identidad 57-66 %) cuando `nr` tiene proteinas identicas en
  *B. subtilis*: exagera la distancia del origen.

El resultado oficial es el de NCBI (`<especie>/exclusivos_verificados.tsv`);
el del EBI queda en `<especie>/ebi/` como comparacion.
