# Pies de figura y de tabla — Fase 2 (OE2)

Textos listos para pegar en la tesis. Ajusta la numeracion segun el capitulo.

---

## Figura 1. Embudo de seleccion de genomas desde los bins crudos hasta los linajes.

**(a)** Numero de genomas que sobrevive a cada etapa del flujo: 97 bins crudos
procedentes de 22 librerias Illumina, 66 con dominio mayoritariamente procariota
segun Tiara, 21 que superan los umbrales de calidad de CheckM2 (completitud
> 70 %, contaminacion < 5 %) y 19 linajes tras la desreplicacion con dRep a
95 % de ANI. Las cifras en blanco indican las perdidas en cada paso.
**(b)** Motivos de descarte de los 76 bins eliminados.
**(c)** Composicion por dominio de cada uno de los 97 bins segun Tiara,
expresada como porcentaje de la longitud del bin clasificada en cada categoria;
los bins se separan segun su dominio mayoritario y se ordenan por la fraccion
procariota.

---

## Figura 2. Calidad y metricas de ensamblaje de los genomas.

**(a)** Completitud frente a contaminacion (CheckM2) de los 66 bins procariotas.
El eje de contaminacion usa una escala logaritmica (log(1+x)) porque unos pocos
bins superan el 60 % y comprimirian la region de interes. Las lineas discontinuas
marcan los umbrales del estudio (completitud >= 70 %, contaminacion < 5 %) y la
region sombreada, los 21 genomas seleccionados.
**(b)** Tamano del genoma frente a contenido GC de los 21 genomas
seleccionados; el tamano del punto es proporcional al N50 de los contigs y el
color indica el filo.
**(c)** Distribucion de N50, numero de contigs, densidad codificante y
completitud por filo. Cada punto es un genoma; las cajas muestran mediana y
cuartiles.

---

## Figura 3. Arbol filogenomico de los 21 genomas seleccionados con sus metadatos.

Arbol de Neighbor-Joining calculado sobre distancias-p a partir del alineamiento
concatenado de los 120 marcadores bacterianos de copia unica (bac120, 5035
posiciones) generado por GTDB-Tk (datos R220), enraizado en el filo
Pseudomonadota. Los numeros sobre las ramas son valores de soporte de bootstrap
(500 replicas); solo se muestran los >= 70 %. Las columnas de la derecha indican,
para cada genoma: la muestra de origen, el linaje asignado por dRep (el asterisco
marca el genoma representante del linaje; el fondo amarillo, los tres linajes que
pasan a la Fase 3) y los valores de completitud y contaminacion de CheckM2.

> **Nota metodologica:** se trata de un arbol de distancias sobre los genomas de
> consulta, no de un arbol de maxima verosimilitud con genomas de referencia de
> GTDB. El arbol de novo focalizado por linaje esta previsto para la Fase 4.

---

## Figura 4. Matriz de similitud genomica entre los 21 genomas seleccionados.

Similitud MASH (100 − distancia MASH) entre todas las parejas de genomas,
ordenadas por filo y genero. Solo dos parejas presentan similitud maxima: los dos
*Telluria timonae* (muestras 50 y 52) y los dos *Corynebacterium* sp. (muestras
68 y 69). El resto de los genomas son mutuamente distantes —incluido el par de
*Cytobacillus* (80 % de similitud MASH), que pese a compartir genero permanece
como dos linajes separados—, lo que explica que 21 genomas den lugar a 19
linajes.

Las dos parejas que dRep somete a fastANI dentro de su cluster primario superan
holgadamente el umbral de especie: *Corynebacterium* sp. (bin-1-68 ↔ bin-2-69)
alcanza 99,96 % de ANI con 97,1 % de cobertura del alineamiento, y *Telluria
timonae* (bin-3-52 ↔ bin-6-50), 99,78 % de ANI con 85,2 % de cobertura. Ambas
colapsan por tanto en un unico linaje. Los valores completos estan en la
Tabla S2.

> **Nota metodologica:** dRep ejecuta fastANI unicamente dentro de cada cluster
> primario definido por MASH, de modo que no existe una matriz de ANI completa
> para los 21 genomas; la matriz de todos contra todos disponible es la de MASH.

---

## Figura 5. Evaluacion de novedad taxonomica frente a la referencia GTDB mas cercana.

ANI frente a fraccion alineada (AF) de cada genoma respecto de su referencia GTDB
mas cercana. Las lineas discontinuas marcan los criterios de confirmacion de
especie (ANI >= 95 % y AF >= 65 %). Los 16 genomas del cuadrante sombreado tienen
especie confirmada; los tres que quedan fuera —los dos *Corynebacterium* sp.
(ANI ~90,5 %, AF ~44 %) y *Sporosarcina* sp. (ANI 91,3 %, AF 60,4 %)— son
candidatos a especies aun no descritas en GTDB R220. Dos genomas adicionales
(bin-1-51 *Pantoea* sp. y bin-9-71 *Aristophania* sp.) fueron clasificados
unicamente por topologia del arbol, sin calculo de ANI.

En conjunto, GTDB-Tk asigno especie a 16 de los 21 genomas (12 Bacillota, 3
Pseudomonadota y 1 Actinomycetota) y se detuvo en el rango de genero en los 5
restantes (1 Bacillota, 2 Pseudomonadota y 2 Actinomycetota). El detalle de esos
cinco genomas figura en la Tabla S1.

---

## Figura 6. Prevalencia espacial de los linajes y priorizacion para la Fase 3.

**(a)** Matriz de presencia de cada linaje (filas, ordenadas por prevalencia
descendente) en cada muestra (columnas). El color indica el filo.
**(b)** Prevalencia de cada linaje, definida como el numero de muestras distintas
en que aparece. En amarillo, los tres linajes priorizados para el analisis
pangenomico de la Fase 3: L1 (*Telluria timonae*), L2 (*Corynebacterium* sp.) y
L3 (*Bacillus_AB infantis*).
**(c)** Numero de genomas de alta calidad recuperados por muestra, coloreado por
filo.

> **Limitacion:** la prevalencia se calcula sobre el numero de muestras y no
> sobre zonas fisicas de la estela, porque no se dispone de un mapeo
> muestra → zona en formato analizable.

---

## Figura 7. Composicion taxonomica del conjunto de genomas seleccionados.

**(a)** Numero de genomas por genero, agrupados por filo.
**(b)** Numero de genomas por orden.
**(c)** Reparto de los 21 genomas entre los tres filos detectados.
La comunidad viable esta dominada por Bacillota (13 genomas, 62 %), un patron
compatible con superficies petreas expuestas, donde la formacion de endosporas
favorece la supervivencia.

---

## Figura S1. Calidad de la asignacion taxonomica.

**(a)** Porcentaje de columnas del alineamiento bac120 retenidas para cada genoma
(`msa_percent`). La linea roja marca el minimo del 50 % que exige GTDB-Tk para
emitir una clasificacion.
**(b)** Numero de marcadores bac120 unicos recuperados y ausentes en cada genoma.
Los genomas con menor cobertura del alineamiento (bin-5-50 y bin-5-63, ~57-58 %)
son tambien los de menor completitud, por lo que su asignacion debe interpretarse
con mayor cautela.

---

## Tablas

- **Tabla 1.** Genomas seleccionados en la Fase 2 (n = 21): metricas de
  ensamblaje, calidad CheckM2, taxonomia GTDB-Tk completa, referencia mas
  cercana, ANI, AF y linaje asignado.
- **Tabla 2.** Linajes definidos por desreplicacion a 95 % de ANI (n = 19),
  ordenados por prevalencia espacial, con su genoma representante y la decision
  de paso a la Fase 3.
- **Tabla 3.** Estadisticas agregadas de calidad y ensamblaje por filo.
- **Tabla S1.** Candidatos a novedad taxonomica: genomas sin especie asignada,
  con su referencia mas cercana y los valores de ANI/AF que sustentan la
  asignacion incompleta.
- **Tabla S2.** ANI y cobertura del alineamiento de las parejas comparadas por
  fastANI dentro de un mismo cluster primario (sustenta la Fig. 4).
- **Tabla S3.** Versiones de software y parametros usados en la Fase 2.
