# Tesis: Diseño de un modelo predictivo de riesgo de deterioro de la Estela de Raimondi basado en genómica comparativa y potencial metabólico

Este repositorio contiene los archivos fuente en **LaTeX** correspondientes a la tesis de licenciatura en **Bioingeniería** (UTEC).

---

## 👥 Autores y Asesores

- **Autores:**
  - Eduardo Esteban Tello Osorio
  - Ana Sofía Barrientos Uriarte
- **Asesores:**
  - Pablo Tsukayama Cisneros
  - Monica Cecilia Santa María Fuster
  - Harry Gustavo Saavedra Espinoza

---

## 📁 Estructura del Proyecto

```text
├── main.tex                 # Archivo raíz de LaTeX a compilar
├── main.pdf                 # Documento final compilado (versión actual)
├── tesisutec.cls            # Clase de estilo de tesis UTEC
├── IEEEtran.bst             # Estilo de citas IEEE
├── referencias.bib          # Base de datos de citas bibliográficas
├── compilar.bat             # Script de Windows para compilar con un clic
├── secciones/               # Capítulos individuales (resumen, metodología, etc.)
├── encabezados/             # Dedicatoria, agradecimientos
└── images/                  # Figuras y diagramas de la tesis
```

---

## ⚙️ Cómo compilar el documento

### Opción 1: Con Tectonic (Recomendado - Sin límites de tiempo)
1. Descarga el ejecutable `tectonic.exe` desde [Tectonic Releases](https://github.com/tectonic-typesetting/tectonic/releases) y colócalo en esta carpeta.
2. Haz doble clic en `compilar.bat` o ejecuta en la consola:
   ```bash
   tectonic main.tex
   ```

### Opción 2: Con TeX Live / MiKTeX / pdflatex
Puedes compilar utilizando la secuencia estándar:
```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```
o con `latexmk`:
```bash
latexmk -pdf main.tex
```
