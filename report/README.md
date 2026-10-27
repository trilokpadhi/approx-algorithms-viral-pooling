# LaTeX progress report

This is the student progress report in the official NeurIPS 2025 LaTeX style,
using the `preprint` option. It is a course progress report, not a conference
submission or an accepted NeurIPS paper.

## Overleaf

Upload `viral-pooling-progress-report.zip` through **New project > Upload project**.
Select `main.tex` as the main document and use pdfLaTeX. No external files,
BibTeX run, shell escape, or data downloads are needed.

- `main.tex`: title, packages, and main document.
- `content.tex`: editable report text, tables, and references.
- `neurips_2025.sty`: unmodified official style file.
- `progress-report.pdf`: compiled, visually checked five-page report.

The ZIP contains the three source files and this README. The report focuses on the public HCV dataset, including the methodology,
algorithm, specimen-level results, limitations, and planned work. It retains
the Markdown report's student voice.
The author block lists Trilok Padhi (Panther ID: 02716700) and Mark Grinshpon
(Panther ID: TBD).

## Local build

With a LaTeX installation:

```sh
latexmk -pdf main.tex
```

The style was obtained from the [official NeurIPS 2025 style archive](https://media.neurips.cc/Conferences/NeurIPS2025/Styles.zip),
linked from the [conference formatting instructions](https://neurips.cc/Conferences/2025/CallForPapers).
The original style-file notices are preserved.
