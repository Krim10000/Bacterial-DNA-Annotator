# Demo incluida

Datos generados con `python tools/build_demo_data.py` a partir de
`webapp/samples/demo.fna` (150 kb de *Chlorobaculum tepidum*).

- `genes.tsv`: predicción ab initio con pyrodigal.
- `operons.tsv`: agrupación **heurística** por proximidad (≤120 pb) y
  misma hebra. No es evidencia experimental; se marca `demo`/`low`.
- `promoters.tsv` y `hairpins.tsv`: vacíos; la demo no inventa esas
  pistas.

La demo sirve para probar el navegador y las fichas al hacer clic. No
se usa para las métricas del proyecto ni para conclusiones biológicas.
