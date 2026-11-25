# Benchmark externo CPU

Runner reproducible para comparar CDS, tRNA y rRNA de **Bakta**, **Prokka** y
**DFAST** contra una referencia GFF3 o GenBank/GBFF.

- `external_benchmark.py`: detector (incluida Bakta DB), ejecutor, comparador,
  validador y CLI. Sin binarios/base deja `not_available` y no inventa métricas.
- `make_gold.py`: construye un gold pequeño desde un GFF3/GBFF local y un
  comando `--run` reproducible por genoma.
- `gold/chlorobaculum_tepidum_60kb/`: gold real de 60 kb (55 CDS + 1 tRNA) con
  `genome.fna`, `reference.gff3`, `gold_manifest.json` y `run_benchmark.sh`.
- `Dockerfile`: imagen CPU con las tres herramientas y versiones fijadas.
- `docker-compose.yml`: ejecución con el repositorio montado en `/work`.
- `prepare_tools.sh`: alternativa local con entornos micromamba separados.

La guía completa, estados, Bakta DB, validación, métricas y límites está en
[`docs/EXTERNAL_BENCHMARK_RUN.md`](../../docs/EXTERNAL_BENCHMARK_RUN.md).

Estado real: en el entorno de desarrollo no hay `bakta`, `prokka`, `dfast` ni
Bakta DB. El informe se genera con las tres filas en `not_available`.

Pruebas (solo datos sintéticos, sin red):

```bash
.venv/bin/python -m pytest tests/test_external_benchmark_real.py -q
.venv/bin/python -m pytest tests/test_external_benchmark.py -q
```
