# Benchmark externo CPU

Runner reproducible para comparar CDS, tRNA y rRNA de **Bakta**, **Prokka** y
**DFAST** contra una referencia GFF3 o GenBank/GBFF.

- `external_benchmark.py`: parser, comparador, ejecutor y CLI.
- `Dockerfile`: imagen CPU con las tres herramientas y versiones fijadas.
- `docker-compose.yml`: ejecución con el repositorio montado en `/work`.
- `prepare_tools.sh`: alternativa local con entornos micromamba separados.

La guía completa, métricas, estados y límites está en
[`docs/EXTERNAL_BENCHMARK_RUN.md`](../../docs/EXTERNAL_BENCHMARK_RUN.md).

Pruebas (solo datos sintéticos, sin red):

```bash
.venv/bin/python -m pytest tests/test_external_benchmark.py -q
```
