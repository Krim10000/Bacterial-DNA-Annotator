# Benchmark externo contra anotadores

**Generado:** 2026-11-24T02:26:42Z

**Alcance:** Benchmark externo con datos locales, sin red ni descargas. La ejecución de terceros solo ocurre con --run.

**Comando de ejemplo:**

```bash
python tools/benchmark_external.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --reference-gff data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --tool all
# añadir --run para ejecutar un binario ya instalado (sin red)
```

Este informe no declara superioridad entre herramientas: sin gold independiente y ejecuciones pareadas solo se describen diferencias.

## Entorno

| herramienta | binario | disponible | PATH / razón |
|---|---|---|---|
| Bakta | `bakta` | no | binario 'bakta' no encontrado en PATH (shutil.which) |
| Prokka | `prokka` | no | binario 'prokka' no encontrado en PATH (shutil.which) |
| DFAST | `dfast` | no | binario 'dfast' no encontrado en PATH (shutil.which) |

- Herramientas disponibles: **ninguna**.
- Herramientas ausentes: **bakta, prokka, dfast**.
- `--run` solicitado: **no**.
- Umbral de solapamiento recíproco: **0,50**.

## Referencia

- GFF: `data/annotations/GCF_000006985.1_ASM698v1.gff.gz`.
- Conteos: CDS **1.989**, tRNA **50**, rRNA **6**.
- Contigs: **1**; CDS pseudo excluidos: **47**.
- CDS con product: **1.989**; hipotéticas: **289** (14,5 %).

## Resultados

| herramienta | binario | ejecución | comparación | precisión | recall | F1 | % hypothetical (query) |
|---|---|---|---|---|---|---|---|
| Bakta | no | no | no evaluable | — | — | — | — |
| Prokka | no | no | no evaluable | — | — | — | — |
| DFAST | no | no | no evaluable | — | — | — | — |

La precisión/recall/F1 se calculan contra el GFF de referencia, no contra una verdad funcional. No son una nota global del anotador.

### Bakta

- Estado: **no evaluable**.
- Razón: binario 'bakta' no encontrado en PATH (shutil.which).
- Comparación: **no evaluable** (ver razones abajo).

### Prokka

- Estado: **no evaluable**.
- Razón: binario 'prokka' no encontrado en PATH (shutil.which).
- Comparación: **no evaluable** (ver razones abajo).

### DFAST

- Estado: **no evaluable**.
- Razón: binario 'dfast' no encontrado en PATH (shutil.which).
- Comparación: **no evaluable** (ver razones abajo).

## No evaluable

- **bakta / ejecución externa**: binario 'bakta' no encontrado en PATH (shutil.which)
- **prokka / ejecución externa**: binario 'prokka' no encontrado en PATH (shutil.which)
- **dfast / ejecución externa**: binario 'dfast' no encontrado en PATH (shutil.which)

## Limitaciones

- No se descargó nada: el script solo lee datos locales y, con --run, ejecuta un binario ya instalado.
- Una herramienta ausente no es un cero ni un fracaso del método: queda como not_evaluable.
- La referencia GFF es el gold de esta comparación; no se evalúa función, rutas ni dominios.
- Los porcentajes de hipotéticas dependen de cómo cada herramienta nombra 'product'.
- No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.

## Reproducción

```bash
python tools/benchmark_external.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --reference-gff data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --tool all \
  --md reports/external_benchmark.md \
  --json reports/external_benchmark.json
```

Para ejecutar de verdad una herramienta instalada, añadir `--run`. La ejecución usa un directorio temporal; el script no descarga bases ni consulta servicios web.
