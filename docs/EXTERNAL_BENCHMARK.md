# Benchmark externo contra anotadores

**Herramienta:** `tools/benchmark_external.py`  
**Salidas:** `reports/external_benchmark.md`, `reports/external_benchmark.json`  
**Alcance:** comparación local, sin red y sin descargas.

## Qué mide y qué no

El benchmark compara un GFF de **referencia** (gold local, por ejemplo el GFF
NCBI de `data/annotations/`) contra el GFF que produce `bakta`, `prokka` o
`dfast`. Compara:

- CDS, tRNA y rRNA;
- coordenadas 1-based inclusive, `seqid` y hebra;
- emparejamiento 1:1 por solapamiento recíproco (umbral por defecto 0,5);
- coincidencia exacta de coordenadas;
- presencia y coincidencia literal de `product`, `gene` y `locus_tag`;
- `% hypothetical` de CDS cuando el GFF trae `product`;
- tiempo de pared y RSS máximo de la ejecución externa, si se usa `--run`.

No mide función, dominios, rutas metabólicas ni expresión. Tampoco declara que
una herramienta sea superior: para eso haría falta un gold independiente y
ejecuciones pareadas.

## Detección de herramientas

La detección usa `shutil.which` para `bakta`, `prokka` y `dfast`. En el entorno
actual **no hay ninguna de las tres instalada**. Por eso el informe generado
marca las tres como `not_evaluable` con la razón:

```text
binario '<tool>' no encontrado en PATH (shutil.which)
```

No se descarga, instala ni empaqueta ninguna herramienta. Una ausencia no es un
cero ni un fracaso del método: simplemente no hay medición posible.

## Uso

Solo detección y referencia (no ejecuta terceros):

```bash
python tools/benchmark_external.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --reference-gff data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --tool all \
  --md reports/external_benchmark.md \
  --json reports/external_benchmark.json
```

Ejecución externa explícita (solo si el binario existe):

```bash
python tools/benchmark_external.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --reference-gff data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --tool prokka --run
```

Bakta necesita una base local. Se puede indicar sin descargarla:

```bash
python tools/benchmark_external.py ... --tool bakta --run \
  --bakta-db /ruta/local/a/bakta_db
```

Comandos que se usan cuando hay binario y `--run`:

| herramienta | comando base |
|---|---|
| Bakta | `bakta --output <tmp>/out --prefix result --threads 1 --force [--db <db>] <genoma>` |
| Prokka | `prokka --outdir <tmp>/out --prefix result --cpus 1 --force <genoma>` |
| DFAST | `dfast --genome <genoma> --out <tmp>/out --cpu 1` |

La ejecución ocurre en `tempfile.TemporaryDirectory`. Si falla, se guarda el
código de salida y la cola de `stderr`; no se inventan métricas de comparación.

## Métricas

Para cada categoría:

- `reference_n`, `query_n`;
- `TP`, `FP`, `FN` según el emparejamiento 1:1;
- `precision = TP/(TP+FP)`, `recall = TP/(TP+FN)`, `F1`;
- `exact_matches` y `mean_reciprocal_overlap`;
- atributos en pares emparejados;
- `% hypothetical` de referencia y de la consulta.

`% hypothetical` se calcula solo sobre CDS con atributo `product` no vacío. Se
consideran hipotéticas las cadenas que contienen `hypothetical`,
`uncharacterized`, `unknown protein`, `protein of unknown function`, `DUF...`
o `UPF...`. `putative` a secas **no** se considera hipotética.

La coincidencia literal de `product` no es accuracy funcional: dos anotadores
pueden describir la misma proteína con nombres distintos. Por eso también se
reporta el Jaccard de tokens de `product`, que sigue siendo descriptivo.

## Tiempo y RAM

Con `--run`, `execute_command` mide:

- `runtime_seconds`: tiempo de pared del proceso hijo;
- `max_rss_kb` / `max_rss_mib`: `ru_maxrss` de `RUSAGE_CHILDREN`;
- `max_rss_delta_kb`: diferencia respecto del valor previo;
- `stdout_tail`, `stderr_tail`, `returncode` y estado (`ok`, `failed`, `timeout`).

`ru_maxrss` es el máximo acumulado de procesos hijos, no un perfil instantáneo.
El límite de tiempo por herramienta es `--timeout` (1800 s por defecto).

## Resultado actual (entorno sin herramientas)

Referencia local usada: `data/annotations/GCF_000006985.1_ASM698v1.gff.gz`.

| métrica | valor |
|---|---:|
| CDS (pseudo excluidos) | 1.989 |
| CDS pseudo excluidos | 47 |
| tRNA | 50 |
| rRNA | 6 |
| CDS con `product` | 1.989 |
| CDS hipotéticas | 289 (14,5 %) |

`bakta`, `prokka` y `dfast` no están instalados; no se ejecutó nada y las tres
comparaciones quedan `not_evaluable`. No hay métricas de precisión/recall/F1 de
herramientas externas que reportar honestamente en este entorno.

## Limitaciones

1. La referencia GFF es el gold de esta comparación; si la referencia es NCBI,
   hereda sus criterios y sus errores.
2. El emparejamiento es voraz y 1:1; no es un alineamiento global ni un óptimo
   bipartito. Un cambio de umbral puede cambiar TP/FP/FN.
3. Los nombres de `product` no son comparables literalmente entre bases; la
   coincidencia de tokens es solo descriptiva.
4. No se comparan ncRNA, tmRNA, pseudogenes (se excluyen si están marcados),
   función, dominios ni rutas.
5. Las herramientas externas pueden requerir bases de datos locales y licencias
   propias; este script no las descarga.
6. Aunque el script no usa red, una herramienta de terceros podría intentar
   acceder a servicios externos. Ejecutar `--run` es decisión del usuario.
7. Sin ejecuciones pareadas y un gold independiente no se puede concluir qué
   herramienta es superior.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_benchmark_external.py -q
```

Las pruebas usan GFF sintéticos, binarios ausentes y un subproceso local de
Python para verificar tiempo/RSS. No usan red ni ejecutan bakta/prokka/dfast.
