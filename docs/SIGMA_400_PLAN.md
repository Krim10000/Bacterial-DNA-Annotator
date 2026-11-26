# Extensión de sigma a todo el catálogo

## Regla

No se aplica el calibrador de una especie a otra. La extensión se hace por
niveles de evidencia.

## Nivel 0 — Inventario

Para cada genoma:

```text
genes.tsv → candidatos sigma → native_name + family + evidence
```

Resultado actual:

```text
442 genomas con genes.tsv
2.314 candidatos sigma
40 genomas sin candidatos anotados
```

Familias amplias:

```text
sigma70_like: 1.140
ecf_like:       873
sigma54_like:   301
```

Esto **no** es una asignación de promotor; es el inventario de partida.

## Nivel 1 — Evidencia de promotor

Para cada genoma:

- ventanas aguas arriba de operones;
- PWM por familia ortóloga, no por nombre global;
- fondo nulo específico de la especie;
- z, p-valor y percentil;
- arquitectura: −35/−10, extended −10, −24/−12;
- salida `candidate`, nunca `confirmed`.

Para especies sin TSS, la ventana es un **proxy de inicio de operón**; debe
declararse.

## Nivel 2 — Calibración

Solo cuando hay gold:

```text
K-12, B. subtilis: calibración específica
Mtb, PAO1, S. aureus, S. coelicolor, C. glutamicum: parcial
```

Las 400 sin gold:

- no se calibran;
- usan prior de familia;
- salida con `confidence=low` o `sin_sospecha`;
- no se comparan con accuracy.

## Nivel 3 — Generalización

Objetivo:

```text
modelo de familia entrenado en varias especies
+ ajuste de especie
+ validación leave-one-species-out
```

No:

```text
k-mer de B. subtilis aplicado a las 400
```

## Salida por genoma

```text
genome
contig
operon / TSS
sigma_family
native_name
score, z, pvalue
architecture
model_source
calibration_source
status: candidate | no_model | sin_sospecha
```

## Métricas para 400

No hay gold para las 400, así que no se reporta accuracy global.

Se reporta:

- cobertura de inventario;
- cobertura de modelos;
- fracción de candidatos con evidencia;
- abstención;
- y en las 7 especies con gold, métricas supervisadas.

## Orden de ejecución

1. Inventario masivo (hecho).
2. Detección por dominios/ortología.
3. PWM por familia.
4. Evidencia de promotor para las 400.
5. Calibración solo donde hay gold.
6. Held-out de especie.
