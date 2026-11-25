# Sigma: por qué cambia el número entre especies y cómo agruparlo

> **Pregunta:** veo que distintas especies tienen distinta cantidad y tipos de
> sigma. ¿Nuestro sistema los anota/agrupa bien?
>
> **Respuesta corta:** el sistema agrupa por **tipo/familia** (`sigma70`,
> `sigma54`, `sigmaE`, …), pero **no por ortólogo**. Los conteos son de genes,
> no de grupos ortólogos. Parte de la variación es biología real; parte es
> artefacto de anotación, sobre todo el cajón `sigma_unresolved`.

---

## 1. Lo que se observa en el lote

401 genomas (400 + K-12), 5.993 genes sigma anotados por el inventario:

| Estadística | Valor |
|---|---:|
| Mediana de sigma por genoma | 7 |
| Percentil 90 | 38 |
| Percentil 95 | 52 |
| Máximo | 157 |
| `sigma_unresolved` | 3.765 (62,8 %) |
| Genomas con >10 `sigma_unresolved` | 88 |
| Genomas con >50 `sigma_unresolved` | 11 |
| Genomas con >100 `sigma_unresolved` | 3 |

Top de genomas con más factores:

| Genoma | Total | `sigma_unresolved` | Tipos presentes |
|---|---:|---:|---|
| GCF_901538265.1 | 157 | 134 | σE 13, σ70 5, σ54 1, σB 1, σF 1, σH 1, σK 1 |
| GCF_013128195.2 | 156 | 131 | σE 16, σ70 4, σH 2, σ54 1, σB 1, σF 1 |
| GCF_025148635.1 | 106 | 102 | σ70 2, σ54 1, σE 1 |
| GCF_008087625.1 | 133 | 86 | σE 33, σ70 11, σ38 2, σ54 1 |
| GCF_026965555.1 | 113 | 78 | σE 21, σ70 5, σ54 2, … |

**Lectura:** los genomas con 100–157 sigma están dominados por
`sigma_unresolved`, no por familias bien asignadas. Eso no es un censo
biológico: es una mezcla de ECF reales, fragmentos y hits parciales.

---

## 2. Por qué cambia de verdad entre especies

La cantidad de sigma **sí** varía biológicamente:

- *E. coli*: ~7 sigma (1 housekeeping + alternativos).
- *Bacillus subtilis*: ~18, incluyendo sigma de esporulación.
- Actinobacterias (*Streptomyces*) y mixobacterias: decenas de sigma ECF.
- Genomas pequeños, simbiontes o parásitos: pocos sigma, a veces uno solo.

Las expansiones de ECF son reales en algunos linajes. Por eso no se puede
asumir que “más sigma = error”.

Pero también hay artefactos:

1. **Hits parciales de dominio.** `Sigma70_r2` solo, sin `r4`, sin producto y
   sin nombre, es un candidato débil, no un factor sigma confirmado.
2. **Fragmentos de proteína.** Genes truncados o mal llamados.
3. **Duplicados de nombre.** Varias copias de `sigA`/`sigE` pueden ser parálogos
   reales o anotaciones colapsadas; el nombre no es un grupo ortólogo.
4. **Panel HMM incompleto.** No cubre todas las familias; lo que no cubre cae en
   `sigma_unresolved` o se pierde.
5. **Anotación de producto.** Un “sigma factor” genérico no dice familia.

---

## 3. Cómo agrupa el sistema hoy

`sigma_factors.tsv` tiene:

| Campo | Qué agrupa |
|---|---|
| `sigma_type` | Familia asignada (`sigma70`, `sigma54`, `sigmaE`, …) |
| `sigma_class` | Clase (`housekeeping`, `alternative`, `ecf`, `sporulation`, `unknown`) |
| `status` | `validated`, `consensus`, `novel`, `unresolved` |
| `confidence` | `high`, `medium`, `low` |

Lo que **no** tiene:

- `ortholog_group`: no hay clustering de proteínas entre genomas.
- `copy_number`: no se distingue una copia única de una expansión de parálogos.
- `countable`: no se separa un factor completo de un fragmento.
- `quality_flag`: no se marca el motivo del dudoso.

Por eso dos genomas con el mismo número de `sigma_unresolved` pueden tener
biologías muy distintas.

---

## 4. Agrupación propuesta: tres niveles

### Nivel 1 — Clase

```text
housekeeping | alternative | ecf | sporulation | unknown
```

Sirve para comparar genomas sin mezclar σ70 con σE o σK.

### Nivel 2 — Familia

```text
sigma70 | sigma54 | sigma32 | sigma38 | sigma28 | sigma24 | sigmaE | sigmaB | ...
```

Es lo que ya existe. No es un ortólogo; es una familia.

### Nivel 3 — Grupo ortólogo

Clusterizar las proteínas sigma de los 400 genomas (DIAMOND/MMseqs2) y asignar:

```text
ortholog_group   # ID estable del grupo
copy_number      # copias de ese grupo en el genoma
is_expansion     # copy_number > 1
```

Así se distingue:

- una expansión real de ECF (varias copias divergentes del mismo grupo),
- una duplicación reciente (copias casi idénticas),
- un artefacto de anotación (fragmentos con el mismo dominio).

### Nivel 4 — Calidad

```text
countable        # true/false
quality_flag     # validated | consensus | partial_domain | no_product | fragment
```

Regla propuesta:

- `countable=true` si tiene nombre/producto y ≥2 dominios, o si está
  `validated`.
- `countable=false` si es `novel` con un solo dominio, sin producto o truncado.
- `sigma_unresolved` no se cuenta como familia; se reporta aparte.

---

## 4.1 Resultado del agrupamiento

Corrida del 2026-11-25 sobre 401 genomas y 5.993 factores:

| Métrica | Valor |
|---|---:|
| Factores totales | 5.993 |
| Factores contables | **2.135** |
| `sigma_unresolved` | **3.765** |
| Grupos ortólogos | **2.432** (663 contables) |
| Expansiones (grupo contable con >1 copia) | **213** |
| Factores con secuencia | 5.993/5.993 |

Los genomas con 157/156 factores bajan a **23/25 contables**, respectivamente.
Eso confirma que el conteo bruto estaba inflado por `sigma_unresolved`.

## 5. Qué cambia en los conteos

Con `countable` y `ortholog_group`, el reporte por genoma debería tener:

| Columna | Ejemplo |
|---|---|
| `n_sigma_genes` | 157 |
| `n_sigma_countable` | 23 |
| `n_sigma_unresolved` | 134 |
| `n_ortholog_groups` | 18 |
| `n_expansions` | 3 |
| `n_fragments` | 5 |

Eso responde “¿cuántos sigma tiene esta especie?” sin mezclar genes, grupos y
basura.

---

## 6. Qué implementar

| Tarea | Salida |
|---|---|
| Añadir `countable`/`quality_flag` al inventario | `sigma_factors.tsv` |
| Clusterizar proteínas sigma de los 400 | `sigma_orthologs.tsv` |
| Añadir `ortholog_group` y `copy_number` | `sigma_factors.tsv` |
| Resumen por genoma con `n_countable` y `n_groups` | `summary_all.tsv` |
| Reporte por filo/clase | `sigma_groups_report.md` |

El clustering se puede hacer con DIAMOND (`diamond cluster`/`linclust`) sobre
las proteínas de `data/annotations_ab_initio/<genoma>/proteins.faa`, sin GPU.

---

## 7. Respuesta en una frase

**Sí, el sistema agrupa por familia y clase, pero no por ortólogo. Las
diferencias entre especies son en parte biología real (expansiones de ECF,
esporulación) y en parte artefacto (`sigma_unresolved` con dominios parciales).
Para comparar especies hay que añadir `countable` y `ortholog_group`, y reportar
genes, grupos y copias por separado.**
