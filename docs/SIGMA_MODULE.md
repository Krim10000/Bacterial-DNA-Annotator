# Módulo Sigma: factores sigma y promotores sigma-dependientes

> **Estado:** v1.0.0 · módulo experimental · CPU · sin GPU · sin descargas.
>
> **Qué hace:** anota, para cada operón de un genoma, el factor sigma más
> probable de su promotor; anota también el inventario de genes sigma del
> genoma. Cada fila conserva la evidencia, la confianza y el paso del pipeline
> que la produjo.
>
> **Qué no hace:** no confirma experimentalmente nada. `sospecha_fuerte` es una
> hipótesis priorizada, no un hecho.
>
> **Agrupación entre especies:** [`SIGMA_GROUPS.md`](SIGMA_GROUPS.md).

---

## 1. Por qué existe

El sistema anterior (`12_promoters_sigma.py`) detectaba genes sigma en un
inventario, pero **no asignaba promotores de sigma alternativos**. Buscaba el
consenso σ70 y, si encontraba un par −35/−10, etiquetaba `sigma70` por defecto.
Las PWMs de σ24/σ28/σ32/σ38/σ54 ya estaban en el repositorio, pero no se usaban.

Este módulo corrige eso:

- escanea cada promotor con las PWMs de σ70/σ54/σ32/σ38/σ28/σ24;
- respeta la arquitectura de cada sigma (σ54 usa −24/−12, no −35/−10);
- calcula un p-valor empírico contra un fondo de ventanas aleatorias;
- compite entre sigmas y reporta alternativas;
- cruza con inventario, función del operón y expresión opcional;
- deja la evidencia en columnas separadas.

---

## 2. Uso

### 2.1 GUI

1. Abre **🧩 Módulos** (o usa el atajo **Más → 🧪 Anotador de sigmas**)
   para **ejecutar** el módulo.
2. Elige **Sigmas y promotores**.
3. Para **explorar** resultados ya calculados, usa
   **🧭 Explorar genoma → 🧪 Sigmas**.
3. Selecciona el genoma.
4. Pulsa **▶️ Anotar sigmas**.
5. La vista muestra una **tabla de evidencia**; selecciona una fila para ver el
   JSON completo.
6. Descarga `sigma_promoters.tsv`, `sigma_factors.tsv` o `sigma_report.md`.

### 2.2 CLI, un genoma

```bash
python scripts/sigma/10_annotate_sigma.py --genome ecoli_k12
```

### 2.3 CLI, los 400

```bash
python scripts/sigma/10_annotate_sigma.py --all --jobs 8
```

Salida agregada: `data/sigma/summary_all.tsv`. `--jobs` usa procesos
independientes; cada genoma escribe su propia carpeta.

### 2.4 Agrupar por ortólogo (global)

Después de anotar los 400:

```bash
python scripts/sigma/20_group_sigma.py --all --threads 8
```

Añade `countable`, `quality_flag`, `ortholog_group`, `copy_number` e
`is_expansion` a cada `sigma_factors.tsv`, y escribe
`data/sigma/sigma_orthologs.tsv`.

### 2.5 Autoprueba

```bash
python -m modules --selftest sigma
python scripts/sigma/10_annotate_sigma.py --selftest
python scripts/sigma/20_group_sigma.py --selftest
```

---

## 3. Entradas

| Entrada | Ruta | Uso |
|---|---|---|
| Genes | `data/regulation/<genoma>/genes.tsv` | Coordenadas, locus, producto |
| Operones | `data/regulation/<genoma>/operons.tsv` | Unidades a anotar |
| Genes del operón | `data/regulation/<genoma>/operon_genes.tsv` | Primer gen y contexto |
| Inventario sigma | `data/regulation/<genoma>/regulators_inventory.tsv` | Factores del genoma |
| PWMs sigma | `data/regulondb/gt/promoter_pwms.tsv` | σ24/28/32/38/54/70 |
| Secuencia | `genes_meta.json` → raw/gbff; fallback `data/annotations/<genoma>_genomic.gbff.gz` | Ventana del promotor |
| Expresión (opcional) | `data/expression/<genoma>.tsv` | Evidencia adicional |

### 3.1 Formato de expresión opcional

TSV con una fila por gen y columnas:

```text
gene_id   tpm
b0001     12.3
b0002     0.0
```

También acepta `locus_tag`/`gene` e `expression`/`value`. Si no existe el
archivo, el módulo no falla: escribe `expression_support=not_available`.

---

## 4. Arquitecturas escaneadas

| Sigma | Arquitectura | Espaciador | Clase |
|---|---|---|---|
| σ70 | −35/−10 | 15–19 | housekeeping |
| σ54 | −24/−12 | 10–11 | alternative |
| σ32 | −35/−10 | 12–16 | stress |
| σ38 | −35/−10 | 15–19 | stationary |
| σ28 | −35/−10 | 15–19 | flagellar |
| σ24/ECF | −35/−10 | 14–18 | ECF |

σ54 **no** se escanea con −35/−10: la fuente etiqueta las cajas como −35/−10,
pero la geometría funcional es −24/−12. El módulo lo declara en `sigma_arch`.

---

## 5. Cómo se calcula la evidencia

Para cada operón:

1. **Ventana:** `window` nt aguas arriba del primer gen, en orientación de
   transcripción.
2. **Score PWM:** suma de log-odds de la mejor pareja de cajas para cada sigma.
3. **Fondo nulo:** se puntúan `null_windows` ventanas aleatorias del genoma con
   el mismo scanner; de ahí salen `zscore` y `pvalue` empíricos.
4. **Competencia:** se ordenan los sigmas por z; `score_margin` es la diferencia
   entre el primero y el segundo.
5. **Contexto:** se buscan palabras clave de función del operón compatibles con
   el tipo de sigma (flagelo, chaperonas, nitrógeno, envelope, etc.).
6. **Expresión:** si hay tabla, se calcula la fracción de genes del operón con
   TPM > 0.
7. **Inventario:** se comprueba que el genoma tenga un factor de ese tipo y se
   reporta el más cercano.
8. **Confianza:** se combinan p-valor, z, margen, inventario, contexto y
   expresión.

---

## 6. Confianza

| Etiqueta | `confidence` | Criterio | Lectura |
|---|---|---|---|
| `sospecha_fuerte` | `high` | p ≤ 0.01, z ≥ 3, margen ≥ 2, inventario, contexto/expresión | Hipótesis prioritaria |
| `sospecha_media` | `medium` | p ≤ 0.05, z ≥ 2, margen ≥ 1, inventario | Candidato con evidencia PWM |
| `sospecha_baja` | `low` | p ≤ 0.10 o z ≥ 1 | Candidato débil |
| `sin_sospecha` | `low` | No supera umbral | Sin asignación |

`sigma_status` es independiente:

- `assigned`: confianza fuerte o media;
- `ambiguous`: hay un top pero no alcanza para asignar;
- `unresolved`: no hay candidato creíble; la columna `sigma` puede conservar
  el top exploratorio, pero `confidence_label=sin_sospecha`.

**Importante:** `sospecha_media` no es validación. Es un candidato que combina
PWM, arquitectura e inventario; el contexto funcional es lo que lo promueve a
`sospecha_fuerte`.

---

## 7. Salidas

### 7.1 `sigma_factors.tsv`

Una fila por factor sigma del inventario:

```text
genome, gene_id, locus_tag, gene_name, product, sigma_type, sigma_class,
status, confidence, confidence_label, confidence_score, evidence,
detected_by, pipeline_step, promoter_hits, source_status, source_confidence,
source_evidence, factor_note, countable, quality_flag, ortholog_group,
copy_number, is_expansion, provenance
```

### 7.2 `sigma_promoters.tsv`

Una fila por operón:

```text
genome, operon_id, operon_label, contig, strand, first_gene, first_gene_name,
sigma, sigma_class, sigma_arch, spacer, score, zscore, pvalue, score_margin,
dist_to_start, upstream_box, downstream_box, sigma_status, confidence,
confidence_label, confidence_score, status, evidence, detected_by,
pipeline_step, alternatives, pathway_support, pathway_hits, expression_support,
expression_detail, inventory_support, inventory_detail, location_support,
provenance
```

Campos clave para auditar:

| Campo | Qué responde |
|---|---|
| `evidence` | Resumen legible de la evidencia |
| `detected_by` | Pasos que la detectaron (`promoter_pwm`, `architecture`, `inventory`, `pathway`, `expression`, `location`) |
| `pipeline_step` | Versión del paso (`sigma.scan.v1`, `sigma.inventory.v1`) |
| `alternatives` | Top 3 sigmas con z, p-valor y spacer |
| `pathway_hits` | Genes/palabras que apoyan el contexto |
| `inventory_detail` | Factor sigma del genoma y distancia |
| `provenance` | Genoma, operón, ventana, semilla y versión |

### 7.3 `sigma_summary.json`

Conteos por tipo y confianza, cobertura, parámetros, `sequence_available` y
hash SHA-256. Si un genoma no tiene secuencia local, el módulo no falla:
reporta `sequence_available=false` y deja los operones como `unresolved`.

### 7.4 `sigma_report.md`

Informe legible con top promotores y límites.

---

## 8. Ejemplo de fila

```text
operon_label: otsBA
first_gene_name: otsB
sigma: sigma38
sigma_arch: -35/-10
spacer: 18
zscore: 4.388
pvalue: 0.00398
score_margin: 2.801
confidence_label: sospecha_fuerte
evidence: PWM sigma38 z=4.388 p=0.00398 spacer=18; margen=2.801;
          inventario=rpoS; contexto=sí; expresión=not_available
detected_by: promoter_pwm:sigma38;architecture:-35/-10;spacer:18;
             inventory:sigma38;pathway:sigma38;location:upstream
pipeline_step: sigma.scan.v1
```

---

## 8.1 Resultado del lote de los 400

Corrida del 2026-11-25: **401 genomas** (400 del proyecto + K-12), 8 procesos,
200 ventanas nulas, ventana 180 nt.

| Métrica | Valor |
|---|---:|
| Operones escaneados | 786.805 |
| Asignados (`assigned`) | 73.202 (9,3 %) |
| `sospecha_fuerte` | 752 (0,1 %) |
| `sospecha_media` | 72.450 |
| `sospecha_baja` | 423.310 |
| `sin_sospecha` | 289.877 |
| Genomas con secuencia | 401/401 |
| Genomas con 0 asignados | 10 |
| Factores contables | 2.135 de 5.993 |
| Factores `sigma_unresolved` | 3.765 |
| Grupos ortólogos | 2.432 (663 contables) |
| Expansiones detectadas | 213 |

Distribución de los promotores con candidato:

| Sigma | Promotores |
|---|---:|
| σ70 | 253.201 |
| σ38 | 161.962 |
| σ24 | 118.603 |
| σ28 | 87.254 |
| σ32 | 84.913 |
| σ54 | 80.456 |

**Lectura:** el módulo es conservador. La mayoría de los operones queda en
`sospecha_baja` o `sin_sospecha`; solo el 9,3 % alcanza asignación media/fuerte.
Eso es intencional: el PWM solo no basta para distinguir σ70 de σ38/σ24. Los
10 genomas con 0 asignados conservan sus candidatos ambiguos y su evidencia.

En el inventario, **2.135 de 5.993** factores son contables; **3.765** son
`sigma_unresolved`. Los 2.432 grupos ortólogos (663 contables) y las 213 expansiones permiten
comparar especies sin confundir genes, grupos y fragmentos. Ver
[`SIGMA_GROUPS.md`](SIGMA_GROUPS.md).

---

## 8.2 Validación K-12 y B. subtilis

> Evaluación fina por capas: [`SIGMA_LAYERED_EVAL.md`](SIGMA_LAYERED_EVAL.md).
> Combinaciones de features: [`SIGMA_FEATURE_COMBINATIONS.md`](SIGMA_FEATURE_COMBINATIONS.md).
> Rutas metabólicas como prior: [`SIGMA_PATHWAY_PRIOR.md`](SIGMA_PATHWAY_PRIOR.md).

Resultado real del benchmark TSS-level ([`SIGMA_BENCHMARK.md`](SIGMA_BENCHMARK.md)):

| Métrica | K-12 | B. subtilis |
|---|---:|---:|
| TSS evaluables | 3.189 | 2.226 |
| Top-1 | **57,0 %** | **44,8 %** |
| Top-3 cualquiera | 90,2 % | 85,5 % |
| Mayoría trivial | 62,6 % (σ70) | 87,6 % (A→σ70) |
| TSS multi-sigma | 124 | 0 |

La confianza ordena algo mejor que el azar, pero no es una probabilidad
calibrada:

| Confianza | Acierto K-12 |
|---|---:|
| `sospecha_media` | 65,83 % |
| `sospecha_baja` | 51,25 % |
| `sin_sospecha` | 36,15 % |

**Lectura honesta:** el top-1 no supera la mayoría trivial; el top-3 está
inflado por σ70. El módulo sirve para **priorizar y explorar**, no para
confirmar. El gold ya trae 124 TSS de K-12 con **más de un sigma**, así que un
top-1 único no puede acertarlos todos. La comparación por operón expone la
limitación de una fila por operón.

---

## 9. Límites

1. **PWM no es verdad.** Un score alto indica parecido al motivo, no actividad.
2. **σ70 y σ38/σ24 se confunden** en promotores degenerados. Por eso existe
   `ambiguous` y `alternatives`.
3. **La precisión por sigma no es uniforme.** En K-12, σ70 es el más frecuente;
   σ38/σ24 tienen más falsos positivos. El módulo lo declara en la confianza.
4. **El inventario puede tener falsos positivos.** Ejemplo: `rseD` (péptido
   líder de `rpoE`) aparece como `sigmaE` con estado `novel`.
5. **σ24/ECF y σ19 no siempre están en el inventario**; el escaneo puede
   encontrar el promotor sin que el gen esté anotado como sigma.
6. **La expresión es opcional.** Sin tabla, `expression_support=not_available`.
7. **No hay RNA-seq en el repositorio.** El contexto metabólico es un prior, no
   una medida de expresión.
8. **No es confirmación experimental.** Para eso hace falta TSS, EMSA, ITC o
   reportero.

---

## 10. Criterios de aceptación

- [ ] Ningún promotor se etiqueta `sigma70` por defecto.
- [ ] σ54 se escanea con −24/−12.
- [ ] Cada fila tiene `evidence`, `detected_by`, `pipeline_step` y
      `provenance`.
- [ ] Cada fila tiene `confidence_label` y `sigma_status`.
- [ ] Las alternativas (top 3) están en la salida.
- [ ] El módulo aparece en **🧩 Módulos** y se ejecuta desde la GUI.
- [ ] La GUI muestra tabla de evidencia y fila completa.
- [ ] `--all` genera `summary_all.tsv`.
- [ ] La autoprueba pasa en CI.
- [ ] Los límites están en el informe.

---

## 11. Siguiente etapa

Con los sigmas marcados, la siguiente etapa puede:

1. cruzar promotores sigma con regulones TF→operón;
2. priorizar TSS candidatos para RNA-seq/EMSA;
3. construir la red sigma → operón → función;
4. comparar la distribución de sigmas entre los 400 genomas;
5. elegir los casos de `sospecha_fuerte` para validación experimental.

---

## 12. Frase final

**El módulo no “decide” el sigma: enumera promotores, puntúa cada arquitectura,
compite entre sigmas y publica la evidencia con confianza. Donde no alcanza,
dice `ambiguous` o `unresolved` en vez de inventar.**
