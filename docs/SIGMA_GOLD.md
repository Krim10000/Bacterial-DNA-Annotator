# Golden set de sigma (TSS-level)

> **Estado:** snapshot congelado y offline. Fecha del snapshot: **2026-11-25**.
> **Pregunta que responde:** ¿cuál es la etiqueta de factor sigma curada para
> cada TSS de *E. coli* K-12 y *B. subtilis* 168, cuando la fuente y el TSS
> tienen coordenada válida y el sigma tiene PWM equivalente?
>
> **Respuesta corta:** **5.540** pares TSS–sigma: **3.314** de K-12
> (RegulonDB) y **2.226** de *B. subtilis* (**BSGatlas v1.1**, no DBTBS puro).
> Corrección 2026-11-26: de las 2.226 filas, **1.777 (79,8 %)** solo citan
> Nicolas et al. (tiling array), **379 (17,0 %)** mencionan DBTBS y solo
> **430** cumplen DBTBS/BsubCyc + `Resolution: 0`. No es un benchmark de los
> 6 mecanismos de [`SIGMA_CASES.md`](SIGMA_CASES.md), no tiene conjunto
> negativo y no reemplaza las predicciones de
> [`SIGMA_BENCHMARK.md`](SIGMA_BENCHMARK.md).

---

## 1. Qué es el gold

Un conjunto de referencia **TSS-level** con una fila por TSS, hebra, sigma y
genoma. Sirve para:

- evaluar un detector de sigma contra etiquetas curadas;
- comparar predicciones top-1/top-k por tipo de sigma;
- congelar un denominador reproducible para reportes y tests;
- auditar la procedencia de cada etiqueta (ID de promotor, evidencia, PMIDs).

Cada fila tiene los 13 campos del esquema
[`data/sigma/gold/schema.json`](../data/sigma/gold/schema.json):

| Columna | Tipo | Nulo | K-12 | *B. subtilis* |
|---|---|---|---|---|
| `species` | string | no | `Escherichia coli K-12 MG1655` | `Bacillus subtilis 168` |
| `genome_id` | string | no | `ecoli_k12` | `subtilis` |
| `accession` | string | no | `GCF_000005845.2` | `GCF_000009045.1` |
| `contig` | string | no | `NC_000913.3` | `basu168` |
| `pos` | int ≥ 1 | no | `PromoterSet.tsv:4)posTSS` | `gt_tss.tsv:pos` |
| `strand` | `+` / `-` | no | normalizada desde `forward`/`reverse` | tal cual |
| `sigma` | enum de 6 | no | `sigmaFactor` sin cambios | letra DBTBS mapeada |
| `gene` | string | **sí** | `firstGeneName` | vacío: la fuente no trae gen |
| `evidence` | string | no | `promoterEvidence` | `evidence` completo |
| `pmid` | string | **sí** | `pmids`, normalizado a `;` | extraído de `evidence` |
| `source` | enum | no | `RegulonDB` | `DBTBS` |
| `confidence` | enum | no | `S`/`W`/`C` o `curated` | `curated` |
| `notes` | string | **sí** | IDs, nombres, cajas, `addEvidence` | `dbtbs_label` y `Resolution` |

La coordenada es **1-based** y `pos` es la primera base transcrita. La clave
primaria es `(genome_id, contig, pos, strand, sigma)`.

### Archivos

| Archivo | Filas | Rol |
|---|---:|---|
| [`data/sigma/gold/schema.json`](../data/sigma/gold/schema.json) | — | contrato, enums y procedencia |
| [`data/sigma/gold/ecoli_k12/tss_sigma.tsv`](../data/sigma/gold/ecoli_k12/tss_sigma.tsv) | 3.314 | gold de K-12 |
| [`data/sigma/gold/subtilis/tss_sigma.tsv`](../data/sigma/gold/subtilis/tss_sigma.tsv) | 2.226 | gold de *B. subtilis* |
| [`scripts/sigma/40_validate_gold.py`](../scripts/sigma/40_validate_gold.py) | — | validador offline |
| [`tests/test_sigma_gold.py`](../tests/test_sigma_gold.py) | — | esquema, conteos y procedencia |

---

## 2. Qué **no** es el gold

- **No es un benchmark de mecanismos.** No tiene columna `case` ni
  `mechanism`. No distingue activador clase I/II, MerR, represor estérico,
  looping ni anti-activación. Para eso ver
  [`SIGMA_CASES.md`](SIGMA_CASES.md).
- **No es un conjunto negativo.** La ausencia de una fila no significa que no
  haya promotor; solo que la fuente no tiene una etiqueta evaluable.
- **No es el universo de TSS.** Solo incluye TSS con sigma etiquetado,
  coordenada válida y PWM equivalente.
- **No es una re-anotación.** Copia la etiqueta de la fuente; no la corrige ni
  la infiere por homología de secuencia.
- **No es una probabilidad.** `confidence` es el nivel de la fuente o
  `curated`; no está calibrado.
- **No es una tabla de genes.** `gene` puede estar vacío y no se imputa desde
  otras tablas.
- **No es la salida del modelo.** Es referencia para evaluar salidas como
  [`data/sigma/benchmark/`](../data/sigma/benchmark/).

---

## 3. Conteos

### 3.1 Embudo

| Genoma | Fuente | Filas fuente | Sigma con PWM (pre-coord.) | Coord. inválida | Duplicado TSS | **Gold** |
|---|---|---:|---:|---:|---:|---:|
| K-12 | RegulonDB | 4.057 | 3.340 | 25 | 1 | **3.314** |
| *B. subtilis* | BSGatlas v1.1 | 3.390 | 2.226 | 0 | 0 | **2.226** (430 estrictos DBTBS) |
| **Total** | | **7.447** | **5.566** | **25** | **1** | **5.540** |

K-12: de 4.057 filas, 716 no tienen `sigmaFactor`, 1 tiene `sigma19` (sin PWM),
25 tienen sigma permitido pero `posTSS` vacío/0 y 1 es un duplicado
TSS+hebra+sigma. La suma de descartes es 743: 716 + 1 + 25 + 1 duplicado
colapsado.

*B. subtilis*: de 3.390 filas, 316 tienen etiqueta `?` y 848 tienen etiquetas
sin equivalente PWM (`E`, `F`, `G`, `H`, `K`, `W`, `X`, `M`, `YlaC`). Total
excluido: 1.164.

### 3.2 Por sigma

| Sigma | K-12 | *B. subtilis* | Total |
|---|---:|---:|---:|
| `sigma70` | 1.997 | 1.949 | 3.946 |
| `sigma24` | 522 | 0 | 522 |
| `sigma32` | 312 | 13 | 325 |
| `sigma38` | 241 | 182 | 423 |
| `sigma28` | 143 | 71 | 214 |
| `sigma54` | 99 | 11 | 110 |
| **Total** | **3.314** | **2.226** | **5.540** |

En K-12 hay **124 posiciones TSS/hebra con más de un sigma** (125 filas
extra): la misma coordenada puede tener, por ejemplo, `sigma70` y `sigma24`.
En *B. subtilis* no hay posiciones TSS/hebra repetidas en el gold.

### 3.3 Confianza, evidencia y PMIDs

| Campo | K-12 | *B. subtilis* | Total |
|---|---:|---:|---:|
| `confidence=S` | 928 | 0 | 928 |
| `confidence=W` | 2.307 | 0 | 2.307 |
| `confidence=C` | 73 | 0 | 73 |
| `confidence=curated` | 6 | 2.226 | 2.232 |
| Con `pmid` | 3.314 (100,0%) | 429 (19,3%) | 3.743 (67,6%) |
| Con `gene` | 1.545 (46,6%) | 0 (0,0%) | 1.545 (27,9%) |
| Con `notes` | 3.314 | 2.226 | 5.540 |

`S`/`W`/`C` son los valores de `confidenceLevel` de RegulonDB copiados
literalmente; **no se reinterpretan**. `curated` significa que la fuente no
aporta confianza (todas las filas de *B. subtilis* y 6 filas de K-12 cuyo
`confidenceLevel` era `None`).

---

## 4. Derivación

### 4.1 *E. coli* K-12

Fuente primaria:
[`data/regulondb/PromoterSet.tsv`](../data/regulondb/PromoterSet.tsv).
Contraste:
[`data/regulondb/gt/gt_promoters.tsv`](../data/regulondb/gt/gt_promoters.tsv).

Reglas:

1. Se conserva `5)sigmaFactor` en `sigma` si está en el set de seis.
2. Se conserva `4)posTSS` en `pos`; se exige entero > 0.
3. Se normaliza `3)strand`: `forward → +`, `reverse → -`.
4. Se conserva `7)firstGeneName` en `gene` (puede quedar vacío).
5. Se conserva `13)promoterEvidence` en `evidence`.
6. Se conserva `16)pmids` en `pmid`, separado por `;`.
7. `confidence` = `15)confidenceLevel` si es `S`, `W` o `C`; si no,
   `curated`.
8. `notes` guarda `regulondb_id`, `name`, `dist_to_first_gene`, cajas
   `minus10`/`minus35` y `addEvidence` cuando existen.
9. Cada `promoter_id` se contrasta con `gt_promoters.tsv`; si `tss`, `strand`,
   `sigma_factor` o `first_gene` no coinciden, la generación falla.

Duplicado colapsado: `RDBECOLIPMC02190` (`aspTp1`) y `RDBECOLIPMC02191`
(`trpTp2`) comparten TSS 3.946.780, hebra `+` y `sigma28`. Se conserva una
fila y ambos IDs/nombres quedan en `notes`. El gold es TSS-level, no
promoter-level.

### 4.2 *Bacillus subtilis* 168

Fuente:
[`data/dbtbs/gt_tss.tsv`](../data/dbtbs/gt_tss.tsv).

Mapeo explícito de etiquetas DBTBS:

| Etiqueta DBTBS | Sigma normalizado |
|---|---|
| `A` | `sigma70` |
| `B` | `sigma38` |
| `D` | `sigma28` |
| `I` | `sigma32` |
| `L` | `sigma54` |

Reglas:

1. Se conserva `pos`; se exige entero > 0.
2. Se conserva `strand` (`+`/`-`).
3. Se mapea `sigma` solo para `A`, `B`, `D`, `I`, `L`.
4. Se conserva `evidence` completo, sin cambios.
5. `pmid` se extrae de `evidence` con `PubMed:` y se normaliza a `;`.
   Si no hay números, queda vacío.
6. `gene` queda vacío: `gt_tss.tsv` no trae gen.
7. `confidence` = `curated`: la fuente no trae confianza.
8. `notes` guarda `dbtbs_label` y `Resolution` cuando existe.

### 4.3 Exclusiones explícitas

No se inventa sigma. Quedan fuera del gold evaluable:

- K-12: 716 filas sin etiqueta, 1 fila `sigma19` (sin PWM), 25 filas con
  coordenada inválida, 1 duplicado TSS-level.
- *B. subtilis*: 316 filas con `?` y 848 filas con etiquetas sin PWM
  equivalente (`E`, `F`, `G`, `H`, `K`, `W`, `X`, `M`, `YlaC`).

No se convierten en negativos: son **no evaluables**, no ausencias.

---

## 5. Diferencias con `SIGMA_CASES`

| | `SIGMA_CASES.md` | Este gold |
|---|---|---|
| Unidad | mecanismo regulatorio | TSS–sigma |
| Etiquetas | 6 casos (clase I/II, MerR, estérico, looping, anti-activación) | 6 tipos de sigma |
| ¿Tiene columna de caso? | no (el documento propone construirla) | no |
| ¿Mide mecanismos? | no | no |
| ¿Sirve para evaluar sigma? | no | sí |
| ¿Tiene negativos? | no | no |

`SIGMA_CASES.md` concluye que **no hay rendimiento por caso** porque faltan
etiquetas de mecanismo y detectores. Este gold **no cambia esa conclusión**:
solo aporta un denominador TSS-level para el **tipo de sigma**.

Tampoco reemplaza [`SIGMA_BENCHMARK.md`](SIGMA_BENCHMARK.md): ese documento
describe predicciones y métricas; el gold es la referencia estática que puede
usarse para recalcularlas o para evaluar otro detector.

---

## 6. Cómo usarlo

### 6.1 Cargar y validar

```python
import csv
from pathlib import Path

gold = Path("data/sigma/gold")
rows = []
for genome in ("ecoli_k12", "subtilis"):
    with (gold / genome / "tss_sigma.tsv").open(newline="", encoding="utf-8") as fh:
        rows.extend(csv.DictReader(fh, delimiter="\t"))
```

Antes de usarlo:

```bash
python scripts/sigma/40_validate_gold.py
```

Sale con código 0 solo si no hay errores. El resumen incluye filas, genomas,
sigma, fuente, confianza, cobertura de PMID/gen y duplicados.

El consumidor existente
[`scripts/sigma/50_eval_layers.py`](../scripts/sigma/50_eval_layers.py) ya
espera este layout (`data/sigma/gold/<genoma>/tss_sigma.tsv`), acepta estas
columnas y agrupa los multi-sigma por TSS. Se puede usar con
`--gold-root data/sigma/gold --all`.

### 6.2 Unir predicciones

La clave natural de un TSS es `(genome_id, contig, pos, strand)`. La clave
completa del gold es `(genome_id, contig, pos, strand, sigma)`.

- Si el detector predice **un sigma por TSS**, unir por
  `(genome_id, contig, pos, strand)` y comparar `sigma`.
- Si predice **varios sigma por TSS**, conservar todas las filas del gold para
  esa coordenada: K-12 tiene 124 posiciones multi-sigma.
- No unir por `gene`: en *B. subtilis* está vacío y en K-12 falta en 1.769
  filas.
- No usar `notes` como etiqueta: es procedencia.

### 6.3 Métricas sugeridas

- accuracy top-1 por TSS y por sigma;
- precision/recall/F1 por sigma (macro y micro);
- top-k cuando el detector ordene varios sigma;
- cobertura: cuántos TSS del gold tienen predicción;
- no mezclar `sigma70` con las alternativas sin reportar el desbalance: el
  71,2% del gold es `sigma70` (3.946/5.540).

### 6.4 Qué no hacer

- No entrenar con el gold y reportar métricas sobre él.
- No usar el gold como conjunto negativo.
- No interpretar `confidence` como probabilidad.
- No imputar `gene` desde otra tabla y presentarlo como gold.
- No reportar performance por caso usando este archivo.

---

## 7. Validación

El validador
[`scripts/sigma/40_validate_gold.py`](../scripts/sigma/40_validate_gold.py)
es **offline** (solo stdlib, sin red) y comprueba:

- cabecera exacta contra `schema.json`;
- columnas obligatorias no vacías;
- `pos` entero > 0;
- `strand` en `+`/`-`;
- `sigma` en el enum de seis;
- `genome_id`, `species`, `accession`, `contig` y `source` coherentes con el
  genoma del archivo;
- duplicados por clave primaria (dentro y entre archivos);
- PMIDs numéricos y separados por `;`/`,`;
- resumen por genoma, sigma, fuente, confianza y cobertura.

Comandos:

```bash
python scripts/sigma/40_validate_gold.py
python scripts/sigma/40_validate_gold.py --json
.venv/bin/python -m pytest tests/test_sigma_gold.py -q
```

Salida esperada del validador (resumen):

```text
filas: 5540
por genoma: ecoli_k12=3314, subtilis=2226
por sigma: sigma24=522, sigma28=214, sigma32=325, sigma38=423, sigma54=110, sigma70=3946
por fuente: BSGatlas=2226 (430 DBTBS estrictos), RegulonDB=3314
por confianza: C=73, S=928, W=2307, curated=2232
PMID: 3743/5540 filas (67.6%)
duplicados: 0
errores: 0
```

Los tests de procedencia
[`tests/test_sigma_gold.py`](../tests/test_sigma_gold.py) re-derivan las claves
y campos desde `PromoterSet.tsv`, `gt_promoters.tsv` y `gt_tss.tsv`, y
comprueban que no aparezcan etiquetas excluidas. Si cambian las fuentes, hay
que regenerar el snapshot, subir la versión del esquema y volver a correr
validador y tests.

---

## 8. Límites

1. **No es exhaustivo.** Solo TSS etiquetados y con coordenada. K-12 deja
   fuera 742 filas no evaluables y colapsa 1 duplicado; *B. subtilis* deja
   fuera 1.164 filas.
2. **Solo seis sigma.** `sigma19` (K-12) y `E`/`F`/`G`/`H`/`K`/`W`/`X`/`M`/
   `YlaC` (*B. subtilis*) no tienen PWM equivalente en el set usado y no son
   evaluables.
3. **No hay negativos.** No se puede medir especificidad ni tasa de falsos
   positivos sin un conjunto negativo adicional.
4. **K-12 es multi-label en 124 TSS.** Una fila por sigma; no asumir un único
   sigma por coordenada.
5. **`gene` incompleto.** 1.769 filas de K-12 y todas las de *B. subtilis* no
   tienen gen. No se imputa.
6. **`pmid` incompleto.** 1.797 filas no tienen PMID (todas de *B. subtilis*
   sin PubMed; K-12 tiene 100% de cobertura).
7. **`confidence` no calibrada.** `S`/`W`/`C` son niveles de RegulonDB;
   `curated` no es una medida de calidad.
8. **Snapshot, no consenso.** Si otra base de datos discrepa, el gold conserva
   la etiqueta de RegulonDB/DBTBS. No se reconciliaron fuentes externas.
9. **El mapeo de *B. subtilis* es a equivalentes PWM**, no una afirmación de
   que la nomenclatura de factores sea idéntica entre especies.
10. **No responde a los 6 casos.** Ver [`SIGMA_CASES.md`](SIGMA_CASES.md).

---

## 9. Procedencia y licencias

| Dato | Ruta | Uso |
|---|---|---|
| Promotores K-12 | [`data/regulondb/PromoterSet.tsv`](../data/regulondb/PromoterSet.tsv) | etiquetas, TSS, evidencia, PMIDs |
| Promotores K-12 (GT) | [`data/regulondb/gt/gt_promoters.tsv`](../data/regulondb/gt/gt_promoters.tsv) | contraste de TSS/hebra/sigma/gen |
| TSS *B. subtilis* | [`data/dbtbs/gt_tss.tsv`](../data/dbtbs/gt_tss.tsv) | etiquetas, TSS, evidencia |
| PWMs | [`data/regulondb/gt/promoter_pwms.tsv`](../data/regulondb/gt/promoter_pwms.tsv) | define qué sigma tienen equivalente |

El gold no descarga nada. Las licencias y condiciones de las fuentes están en
[`DATA_NOTICE.md`](../DATA_NOTICE.md) y
[`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md).
