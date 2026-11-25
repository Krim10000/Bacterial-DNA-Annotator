# Combinaciones de features para refinar sigma

> **Qué es:** un análisis descriptivo de reglas interpretables sobre
> `layer_sigma_features.tsv` (una fila por TSS × sigma). No es un modelo
> entrenado ni una validación *held-out*.
>
> **Para qué:** ver qué combinaciones de señal, arquitectura y contexto
> mejoran la precisión de un tipo de sigma, aunque bajen el recall.

---

## 1. Features

| Feature | Origen | Uso |
|---|---|---|
| `z` | score del PWM contra señuelos | señal |
| `pvalue` | empírico contra señuelos | señal |
| `margin` | z del sigma − mejor z rival | competición |
| `spacer_dev` | desviación del espaciador canónico | arquitectura |
| `dist_to_start` | distancia del mejor hit al inicio | localización |
| `inventory` | factor del genoma presente | contexto |

No se usan `n_true`, `gold_confidence` ni `source`: son metadatos del gold, no
señales del pipeline.

---

## 2. Resultados K-12

Reglas evaluadas sobre 3.189 TSS. Se reporta la mejor regla por F1 con
`min_support=10`.

| Sigma | Regla | Precision | Recall | F1 | Lift vs top-1 |
|---|---|---|---:|---:|---:|---:|
| σ54 | p ≤ 0,05 **y** margen ≥ 1 **y** distancia ≤ 60 | **0,600** | 0,212 | 0,313 | **3,30×** |
| σ28 | margen ≥ 1 **y** distancia ≤ 60 | **0,500** | 0,259 | 0,341 | **1,75×** |
| σ38 | margen ≥ 1 **y** distancia ≤ 60 | **0,314** | 0,249 | 0,278 | **1,67×** |
| σ24 | z ≥ 2 **y** distancia ≤ 60 | **0,734** | 0,270 | 0,395 | 1,20× |
| σ32 | z ≥ 2 **y** distancia ≤ 60 | **0,733** | 0,353 | 0,476 | 1,22× |
| σ70 | z ≥ 2 | **0,869** | 0,329 | 0,477 | 1,01× |

Lectura:

- **σ54 y σ38** son los que más se benefician de exigir **margen + distancia**:
  suben precisión, pero pierden recall. Sirven como filtro de alta confianza.
- **σ24 y σ32** mejoran precisión con `z ≥ 2` y distancia ≤ 60, pero su F1
  global baja porque el top-1 ya tenía buen recall.
- **σ70** no mejora F1 con estas reglas: ya es la clase mayoritaria; filtrar
  solo quita recall.

Ninguna regla debe presentarse como “la regla del sigma”. Son **filtros
candidatos** para una revisión contextual.

---

## 3. Resultados *B. subtilis*

Con el gold actual, las reglas no mejoran el top-1:

- σ28: `distancia ≤ 60` → precisión 0,037, recall 0,972, F1 0,071.
- σ38: `spacer en rango` → precisión 0,082, recall 1,000, F1 0,151.
- σ54: `distancia ≤ 60` → precisión 0,006, recall 0,909, F1 0,012.

**Interpretación:** no es que las features sean malas; es que las PWMs son de
*E. coli* y la transferencia a *B. subtilis* falla. Ajustar reglas sobre este
gold sería overfitting a una transferencia mala.

---

## 4. Cómo usar estas reglas

1. **No reemplazan el ranking.** Se aplican después del top-1/top-3.
2. **Suben precisión a costa de recall.** Úsalas para `sospecha_fuerte` o para
   priorizar revisión, no para “corregir” el top-1.
3. **Dependen de la especie.** Las reglas de K-12 no se copian a Subtilis.
4. **Necesitan validación *held-out*.** Elegir la regla sobre el mismo gold no
   es validación.
5. **El contexto manda.** Una regla sin `inventory`/ruta/expresión no distingue
   un falso positivo de un hallazgo nuevo.

---

## 5. Ejemplo de uso

```python
# Filtro de alta precisión para σ54, solo como candidato
if row["sigma"] == "sigma54":
    if row["pvalue"] <= 0.05 and row["margin"] >= 1 and row["dist_to_start"] <= 60:
        row["confidence_label"] = "sospecha_fuerte"
```

Esto **no** cambia el `sigma` asignado: cambia la confianza de la fila. Es la
forma correcta de incorporar una regla de features.

---

## 6. Reproducir

```bash
python scripts/sigma/50_eval_layers.py --gold-root data/sigma/gold --all
python scripts/sigma/60_feature_combinations.py --all --min-support 10 --top 8
```

Salidas:

```text
data/sigma/feature_combinations/
  feature_combinations.tsv
  feature_combinations.json
  feature_combinations.md
```

---

## 7. ¿Y las rutas metabólicas?

Se probaron como prior de desempate: **+0,35 pp** en K-12 y **+0,40 pp** en
margen bajo; con peso alto empeoran. No se incorporan como feature fuerte.
Detalle: [`SIGMA_PATHWAY_PRIOR.md`](SIGMA_PATHWAY_PRIOR.md).

## 8. Límite

**Una regla con lift alto puede ser un artefacto del gold.** Por eso el
siguiente paso no es aplicarla a los 400, sino probarla en un *held-out*
(otra especie, otro conjunto de TSS) y ver si mantiene precisión y recall.
