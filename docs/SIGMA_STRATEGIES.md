# Estrategias probadas para sigma

Evaluación GroupKFold (5 folds, grupos de 10 kb) sobre el gold curado. Son
resultados internos, no validación externa.

## K-12

| Estrategia | Top-1 | Macro-F1 | Balanced acc |
|---|---:|---:|---:|
| Mayoría σ70 | 0,6127 | 0,1284 | 0,1667 |
| argmax z (actual) | 0,5667 | 0,4645 | 0,5560 |
| argmax score bruto | 0,5925 | 0,4799 | 0,5808 |
| two-stage margen 0,5 | **0,6388** | 0,4866 | 0,5384 |
| two-stage margen 1,0 | **0,6669** | 0,4781 | 0,4872 |
| jerárquico logístico | **0,7057** | **0,5099** | 0,5165 |
| jerárquico + contexto | 0,6551 | 0,4949 | 0,5452 |
| k-mer 4/5 | 0,476 | 0,296 | 0,313 |

La mejor opción K-12 es el modelo jerárquico, pero usa entrenamiento sobre el
mismo gold y PWMs de RegulonDB: **no es generalización independiente**.

## B. subtilis

| Estrategia | Top-1 | Macro-F1 | Balanced acc |
|---|---:|---:|---:|
| Mayoría σ70 | 0,8756 | 0,1867 | 0,2000 |
| PWM propia argmax | 0,2965 | 0,1536 | 0,3889 |
| two-stage margen 2 | 0,8585 | 0,2549 | 0,3058 |
| two-stage margen 3 | 0,8751 | 0,2153 | 0,2180 |
| k-mer 4-mer | 0,8531 | **0,4191** | **0,4162** |
| híbrido k-mer + two-stage, umbral 0,75 | **0,8756** | **0,4220** | 0,3872 |

El híbrido **iguala** la mayoría trivial en top-1 y mejora mucho el macro-F1:
detecta algunas alternativas sin destruir σ70. No supera la mayoría en top-1.

## Conclusión

- No hay una estrategia que gane en todo.
- K-12: jerárquico local, top-1 70,6 %.
- B. subtilis: híbrido k-mer + two-stage, top-1 87,6 %, macro-F1 0,422.
- La PWM propia sola es peor que el fallback σ70 en top-1.
- El contexto ayuda a B. subtilis, pero no a K-12.
- El k-mer aporta a B. subtilis, no a K-12.
- Falta validación en especie held-out y negativos reales.

## Políticas disponibles

`webapp/sigma_policy.py`:

```python
argmax_order(scores)
two_stage_order(scores, delta=1.0)
hybrid_order(scores, external, threshold=0.75, delta=2.0)
```

## Implementación

`sigma_annotator.annotate_genome(..., policy="two_stage", delta=1.0)` y el
módulo sigma (`policy`, `policy_delta`) permiten activar la política sin
cambiar el escáner. El valor por defecto sigue siendo `argmax` para no alterar
resultados históricos sin una decisión explícita.

CLI:

```bash
python scripts/sigma/10_annotate_sigma.py --genome subtilis \
  --policy two_stage --delta 1.0
```

## Mejor ajuste encontrado

Protocolo: GroupKFold 5, grupos de 10–100 kb; sin gold de otra especie.

| Genoma | Estrategia | Top-1 | Macro-F1 |
|---|---|---:|---:|
| K-12 | jerárquico logístico | **0,6972** | **0,5107** |
| B. subtilis | híbrido k-mer 4/5 (−40/+5) + two-stage, umbral 0,6 | **0,9003** | **0,5184** |
| B. subtilis | mayoría σ70 | 0,8756 | 0,1867 |

El híbrido de B. subtilis usa una probabilidad k-mer calculada out-of-fold;
el two-stage solo decide cuando la señal k-mer no es concluyente.

No se declara modelo multi-especie: cada especie tiene su propio calibrador.
