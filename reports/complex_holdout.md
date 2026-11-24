# Holdout de complejos contra catálogos ComplexPortal locales

> **Veredicto: `not_evaluable` para un holdout independiente.** Los únicos catálogos locales son de *E. coli* K-12 y el modelo se entrenó en K-12 con 273 complejos gold. No hay catálogo de otra especie ni validación experimental. Lo que sigue separa lo **in-sample** de lo **exploratorio**; nada de esto es precisión de generalización.

Generado: `2026-11-24T02:28:47+00:00` · herramienta: `tools/complex_holdout.py`

## 1. Veredicto del holdout

| Estrategia | Estado | Detalle |
|---|---|---|
| Por genoma | not_evaluable | No hay catálogo ComplexPortal para un genoma fuera del entrenamiento; el único catálogo local es K-12. |
| Por complejo | exploratorio_no_independiente | Complejos raw-only ausentes del gold de entrenamiento |

- Genomas de entrenamiento declarados: `ecoli_k12`.
- Catálogos independientes encontrados: 0.
- El holdout por genoma es **not_evaluable** porque no existe catálogo de un genoma no entrenado.
- El holdout por complejo es **exploratorio**: los complejos raw-only no estaban en el gold, pero el genoma K-12 sí es el genoma de entrenamiento.
- **No se puede concluir** precisión/recall de generalización. **Sí se puede** reportar cobertura de catálogo, recall exploratorio de los complejos raw-only y contexto in-sample.

## 2. Procedencia del modelo y catálogos

| Dato | Valor |
|---|---|
| Genoma de entrenamiento | ecoli_k12 |
| Gold de entrenamiento (declarado) | 273 |
| Gold de entrenamiento (encontrado) | 273 |
| Pares positivos (`n_pos`) | 1.773 |
| Pares negativos (`n_neg`) | 265.042 |
| Fuente de procedencia | model_json |

| Catálogo | Genoma | Complejos | ≥2 genes | <2 genes | En gold |
|---|---|---|---|---|---|
| complexportal_83333 | ecoli_k12 | 324 | 283 | 41 | 273 |
| gold_k12 | ecoli_k12 | 273 | 273 | 0 | 273 |

El catálogo crudo ComplexPortal contiene 324 complejos; el gold de entrenamiento tiene 273. Los 51 ausentes del gold son la base del holdout por complejo, pero solo 10 tienen ≥2 genes mapeados.

## 3. Holdout por genoma

**No evaluable.** No hay `gold_<genoma>.tsv` ni `raw/complexportal_<taxid>.tsv` para un genoma fuera de `ecoli_k12`. No se puede estimar precisión/recall de generalización con los datos locales.

## 4. Holdout por complejo (exploratorio)

Se usan 10 complejos de ComplexPortal crudo que no están en el gold de entrenamiento y tienen ≥2 genes mapeados. Se excluyen 41 por tener <2 genes mapeados (homómeros o subunidades no mapeadas). El genoma sigue siendo K-12; esto **no** es un holdout independiente.

| Umbral Jaccard | Recuperados | Catálogo | Recall |
|---|---|---|---|
| 0.3 | 8 | 10 | 0,800 |
| 0.5 | 6 | 10 | 0,600 |
| 0.7 | 4 | 10 | 0,400 |
| 0.9 | 3 | 10 | 0,300 |

| Complejo | Genes | Mejor predicción | Jaccard | Recuperado @ 0,5 |
|---|---|---|---|---|
| CPX-4389 | 3 | ecoli_k12:cx373 | 1,000 | sí |
| CPX-4450 | 3 | ecoli_k12:cx509 | 1,000 | sí |
| CPX-5125 | 3 | ecoli_k12:cx424 | 1,000 | sí |
| CPX-4449 | 3 | ecoli_k12:cx185 | 0,750 | sí |
| CPX-282 | 4 | ecoli_k12:cx42 | 0,667 | sí |
| CPX-5942 | 2 | ecoli_k12:cx182 | 0,500 | sí |
| CPX-1924 | 8 | ecoli_k12:cx51 | 0,444 | no |
| CPX-5780 | 8 | ecoli_k12:cx511 | 0,375 | no |
| CPX-1926 | 4 | ecoli_k12:cx51 | 0,286 | no |
| CPX-2244 | 2 | — | 0,000 | no |

**Precisión: `not_evaluable`.** La precisión requeriría puntuar contra todo el catálogo, que incluye complejos de entrenamiento; se reporta solo recall del holdout.

_La selección no es aleatoria: muchos complejos raw-only son homómeros/unidades sin mapeo. No usar estos recall como generalización._

## 5. Contexto in-sample K-12 (no es holdout)

Predicciones `clusters_ecoli_k12_evidencia.tsv` contra 273 complejos gold. El modelo vio estos complejos al entrenar. **No usar como precisión de generalización.**

### 5.1 Barrido de umbral

| Umbral | Pred. TP | Catálogo recuperado | Precisión | Recall | F1 |
|---|---|---|---|---|---|
| 0.3 | 189 | 232 | 0,197 | 0,850 | 0,320 |
| 0.5 | 167 | 198 | 0,174 | 0,725 | 0,281 |
| 0.7 | 103 | 106 | 0,107 | 0,388 | 0,168 |
| 0.9 | 80 | 81 | 0,083 | 0,297 | 0,130 |

En el umbral principal 0.5: precisión **0,174**, recall **0,725**, F1 **0,281**.

### 5.2 Por confianza

| Confianza | Predicciones | Precisión | Recall del catálogo | F1 |
|---|---|---|---|---|
| alta | 498 | 0,261 | 0,546 | 0,353 |
| baja | 233 | 0,056 | 0,055 | 0,055 |
| media | 229 | 0,105 | 0,124 | 0,114 |

### 5.3 Por evidencia

| Evidencia (clase) | Predicciones | Precisión | Recall del catálogo | F1 |
|---|---|---|---|---|
| documentado | 86 | 1,000 | 0,374 | 0,544 |
| nuevo | 757 | 0,000 | 0,000 | 0,000 |
| parcial | 84 | 0,691 | 0,242 | 0,358 |
| subcomplejo | 33 | 0,697 | 0,110 | 0,190 |

`nuevo` significa no documentado en el gold; su precisión 0 es por construcción y no un error del modelo.

_K-12 es el genoma de entrenamiento y el gold de 273 complejos se usó para entrenar. Estas cifras son in-sample; no son precisión de generalización._

## 6. Cobertura de transferencias K-12 (no es precisión)

| Métrica | Valor |
|---|---|
| Filas de transferencia | 19.085 |
| Genomas con transferencia | 400 |
| Complejos transferidos | 258 |
| Transferidos en gold | 258 |
| Cobertura de complejos gold | 94,5 % |
| Genes K-12 en transferencias | 701 |
| Cobertura de genes gold | 95,9 % |
| Complejos gold completos | 255 |
| Complejos gold parciales | 3 |
| Complejos gold sin cobertura | 15 |
| Fracción media por complejo | 0,942 |

_Las transferencias se derivan del mismo catálogo K-12; esto mide cobertura, no precisión ni generalización._

## 7. Limitaciones

- Hay 273 complejos gold + 324 crudos locales (597 filas de catálogo) y todos pertenecen a K-12 (ComplexPortal taxón 83333). No hay catálogo curado de otra especie; por eso el holdout independiente es `not_evaluable`.
- El modelo declara `trained_genome=ecoli_k12` y 273 complejos gold usados en entrenamiento. Cualquier evaluación sobre K-12 es in-sample y no mide generalización.
- El holdout por complejo usa 10 complejos raw-only con >=2 genes mapeados; excluye 41 por falta de genes/estequiometría. No es una muestra aleatoria: los homómeros y complejos sin mapeo quedan fuera, y el genoma sigue siendo K-12.
- El emparejamiento es por Jaccard de conjuntos de genes. No usa `curado_id`, no evalúa estequiometría, isoformas ni interfaces, y puede fallar con complejos que comparten subunidades.
- Las transferencias cubren 258 complejos y se derivan del mismo catálogo K-12; su cobertura no es precisión. RBH no trae cobertura de alineamiento ni demuestra interacción física.
- Los clusters `model_<genoma>.complexes.tsv` de los 400 genomas no se evalúan aquí porque no hay catálogo ComplexPortal de sus especies; solo K-12 tiene gold local.
- No hay validación experimental en otras especies, ni gold independiente para los 400 genomas. El informe no cierra brechas ni cambia el pipeline, la app/GUI o los datos.

## 8. Reproducir

```bash
.venv/bin/python -m pytest tests/test_complex_holdout.py -q
.venv/bin/python tools/complex_holdout.py \
    --out-json reports/complex_holdout.json \
    --out-md reports/complex_holdout.md
```

El script no usa red y solo escribe los reportes indicados.
