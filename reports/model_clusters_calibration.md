# Calibración de clusters del modelo K-12 en los 400 genomas

> **Hipótesis, no validación.** Los clusters son candidatos del modelo K-12. Este informe mide solapamiento descriptivo y prioriza revisiones; no estima precisión, recall ni interacción física. Sin red.

Generado: `2026-11-24T01:27:22+00:00` · herramienta: `tools/calibrate_model_clusters.py`

## 1. Insumos y alcance

| Insumo | Valor |
|---|---|
| Genomas con `model_<genoma>.complexes.tsv` | 400 |
| Filas de clusters | 281.090 |
| Genomas con `genes.tsv` | 400 |
| Genomas con `transfer_<genoma>.tsv` | 400 |
| Filas de transferencia | 19.085 |
| Complejos K-12 transferidos (únicos) | 258 |
| Firmas K-12 | 960 |
| Complejos gold K-12 | 273 |
| Genes K-12 con coordenadas | 4.308 |
| Umbral de intervalo colocalizado | 50.000 pb |
| Transferencias con intervalo colocalizado | 7.341 |
| Firmas con intervalo colocalizado | 722 |

Los 281.090 clusters contienen 1.210.835 subunidades listadas y 1.210.739 genes resueltos con coordenadas. Hay 0 clusters sin ningún gen resuelto. El estado es siempre `complex_candidate`.

**Hallazgos principales**

- Confianza: 37,9% alta, 3,2% media y 58,8% baja. La confianza alta es mayoritariamente de tamaño 2 (79,0% de los clusters alta).
- Transferencias: 19.580 clusters (6,97%) comparten al menos un gen y 3.306 (1,18%) tienen cobertura ≥ 0,5.
- Intervalo: 93,72% solapa algún intervalo amplio, pero solo 2,67% uno colocalizado (≤ 50.000 pb). El intervalo amplio no es evidencia.
- Firmas: 19.476 clusters (6,93 %) tienen match por gen; se alcanzan 204/960 firmas.
- Priorización: 207 P1 y 1.350 P2; el primer cluster es `GCF_000006985.1_ASM698v1:cx346` (GCF_000006985.1_ASM698v1, score 5,12).

## 2. Agregación

### 2.1 Por confianza

| Confianza | Clusters | % del total | Tamaño medio | Tamaño mediana | Clusters de 2 | % de 2 | Genes resueltos (media) |
|---|---|---|---|---|---|---|---|
| alta | 106.547 | 37,9 % | 2,25 | 2,0 | 84.204 | 79,0 % | 2,25 |
| media | 9.126 | 3,2 % | 4,60 | 4,0 | 2.906 | 31,8 % | 4,60 |
| baja | 165.417 | 58,8 % | 5,62 | 4,0 | 0 | 0,0 % | 5,62 |

La confianza alta está dominada por clusters de dos subunidades: el tamaño no crece con la confianza. Esto es una propiedad del modelo y debe considerarse al leer cualquier solapamiento. La confianza alta tampoco implica más anclaje: 1,8% de los clusters alta comparte un gen con transferencia frente a 10,7% de los baja.

### 2.2 Por tamaño del cluster

| Subunidades | Clusters | % del total | Media | Mediana | Cobertura gen transf. (media) | Intervalo coloc. (media) |
|---|---|---|---|---|---|---|
| 2 | 87.110 | 31,0 % | 2,00 | 2,0 | 0,012 | 0,008 |
| 3 | 70.057 | 24,9 % | 3,00 | 3,0 | 0,013 | 0,010 |
| 4 | 40.708 | 14,5 % | 4,00 | 4,0 | 0,017 | 0,014 |
| 5-6 | 41.728 | 14,8 % | 5,40 | 5,0 | 0,024 | 0,019 |
| 7-9 | 24.504 | 8,7 % | 7,77 | 8,0 | 0,032 | 0,028 |
| 10 o más | 16.983 | 6,0 % | 14,59 | 12,0 | 0,045 | 0,045 |

El solapamiento genético y el intervalo colocalizado crecen con el tamaño del cluster (cobertura media 0,012 en tamaño 2 frente a 0,045 en 10+). Es consistente con que un cluster más grande tiene más oportunidades de contener un gen transferido; no implica mejor calidad.

### 2.3 Por genoma

Hay 400 genomas en el JSON. Clusters por genoma: mínimo 66, mediana 575,0, máximo 2.750. Tamaño medio por genoma: 3,09–8,35 subunidades.

Genomas con más clusters:

| Genoma | Clusters | Confianza alta | Tamaño medio | Con gen transf. |
|---|---|---|---|---|
| GCF_060030665.1_ASM6003066v1 | 2.750 | 1.258 | 3,59 | 2,9 % |
| GCF_060026435.1_ASM6002643v1 | 2.557 | 1.087 | 3,56 | 3,2 % |
| GCF_026965555.1_ASM2696555v1 | 2.537 | 1.135 | 3,41 | 2,2 % |
| GCF_037157775.1_ASM3715777v1 | 2.532 | 1.102 | 3,64 | 3,6 % |
| GCF_033547035.1_ASM3354703v1 | 2.429 | 1.065 | 3,61 | 3,1 % |
| GCF_001027285.1_ASM102728v1 | 2.372 | 1.083 | 3,57 | 3,0 % |
| GCF_002813575.1_ASM281357v1 | 2.241 | 874 | 3,76 | 3,1 % |
| GCF_060030825.1_ASM6003082v1 | 2.132 | 898 | 3,60 | 3,8 % |

La tabla completa por genoma (confianza, tamaño, evidencia y solapamiento) está en el JSON, campo `por_genoma`.

### 2.4 Evidencia

| Token de evidencia | Clusters | % del total |
|---|---|---|
| modelo K-12 | 281.090 | 100,0 % |
| score medio | 281.090 | 100,0 % |
| score mínimo | 281.090 | 100,0 % |
| adyacente | 281.043 | 100,0 % |
| operón | 101.015 | 35,9 % |

La evidencia del modelo contiene el prior K-12 y señales genómicas (`operón`, `adyacente`, scores). No contiene STRING ni sintenia de K-12 aplicada a otros genomas.

## 3. Solapamiento

**Gen**: Jaccard e intersección de genes. **Intervalo**: fracción del intervalo del cluster cubierta por el intervalo de referencia **colocalizado** (mismo contig y span <= 50.000 pb). Son planos independientes: un intervalo puede solaparse sin compartir genes. No se afirma precisión.

### 3.1 Con transferencias K-12 (mismo genoma)

| Métrica | Clusters | % del total |
|---|---|---|
| Con al menos un gen compartido | 19.580 | 7,0 % |
| Cobertura de gen ≥ 0,5 | 3.306 | 1,2 % |
| Con intervalo amplio (diagnóstico) | 263.439 | 93,7 % |
| Con intervalo colocalizado | 7.501 | 2,7 % |
| Intervalo colocalizado ≥ 0,5 | 3.831 | 1,4 % |
| Con puente K-12 (≥1 fila) | 19.580 | 7,0 % |

| Confianza | Clusters | Con gen | % | Int. amplio | % | Int. coloc. | % | Cob. gen mediana | Jaccard mediano |
|---|---|---|---|---|---|---|---|---|---|
| alta | 106.547 | 1.876 | 1,8 % | 99.664 | 93,5 % | 892 | 0,8 % | 0,000 | 0,000 |
| media | 9.126 | 84 | 0,9 % | 8.185 | 89,7 % | 58 | 0,6 % | 0,000 | 0,000 |
| baja | 165.417 | 17.620 | 10,7 % | 155.590 | 94,1 % | 6.551 | 4,0 % | 0,000 | 0,000 |

Las transferencias cubren 19.085 filas y 258 complejos K-12. 7.341 filas tienen intervalo colocalizado (span <= 50.000 pb); el span mediano de todas las transferencias es 364.935 pb. 251 complejos tienen al menos un cluster con gen compartido; 213 tienen al menos un cluster con intervalo colocalizado.

**Lectura:** el intervalo amplio es casi universal porque las subunidades transferidas suelen estar dispersas; por eso se reporta solo como diagnóstico. El intervalo colocalizado es la lectura conservadora y aun así más permisiva que el gen. Cualquier conclusión debe apoyarse en el gen.

### 3.2 Con firmas K-12

| Métrica | Valor |
|---|---|
| Firmas K-12 | 960 |
| Firmas con match por gen | 204 |
| Firmas con match por intervalo colocalizado | 142 |
| Firmas con intervalo colocalizado | 722 |
| Firmas sin intervalo colocalizado | 238 |
| Clusters con puente K-12 | 19.580 |
| Clusters con puente colocalizado | 8.906 |
| Clusters con puente a complejo con `curado_id` | 17.718 |
| Clusters sin puente K-12 | 261.510 |
| Clusters con match de firma por gen | 19.476 |
| Clusters con match de firma por intervalo | 7.284 |

| Confianza | Clusters con puente | % | Puente coloc. | % | Con gen firma | % | Con int. firma | % |
|---|---|---|---|---|---|---|---|---|
| alta | 1.876 | 1,8 % | 642 | 0,6 % | 1.870 | 1,8 % | 507 | 0,5 % |
| media | 84 | 0,9 % | 23 | 0,2 % | 84 | 0,9 % | 18 | 0,2 % |
| baja | 17.620 | 10,7 % | 8.241 | 5,0 % | 17.522 | 10,6 % | 6.759 | 4,1 % |

El match con firmas exige un puente de transferencia; no es una comparación directa de secuencias. El match por intervalo exige además que la firma y la fila-puente tengan intervalo colocalizado, por eso es menor que el match por gen. Las firmas no cubren todos los clusters ni todos los complejos.

### 3.3 Cobertura de gold K-12

| Referencia | Valor |
|---|---|
| Complejos gold K-12 | 273 |
| Gold con al menos una transferencia | 258 |
| Complejos transferidos | 258 |
| Transferidos que están en gold | 258 |
| Firmas con `curado_id` | 239 |
| `curado_id` únicos | 213 |
| `curado_id` presentes en gold | 213 |
| Firmas sin `curado_id` | 721 |

Gold K-12 es la única referencia curada; no hay gold para las otras especies. La cobertura de gold mide alcance del catálogo, no acierto de los clusters.

## 4. Priorización para revisión

**Fórmula:** `score = peso_confianza + (cobertura_gen + jaccard_gen) + 0,5*jaccard_firma + 0,5*fraccion_intervalo_colocalizado + 0,5*min(tamano,12)/12`. Es una heurística: suma confianza (alta 2,0, media 1,0, baja 0,0), solapamiento genético con transferencia, y componentes menores de firma, intervalo y tamaño. **No es una probabilidad.**

| Tier | Criterio | Clusters |
|---|---|---|
| P1 | score >= 4,0 | 207 |
| P2 | 3,0 <= score < 4,0 | 1.350 |
| P3 | score < 3,0 | 279.533 |

Top 30 de 200 candidatos en el JSON (filtrado a clusters con cobertura de gen de transferencia > 0). La tabla completa está en `priorizacion.top`.

| # | Genoma | Cluster | Conf. | Tam. | Cob. gen | Jacc. gen | Int. coloc. | Jacc. firma | Transferencia | Firma | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | GCF_000006985.1_ASM698v1 | cx346 | alta | 3 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2197 | ecoli_k12:cx485 | 5,12 |
| 2 | GCF_000024985.1_ASM2498v1 | cx671 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 3 | GCF_018771645.1_ASM1877164v1 | cx540 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 4 | GCF_060026435.1_ASM6002643v1 | cx2520 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5222 | ecoli_k12:cx750 | 5,08 |
| 5 | GCF_044336245.1_ASM4433624v1 | cx1155 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 6 | GCF_060030825.1_ASM6003082v1 | cx1784 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5222 | ecoli_k12:cx750 | 5,08 |
| 7 | GCF_025558325.1_ASM2555832v1 | cx329 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-4862 | ecoli_k12:cx651 | 5,08 |
| 8 | GCF_000020505.1_ASM2050v1 | cx345 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 9 | GCF_015277775.1_ASM1527777v1 | cx1372 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5222 | ecoli_k12:cx750 | 5,08 |
| 10 | GCF_019048065.1_ASM1904806v1 | cx783 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-4021 | ecoli_k12:cx630 | 5,08 |
| 11 | GCF_033971325.1_ASM3397132v1 | cx404 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 12 | GCF_000008625.1_ASM862v1 | cx340 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 13 | GCF_019048105.1_ASM1904810v1 | cx714 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5623 | ecoli_k12:cx742 | 5,08 |
| 14 | GCF_000800395.1_ASM80039v1 | cx499 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 15 | GCF_055384745.1_ASM5538474v1 | cx300 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 16 | GCF_019048105.1_ASM1904810v1 | cx601 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-4021 | ecoli_k12:cx630 | 5,08 |
| 17 | GCF_002302565.1_ASM230256v1 | cx361 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-4021 | ecoli_k12:cx630 | 5,08 |
| 18 | GCF_001953175.1_ASM195317v1 | cx357 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5542 | ecoli_k12:cx837 | 5,08 |
| 19 | GCF_000317575.1_ASM31757v1 | cx984 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5623 | ecoli_k12:cx742 | 5,08 |
| 20 | GCF_000010305.1_ASM1030v1 | cx703 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5884 | ecoli_k12:cx770 | 5,08 |
| 21 | GCF_000185805.1_ASM18580v1 | cx292 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 22 | GCF_007993755.1_ASM799375v1 | cx904 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5623 | ecoli_k12:cx742 | 5,08 |
| 23 | GCF_014672715.2_ASM1467271v2 | cx346 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5542 | ecoli_k12:cx837 | 5,08 |
| 24 | GCF_000011905.1_ASM1190v1 | cx248 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5542 | ecoli_k12:cx837 | 5,08 |
| 25 | GCF_001953175.1_ASM195317v1 | cx370 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5222 | ecoli_k12:cx750 | 5,08 |
| 26 | GCF_918320885.1_Harlequin_e… | cx340 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5201 | ecoli_k12:cx897 | 5,08 |
| 27 | GCF_001747405.1_ASM174740v1 | cx438 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-2113 | ecoli_k12:cx940 | 5,08 |
| 28 | GCF_050039385.1_ASM5003938v1 | cx1417 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5245 | ecoli_k12:cx759 | 5,08 |
| 29 | GCF_032102965.1_ASM3210296v1 | cx601 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5623 | ecoli_k12:cx742 | 5,08 |
| 30 | GCF_013004105.1_ASM1300410v1 | cx751 | alta | 2 | 1,00 | 1,00 | 1,00 | 1,00 | CPX-5884 | ecoli_k12:cx770 | 5,08 |

`Transferencia` es el complejo K-12 de la mejor fila de transferencia por gen; `Firma`, la firma K-12 con mayor Jaccard. Ordenar por `priority_score` prioriza clusters con anclaje de transferencia y confianza alta, pero no garantiza que el cluster sea correcto. Un cluster sin transferencia puede seguir siendo una hipótesis válida; simplemente no tiene anclaje en este análisis.

## 5. Limitaciones

- El modelo es un prior entrenado/derivado de K-12: los 281.090 clusters de los 400 genomas heredan esa distribución. No hay gold independiente para las otras especies; gold_k12.tsv cubre solo 273 complejos de K-12.
- Los clusters son candidatos (`complex_candidate`): no demuestran interacción física, estequiometría ni función. La confianza alta/media/baja no es una probabilidad calibrada.
- Las 19.085 filas de transferencia cubren 258 complejos K-12 y usan RBH; no hay cobertura de alineamiento en los TSV actuales. La evidencia transferida es `documentado_inferido`/ECO:0005544, no experimental.
- Las 960 firmas K-12 solo cubren 239 filas con `curado_id`; el resto no tiene ancla curada. Solo 204 firmas reciben al menos un cluster por gen y 142 por intervalo colocalizado. La comparación con firmas es indirecta (puente de transferencia) y puede inflarse o perderse por parálogos y mapeos muchos-a-muchos.
- El intervalo amplio es engañoso: los genes de una transferencia tienen span mediano de 364.935 pb y los de una firma 3.383 pb; por eso 93,7% de los clusters tocan algún intervalo amplio frente a 2,7% con intervalo colocalizado (<= 50.000 pb) y 7,0% con al menos un gen compartido. Solo 7.341/19.085 transferencias y 722/960 firmas tienen intervalo colocalizado. Un solapamiento de intervalo por sí solo no es evidencia.
- La priorización usa pesos elegidos a mano (confianza, gen, intervalo colocalizado y tamaño). Cambiar los pesos cambia el orden. Los tiers P1/P2/P3 son cortes operativos, no clases biológicas.
- Los intervalos se derivan del mínimo/máximo de genes mapeados; en transferencias multi-contig se elige el contig con más genes y se pierde el resto. Los alias `gene_id`/`locus_tag` pueden colapsar duplicados.
- El informe es descriptivo y sin red. No se editó el pipeline, la GUI ni los datos; no se hizo commit.

## 6. Reproducir

```bash
.venv/bin/python -m pytest tests/test_calibrate_model_clusters.py -q
.venv/bin/python tools/calibrate_model_clusters.py \
    --out-json reports/model_clusters_calibration.json \
    --out-md reports/model_clusters_calibration.md
```

El script no usa red y solo escribe los dos reportes indicados.
