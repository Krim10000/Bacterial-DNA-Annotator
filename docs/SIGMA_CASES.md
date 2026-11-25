# Los 6 casos: qué son y qué rendimiento tenemos

> **Pregunta:** ¿tenemos el rendimiento de cómo el sistema detecta cada uno de
> los 6 casos?
>
> **Respuesta corta:** **no**. Los 6 casos son una taxonomía de mecanismos
> propuesta en el diseño; el sistema actual detecta **tipos de sigma** y
> **promotores sigma-dependientes**, no los 6 mecanismos. No hay detector por
> caso ni benchmark por caso.
>
> El gold nuevo (`data/sigma/gold/`) tiene TSS + sigma, **no** etiqueta de
> mecanismo. Por eso no cambia esta respuesta; sí mejora la evaluación por
> capas: [`SIGMA_LAYERED_EVAL.md`](SIGMA_LAYERED_EVAL.md).

---

## 1. Qué son los 6 casos

| Caso | Mecanismo | Qué habría que detectar |
|---|---|---|
| 1 | Activador clase I | Sitio upstream de −35, contacto con α-CTD de RNAP |
| 2 | Activador clase II | Sitio que solapa −35, contacto con σ |
| 3 | Activador MerR | Motivo MerR, spacer 19–20, torsión del DNA |
| 4 | Represor estérico | Sitio en el promotor que bloquea la RNAP |
| 5 | DNA looping | Dos sitios, fase helicoidal, oligómero |
| 6 | Anti-activación / trampa | Sitio solapado con activador, bloqueo de la activación |

Son **mecanismos**, no tipos de sigma. Un mismo σ70 puede ser caso 4 o caso 5;
MerR es una familia de activadores, no un sigma. Por eso el benchmark de sigma
no responde esta pregunta.

---

## 2. Qué detecta el sistema hoy

| Pieza | Qué hace | ¿Responde a los 6 casos? |
|---|---|---|
| `10_annotate_sigma.py` | Asigna σ70/σ54/σ32/σ38/σ28/σ24 a promotores | **No**; solo tipo de sigma |
| `20_group_sigma.py` | Agrupa factores por ortólogo | **No** |
| `40_build_sigma_regulon.py` | Une sigma → operón → TF del regulón | **No**; no infiere mecanismo |
| `browser.py` / `sigma_view.py` | Muestra sigma, confianza y evidencia | **No** |
| `30_validate_sigma.py` | Benchmark por tipo de sigma | **No**; no hay etiqueta de caso |

No hay un campo `case`/`mechanism` en `sigma_promoters.tsv`. Tampoco hay
estructuras 3D, superposición con RNAP, looping ni anti-activación en el
módulo actual.

---

## 3. Lo que sí tenemos: rendimiento por tipo de sigma

K-12, benchmark TSS-level (3.315 TSS evaluables):

| Sigma | Support | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
| σ70 | 1.997 | 0,857 | 0,480 | 0,616 |
| σ24 | 522 | 0,528 | 0,676 | 0,593 |
| σ32 | 312 | 0,523 | 0,615 | 0,566 |
| σ28 | 144 | 0,253 | 0,458 | 0,326 |
| σ38 | 241 | 0,181 | 0,523 | 0,269 |
| σ54 | 99 | 0,172 | 0,354 | 0,231 |

Esto **no** es rendimiento por caso:

- no distingue activador vs represor;
- no distingue clase I vs clase II;
- no distingue MerR de otro activador;
- no distingue looping de un sitio único;
- no distingue anti-activación de represión estérica.

---

## 4. Por qué no se puede derivar de los datos actuales

Para medir un caso necesitás:

1. **Etiqueta de mecanismo** por sitio/operón.
2. **Evidencia independiente** (estructura, mutantes, TSS, RNA-seq, literatura).
3. **Un detector implementado** para ese mecanismo.
4. **Controles negativos** del mismo tipo.

Hoy hay etiqueta de **sigma** (RegulonDB/DBTBS), no de mecanismo. Intentar
inferir “caso 4” desde el score PWM sería circular: el PWM no codifica el
mecanismo.

---

## 5. Qué se necesitaría para medir cada caso

| Caso | Detector mínimo | Benchmark mínimo |
|---|---|---|
| 1 | Posición upstream + motivo + contacto α-CTD | 10–20 activadores clase I curados |
| 2 | Solapamiento con −35 + motivo + contacto σ | 10–20 activadores clase II |
| 3 | Motivo MerR + spacer 19–20 + geometría | 10 MerR con estructura/estructura modelada |
| 4 | Sitio que solapa promotor/TSS + bloqueo RNAP | 20 represores con sitio y TSS |
| 5 | Dos sitios + distancia/fase + oligómero | 10 operones con looping documentado |
| 6 | Sitio solapado con activador + efecto | 10 casos de anti-activación |

Sin ese set, cualquier “rendimiento por caso” sería inventado.

---

## 6. Qué haría primero

1. **No** presentar los 6 casos como detectados.
2. Construir un `case_benchmark.tsv` con casos curados de K-12 y B. subtilis
   (aunque sean 5–10 por caso).
3. Implementar primero los casos más decidibles:
   - **caso 4** (represor estérico) por solapamiento con promotor/TSS;
   - **caso 3** (MerR) por motivo + spacer;
   - **caso 1/2** por posición relativa a −35 y TSS.
4. Dejar **caso 5/6** para cuando haya estructura/RNAP o datos de mutantes.
5. Reportar precisión/recall por caso, no un accuracy global.

---

## 7. Respuesta en una frase

**No tenemos rendimiento por los 6 casos: tenemos rendimiento por tipo de sigma
(σ70, σ24, σ32, σ28, σ38, σ54) y una capa sigma→regulón. Los 6 casos son
mecanismos que aún no están implementados ni etiquetados; para medirlos hace
falta un benchmark de mecanismos, no reutilizar el de sigma.**
