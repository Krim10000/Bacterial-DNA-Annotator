# Resultados negativos

Lo que probamos y **no funcionó**, con los números. Sirve para que nadie repita
el experimento ni lo reintente con la misma configuración creyendo que es nuevo.

Regla de lectura: un resultado negativo no es un fallo; es información. Lo que
sí sería un fallo es volver a intentarlo sin cambiar la premisa.

---

## 1 · Boltz-2 no predice la dirección del efector (alostería)

**Qué se probó.** Co-folding de 6 represores con su ligando, con y sin operador
de ADN, comparando afinidades. Hipótesis: `Δ = afinidad(DNA) − afinidad(libre)`
< 0 ⇒ correpresor; > 0 ⇒ inductor.

**Resultado.** 2 de 6 aciertos de signo (33%), **p = 0,69**. Indistinguible del
azar. El diseño "reparado" (cognato vs barajado, mismo número de cadenas)
tampoco lo rescata: TrpR, que es correpresor, da Δ = +0,27 (inductor).

**Coste.** ~0,5 horas-GPU en Modal (3 × A10G, ~10 min).

**Qué haría falta para revisitarlo.** Otra fuente de señal: comparación de
bolsillos, dinámica molecular, o datos experimentales de conformación. No más
co-folding.

---

## 2 · La sintenia densa (ortología recíproca) no aporta

**Qué se probó.** RBH entre K-12 y 402 genomas (DIAMOND bidireccional,
1.337.745 proteínas). 175.910 pares recíprocos, 210 pares de K-12 conservados en
≥5 genomas.

**Resultado.**

| configuración | F1 |
|---|---:|
| unión de TSS | 0,798 |
| unión + fusión por sintenia dispersa | 0,801 |
| unión + fusión por sintenia densa | 0,801 |
| clasificador con sintenia dispersa | 0,804 |
| clasificador con sintenia densa | 0,804 |

Densa y dispersa dan **exactamente lo mismo**. El peso del rasgo en el modelo es
−0,03.

**Por qué.** Los 400 genomas son filogenéticamente diversos: el orden génico se
conserva poco entre especies lejanas. La cobertura sigue siendo baja (210 pares
de 3.025).

**Qué haría falta.** Genomas **cercanos** (cepas de la misma especie o género),
donde el orden sí se conserva. No más ortología entre especies lejanas.

---

## 3 · El detector de terminadores dependientes de Rho por secuencia

**Qué se probó.** Clasificador con rasgos de composición y motivos *rut*
(fracción de C, G, C−G, rachas, `CC`, `CAA`/`CAC`, ΔG) sobre la ventana de 100 nt
aguas arriba. Positivos: 37 terminadores dependientes de Rho; negativos: 338
independientes + 5.000 posiciones intergénicas.

**Resultado.** **AUC 0,673** (validación cruzada, ±0,075); sensibilidad 19–22%
con 5% de falsos positivos. Añadir los motivos *rut* no movió el AUC.

**Por qué.** Rho no reconoce una secuencia: es un motor de ARN que se engancha
en regiones C-ricas y pobres en G y transloca hasta la polimerasa. Su
especificidad es contextual, no local.

**Qué haría falta.** Datos experimentales de Rho (ChIP, Term-seq con y sin Rho)
o un modelo que incluya estructura del transcrito y velocidad de traducción.

---

## 4 · Unir TSS no es universal: depende de la independencia de las fuentes

**Qué se probó.** Unión de dos fuentes de TSS en las dos especies.

| especie | fuentes | ¿independientes? | F1 solo | F1 unión |
|---|---|---|---:|---:|
| K-12 | curados + dRNA-seq (15 estudios) | **sí** | 0,785 | **0,798** |
| *B. subtilis* | BSGatlas + DBTBS | **no** (43% repetido) | 0,799 | 0,795 |

**Resultado.** La unión mejora cuando las fuentes son independientes y empeora
cuando una contiene a la otra.

**La regla que queda.** Antes de unir, medir solapamiento: si >40% de la fuente
nueva ya está presente, aporta ruido.

---

## 5 · Clasificar "¿termina aquí?" está dominado por la tasa base

**Qué se probó.** Clasificador por frontera con evidencia de terminación
(horquilla + poli-U, Rho, composición).

**Resultado.** F1 0,473, **por debajo** del predictor constante "siempre
termina" (0,785), porque el 65% de los genes son ya final de TU. Con prior
calibrado y evidencia de read-through sube a 0,804, pero el margen sobre el
constante es de 0,019.

**Por qué.** La pregunta útil no es "¿termina?" (mayoritario) sino "¿hay
read-through?" (minoritario). Las señales de terminación pesan casi cero en esa
decisión: la horquilla tiene peso 0,02.

**Consecuencia.** No construir detectores de TTS para decidir fronteras de
operón. Son dos tareas distintas: anotar TTS (secuencia local) y decidir
fronteras (contexto del vecino).

---

## 6 · BATTER-TPE como clasificador de fronteras

**Qué se probó.** Las 1.389 predicciones de BATTER sobre las regiones candidatas,
usadas para decidir si un gen termina su TU.

**Resultado.** P 0,770, R 0,206, **F1 0,325**. Más preciso que las horquillas
pero encuentra mucho menos, porque su umbral está calibrado para alta precisión.

**Para qué sirve entonces.** Como **detector de terminadores** (anotación), no
como clasificador de fronteras. Ahí sí aporta: complementa a las horquillas y
entre ambos cubren familias distintas.

---

## 7 · Transferir TSS por ortología no funciona entre filo distintos

**Qué se probó.** Para extender lo documentado a genomas sin TSS: tomar el TSS
del ortólogo en K-12 y colocarlo a la misma distancia del inicio del gen.

**Validación en *B. subtilis* (donde sí hay verdad).**

| | |
|---|---:|
| proteínas traducidas | 4.330 |
| genes con ortólogo en K-12 | **430 (10%)** |
| TSS transferidos | 277 |
| **F1 contra operones BSGatlas** | **0,503** |
| referencia: TSS propios | 0,799 |
| referencia: motor de distancia | 0,752 |

**Por qué.** K-12 (proteobacteria) y *B. subtilis* (firmicute) son filo distintos:
solo el 10% de los genes tiene ortólogo. Con esa cobertura, la transferencia deja
la mayoría de los genes sin TSS y los pocos que lo reciben arrastran desfases.

**Consecuencia para los 400 genomas.** No se puede poblar su anotación de
transcripción por transferencia desde K-12. Opciones reales:

1. **Motor de distancia** (0,752), con confianza declarada — lo que ya hacemos.
2. TSS experimentales propios de cada genoma (no los hay).
3. Transferir desde un **pariente cercano** con TSS medidos (misma especie o
   género). No probado, pero es donde la sintenia y la ortología sí funcionan.
