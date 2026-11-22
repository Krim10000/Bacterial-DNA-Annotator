# Proteínas reguladoras: lo que aprendimos

> Síntesis de la línea de trabajo sobre represores bacterianos y sus efectores
> (AlloReg / módulo `allostery`). Documento de traspaso: qué se estableció, qué se
> midió, qué falló y por qué, y qué haría falta para cerrar cada pregunta.
>
> Fecha: 2026-11-21. Todo lo cuantitativo está enlazado a su archivo de salida.

---

## 1 · El problema, en una frase

Un represor tiene al menos dos conformaciones: una que une su operador (reprime) y otra
que no. Un metabolito puede estabilizar una de las dos. Saber **cuál** y **cuánto** es lo
que determina si el sistema responde, y es justo lo que no está en ninguna base de datos.

La pregunta se descompone en tres, y solo la primera es decidible hoy con computación:

| Pregunta | Naturaleza | Estado |
|---|---|---|
| ¿Tiene el represor una bisagra entre el dominio de unión a ADN y el cuerpo? | estructural | **decidible, con número** |
| ¿Qué ligando ocupa el bolsillo del efector? | termodinámica | sin resolver |
| ¿Es inductor o correpresor? | termodinámica | sin resolver |

---

## 2 · Hechos establecidos sobre la biología

### 2.1 "Activa" no es un estado de la proteína: es una relación con el operador

RegulonDB no guarda conformaciones estructurales. Guarda **especies moleculares**
(`LacI`, `LacI-allolactosa`, `CueR-Cu⁺`) y las clasifica en dos conjuntos: las que unen
el operador y las que no. "Activa" significa **competente para unir ADN**, no "estado A".

Consecuencia práctica: no se puede validar una conformación mirando la proteína sola. Se
valida contra el ADN (afinidad por el operador) o contra la geometría del motivo HTH.

### 2.2 Inductor y correpresor son propiedades del par (TF, rol), no de la proteína

Para un represor, un efector que estabiliza la forma unida a ADN es un **correpresor**;
para un activador, el mismo efector es un **coactivador**. Hay TFs duales (XylR, AraC) y
casos con el mismo efector anotado en ambas conformaciones (`MetJ`, `IclR`, `ExuR`,
`ZntR`, `NikR`, `AcrR`, `AllR`, `CueR`) — 11 en total. Esos no se fuerzan a una etiqueta.

### 2.3 El efector no es el sustrato

El caso canónico: el efector de LacI es **alolactosa**, no lactosa. Lo mismo con SAM
(no metionina), PRPP (no purina), trehalosa-6-P (no trehalosa). Cualquier modelo que
busque "el sustrato de la ruta" falla sistemáticamente por esto.

### 2.4 Hay dos movimientos distintos, y confundirlos arruina el análisis

| Movimiento | Qué cambia | Qué controla | Cómo se detecta |
|---|---|---|---|
| **Bisagra DBD–core** | el DBD gira respecto al cuerpo | unión al operador | ANM sobre el oligómero |
| **Gating del bolsillo** | movimiento interno del core | entrada/salida del efector | requiere estructura holo o co-folding |

Medimos el primero y **no** el segundo. A lo largo del modo bisagra, el bolsillo anotado
del efector no se deforma (z ≈ −0,3 en PurR y TreR). Son grados de libertad separados.

### 2.5 La unidad funcional es el oligómero, no el monómero

LacI es tetrámero. Su transición es cuaternaria: los dos DBD se abren y giran en
antiparalelo (modo "mariposa"). Con el monómero de AlphaFold, el modo más blando no mueve
el DBD preferentemente (enriquecimiento 0,50). Con el tetrámero experimental (1LBG),
2,35 y coseno antiparalelo −0,89.

### 2.6 No todos los represores son sensores metabólicos

Una fracción importante responde a iones o xenobióticos: MarR (salicilato), MprA
(salicilato, 2,4-dinitrofenol), ArsR (arsenito), RcnR (Co/Ni), Zur (Zn), CueR (Cu).
Sus regulones no tienen genes metabólicos, así que cualquier método que asuma "el efector
está en la ruta" es ciego por construcción en estos casos.

### 2.7 El bolsillo del efector está anotado en pocos casos

De los 35 represores decisivos, UniProt anota sitios de unión en **5**, y solo **4**
fuera del DBD (PurR/hipoxantina, TreR/trehalosa-6-P, BirA/biotina, NemR). Es la única
fuente que da coordenadas del bolsillo sin inferencia.

---

## 3 · Lo que medimos (y con qué números)

### 3.1 Base de datos disponible

| Dato | Valor | Archivo |
|---|---:|---|
| TFs en RegulonDB | 243 | `data/regulondb/TFSet.tsv` |
| Con efector anotado | 93 (38 %) | ídem |
| Efectores distintos | 109 | ídem |
| TFs con PMID de conformación | 230 | ídem |
| Interacciones TF–promotor con secuencia | 3.444 (216 TFs) | `data/regulondb/TF-RISet.tsv` |
| Pares represor + confianza S/C + evidencia experimental | 1.485 (148 TFs) | ídem |
| **Represores con efector en una sola conformación (decisivos)** | **36** | 27 inductor / 9 correpresor |
| Represores anotados "repressor" en 400 genomas | 1.555 | `data/parsed_annotations/` |
| TFs con familia declarada en 400 genomas | 14.850 | ídem |

El cuello no son los represores: es la **etiqueta de efector**. Con 36 casos decisivos, no
hay forma de entrenar un modelo desde cero; hay que transferir y validar.

### 3.2 Guilt-by-association: funciona, pero solo a distancia cero

Hipótesis: el efector es un metabolito de la ruta que el represor controla.
Evaluación retrospectiva sobre 40 casos evaluables (`00_c0_gba.py`):

| Estadístico | Real | Regulón aleatorio del mismo tamaño |
|---|---:|---:|
| Efector es participante directo de una reacción del regulón (d=0) | **0,80** | 0,053 (2,1 esperados) |
| Efector a un salto de distancia (d≤1) | 1,00 | 0,93 |
| recall@5, ranking dentro del conjunto | 0,425 | — |
| rango mediano | 7 | sobre ~10.700 candidatos |

**15× de enriquecimiento a distancia 0.** La expansión por saltos no aporta nada: cualquier
regulón aleatorio cubre el 90 % del metabolismo a un salto. Conclusión: GBA sirve como
**filtro de candidatos**, nunca como predictor.

### 3.3 Los baselines de secuencia fallan

| Baseline | balanced acc | MCC | p vs nulo |
|---|---:|---:|---:|
| Familia Pfam del DBD (LOO dentro de familia) | 0,423 | −0,211 | — (peor que la mayoría global) |
| ESM-2 8M + PCA10 + LogReg, **leave-one-family-out** | 0,585 | 0,151 | 0,239 |
| 5-fold aleatorio (con fuga de familia) | 0,513 | 0,024 | — |

Dos lecciones: (1) **LOFO con solo la familia es indefinido**, no difícil — al sacar la
familia no queda entrenamiento; (2) el split aleatorio rinde *peor* que el honesto, señal
de que la varianza la domina n, no la fuga. Estamos en el piso de ruido.

### 3.4 Los métodos estructurales: qué funcionó y qué no

| Método | Resultado | Veredicto |
|---|---|---|
| Docking ciego (Vina) sobre AFDB, 13 casos | top-1 0,308 vs 0,25 azar; Δ +0,07 kcal/mol | sin señal |
| Co-folding Boltz-2, 6 TFs, estado con ADN vs libre | 1/6 signos correctos (p=0,22) | ensayo inválido (composición distinta) |
| Boltz-2 reparado: ADN cognato vs barajado, MSA real | 2/6 (p=0,69); absoluto 7/12 | sin señal, ensayo válido |
| ANM sobre monómero (35 TFs) | LacI 0,50 de enriquecimiento DBD | mide la bisagra equivocada |
| ANM sobre tetrámero (1LBG) | 2,35; antiparalelo −0,89; subespacio 0,453 | **identifica el modo** |
| Detector de bisagra, 8 negativos | sensibilidad 1/1, especificidad 5/8 | **útil, acotado** |
| Gating de cavidad con P2Rank | LacI: 0 bolsillos ligandables en conformeros | inválido |
| Deformación del bolsillo anotado | z ≈ −0,3 (PurR, TreR, BirA) | el core no se deforma |

### 3.5 El detector de bisagra

Regla: algún modo con la varianza concentrada en el DBD (z > 2 contra ventanas aleatorias)
y movimiento antiparalelo entre cadenas (< −0,5).

| Proteína | Rol | z_dbd | anti | Bisagra |
|---|---|---:|---:|---|
| LacI | positivo | 3,49 | −0,89 | sí |
| λ cI | negativo | 3,98 | −0,97 | **sí (falso positivo)** |
| LexA | negativo | 2,70 | −0,90 | **sí (falso positivo)** |
| 434 Cro | negativo | 2,02 | −0,43 | borde |
| 434 cI, Fis, IHF, SinR, λ Cro | negativos | < 2 | — | no |

Los falsos positivos son interpretables: λ cI y LexA **tienen** bisagra (linker flexible,
inactivación por proteólisis), pero no bolsillo de efector. El detector mide bisagra
DBD–core, no interruptor por metabolito.

---

## 4 · Errores que costaron tiempo (y cómo evitarlos)

1. **AlphaFold DB es una estructura estática.** No contiene conformaciones alternativas.
   El ensamble se obtiene por otros medios (AF2 con submuestreo de MSA, AF-Cluster,
   BioEmu, AlphaFlow), y el ANM solo da direcciones armónicas, no poblaciones.
2. **El API de AFDB no sirve las anotaciones funcionales** (`/api/annotations/<acc>` → 404).
   Los sitios de unión están en UniProt; hay que pedirlos allí.
3. **Comparar afinidades entre complejos de composición distinta no es válido.** Añadir ADN
   sube `ligand_iptm` también para el ligando señuelo. La corrección es mantener la misma
   composición (ADN cognato vs barajado).
4. **Los SMILES no se comparan por cadena.** Con coincidencia exacta se perdieron 17 de 44
   efectores. Con InChIKey (bloque esqueleto) se recuperaron todos.
5. **La extrapolación lineal de un modo produce geometrías no físicas.** P2Rank no
   encontró ni un bolsillo ligandable en LacI desplazado.
6. **El control positivo es obligatorio.** LacI salvó dos conclusiones falsas: el ranking
   del ANM de monómero y el gating de cavidad.
7. **Optimismo con n pequeño.** Un piloto con n=4 daba top-1 0,75; el run completo (n=13)
   bajó a 0,31.
8. **Desbalance de clases.** 27 inductores y 9 correpresores: predecir siempre "inductor"
   acierta el 75 %. Hay que reportar balanced accuracy, MCC y recall por clase.

---

## 5 · Recursos: qué hay y qué vale la pena

| Fuente | Qué aporta | Calidad |
|---|---|---|
| **RegulonDB** | sitios con secuencia, evidencia ECO, PMID, efector y conformación | A (E. coli) |
| **UniProt** | sitios de unión, metal, activos | A, pero escaso (5/35) |
| **AlphaFold DB** | estructura única, pLDDT, PAE | B (sin ensamble) |
| **PDB** | las tres cosas a la vez: proteína + operador + ligando | A (decenas de complejos) |
| **AlloRep** | >5.750 variantes fenotípicas de LacI/GalR, 65 estructuras | A (familia única) |
| **RegPrecise** | efectores por regulón, 640 genomas | B (inferencia comparativa) |
| **CollecTF / DBTBS** | sitios fuera de E. coli, filtrables por técnica | A/B |
| **Reg-Seq** | matrices de energía con ΔG medido | A (E. coli) |
| **P2TF** | familia y arquitectura de dominios | B |
| **Literatura** | solo 7/35 TFs tienen datos cuantitativos recuperables | D (dispersa) |

---

## 6 · Implicaciones de diseño

- **No entrenar desde cero con 36 casos.** Transferir (ESM-2, ProkBERT) y validar; el
  modelo grande con pocos datos solo produce sobreajuste.
- **Partir por familia, no por TF.** La química del bolsillo la fija la familia; si no,
  el test miente.
- **Normalizar por señuelo.** Cualquier Δ entre condiciones debe compararse contra el mismo
  Δ con un ligando señuelo (diferencia en diferencias).
- **Separar las tres preguntas.** Un modelo que intente responder las tres a la vez no
  responde ninguna.
- **Documentar el modo de fallo.** λ cI y LexA fallan el detector por una razón concreta
  (bisagra sin bolsillo), lo que permite acotar el alcance en vez de descartar la métrica.

---

## 7 · Qué haría falta para cerrar cada pregunta

| Pregunta | Qué se necesita | Coste |
|---|---|---|
| ¿Bisagra DBD–core? | oligómero por plantilla PDB + ANM + nulo de ventanas | hecho; extender a 24 TFs, ~30 min GPU |
| ¿Qué ligando? | ΔΔG medido o co-folding con transición inducida | GPU + metodología que hoy no funciona |
| ¿Inductor o correpresor? | `K_d(operador, ± efector)` medido (EMSA/ITC) | banco |
| Validación de todo | >5.750 variantes de AlloRep como ground truth de mutantes | CPU |

La vía más rentable sin banco es **AlloRep**: es la única fuente con etiquetas de fenotipo
de mutantes para una familia entera, y permite validar cualquier score estructural.

---

## 8 · Preguntas abiertas

1. ¿La bisagra DBD–core predice algo funcional, o es genérica en los pliegues HTH de dos
   dominios? Hace falta un set negativo más grande y diverso.
2. ¿Cuántos represores de los 400 genomas anotados tienen bisagra? Es una pregunta de
   census, no de mecanismo, y es respondible con el detector actual.
3. ¿El gating del bolsillo es intra-core? La evidencia de PurR/TreR (z ≈ 0 con el modo
   DBD) sugiere que sí, pero faltan estructuras holo para medirlo.
4. ¿La mayoría de los represores con efector son realmente alostéricos por un metabolito,
   o hay una fracción grande con mecanismos distintos (iones, xenobióticos, proteólisis)?
   El 20 % de los 36 no son metabólicos, y esa fracción puede ser mayor en otros filos.

---

## 9 · Archivos de esta línea

```
scripts/allostery/
  00_c0_gba.py               evaluador GBA (C0)
  01_family_baseline.py      baseline de familia
  02_esm_baseline.py         baseline ESM-2 con LOFO
  03_pocket_specificity.py   docking ciego
  04_boltz_modal.py          co-folding Boltz-2 (Modal)
  05_boltz_analyze.py        análisis de Δ entre estados
  06_anm_modal.py            ANM sobre 35 estructuras AF
  07_laci_oligomer.py        prueba de oligómero
  08_hinge_specificity.py    control negativo del detector
  09_cavity_gating.py        gating de cavidad (descartado)
  10_pocket_patch.py         deformación del bolsillo anotado
  modal_guard.sh             control de GPU en Modal
data/allostery/
  c0_gba.tsv / .json / .md   GBA por caso
  baseline_family.*  baseline_esm.*
  pocket_specificity.*  boltz_analysis*.*
  hinge_specificity.json  cavity_gating.json  pocket_patch.json
  laci_oligomer.json  anm_results.json
docs/
  PLAN_ALOSTERIA.md          plan y decisiones fase a fase
  MODAL.md                   regla de GPU en Modal
```

Los cuatro negativos (secuencia, docking, co-folding, gating) y el positivo acotado
(bisagra) están documentados con su métrica en `PLAN_ALOSTERIA.md` §8.
