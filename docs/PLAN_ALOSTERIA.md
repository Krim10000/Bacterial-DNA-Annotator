# Plan de trabajo — módulo de alosteria y especificidad de efector

> **ESTADO v1.0 (2026-11-22): línea cerrada.** La energía libre y la magnitud del
> efector quedan fuera de v1.0 por falta de recursos de cómputo (GPU). Lo que entra es
> el sistema de detección y anotación de represores, con preparación local y puntuación
> de contacto. Ver [`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md).

> Objetivo: dado un represor, decidir **si responde a un efector**, **en qué dirección**
> (inductor / correpresor) y **a cuál molécula**. Es el módulo que falta para responder
> la pregunta biológica del proyecto; no requiere las tres torres ni el DNA-LLM propio.
>
> Formulación textual del objetivo y criterios de éxito: ver
> [`ALOSTERIA_V1_CIERRE.md` §0](ALOSTERIA_V1_CIERRE.md).
>
> Reglas heredadas del repo: un dueño por archivo, salidas con `evidence`/`confidence`/
> `status`, nada `confirmed` sin experimento, todo idempotente y por etapa de CLI.
>
> **Síntesis de hallazgos y traspaso:** [`PROTEINAS_REGULADORAS.md`](PROTEINAS_REGULADORAS.md).
> **No toca** `scripts/regulation/` ni `data/regulation/` (dueño: `regulon-builder`).

## 0 · Estado de partida (verificado, no estimado)

Fuente: `data/regulondb/` (RegulonDB, *E. coli* K-12 NC_000913.3, licencia académica).

| Dato | Valor | Archivo |
|---|---:|---|
| TFs totales | 243 | `TFSet.tsv` |
| TFs con efector anotado | 93 (38%) | `TFSet.tsv` |
| Efectores distintos | 109 | `TFSet.tsv` |
| TFs con PMID de conformación | 230 | `TFSet.tsv` |
| Interacciones TF–promotor con secuencia | 3.444 (216 TFs) | `TF-RISet.tsv` |
| Pares represor + `S`/`C` + `EXP-*` | 1.485 (148 TFs) | `TF-RISet.tsv` |

**Conjunto de evaluación decisivo** (represores con efector en una sola conformación):

| Rol | n | confianza `S` |
|---|---:|---:|
| Inductor (efector en conformación inactiva) | 27 | 11 |
| Correpresor (efector en conformación activa) | 9 | 7 |
| **Total decisivo** | **36** | **18** |
| Ambiguos (efector en ambas conformaciones) | 7 | — |

Los 7 ambiguos (`MetJ`, `IclR`, `ExuR`, `ZntR`, `NikR`, `AcrR`, `AllR`) se excluyen del test.

**Consecuencia dura:** n=36 no es un set de entrenamiento, es un **benchmark**. Y el
desbalance 27/9 hace que un modelo que siempre diga "inductor" acierte 75%. Se reporta
balanced accuracy, MCC y recall por clase; **nunca** accuracy agregada.

## 1 · Target operativo (evita la circularidad)

No se predice "cuál conformación es activa" (eso es definicional). Se predice:

```
ΔΔG_efector = ΔG_unir(ligando | TF unido a ADN) − ΔG_unir(ligando | TF libre)
```

| Signo (represor) | Lectura |
|---|---|
| negativo | prefiere la forma unida a ADN ⇒ **correpresor** |
| positivo | prefiere la forma libre ⇒ **inductor** |
| ≈ 0 | no es efector de unión (puede actuar en otro paso) |

Latencia 2: `|ΔΔG| ≳ 1–2 kcal/mol` **y** `[S]` dentro de ~10× del K_d del estado
relevante. Sin eso, el interruptor no conmuta in vivo y no se reporta como efector.

## 2 · Fases, entregables y puertas

### F0 · Verdad de referencia — 3 días

`00_build_ground_truth.py` (planificado)

- Join `TFSet` × `NetworkRegulatorGene` (rol) × `TF-RISet` (sitios) × PDB (holo/apo)
  × AlloRep (variantes LacI/GalR).
- Salida `data/allostery/ground_truth.tsv`: `tf`, `bnumber`, `familia`, `rol`,
  `efector`, `estado` (`A`/`I`/`ambiguo`), `modo` (`inductor`/`correpresor`/`ambiguo`),
  `evidencia` (código ECO), `pmid`, `confianza`, `fuente`.
- Aceptación: 93 filas, 36 decisivas, 7 ambiguas excluidas, 100% con PMID o marcadas `assumed`.
- **Puerta G0:** ninguna fila `documented` sin PMID de experimento de unión.

### F1 · Baselines sin GPU — 1 semana

`01_family_baseline.py` (familia + arquitectura de dominios, P2TF/Pfam) y
`02_esm_baseline.py` (ESM-2 congelado + regresión logística, split **leave-one-family-out**).

- Salida: `data/allostery/baseline.json` con balanced acc, MCC, recall por clase y
  desglose por familia (LacI/GalR, TetR, LysR, GntR, MarR, AraC, DtxR…).
- **Puerta G1 (criterio de muerte):** si el modo se predice con balanced acc ≥ 0,80 en
  LOFO, **la Fase 2 no se justifica para "modo"**. Cae cerca de 0,50 ⇒ hay señal en el
  bolsillo y F2 vale la pena.

**Resultado de F1 (2026-11-20).** Dos baselines, ninguno pasa.

| baseline | balanced acc | MCC | p vs nulo | veredicto |
|---|---:|---:|---:|---|
| Familia (Pfam DBD), LOO dentro de familia | 0,423 | −0,211 | — | peor que la mayoría global (0,500) |
| ESM-2 8M + PCA10 + LogReg, **LOFO** | 0,585 | 0,151 | **0,239** | indistinguible del azar |

Detalles que importan: (1) **LOFO con solo la familia es indefinido**, no difícil: 6 de 14
familias tienen un único TF decisivo, así que al sacar la familia no queda entrenamiento;
la medida honesta es LOO dentro de familia y da **0,423**, por debajo de predecir siempre
"inductor". (2) El 5-fold aleatorio (con fuga de familia) da **0,513**, *peor* que LOFO:
la varianza la domina n=35, no la fuga — estamos en el piso de ruido. (3) Nulo por
permutación: 0,502 ± 0,109, p95 = 0,682, p = 0,239.

**Lectura:** el modo no está escrito en la secuencia de forma que un LM pequeño pueda leer
con n=35. Es consistente con que lo determine la química del bolsillo frente a un ligando
concreto — exactamente lo que F2 mide. Los dos baselines fallan, así que **F2 queda
justificada**, con la advertencia de que 8M es el ESM-2 más pequeño: repetir con 650M es
una mejora barata antes de dar el veredicto por definitivo.

### F2 · Bolsillo y ΔΔG — 3 semanas

- `03_pockets.py`: fpocket/P2Rank sobre AlphaFold/PDB; features del bolsillo
  (volumen, hidrofobicidad, residuos del sitio, distancia al HTH).
- `04_state_models.py`: dos referencias — complejo TF–operador de PDB (estado unido) y
  apo (AF3, o estructura PDB sin operador). Co-folding Boltz-2 del ligando en ambas.
  **Documentar que el score no es kcal/mol**: solo diferencia relativa.
- `05_rank_ligands.py`: retrieval del bolsillo contra los 109 efectores + ChEBI;
  ranking por `ΔΔG × [S]` fisiológica (mM en *E. coli*).
- Aceptación: top-5 dentro de familia; IPTG y tetraciclina deben salir como ligandos de
  LacI/TetR **sin** ser metabolitos naturales (inductores gratuito: test de fuga).
- **Puerta G2:** Spearman ρ > 0,5 frente a K_d experimentales en LacI/GalR, y acierto de
  signo ≥ 70% en los 36 decisivos. Si no: entregable = modo de fallo documentado.

**Resultado F2-reducido (2026-11-21).** Sin co-folding local se usó un sustituto
decidible: docking ciego con AutoDock Vina sobre estructuras AlphaFold DB, con diseño
**específico por ligando** (el mismo ligando contra la proteína propia y contra 3
proteínas sorteadas; así se evita el sesgo de tamaño del ligando).

| métrica (n=13, 3 controles) | observado | azar |
|---|---:|---:|
| top-1 | **0,308** | 0,250 |
| top-2 | 0,462 | 0,500 |
| rango normalizado medio | 0,654 | 0,500 |
| Δ medio (propio − controles) | **+0,07 kcal/mol** | 0 |
| fracción con Δ<0 | 0,462 | 0,500 |

**Veredicto: no se detecta especificidad de bolsillo con este ensayo.** 4/13 top-1 contra
3,25 esperados por azar (p ≈ 0,65). Lección metodológica: el piloto de n=4 daba top-1
0,75 y el run completo lo bajó a 0,31 — con n pequeño la señal es ruido.

Esto **no** resuelve G2 ni prueba que el bolsillo no importe: prueba que el *instrumento*
(docking ciego, exhaustiveness 2, caja de 22 Å, n=13) no tiene sensibilidad para la
hipótesis. G2 sigue pendiente del método previsto (co-folding Boltz-2/AF3 en GPU).
Excluido `PaaX` (fenilacetil-CoA, >12 enlaces rotables: no converge). El ensayo es
reproducible con `03_pocket_specificity.py --serial`; con más núcleos libres conviene
subir `--exhaustiveness` y `--per-ligand`.

**Resultado co-folding Boltz-2 en Modal (2026-11-21).** 24 predicciones (6 TFs ×
ligando verdadero/señuelo × estado con ADN/libre) en A10G, ~25 min, ~0,5 USD.
`msa: empty`, operador desde `TF-RISet`, afinidad del módulo de Boltz-2.

| TF | modo | Δ verdadero | Δ señuelo | esperado | observado |
|---|---|---:|---:|---|---|
| BetI | inductor | −0,313 | −0,223 | positivo | negativo ✗ |
| CytR | inductor | −0,271 | −0,150 | positivo | negativo ✗ |
| GalR | correpresor | +0,033 | +0,001 | negativo | positivo ✗ |
| MurR | inductor | −0,173 | +0,184 | positivo | negativo ✗ |
| PurR | correpresor | −0,095 | −0,117 | negativo | negativo ✓ |
| TrpR | correpresor | +0,214 | +0,005 | negativo | positivo ✗ |

Δ = afinidad(con ADN) − afinidad(libre), en log10(IC50 µM).

**Veredicto: 1/6 signos correctos (p=0,22), G2 no se cumple.** Las magnitudes
(0,03–0,31 log) son ruido: 1–2× en IC50. La discriminación absoluta verdadero vs
señuelo es 8/12 (gap medio −0,65 log), tampoco significativa.

**Confusión diagnosticada:** `ligand_iptm` sube al añadir ADN *también para el
señuelo* (BetI 0,36→0,73; TrpR 0,13→0,76). Es un artefacto global del complejo, no
biología. Comparar la cabeza de afinidad entre composiciones distintas de complejo
no es una medida válida de ΔΔG.

**Ensayo reparado (2026-11-21, misma noche).** Se corrigieron las dos objeciones
principales: MSA real (`--use_msa_server`, no `empty`) y comparación **sin cambiar la
composición del complejo** — ambos estados son complejos de 4 cadenas del mismo tamaño;
solo cambia la identidad del ADN (operador cognato vs el mismo operador barajado). La
normalización por señuelo (DiD) sigue aplicada.

| TF | modo | Δ verdadero | Δ señuelo | DiD | esperado | observado |
|---|---|---:|---:|---:|---|---|
| BetI | inductor | −0,068 | −0,013 | −0,056 | positivo | negativo ✗ |
| CytR | inductor | −0,111 | −0,068 | −0,043 | positivo | negativo ✗ |
| GalR | correpresor | −0,012 | +0,128 | −0,139 | negativo | negativo ✓ |
| MurR | inductor | +1,259 | −0,076 | +1,335 | positivo | positivo ✓ |
| PurR | correpresor | +0,012 | −0,141 | +0,153 | negativo | positivo ✗ |
| TrpR | correpresor | +0,269 | −0,190 | +0,459 | negativo | positivo ✗ |

Δ = afinidad(cognato) − afinidad(barajado). **2/6 signos correctos (p=0,69).**
Discriminación absoluta verdadero vs señuelo: 7/12.

**Veredicto: reparar el ensayo no rescató la señal.** El único efecto grande y correcto
es `MurR` (DiD +1,34); el resto está en 0,04–0,46, dentro del ruido. Con la confusión de
composición eliminada, la conclusión es más fuerte que antes: la cabeza de afinidad de
Boltz-2 no recupera el modo inductor/correpresor en este conjunto, y no hay evidencia de
que un ΔΔG predicho entre estados sea medible con esta herramienta.

**Qué haría falta para un test válido:** (1) MSA real (`--use_msa_server`), no
`empty`; (2) fijar el estado del receptor con plantillas PDB (TF unido a operador)
en lugar de dejar que Boltz lo prediga; (3) no comparar afinidades entre complejos
de composición distinta: usar la confianza de los residuos del bolsillo con el ADN
como cuerpo rígido, o ΔΔG de mutaciones del bolsillo; (4) réplicas con varias
semillas y promedio del ensemble.

### F3 · Falsación con mutantes — 1,5 semanas

`06_mutant_phenotype.py`: predecir fenotipo de mutantes puntuales de LacI/GalR en 3 clases
(constitutivo / super-represor / no une ADN) y validar contra AlloRep (>5.750 variantes).

- **Puerta G3:** balanced acc ≥ 0,65, estratificado por posición y por familia, contra
  baseline de "posición en el bolsillo". Si no pasa, el modelo de bolsillo no captura
  alosteria y el módulo **no** se publica como predictivo.

### F4 · Integración — 1 semana

- `07_report.py` → `data/allostery/report.md` + `predictions.tsv` con
  `evidence`/`confidence`/`status`.
- Pestaña en la GUI **coordinada con `lead`** (no editar `webapp/app.py` sin acuerdo).
- Lista de experimentos propuestos (EMSA ± ligando, reportero dosis-respuesta) para los
  casos de mayor incertidumbre: son las hipótesis que valen un paper.

## 3 · Estructura del módulo

```
scripts/allostery/
  00_build_ground_truth.py  01_family_baseline.py  02_esm_baseline.py
  03_pockets.py             04_state_models.py     05_rank_ligands.py
  06_mutant_phenotype.py    07_report.py           _common.py
data/allostery/
  ground_truth.tsv  baseline.json  pockets/  predictions.tsv  report.md
```

`run_pipeline.sh`: añadir una etapa `allostery` con el mismo patrón que `components`/`pwms`
(una línea en el `case`, sin tocar el resto).

## 4 · Datos externos a incorporar

| Fuente | Uso | Nota |
|---|---|---|
| PDB / BioLiP | pares holo/apo, complejos TF–operador–ligando | única fuente con las tres cosas a la vez |
| [P2TF](http://www.p2tf.org) | familia y arquitectura de dominios | lado proteína, no sitios |
| [AlloRep](https://www.sciencedirect.com/science/article/abs/pii/S0022283615005380) | verdad de mutantes LacI/GalR | >5.750 variantes, 65 estructuras |
| [CollecTF](https://collectf.org/) / [DBTBS](https://dbtbs.hgc.jp/) | sitios fuera de *E. coli*, filtrables por técnica | Tier A |
| [RegPrecise](https://regprecise.lbl.gov/) | efector por regulón, 640 genomas | inferencia comparativa → Tier B |
| [Reg-Seq](https://github.com/RPGroup-PBoC/RegSeq) | matrices de energía con ΔG medido | calibración de F2 |
| ChEBI | identidad química y concentraciones | normaliza nombres del `TFSet` |

Licencias: RegulonDB académico/no comercial; verificar AlloRep y PRODORIC antes de publicar.

## 5 · Riesgos y mitigación

| Riesgo | Mitigación |
|---|---|
| 9 correpresores, desbalance 3:1 | MCC + recall por clase; no ajustar umbrales sobre el test |
| Apo desordenada → AF3 poco fiable | usar pares PDB holo/apo; marcar confianza baja y excluir de G2 |
| Boltz-2 sesgado a fármacos (nM–µM) vs efectores fisiológicos (mM) | solo ranking relativo; nunca ΔG absoluto |
| Fuga por familia | LOFO obligatorio en F1 y G2 |
| Conflicto con `regulon-builder` | módulo y `data/` propios; solo lectura de `data/regulondb/` |
| n=36 y sobreajuste del plan | F1 y G1 existen para decidir antes de gastar GPU |

## 6 · Cronograma y coste

| Fase | Duración | GPU |
|---|---|---|
| F0 | 3 días | — |
| F1 | 1 semana | — |
| F2 | 3 semanas | ~900 co-folds ⇒ ~20 h A10G ⇒ **$25–40** |
| F3 | 1,5 semanas | — |
| F4 | 1 semana | — |
| **Total** | **~7 semanas** (tiempo parcial) | **< $50** |

## 7 · Criterio de éxito global

1. `ground_truth.tsv` reproducible y con evidencia por fila.
2. Un número honesto por subtarea: ¿tiene efector? / ¿inductor o correpresor? / ¿qué ligando?
3. G2 y G3 pasan, **o** queda documentado dónde falla y por qué (eso también es un resultado).
4. Tres hipótesis experimentales concretas, priorizadas, para validar en banco.

## 8 · Cómo crecer el set de represores

Los represores **no** son el recurso escaso. En los 400 genomas del proyecto hay
**1.555 CDS anotados "repressor"**, **14.850** con familia declarada en el producto
(LysR 2.249, MarR 1.863, GntR 1.802, LacI 1.562, FNR 1.143, AraC 1.122, MerR 1.043,
PadR 907, Fur 772, Rrf2 690, IclR 531, AsnC 493…) y 30.069 "transcriptional regulator"
genéricos. Ninguno tiene efector anotado. El cuello es la etiqueta, no la proteína.

| Eje | Qué añade | Orden | Coste | Calidad |
|---|---|---|---|---|
| Especies | CollecTF, PRODORIC, DBTBS, RegPrecise | ~10³ TFs con sitios, efectores por regulón | 1–2 sem | B |
| Familias | revisiones TetR/LysR/GntR/MarR/AraC + AlloRep | +200–400 pares (represor, efector) | 2–3 sem | A |
| Estructura | PDB + AlphaFold DB (millones) | bolsillos sin etiqueta | 1 sem | preentrenamiento |
| Banco | síntesis + EMSA + reportero | 20 TFs elegidos por información mutua | 3–4 sem | A, y es tuyo |

Reglas: (1) el producto de NCBI subestima el censo — TetR aparece 157 veces donde
deberían ser miles; para contar familias hace falta `hmmsearch` sobre Pfam, no el
nombre. (2) Toda fila nueva entra como `status=hypothesis` hasta tener PMID de unión.

**La vía que escala sin descargar nada: guilt-by-association.** Represor divergente
junto a su ruta ⇒ el efector candidato está entre los metabolitos de esa ruta. Con las
reconstrucciones Rhea/EC ya hechas (`data/metabolism/`) genera cientos de hipótesis
priorizadas — lo que RegPrecise hace a mano. Se refuerza con conservación: si operador
y gen se conservan juntos en los 400 genomas, la hipótesis sube.

Aunque el set pase de 36 a 500, sigue siendo pequeño para entrenar desde cero: el plan
no cambia, mejora el benchmark (G2/G3 más robustos). Nadie tiene 10⁴ pares
(represor, efector) curados.

### 8.1 · Cómo se comprueba el guilt-by-association

La afirmación no es "el efector es X", es "el efector está en los metabolitos de la ruta
vecina **más de lo que estaría por azar**". Se verifica con ranking contra nulos, no con
aciertos. Cuatro controles, del barato al caro:

| Control | Qué mide | Criterio de paso |
|---|---|---|
| **C0 retrospectivo** (LOOCV sobre los 43 con efector) | recall@1/@3/@5 y rango mediano, ciego | recall@5 ≥ 50% |
| **C1 sobre la familia** | ¿aporta sobre el prior trivial (LacI→azúcar, Fur→metal)? | ≥ 15 puntos de mejora |
| **C2 especie held-out** (construir en *E. coli*, probar en DBTBS/RegPrecise) | generalización, no memorización del genoma | recall@5 ≥ 40% |
| **C3 prospectivo pre-registrado** (10–20 predicciones congeladas, commit + hash) | la única prueba real | ≥ 50% confirmadas; < 30% ⇒ no supera al azar |

Nulos obligatorios de C0: azar sobre la red metabólica, azar sobre el genoma, prior de
clase de metabolito, y familia sola. Sin los cuatro, un recall alto no significa nada.

**Seis formas en que GBA miente:** (1) circularidad — el nombre de producto ya viene del
contexto génico, así que no se valida contra la anotación sino contra el PMID; (2) divergón
≠ co-regulación: exigir operador en el intergénico, no solo distancia; (3) efector ≠
sustrato (alolactosa ≠ lactosa, PRPP ≠ purina, SAM ≠ metionina) ⇒ candidatos a ±1 reacción
(Rhea); (4) cofactores ubicuos (ATP, NAD, H₂O, H⁺) inundan el ranking ⇒ excluirlos como
hace RegPrecise; (5) concentración fisiológica; (6) TFs duales con 2+ efectores.

**Escalera de estado** (con la convención del repo):

```
hypothesis   ← solo vecindad + ruta
corroborated ← + conservación del operador y del gen en los 400 genomas, o estructura del bolsillo
documented   ← + PMID de unión (EMSA/ITC/reportero)
```

Base para probarlo ya, en *E. coli*: **142** represores con sitio `S`/`C`, **126** con sitio a
≤100 pb del primer gen; 43 tienen efector (retrospectivo) y ~99 no (cola prospectiva).

**Resultado de C0 (2026-11-20, `scripts/allostery/00_c0_gba.py`, matching por InChIKey).**
40 casos evaluables de 44 (4 sin ruta: `MprA`, `MarR`, `YegW`); 0 artefactos de SMILES.
Veredicto: **GBA pasa en su forma estrecha y falla en la amplia.**

| estadístico (maxdist=1) | regulón real | regulón aleatorio | lectura |
|---|---:|---:|---|
| efector a distancia 0 (participante directo) | **0,800** (32/40) | 0,053 (2,1 esperados) | **15× — señal real** |
| efector a distancia ≤1 | 1,000 | 0,929 | sin señal: la expansión destruye la información |
| recall@5 dentro del conjunto | 0,425 | 0,0005 al azar | filtro, no oráculo |
| rango mediano | 7 | — | 10.749 candidatos |

Por modo: inductor recall@5 0,476 (n=21), correpresor 0,667 (n=6).

Consecuencia: GBA **no** se usa como predictor de efector, solo como **filtro de conjuntos
de candidatos** (regla d=0) para F2: reduce ~10⁴ compuestos a ~10. La expansión por saltos
queda prohibida como criterio. Los priors de familia y de frecuencia **no fueron testables**
(0 casos con prior): no hay dos TFs de la misma familia con efector resuelto.

Modos de fallo identificados: (a) d=0 con ranking bajo cuando el efector es un metabolito
compartido por muchas rutas (`CytR`/citidina 32, `MetJ`/SAM 26, `TrpR`/Trp 17); (b) d=1,
el efector está a un paso de la química regulada (`PurR`/hipoxantina, `FrmR`/formaldehído);
(c) sin ruta: regulones sin genes metabólicos (`MarR`, `MprA`, sensores de xenobióticos),
donde GBA es ciego por construcción.
## 8.2 · Restricción real: no hay laboratorio

Sin banco, la pregunta "¿inductor o correpresor para este TF?" **no es decidible hoy**.
Se sondeó si la literatura puede sustituir la medida: consulta a Europe PMC por TF y
efector con términos cuantitativos (ITC, Kd, SPR, afinidad).

| resultado | valor |
|---|---:|
| TFs con ≥1 paper cuantitativo | **7/35** |
| TFs con ≥2 papers | 5/35 |

La medida que falta **tampoco está en la literatura** de forma recuperable a escala: es
cualitativa y dispersa. Por eso el campo depende de etiquetas curadas, y por eso las
cuatro vías computacionales fallaron igual.

**Qué queda haciendo solo simulación, por valor esperado:**

1. **Paquete de resultado negativo + benchmark** (0 GPU): los 36 decisivos con PMID,
   el evaluador de GBA, los cuatro baselines y los dos protocolos de co-folding. Es el
   entregable que sobrevive y el que nadie tiene.
2. **Generador de candidatos a escala** (reusa C0, sin GPU): aplicar el filtro GBA d=0
   a los ~99 represores con operador y sin efector conocido → hipótesis ordenadas con
   `status=hypothesis`. Útil para la anotación de los 400 genomas aunque no se puedan
   validar localmente.
3. **Ampliar etiquetas por curación** (sin GPU): revisiones de familia + CollecTF +
   RegPrecise para llevar de 36 a ~200. Hace más fuerte cualquier benchmark futuro.
4. **Cambiar el target a algo decidible in silico**: presencia de bolsillo de efector
   (binario), no su modo. Etiqueta parcial: estructuras PDB holo/apo y variantes AlloRep.

Lo que **no** se sostiene sin banco: seguir buscando un ΔΔG predicho que la herramienta
no puede dar. Cuatro instrumentos ya dieron el mismo resultado.

## 8.3 · Modos normales sobre AlphaFold (en Modal, GPU)

`06_anm_modal.py`: ANM (ProDy) sobre las 35 estructuras AFDB en un contenedor A10G
(~2 min, coste despreciable). Métrica: fracción de la varianza del modo más blando que
cae en el DBD, normalizada por la fracción de residuos que el DBD representa
(enriquecimiento; 1,0 = el modo mueve el DBD como al resto de la proteína).

| enriquecimiento | TF | familia | | enriquecimiento | TF |
|---:|---|---|---|---:|---|
| 3,50 | PurR | LacI/GalR | | 0,55 | TrpR |
| 3,46 | GalR | LacI/GalR | | 0,50 | **LacI** |
| 3,27 | TreR | LacI/GalR | | 0,45 | CytR |
| 3,12 | GntR | GntR | | 0,22 | GalS |
| 2,90 | BirA | HTH_11 | | 0,14 | AlsR |
| 2,17 | UlaR | DeoR | | 0,03 | FrlR |

**Autochequeo que falla:** `LacI`, el interruptor alostérico de manual, sale con
enriquecimiento 0,50 — el modo más blando **no** mueve preferentemente su DBD. O el
monómero de AFDB no es la unidad funcional (LacI es tetrámero y la transición es
cuaternaria), o el modo 1 no es el funcional. En cualquier caso, la métrica tal como
está no recupera el caso canónico y **no se puede usar como clasificador** hasta
resolver el oligómero y la selección de modo.

Siguiente paso obligatorio: montar el oligómero (plantilla PDB o AF-Multimer) y
repetir. Sin eso, el ANM del monómero mide la bisagra equivocada.

**Prueba de oligómero (2026-11-21).** `07_laci_oligomer.py` en Modal: tetrámero apo
1LBG vs unido a ADN 1LBI (1.184 CA comunes, RMSD del core 2,9 Å).

| criterio | monómero AF | tetrámero 1LBG | umbral |
|---|---:|---:|---|
| enriquecimiento DBD (modo 1) | 0,50 | **2,35** (2,1–4,1 en modos 1–5) | > 2 ✓ |
| coseno DBD entre cadenas | — | **−0,89 / −0,94** (modos 1–2) | < 0 ✗→✓ |
| solapamiento con la transición observada | — | 0,27 (mejor modo) | > 0,5 ✗ |
| solapamiento de subespacio (modos 1–5) | — | **0,453** | > 0,5 ✗ (borde) |

**Lectura:** el oligómero confirma la premisa mecanística — el modo blando del tetrámero
es DBD-dominado y antiparalelo (el "mariposa" clásico), y el monómero fallaba por eso.
Pero el ANM solo captura ~45 % de la dirección de la transición real: sirve para
**identificar el modo funcional**, no para generar conformaciones cuantitativas.

Consecuencia: usar el ANM como **clasificador de "¿tiene bisagra alostérica?"** (modo
DBD-dominado + antiparalelo + blando), no como generador de conformeros para comparar
bolsillos. Para conformeros cuantitativos haría falta MD con muestreo mejorado o un
ensamble aprendido (BioEmu/AlphaFlow), no NMA.

## 8.4 · Especificidad del detector de bisagra (control negativo)

`08_hinge_specificity.py`. Nulo interno por proteína: la ventana del DBD se compara
contra 200 ventanas aleatorias del mismo tamaño (z_dbd). Regla: bisagra si algún modo
tiene z_dbd > 2 **y** coseno antiparalelo entre cadenas < −0,5.

| proteína | rol | n | bisagra | z_dbd | anti (modo 1) |
|---|---|---:|---|---:|---:|
| LacI (1LBG) | positivo | 1428 | **sí** | 3,49 | −0,89 |
| λ cI (1LMB) | negativo | 179 | **sí** | 3,98 | −0,97 |
| LexA (1JHH) | negativo | 316 | **sí** | 2,70 | −0,90 |
| 434 Cro (3CRO) | negativo | 130 | sí (borde) | 2,02 | −0,43 |
| 434 cI (2OR1) | negativo | 126 | no | 1,12 | −0,82 |
| Fis (1F36) | negativo | 178 | no | 1,67 | −0,97 |
| IHF (1IHF) | negativo | 190 | no | −2,53 | +0,99 |
| SinR (1B0N) | negativo | 134 | no | 0,98 | +0,14 |
| λ Cro (5CRO) | negativo | 243 | no | −9 | +0,59 |

**Sensibilidad 1/1 · especificidad 5/8.** El detector **no** mide "interruptor por
efector": mide **bisagra entre el DBD y el cuerpo**. λ cI y LexA son de dos dominios con
bisagra (linker flexible, inactivación por proteólisis) y disparan correctamente como
bisagra, pero no son interruptores de metabolito. 434 Cro (z=2,02) está en el borde.

**Consecuencia:** renombrar la salida a `bisagra_DBD_core` (no "alostérica"). Es condición
necesaria pero no suficiente para regulación por efector. Para discriminar hace falta un
criterio adicional: que la bisagra **abra/cierre una cavidad** (bolsillo diferencial a lo
largo del modo). λ cI y LexA tienen bisagra sin bolsillo de efector; LacI tiene bisagra
que abre el bolsillo de alolactosa. Ese es el siguiente test.

## 8.5 · Uso de Modal (regla operativa)

**Modal admite 10 GPU simultáneas por workspace.** Este módulo corre todo en A10G con
concurrencia ≤3. Regla: avisar antes de usarlo y cerrar al terminar.

```bash
bash scripts/allostery/modal_guard.sh status     # 0 activas, exit 0
bash scripts/allostery/modal_guard.sh stop-all   # cerrar todo
```

Detalle, convenciones e inventario de scripts: [`MODAL.md`](MODAL.md).

## 8.6 · Gating de cavidad: descartado (control positivo inválido)

`09_cavity_gating.py`: se generan dos conformeros a lo largo del modo bisagra
(±2,5 Å de desplazamiento máximo, extrapolación lineal) y se corre P2Rank en cada uno.

| proteína | rol | P2Rank en estructura original | bolsillos en ±2,5 Å | veredicto |
|---|---|---|---|---|
| **LacI** | positivo | 1 CSV, bolsillos encontrados | **0 ligandables** | **control inválido** |
| λ cI | falso positivo | — | 1 nuevo | no interpretable |
| LexA | falso positivo | — | 3 nuevos | no interpretable |
| 434 Cro | falso positivo | — | 0 | — |
| Fis | negativo | — | 0 | — |

**El control positivo falla por la razón equivocada:** P2Rank sí encuentra bolsillos en la
estructura original de LacI (22.916 puntos SAS), pero **cero ligandables** en los
conformeros desplazados. La extrapolación lineal a lo largo de un modo armónico produce
geometrías no físicas (colisiones locales y superficie colapsada), así que la comparación
de bolsillos no mide gating: mide el artefacto del desplazamiento.

**Decisión: se descarta el criterio de gating.** Para comparar bolsillos entre estados
hacen falta conformeros físicamente válidos (MD con muestreo mejorado o ensamble
aprendido tipo BioEmu/AlphaFlow), no extrapolación de NMA. El detector de bisagra se queda
como está: `bisagra_DBD_core`, sensibilidad 1/1, especificidad 5/8, con su modo de fallo
documentado (bisagras de linker sin bolsillo).

## 8.7 · ¿El modo bisagra deforma el bolsillo anotado? No

`10_pocket_patch.py` usa los sitios de unión de UniProt (no P2Rank) y mide la deformación
interna del parche del bolsillo (cambio de la matriz de distancias) a lo largo del modo
bisagra, contra 200 parches aleatorios compactados igual.

| TF | modo | n parche | deformación | azar | z | pasa |
|---|---|---:|---:|---:|---:|---|
| FrmR | 1 (sin DBD Pfam) | 11 | 0,1819 | 0,0371 | **11,6** | artefacto |
| BirA | 3 | 23 | 0,0791 | 0,0588 | 0,6 | no |
| PurR | 3 | 31 | 0,0033 | 0,0081 | −0,3 | no |
| TreR | 3 | 48 | 0,0054 | 0,0155 | −0,3 | no |
| NemR | 5 | 5 | 0,0393 | 0,0522 | −0,3 | no |
| DgoR | 1 | 18 | 0,0012 | 0,0083 | −0,4 | no |

**Interpretación mecanística:** el bolsillo del efector está en el dominio regulador
(core); el modo bisagra mueve el DBD **respecto** al core y deja el core prácticamente
rígido (z ≈ 0 en PurR y TreR, los mejor anotados). Son **dos movimientos distintos**:

- bisagra DBD–core → transición de unión a ADN;
- gating del bolsillo → movimiento interno del core, inducido por el ligando.

`FrmR` es artefacto: 91 aa, sin dominio Pfam asignado, así que el "modo bisagra" es
simplemente el modo más blando de una proteína pequeña. `AlaS` se descarta (es la
alanil-tRNA sintetasa, no un TF).

**Consecuencia:** el detector de bisagra se queda como está y con su alcance: detecta la
bisagra DBD–core. No localiza ni gatea el bolsillo del efector. Para eso haría falta la
transición inducida por ligando (estructura holo o co-folding), no los modos apo.

## 8.8 · Capa 1 (afinidad por secuencia): ingesta hecha, dos bloqueos semánticos

`11_regseq_ingest.py` descarga y parsea las matrices de energía de Reg-Seq.

| Dato | Valor |
|---|---:|
| Matrices | **996** (de 1.000 listadas) |
| Promotores distintos | 79 |
| Tamaño en disco | **986 KB** |
| Tiempo de cómputo | ~2 min (descarga), segundos de parseo |
| Amplitud máxima (ΔΔG) | **0,43** |

**Bloqueo 1 — unidades.** El README dice "k_BT units", pero la amplitud máxima de las 996
matrices es 0,43. Un sitio de TF real tiene ΔΔG de varios kT (10-20 entre específico y
no específico). Estos archivos contienen los **coeficientes de regresión sobre expresión**
(efecto por base), no ΔΔG absolutos. Para kcal/mol hace falta la conversión termodinámica
de *mpathic* (modelo de promotor), que no está implementada.

**Bloqueo 2 — identidad del TF.** No hay tabla legible promotor→TF en el repositorio
(las figuras son gráficos). El mapeo por RegulonDB da 177 TFs únicos, 278 ambiguos y
**541 sin TF** (54 %). El mapeo por motivo falla porque los PWM de RegulonDB miden 30-40 pb
y las ventanas de matriz 8-50 pb, y el "consenso" de una ventana corta es ruido.

**Resuelto el mapeo (2026-11-21).** El suplemento 3 del paper
(`elife-55308-supp3-v2.xlsx`) trae los sitios de TF (`gene, start, end, type, identity`).
El fallo era de **sistema de coordenadas**: las matrices usan coordenadas del amplicón y
el suplemento relativas al TSS. Alineando `prior_designs/wtsequences.csv` contra el genoma
se obtiene el desplazamiento por promotor: **−113 pb (hebra +) y −116 pb (hebra −)**,
consistente en 107 promotores.

| Resultado del mapeo | n |
|---|---:|
| Matrices con sitio solapante | 207 |
| Con TF único | **174** |
| Ambiguas (2+ TFs) | 33 |
| TFs distintos | 14 |

Top: GlpR 46, YgbI 26, NsrR 22, DicA 16, FNR 15, YdfH 14, IlvY 11, XylR 8.

**Mapeo por motivo: rechazado.** Recortando cada PWM a su núcleo de 14 pb y probando ambas
hebras, se obtienen 501 asignaciones "confiables" pero solo **5/291 coinciden** con
RegulonDB, con NarP asignado a 316 matrices. El consenso derivado de la matriz de energía
no es un motivo usable.

**Consecuencia:** la capa 1 no es un problema de cómputo (1 MB, segundos), es un problema de
**semántica de los datos**. Antes de escribir el modelo hay que: (a) fijar la convención de
unidades con el método del paper; (b) construir la tabla promotor→TF de los 79 promotores
(curación desde las figuras/suplemento) o recortar cada PWM a su núcleo de alta información
para que el match por motivo funcione.

## 8.9 · Capa 1, validación: Reg-Seq no sirve; PWM de RegulonDB sí (relativo)

`16_regseq_energy_model.py`: para los 173 pares matriz→TF se puntúa la secuencia wild-type
de la ventana contra 500 permutaciones de la misma composición.

| hebra | n | z mediano | frac. significativa (WT mejor) | frac. significativa (WT peor) |
|---|---:|---:|---:|---:|
| directa | 173 | **+0,23** | 0,26 | 0,30 |
| reverso-complementaria | 173 | +0,08 | 0,05 | 0,05 |

**La matriz no distingue la secuencia real de una barajada.** Con amplitudes de 0,21
(mediana) y 0,43 (máxima), coherente con coeficientes de expresión y no con ΔΔG. Reg-Seq
queda **descartado como fuente de K_O** salvo que se rehaga con su pipeline original
(`mpathic`, py3.6, datos crudos de secuenciación).

**Alternativa validada: PWM de RegulonDB recortados a núcleo de 14 pb.**

| métrica | valor |
|---|---:|
| TFs evaluados | 62 |
| z mediano (sitio real vs barajado) | **1,8** |
| fracción con z > 2 | **0,42** |
| mejores | EvgA 5,8 · MurR 5,5 · GntR 5,3 · RstA 4,9 · PurR 4,8 |

Discrimina, aunque con margen modesto. Sirve como **capa de afinidad relativa** (percentil
frente al modelo PWM), no como kcal/mol.

**Decisión:** la capa 1 se construye sobre PWM de RegulonDB (82 TFs), no sobre Reg-Seq. La
unidad será relativa y así se reportará; la conversión a kcal/mol queda fuera de alcance.

## 8.10 · Capa 1 construida: tabla de afinidad relativa (955 pares)

`17_affinity_table.py`. Para cada par (represor, sitio) con confianza S/C y evidencia
experimental: score PWM del núcleo de 14 pb, percentil frente a 10.000 ventanas del genoma
y z frente a 200 permutaciones del propio sitio.

| Métrica | Valor |
|---|---:|
| Pares evaluados | **955** |
| TFs cubiertos | **55** |
| **Percentil mediano frente al genoma** | **0,987** |
| Fracción en el top 5 % del genoma | 0,628 |
| z mediano frente a barajado | 1,6 |
| Fracción con z > 2 | 0,43 |
| TFs con percentil mediano ≥ 0,95 | **32 / 55** |

Mejores: GalS, MurR, Cra, MetR, PdhR, ArgR (percentil 1,0).

**Fallo documentado:** `AscG` (percentil 0,04). Su PWM consenso
`TTCACGTGAAATCACATTCAAAAAAA` tiene una cola poli-A; el recorte automático a "ventana de
máxima información" eligió la cola en lugar del núcleo `TGAAAC`, y los sitios reales
puntúan bajo. Afecta a 1 de 55 TFs; la corrección es penalizar cores de baja complejidad.

**Lectura:** la primera capa del objetivo queda construida y validada — *"este represor
reconoce esta secuencia"* con un número (percentil) y un nulo. La unidad es relativa; no es
kcal/mol y así se reporta.

## 8.11 · Prototipo de pipeline y GUI (genoma crudo → represores)

`18_genome_pipeline.py` + `app_alloreg.py`. Flujo completo sobre un FASTA sin anotar:

| Paso | Método | Salida |
|---|---|---|
| 1 | pyrodigal (meta) | genes + proteínas |
| 2 | hmmsearch, 299 HMM de dominio | candidatos TF + familia |
| 3 | Jaccard de 4-mers vs 243 TFs de referencia | homólogo + efector heredado (≥0,35) |
| 4 | PWM de 14 pb de la familia, región −250 pb | sitios + percentil |
| 5 | informe | TSV + Markdown + JSON |

Probado sobre 200 kb de *E. coli*: 170 genes, 8 candidatos, 5 con efector heredado
(araC→D-fucosa, pdhR→piruvato, cra→fructosa-1-P, leuO→2-isopropilmalato, nhaR→Na⁺),
3 `unknown`. Sitios restringidos a la familia (AraC→AraC 0,995; GntR→PdhR 1,0;
LacI/GalR→Cra 0,992).

Estados: `inherited_homolog`, `family_prior`, `unknown`. Ningún `confirmed`.

Documentación, limitaciones y qué haría falta para producción: [`PIPELINE_GENOMA.md`](PIPELINE_GENOMA.md).

## 9 · Fuera de alcance

- Las tres torres y el DNA-LLM desde cero (`G1` del BACKLOG sigue pausado; ProkBERT ya
  mide mejor, ver `NEW_MODEL_DESIGN.md` §0).
- Diseño de novo (RFdiffusion/ProteinMPNN): anterior a tener el módulo de efector validado.
- FBA y metabolismo: sin cambios.
