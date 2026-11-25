# Framework Represor–Operador–Efector: qué falta para aplicarlo

> **Documento de evaluación** del marco descrito en
> [`FRAMEWORK_INFERENCIA_ACOTADA.md`](FRAMEWORK_INFERENCIA_ACOTADA.md).
>
> **Fecha:** 2026-11-24 · **Ámbito:** repositorio `DNA_TOKEN_BAC` / `dna_poc`.
>
> **Actualización (2026-11-24):** la dirección del proyecto es un **censo de
> todos los represores**, no la revisión de un caso. Ver
> [`CENSUS_REPRESSORES.md`](CENSUS_REPRESSORES.md). El piloto de 20 casos de este
> documento se lee como **calibración**, no como alcance final. No se requiere
> cristalizar proteínas; las estructuras existentes son evidencia auxiliar.
>
> **Veredicto corto:** el marco es correcto como **programa de priorización**,
> pero no es un pipeline aplicable tal como está. De sus cuatro módulos, **M2 y
> M4 se pueden aplicar ya** con datos existentes; **M1 requiere una pieza que no
> existe** (estructura → motivo validado); **M3 no es aplicable con los recursos
> actuales** y ya fue cerrado en v1.0 tras cuatro intentos negativos. La línea
> correcta es: aplicar M2 + M4 en modo hipótesis, construir M1 con un piloto
> acotado, y no prometer dirección/magnitud del efector sin GPU o banco.

---

## 0. Resumen ejecutivo

### 0.1 Qué se puede afirmar hoy

El proyecto ya tiene la **mitad barata** del marco:

- 400 genomas anotados y una capa de regulación L2/L3 en *E. coli* K-12.
- 82 PWM de *E. coli*, 2.601 operones predichos, 65 aristas TF→operón (8 TF modelados; 7 con aristas).
- Un prototipo genoma → represores con **92 % de recall** contra los 243 TF de
  RegulonDB y ~22 s por genoma de 2 Mb.
- Una tabla de **955 pares (TF, sitio)** con afinidad relativa por PWM.
- Un *benchmark* de **36 represores decisivos** y **43 con efector**.
- 35 modelos AlphaFold de los TF del *benchmark* y **dos complejos TF–ADN**
  experimentales en el piloto.

### 0.2 Qué no se puede afirmar hoy

- Que una estructura 3D predicha permita derivar un PWM específico para un
  represor bacteriano cualquiera. **No hay motor integrado ni benchmark propio.**
- Que el operador top-1 sea el funcional *in vivo*. Solo hay ranking contra un
  PWM y un percentil genómico.
- Que el efector top-1 sea correcto. En los 5 genomas no modelo probados la
  cobertura de efector es **0**; la herencia por homología no escala.
- Que se pueda predecir **dirección** (inductor/correpresor) o **magnitud**
  (ΔΔG). Cuatro métodos fallaron y el criterio de éxito no se cumplió.
- Que el pipeline escale a 400 genomas con estructuras y DM. El coste crece
  varias órdenes de magnitud.

### 0.3 Decisión recomendada

| Módulo | ¿Aplicar ahora? | Condición |
|---|---|---|
| **M2** Motivo → operador | **Sí, con correcciones** | Arreglar el escaneo, unirlo a operones y calibrar en K-12 |
| **M4** Integración / output | **Sí** | Esquema único, procedencia y confianza |
| **M1** Estructura → motivo | **Opcional, no requisito** | El censo usa motivos de familia; estructura solo para casos sin motivo |
| **M3** Bolsillo → efector | **No como predictor** | Solo generador de candidatos; dirección/magnitud requieren GPU/banco |

---

## 1. Inventario verificado: qué existe ya

### 1.1 Datos y cobertura

| Pieza | Dato actual | Ruta / fuente | Lectura |
|---|---:|---|---|
| Genomas anotados | 400 | `data/parsed_annotations/`, `data/regulation/` | Base del proyecto |
| TF con familia declarada | 14.850 | `data/parsed_annotations/` | Censo amplio, no curado |
| CDS anotados "repressor" | 1.555 | `data/parsed_annotations/` | Etiqueta de producto, no función probada |
| TF de RegulonDB | 243 | `data/regulondb/TFSet.tsv` | Verdad de referencia K-12 |
| TF con efector anotado | 93 (38 %) | `data/regulondb/TFSet.tsv` | Etiqueta escasa |
| Represores decisivos | 36 (27 inductor / 9 correpresor) | `data/allostery/` | *Benchmark*, no entrenamiento |
| Represores con efector | 43 | `data/allostery/c0_gba.md` | Base del GBA |
| Interacciones TF–promotor con secuencia | 3.444 (216 TF) | `data/regulondb/TF-RISet.tsv` | Solo K-12 |
| PWM de TF | 82 | `data/regulondb/gt/tf_pwms.tsv` | Solo *E. coli* |
| Pares (TF, sitio) con afinidad relativa | 955 (55 TF) | `data/allostery/affinity_table.tsv` | Percentil, no kcal/mol |
| Operones K-12 | 2.601 | `data/regulation/ecoli_k12/operons.tsv` | Predichos |
| Aristas TF→operón | 65 (7 TF con aristas; 8 modelados) | `data/regulation/ecoli_k12/regulon.json` | L3 parcial |
| Reglas lógicas | 63 | `data/regulation/ecoli_k12/expression_rules.tsv` | Muchas `assumed` |
| Genomas no modelo con represores | 5 | `data/allostery/campaign/` | 85 represores |
| Efectores en esos 5 genomas | **0** | `data/allostery/campaign/agregado.tsv` | El cuello real |

### 1.2 Estructuras

| Pieza | Dato actual | Ruta | Limitación |
|---|---:|---|---|
| Modelos AlphaFold de TF del *benchmark* | 35 | `data/allostery/work/af/*.pdb` | Monómero, una conformación |
| Complejos TF–ADN experimentales | 2 | `data/allostery/pilot/PurR_hipoxantina/original_1pnr.pdb`, `.../TrpR_Trp/original_1tro.pdb` | Solo dos sistemas |
| Holo TF–ligando sin ADN | 1 | `data/allostery/pilot/LacI_IPTG/original_1lbh.pdb` | No es complejo con operador |
| Estructuras AFDB de otros fines | 1.670 | `data/metabolism/targets/afdb25k/` | No son represores; no sirven para M1 |
| Mapeo `locus_tag` → UniProt → PDB/AFDB | **No existe** como tabla | — | Bloquea el uso masivo |

**Conclusión:** no hay un *structure store* de represores. Hay 35 modelos para el
*benchmark* y dos complejos con ADN. El marco necesita miles de estructuras y,
sobre todo, **complejos proteína–ADN** o plantillas homólogas.

### 1.3 Pipeline y código

| Pieza | Estado | Ruta |
|---|---|---|
| Detección de DBD por HMM | Funciona | `scripts/allostery/18_genome_pipeline.py` |
| Familia por Pfam → familia canónica | Funciona | `_common.py`, `FAMILY_MAP` |
| Homología DIAMOND vs RegulonDB/Swiss-Prot | Funciona | `18_genome_pipeline.py` |
| Herencia de efector | Funciona, cobertura baja | `18_genome_pipeline.py` |
| Escaneo PWM | **Prototipo con un bug** | `18_genome_pipeline.py:545-561` |
| Operones L2 | Funciona | `data/regulation/`, `scripts/regulation/` |
| Regulón L3 | Parcial (8 TF modelados; 7 con aristas) | `scripts/regulation/13_build_regulon.py` |
| Docking Vina | Existe como *benchmark* | `scripts/allostery/03_pocket_specificity.py` |
| ANM / bisagra | Experimental | `scripts/allostery/06_anm_modal.py` |
| GUI de alostería | Prototipo | `scripts/allostery/app_alloreg.py` |
| Tests del pipeline 18 | **No hay** | `tests/` no lo cubre |

**Bug concreto:** en `18_genome_pipeline.py`, el escaneo de sitios está dentro de
`if len(sub_pwms) > 8`, así que **no se ejecuta** cuando la familia tiene 8 PWM o
menos. Es un fallo de indentación, no una decisión de diseño. Antes de cualquier
benchmark hay que corregirlo y añadir un test que verifique que un TF con una
sola PWM produce sitios.

---

## 2. Qué falta, módulo por módulo

### 2.1 MÓDULO 1 — Estructura → motivo

| Brecha | Por qué bloquea | Qué haría falta | Esfuerzo |
|---|---|---|---|
| **Manifest de estructuras** | No se sabe qué represor tiene PDB/AFDB, ni su cobertura, oligómero o pLDDT | Tabla `locus_tag → UniProt → PDB/AFDB → rango DBD → oligómero → calidad` | M |
| **Complejos proteína–ADN** | Sin complejo no hay geometría de interfaz; AF3/RF2NA no están integrados | PDB experimentales + modelado por homología + AF3/RF2NA para los que no tengan plantilla | L |
| **Motor estructura→secuencia** | DeepPBS/NA-MPNN/IDEA no están integrados en el pipeline | Contenedor con el motor, entrada PDB, salida PWM, versión fijada | L |
| **Benchmark propio** | No se sabe si DeepPBS supera al PWM de familia en bacterias | 20–30 TF de K-12 con sitio curado; comparar PWM predicho vs RegulonDB | M |
| **Rosetta / ΔΔG** | No instalado; licencia académica | Protocolo DNA-interface solo para top candidatos | M |
| **Incertidumbre** | Los motores devuelven un score, no una calibración | Calibrar contra sitios reales vs barajados; reportar intervalo | M |

**Realidad:** DeepPBS y NA-MPNN son herramientas de investigación, entrenadas
mayoritariamente con estructuras del PDB. Su transferencia a represores
bacterianos diversos **no está demostrada**. El proyecto ya tiene un baseline
honesto: PWM de RegulonDB con z mediano 1,8 y solo 42 % de TF con z > 2. Cualquier
motor nuevo debe **superar ese baseline** en K-12 antes de usarse en los 400.

### 2.2 MÓDULO 2 — Motivo → operador

| Brecha | Por qué bloquea | Qué haría falta | Esfuerzo |
|---|---|---|---|
| **Escaneo en operones diana** | El código actual escanea −250 pb del **gen del TF**, no los promotores de los operones que regula | Unir `regulon.json`/`operons.tsv` con el escaneo; una fila por (TF, operón diana) | M |
| **FIMO / p-valor** | El score actual es un percentil de un núcleo de 14 pb; no hay p-valor ni corrección por múltiples pruebas | FIMO o implementación equivalente con fondo genómico y FDR | M |
| **Umbrales** | No hay umbral calibrado por familia | Curva precisión–recall contra RegulonDB/CollecTF; umbral por familia | M |
| **Contexto COMMBAT** | No hay *target score* ni coherencia funcional integrada | Reglas de contexto: distancia, orientación, sintenia, función del operón | M |
| **Sintenia** | Hay 39.407 pares de sintenia, pero no se usan para priorizar operadores | Cruzar sitio con conservación del par TF–operón en genomas cercanos | M |
| **Negativos** | No hay controles de motivo barajado, PWM de otra familia ni región no promotora | Incluir controles en el *benchmark* y en la salida | S |
| **Bug de escaneo** | Suprime sitios en familias con ≤8 PWM | Corregir indentación y test de regresión | S |

**Lo que sí se puede hacer ya:** con los 82 PWM y los 2.601 operones de K-12, se
puede construir un escaneo operón-diana y medir recall@k contra los 3.444 sitios
de RegulonDB. Es CPU pura, sin GPU, y es el camino más corto a un resultado
aplicable.

### 2.3 MÓDULO 3 — Bolsillo → efector

| Brecha | Por qué bloquea | Qué haría falta | Esfuerzo |
|---|---|---|---|
| **Bolsillo a escala** | fpocket/SiteMap no están integrados; solo hay un *benchmark* de docking | Pipeline de detección de bolsillo + filtro por distancia al DBD | M |
| **Biblioteca de ligandos** | No hay biblioteca contextual por genoma con ChEBI/InChIKey | Metabolitos de Rhea/KEGG + iones + xenobióticos, deduplicados por InChIKey | M |
| **Docking** | Vina está en `.venv`, pero no en el pipeline; el *benchmark* no mostró señal (top-1 0,31 vs 0,25) | Docking contextual + controles por ligando; no usar como afinidad absoluta | M |
| **Acoplamiento alostérico** | ANM mide bisagra, no gating; el bolsillo no se deforma en el modo bisagra | DM comparativa apo/holo o red de correlación; requiere GPU | L |
| **Dirección** | 1/6 y 2/6 signos correctos en Boltz-2; sin banco no es decidible | Kd(operador, ± efector) medido, o un *benchmark* externo con ΔΔG | XL |
| **Magnitud** | No hay ΔΔG medido; Reg-Seq no sirve como ΔG | Datos experimentales o GPU masiva | XL |

**Realidad:** la línea se cerró en v1.0 con un criterio explícito:
`ALOSTERIA_V1_CIERRE.md` documenta que los criterios 1 y 2 (sitio y ranking)
se cumplieron, y los criterios 3 y 4 (dirección y magnitud) no. Cuatro métodos
—secuencia, docking ciego, co-folding Boltz-2 y gating de cavidad— dieron
resultados negativos o inválidos. Volver a intentarlo sin datos nuevos es
repetir el mismo experimento.

### 2.4 MÓDULO 4 — Integración

| Brecha | Por qué bloquea | Qué haría falta | Esfuerzo |
|---|---|---|---|
| **Esquema único** | Cada script tiene su propio TSV; no hay una tabla (TF, operador, efector) con procedencia | Contrato de columnas + validación + hash | S |
| **Pesos** | No hay pesos calibrados para \(S_{\text{sitio}}\) ni \(S_{\text{efector}}\) | Ajustar en K-12 con validación cruzada; reportar sensibilidad | M |
| **Confianza** | No hay regla explícita Alta/Media/Baja | Definir por evidencia: estructura experimental > homología > familia; acuerdo entre métodos | S |
| **Incertidumbre** | No se reporta rango ni sensibilidad | Bootstrap o barrido de pesos; intervalos por candidato | M |
| **GUI** | La app de alostería es un prototipo separado | Pestaña en la app principal, con procedencia y límites | M |
| **Orquestación** | `run_pipeline.sh` no tiene etapa de inferencia | Etapa `inferencia` idempotente, con caché y manifiesto | S |
| **Tests** | No hay tests del pipeline 18 | `unittest` + casos dorados + controles negativos | S |

### 2.5 Brechas transversales

| Brecha | Detalle | Impacto |
|---|---|---|
| **Ground truth insuficiente** | 36 decisivos, 43 con efector, 955 pares relativos; todo K-12 | Sin *held-out* no se puede afirmar generalización |
| **Cómputo** | 8 núcleos, 23 GB RAM, sin GPU; MD local 0,25 ns/día | M1 a escala y M3 son inviables localmente |
| **Licencias** | RegulonDB académico; AF3/Rosetta/PREDetector/COMMBAT con términos propios | Verificar antes de distribuir |
| **Reproducibilidad** | Versiones de HMM, PWM, estructuras y motores no están fijadas en un manifiesto | Un resultado no es reproducible si cambia la base |
| **Sesgo de familia** | 27 inductores / 9 correpresores; 6/14 familias con un solo caso | Cualquier clasificador necesita LOFO y métricas por clase |
| **Anotación de producto** | 1.555 "repressor" es etiqueta de nombre, no función | El censo no es una verdad biológica |

---

## 3. Las cuatro brechas críticas

Si solo se pudieran cerrar cuatro cosas, serían estas, en orden:

### Brecha 1 — No hay puente estructura → motivo validado

El marco asume que DeepPBS/NA-MPNN convierten un complejo en un PWM. Ese puente
no está construido ni validado en este proyecto. Sin él, M1 no existe y M2 se
queda con PWM de familia.

**Cierre mínimo:** 20 represores de K-12 con estructura y sitio curado; correr
DeepPBS/NA-MPNN; comparar PWM predicho vs PWM de RegulonDB; medir recall@5 y
correlación de scores. Si no supera al PWM de familia, M1 no se justifica.

### Brecha 2 — El escaneo no está unido a los operones diana

El prototipo escanea el promotor del propio gen del TF. El marco pide escanear
los promotores de los operones regulados. Son cosas distintas y la segunda no
está implementada.

**Cierre mínimo:** join `regulon.json` + `operons.tsv` + `promoters.tsv`; una
tabla (TF, operón, sitio, score, distancia, orientación); benchmark contra
RegulonDB.

### Brecha 3 — El efector no tiene termodinámica

Hay candidatos por homología y por GBA, pero no hay ΔG, ni dirección, ni
magnitud. El docking no mostró señal; el co-folding no recuperó el signo; la DM
es inviable localmente.

**Cierre mínimo:** no intentar dirección/magnitud. Entregar **conjunto de
candidatos priorizados** con GBA d=0 + bolsillo + docking, y declarar
`direction=unknown`. Si se quiere dirección, presupuestar GPU/banco.

### Brecha 4 — No hay calibración ni *held-out*

Los scores no están calibrados contra verdad externa. El 92 % de recall es de
detección de TF, no de operador ni de efector. Sin *held-out*, el pipeline puede
estar memorizando K-12.

**Cierre mínimo:** dividir por familia y por genoma; validar en K-12 y en al
menos 3 genomas con datos de CollecTF/DBTBS/RegPrecise; reportar precisión,
recall y calibración.

---

## 4. Lo que se puede aplicar ya (sin GPU)

Estas piezas usan datos y cómputo disponibles. No requieren AF3, Rosetta ni DM.

### 4.1 Escaneo operón-diana con PWM existentes

- Entrada: `data/regulation/<genoma>/operons.tsv`, `promoters.tsv`,
  `data/regulondb/gt/tf_pwms.tsv`.
- Salida: `tf_operator_candidates.tsv` con (TF, operón, coordenadas, hebra,
  secuencia, score PWM, percentil, distancia al TSS, orientación).
- Método: PWM scan con FIMO o implementación equivalente; fondo genómico;
  FDR; controles de motivo barajado y PWM de otra familia.
- Validación: K-12 contra `TF-RISet.tsv` (3.444 sitios).
- Coste: CPU, minutos por genoma.

### 4.2 Generador de candidatos a efector por GBA

- Entrada: regulón del TF + red metabólica Rhea/EC.
- Método: metabolitos a distancia 0 del regulón; excluir cofactores ubicuos;
  deduplicar por InChIKey; ranking por frecuencia y coherencia.
- Salida: top 10 candidatos por TF con `status=hypothesis`.
- Validación: los 43 casos con efector; recall@5 y enriquecimiento contra nulo.
- Coste: CPU, segundos por TF.

### 4.3 Bolsillo y docking para un puñado de casos

- fpocket/P2Rank sobre 35 estructuras AF; Vina sobre 10–20 ligandos por TF.
- No usar la energía como afinidad; usarla solo para ordenar y para descartar
  poses con choques.
- Validación: pose cristalográfica vs *decoy* (ya hecho para PurR/LacI).
- Coste: CPU, horas para decenas de casos.

### 4.4 Integración y salida

- Esquema único, confianza, procedencia, `status`.
- GUI de solo lectura con la tabla y los límites.
- Coste: días.

---

## 5. Lo que no se puede aplicar con los recursos actuales

| Pieza | Por qué | Qué haría falta |
|---|---|---|
| AF3/RF2NA para miles de complejos | Sin GPU local; pesos y términos de uso | Presupuesto GPU + licencia/acuerdo |
| DM 100–500 ns por sistema | 0,25 ns/día local; 100 ns = 400 días | GPU (A10G/H100) o créditos |
| Dirección inductor/correpresor | Cuatro métodos fallaron; sin banco no es decidible | Kd(operador, ± efector) o *benchmark* externo |
| Magnitud ΔΔG | No hay datos medidos; Reg-Seq no es ΔG | Datos experimentales o GPU masiva |
| 400 genomas con M1+M3 | Coste y cobertura de estructuras | Censo: todos en CPU; M1+M3: K-12 y 5–10 genomas, no 400 con estructura/DM |
| Rosetta ΔΔG | No instalado; licencia académica | Instalación + validación |

---

## 6. Piloto mínimo viable

### 6.1 Alcance

> **Esto es calibración, no alcance.** El censo se ejecuta sobre **todos** los
> represores y **todos** los genomas; este piloto solo mide si el censo acierta.

**Fase A — K-12, 20 represores con verdad conocida.** Elegir aquellos con:
- sitio curado en RegulonDB (confianza S/C);
- efector anotado en una sola conformación;
- estructura PDB/AFDB disponible;
- familia diversa (LacI/GalR, TetR, GntR, MarR, Fur, LysR, AraC…).

**Fase B — 3 genomas no modelo.** Elegir genomas con datos de CollecTF,
DBTBS o RegPrecise; no los 5 de la campaña, que tienen 0 efectores.

**Fase C — solo si A y B pasan:** extender a 50–100 represores y luego decidir
si los 400 tienen sentido.

### 6.2 Fases y entregables

| Fase | Qué se hace | Entregable | Gate |
|---|---|---|---|
| **P0. Datos** | Manifest de estructuras, oligómeros, ligandos, sitios y efectores | `tf_structure_manifest.tsv`, `ligand_library.tsv`, `ground_truth.tsv` | Cobertura ≥80 % de los 20 |
| **P1. Operador** | Escaneo PWM operón-diana + FIMO + controles | `tf_operator_candidates.tsv` | Recall@5 ≥0,70; precisión@5 ≥0,50 |
| **P2. Motivo** | DeepPBS/NA-MPNN en los 20; comparar con PWM | `pwm_comparison.tsv` | Supera PWM de familia en ≥60 % de casos |
| **P3. Efector** | GBA d=0 + bolsillo + docking; top 10 | `tf_effector_candidates.tsv` | Recall@5 ≥0,50; enriquecimiento ≥5× |
| **P4. Integración** | Score compuesto, confianza, GUI, informe | `inference_report.md`, `predictions.tsv` | Reproducible; hash; sin `confirmed` |
| **P5. Generalización** | 3 genomas no modelo | `heldout_report.md` | No caer >20 puntos vs K-12 |

### 6.3 Criterios de aceptación

- [ ] El pipeline corre de punta a punta con un comando y sin red.
- [ ] Cada fila tiene `evidence`, `confidence`, `status` y `provenance`.
- [ ] Los controles negativos están incluidos (motivo barajado, PWM de otra
      familia, ligando *decoy*, región no promotora).
- [ ] El benchmark es *held-out* por familia y por genoma.
- [ ] El informe declara lo que falla, no solo lo que acierta.
- [ ] Ningún resultado se etiqueta `confirmed` sin experimento.

### 6.4 Esfuerzo orientativo

| Fase | Esfuerzo | Dependencias |
|---|---|---|
| P0 | 1–2 semanas | UniProt/PDB/AFDB, ChEBI |
| P1 | 1–2 semanas | FIMO o equivalente, operones |
| P2 | 2–4 semanas | DeepPBS/NA-MPNN, estructuras |
| P3 | 2–3 semanas | fpocket/Vina, Rhea |
| P4 | 1–2 semanas | esquema, GUI |
| P5 | 1–2 semanas | datos externos |

**Total realista:** 8–12 semanas a tiempo parcial, sin contar la obtención de
licencias o GPU. La fase P2 es la que puede fallar; por eso va después de P1 y
con un gate de “supera al PWM”.

---

## 7. Backlog priorizado

### P0 — Ahora (sin GPU, alto valor)

| # | Tarea | Salida | Depende de |
|---|---|---|---|
| P0.1 | Corregir el bug de escaneo de `18_genome_pipeline.py` | Sitios para familias con ≤8 PWM | — |
| P0.2 | Añadir test de regresión del escaneo | `tests/test_genome_pipeline_scan.py` | P0.1 |
| P0.3 | Manifest `locus_tag → UniProt → PDB/AFDB` | `data/regulation/tf_structure_manifest.tsv` | UniProt/AFDB |
| P0.4 | Join operón-diana + PWM | `tf_operator_candidates.tsv` | operones, PWM |
| P0.5 | Benchmark operador K-12 | `operator_benchmark.md` | P0.4, `TF-RISet.tsv` |
| P0.6 | Generador GBA de efectores | `tf_effector_candidates.tsv` | regulón, Rhea |
| P0.7 | Esquema único + validador | `predictions_schema.md`, test | — |

### P1 — Después del piloto (CPU/GPU moderada)

| # | Tarea | Salida | Depende de |
|---|---|---|---|
| P1.1 | Integrar DeepPBS o NA-MPNN | `pwm_predicted.tsv` | estructuras, motor |
| P1.2 | Benchmark PWM predicho vs RegulonDB | `pwm_benchmark.md` | P1.1 |
| P1.3 | fpocket + Vina contextual | `pocket_candidates.tsv` | estructuras, ligandos |
| P1.4 | Controles negativos completos | `negative_controls.md` | P0.5, P1.3 |
| P1.5 | GUI de inferencia | pestaña en la app | P0.7 |

### P2 — Condicionado a GPU/banco

| # | Tarea | Salida | Depende de |
|---|---|---|---|
| P2.1 | AF3/RF2NA de complejos | `complex_models/` | GPU, licencia |
| P2.2 | DM apo/holo | `md_coupling/` | GPU |
| P2.3 | Dirección/magnitud | `direction_report.md` | datos experimentales |
| P2.4 | Extensión a 50–100 represores | `scale_report.md` | P1 completo |

---

## 8. Contrato de datos propuesto

Una tabla única, una fila por (TF, operador candidato, efector candidato), con
columnas mínimas:

```text
genome
tf_locus_tag
tf_name
tf_family
tf_uniprot
tf_structure_id
tf_structure_source        # pdb | afdb | model | none
tf_structure_plddt
tf_oligomer                # monomer | dimer | tetramer | unknown

operator_contig
operator_start
operator_end
operator_strand
operator_sequence
operator_score
operator_percentile
operator_pvalue
operator_target_operon
operator_target_genes
operator_distance_to_tss
operator_orientation      # divergon | same | opposite

effector_name
effector_chebi
effector_inchikey
effector_score
effector_evidence         # gba | pocket | docking | homology | prior
effector_direction        # unknown | inducer | corepressor

confidence                # high | medium | low
status                    # hypothesis | corroborated | documented | unknown
provenance_json
```

Reglas:

- `status=documented` solo con PMID de unión o estructura del complejo.
- `effector_direction=unknown` por defecto; no inferir signo sin dato.
- `confidence` no puede ser `high` si `tf_structure_source=none`.
- Toda fila debe poder rastrearse a un archivo y a una versión.

---

## 9. Riesgos y mitigaciones

| Riesgo | Probabilidad | Impacto | Mitigación |
|---|---|---|---|
| DeepPBS no generaliza a bacterias | Alta | M1 no sirve | Gate P2 contra PWM de familia |
| El operador top-1 no es funcional | Alta | Falso positivo | Benchmark K-12 + controles |
| El efector top-1 es un metabolito ubicuo | Alta | Ranking inútil | Excluir cofactores; GBA d=0 |
| El docking no tiene señal | Ya observado | M3 no sirve | No usar energía como afinidad |
| Sesgo de familia | Alta | Métricas infladas | LOFO + métricas por clase |
| Datos insuficientes | Alta | No hay *held-out* | Ampliar con CollecTF/DBTBS/RegPrecise |
| Cómputo | Alta | M3 inviable | GPU o limitar a 10–20 casos |
| Licencias | Media | No publicable | Revisar términos antes de distribuir |

---

## 10. Decisiones que hay que tomar

1. **¿Se quiere afirmar dirección?** Si sí, el proyecto necesita un banco
   (EMSA/ITC/reportero) o un *benchmark* externo con ΔΔG. Sin eso, la respuesta
   honesta es `direction=unknown`.
2. **¿Se escala a 400?** Sí: el censo cubre todos los genomas y todos los
   represores. Lo que no se escala sin calibrar es la **validación**; K-12 y los
   genomas *held-out* siguen siendo los gates.
3. **¿DeepPBS o PWM de familia?** Solo cambiar si el benchmark demuestra mejora.
   El PWM de familia es barato, interpretable y ya está.
4. **¿Se integra en la app principal o se queda como script?** La app es el
   producto; un script aislado no se usa.
5. **¿Qué se hace con los 1.555 "represores"?** Son etiquetas de producto. El
   censo real requiere HMM de familias, no nombres.

---

## 11. Respuesta corta

**Nos falta, en orden:**

1. **Un puente estructura → motivo validado** (manifest + complejos + DeepPBS/
   NA-MPNN + benchmark). Hoy no existe.
2. **Un escaneo de operadores diana** unido a los operones y calibrado contra
   RegulonDB. Hoy el prototipo escanea el gen del TF y tiene un bug.
3. **Termodinámica del efector** (ΔG, dirección, magnitud). Hoy no hay datos ni
   cómputo; cuatro métodos fallaron.
4. **Calibración y *held-out*** para no confundir recall de detección de TF con
   acierto de operador/efector.
5. **Infraestructura de producción**: esquema único, procedencia, tests,
   orquestación, licencias y GPU.

**Lo que sí se puede hacer ya:** M2 + M4 en K-12 con los 82 PWM y los 2.601
operones, y extender el escaneo a los 400 genomas con motivos de familia;
GBA d=0 para candidatos de efector; bolsillo/docking solo para 10–20 casos
difíciles; y un informe honesto con `status=hypothesis`. Eso no es el marco completo, pero
es la parte que hoy produce valor sin prometer lo que no se puede medir.
