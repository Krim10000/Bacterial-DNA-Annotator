# Revisión del pipeline maestro propuesto

> **Veredicto:** la arquitectura es buena como **generador de hipótesis
> jerárquico**: filtra barato primero, reserva el 3D para los pocos pares que
> pasan, y tiene rutas de rescate para los casos raros. Pero no es
> "determinista" ni "cero pérdida": varias reglas son heurísticas fuertes que
> necesitan calibración, y el 3D/RNAP no puede usarse como juez binario.
>
> La idea central se mantiene. Lo que hay que corregir es la **falsa
> precisión**: bits >9 = encendido, spacer 19–20 = MerR, choque 3D = represor,
> Fisher = anti-alucinación.

---

## 1. Lo que está bien

1. **Filtrado por capas.** Promotor → motivos → 3D → RNAP → metabolismo. Es la
   dirección correcta: lo barato filtra antes de lo caro.
2. **Arquitectura del promotor.** Separar −35/−10, extensión TG, UP element,
   spacer y sigma es más informativo que un único PWM.
3. **Cruce estructural.** Usar la geometría del complejo TF–ADN y la RNAP es
   una evidencia que el PWM solo no da.
4. **Rutas de rescate.** Looping, sigma alternativa, anti-activación,
   riboswitch/atenuador son mecanismos reales que un pipeline lineal pierde.
5. **Coherencia metabólica.** Cruzar el regulón con función es un prior útil.
6. **Salida como grafo.** Regulador → operón → rol → mecanismo es la forma
   correcta de representar la red.

---

## 2. Correcciones críticas

### 2.1 "Promotor >9 bits = ENCENDIDO" no es una verdad

El score en bits es una medida de **parecido al motivo**, no de actividad
transcripcional. La fuerza de un promotor depende de:

- disponibilidad de sigma;
- superenrollamiento;
- proteínas nucleoides (H-NS, Fis, IHF);
- competencia entre promotores;
- secuencia UP y α-CTD;
- contexto de crecimiento.

Un promotor "fuerte" puede estar apagado por un represor, por falta de sigma o
porque el gen no se expresa en esa condición. Un promotor "débil" puede estar
encendido por un activador o por un sigma alternativo.

**Corrección:**

- calcular log-odds con fondo genómico, no solo "bits";
- calibrar los umbrales contra TSS/RNA-seq de K-12;
- llamar al resultado `promoter_prior`, no `ON/OFF`;
- reportar percentil y rango, no una etiqueta binaria.

### 2.2 Operones por "distancia <45 pb" es demasiado rígido

La distancia intergénica no define un operón. Hay operones con gaps mayores,
genes solapados, ARN pequeños y promotores internos. El proyecto ya midió esto:

| Métrica L2 | Valor |
|---|---:|
| Pares génicos F1 | 0,808 |
| Conjunto exacto F1 | 0,629 |
| Operones 7+ F1 | 0,442 |
| Cobertura de promotores | 65,9 % |

Además, [`PWM_RECALIBRATION.md`](PWM_RECALIBRATION.md) muestra que la regla de
promotor **no mejora las fronteras de operón** ni siquiera con la matriz
correcta: el mejor resultado se obtiene cuando la regla deja de disparar. La
distancia hace casi todo el trabajo.

**Corrección:**

- usar `operons.tsv` con su `confidence`, no re-delimitar con 45 pb;
- si se proponen sub-operones, hacerlo solo con TSS interno, RBS, terminador o
  RNA-seq, y marcarlo `hypothesis`;
- no usar el promotor como frontera dura.

### 2.3 "Spacer 19–20 pb = MerR" es una pista, no una identidad

MerR sí tiene un spacer largo y tuerce el ADN, pero:

- no todos los promotores con spacer 19–20 son MerR;
- no todos los MerR tienen exactamente 19–20;
- la torsión depende de la secuencia y del complejo, no solo del largo.

**Corrección:**

- usar el spacer como **feature**;
- exigir motivo MerR, homología de familia y geometría 3D;
- reportar `merR_like`, no `MerR`.

### 2.4 Ventanas rígidas represor/activador

"Represor en [−45, +60]" y "activador en [−100, −35]" es una simplificación.
Los represores pueden unir aguas arriba, aguas abajo, en el espaciador o por
looping; los activadores clase I y II tienen posiciones distintas.

**Corrección:**

- escanear una ventana amplia (p. ej. −300 a +100) y puntuar la posición;
- usar la posición como evidencia, no como filtro binario;
- permitir que un TF sea represor y activador (dual).

### 2.5 "Palíndromos 14–22 pb" sesga a dímeros

Muchos represores unen motivos asimétricos, repeticiones directas o sitios
espaciados. Buscar solo palíndromos:

- pierde motivos de monómeros;
- pierde repeticiones directas;
- pierde sitios con espaciador variable;
- sobreajusta a familias tipo LacI/TetR.

**Corrección:**

- usar PWM/PFM de familia, no solo palíndromos;
- permitir motivos asimétricos y repeticiones directas;
- si se hace descubrimiento de novo, usar MEME/DREME con conservación y
  controles, no solo simetría.

### 2.6 El 3D y la RNAP no son un juez binario

Un choque estérico en una estructura estática no prueba represión. Los modelos
AF son una conformación, la DNA se deforma, la α-CTD es flexible y la RNAP
tiene múltiples estados. "Toca α-CTD" o "tuerce 30°" son continuos, no etiquetas.

**Corrección:**

- usar el 3D como **score**, no como veredicto:
  - `steric_clash_score`
  - `rnap_contact_score`
  - `dna_bend_angle`
  - `interface_confidence`
- calibrar contra casos conocidos de K-12 (LacI, GalR, CRP, AraC, MerR);
- usar ensemble/soft clash, no una sola pose;
- reportar mecanismo probable con incertidumbre, no una etiqueta única.

### 2.7 "Choque con RNAP = represor estérico" ignora otros mecanismos

Un represor puede actuar por:

- bloqueo estérico de la RNAP;
- bloqueo del promotor;
- looping;
- anti-activación;
- competencia con activador;
- remodelación del ADN.

El choque es una de varias explicaciones. Etiquetar `REPRESOR ESTÉRICO` sin
comparar alternativas es sobreinterpretar.

**Corrección:**

- generar varias hipótesis de mecanismo;
- puntuar cada una;
- si no se puede distinguir, `mechanism=unknown`.

### 2.8 Las rutas de rescate son buenas, pero no son comodines

- **Looping:** múltiplos de ~10,5 pb es una heurística; la fase helicoidal, la
  longitud del loop, la geometría del tetrámero y el superenrollamiento mandan.
- **Anti-activación/trampa:** "TF en −60 correlaciona negativo" no es evidencia
  suficiente; hay que mostrar el sitio solapado y el mecanismo.
- **Sigma alternativa:** requiere PWMs de σ54/σ32/σS/σ28 y, para σ54, bEBP e
  IHF; no basta con escanear.
- **Riboswitch/atenuador:** requiere UTR 5′, Rfam/Infernal y predicción de
  terminador; no es un fallback para cualquier promotor débil.

**Corrección:**

- implementar cada rescate como módulo opcional con evidencia explícita;
- no usarlo para forzar una explicación;
- si no hay evidencia, `unresolved`.

### 2.9 RNA-seq y lógica AND/OR

La lógica combinatoria no se deduce de una foto estática. Para decir "AND" o
"OR" hace falta:

- múltiples condiciones;
- perturbaciones (knockout/overexpression);
- réplicas;
- modelo de expresión.

Con RNA-seq de una condición, solo se puede decir "coherente/incoherente en esa
condición".

**Corrección:**

- usar RNA-seq como evidencia condicional;
- no inferir compuertas lógicas sin perturbación;
- marcar `logic_status=assumed` cuando corresponda.

### 2.10 Fisher/hipergeométrico no es un "filtro anti-alucinaciones"

El enriquecimiento de rutas es un prior, no una verdad. Tiene problemas de:

- fondo genómico;
- múltiples pruebas;
- regulones pequeños;
- anotación incompleta;
- regulones no metabólicos (LexA, MarR) que fallan por diseño.

**Corrección:**

- reportar p-valor y FDR;
- excluir rutas ubicuas;
- no descartar regulones no metabólicos;
- usar el enriquecimiento como una capa, no como juez final.

### 2.11 "Cero pérdida de información" es imposible

Todo umbral descarta. Lo que sí se puede exigir es **ninguna pérdida
silenciosa**: cada descarte debe tener `reason` y quedar en un bucket
`unresolved`/`discarded`.

**Corrección:**

- registrar cada par descartado y por qué;
- permitir auditoría;
- reportar cobertura y sesgo de los filtros.

### 2.12 Los "casos 1–6" no están definidos en el diagrama

El diagrama usa "Caso 1, 2, 3, 4, 5, 6" sin una tabla que los defina. Eso hace
imposible auditar la lógica. Propuesta de taxonomía, explícita y multi-etiqueta:

| Caso | Mecanismo | Evidencia esperada |
|---|---|---|
| 1 | Activador clase I (upstream, contacta α-CTD) | Sitio upstream + α-CTD + RNAP |
| 2 | Activador clase II (solapa −35, contacta σ) | Sitio en −35 + σ + RNAP |
| 3 | Activador MerR (torsión del spacer) | Motivo MerR + spacer 19–20 + geometría |
| 4 | Represor estérico (bloquea RNAP) | Sitio en promotor + choque con RNAP |
| 5 | DNA looping | Dos sitios + fase helicoidal + oligómero |
| 6 | Anti-activación / trampa | Sitio solapado con activador + efecto |
| σ | Sigma alternativa | PWM de σ + arquitectura de promotor |
| RNA | Riboswitch/atenuador | UTR 5′ + Rfam + terminador |
| Arquitectónico | H-NS/Fis/IHF | ADN curvo/AT-rico, sin motivo discreto |
| Desconocido | No se puede distinguir | `mechanism=unknown` |

Un TF puede tener más de un mecanismo. La salida debe ser un vector de scores,
no una etiqueta única.

### 2.13 "Fase 0/1 en segundos" y "3D solo para el 10 %"

- Puntuar promotores: segundos.
- Escanear motivos con FDR en todos los promotores: minutos.
- 3D de TF+ADN: minutos-horas por par.
- 3D + RNAP + ensemble: horas-días por par.

Aunque el 3D se aplique al 10 %, si hay cientos de TF, siguen siendo cientos de
complejos. Hay que presupuestar GPU y definir un tope.

**Corrección:**

- separar tiers: secuencia → conservación → 3D top N;
- fijar N por genoma (p. ej. 20–50);
- reportar costo por tier.

---

## 3. Tabla de afirmaciones vs correcciones

| Afirmación del diseño | Problema | Corrección |
|---|---|---|
| Promotor >9 bits = ON | No es actividad, es parecido al motivo | `promoter_prior` calibrado |
| Promotor <6 bits = OFF | Puede estar ON por activador/sigma | No binario |
| Spacer 19–20 = MerR | No es único | Feature + motivo + 3D |
| Represor en [−45,+60] | Ventana rígida | Ventana amplia + score |
| Activador en [−100,−35] | Clase I/II difieren | Posición como evidencia |
| Palíndromo 14–22 | Sesgo a dímeros | PWM + motivos asimétricos |
| Choque RNAP = represor | Modelo estático | Score + alternativas |
| Toca α-CTD = activador I | α-CTD flexible | Ensemble + calibración |
| Tuerce 30° = MerR | Continuo | Ángulo + estructura |
| Looping por 10,5 pb | Heurística | Fase + geometría + energía |
| RNA-seq = lógica | Sin perturbación no se deduce | `logic_status=assumed` |
| Fisher = anti-alucinación | Prior, no verdad | FDR + clase |
| Cero pérdida | Imposible | Cero pérdida silenciosa |
| 90 % filtrado en segundos | Depende del tier | Costo por tier |

---

## 4. Flujo revisado: cuatro tiers

```text
TIER 0 — Secuencia y contexto (CPU, minutos)
  genoma + operones + promotores
  → PWM de promotor (σ70/σ alternativos)
  → biblioteca de motivos de familia
  → escaneo exhaustivo con FDR
  → contexto: operón, orientación, distancia, función

TIER 1 — Evolución (CPU, horas)
  → ortólogos de los 400 genomas
  → conservación del sitio y sintenía
  → ranking

TIER 2 — Estructura (GPU, horas-días, top N)
  → PDB/AFDB existentes; sin cristalizar
  → complejo TF–ADN; oligómero si hay plantilla
  → scores de interfaz, curvatura, choque con RNAP
  → mecanismo probable, no etiqueta

TIER 3 — Rescates (opcional)
  → looping, anti-activación, sigma alternativa, riboswitch
  → solo si Tier 0–2 no explican
```

La salida de cada tier es una tabla con `evidence`, `confidence`, `status` y
`provenance`. El tier 3 no reemplaza a los anteriores; añade hipótesis.

### 4.1 Integración con lo que ya existe

No hay que reconstruir L2/L3 ni el pipeline de represores:

| Pieza existente | Uso en el pipeline |
|---|---|
| `data/regulation/<genoma>/operons.tsv` | Unidades transcripcionales y fronteras |
| `data/regulation/<genoma>/promoters.tsv` | Promotores y sigma asignado |
| `data/regulondb/gt/promoter_pwms.tsv` | PWMs de σ24/28/32/38/54/70 |
| `data/regulondb/gt/tf_pwms.tsv` | 82 PWM de TF |
| `scripts/allostery/18_genome_pipeline.py` | Detección de DBD, familia, homología |
| `docs/K12_BACKTEST.md` | Backtest y comparación contra RegulonDB |
| `data/allostery/affinity_table.tsv` | Baseline de sitios conocidos |

Los módulos nuevos van en `scripts/census/`, no dentro de `scripts/regulation/`
ni de `webapp/`. Así se evita pisar archivos compartidos y se mantiene el
contrato de datos.

---

## 5. Qué implementar primero

| Prioridad | Tarea | Por qué |
|---|---|---|
| P0 | Backtest K-12 con PWM y operones existentes | Mide qué calza antes de añadir 3D |
| P0 | Corregir el bug de escaneo de `18_genome_pipeline.py` | Sin sitios no hay backtest |
| P0 | Usar `promoter_pwms.tsv` (σ24/28/32/38/54/70) | Ya existe; no reinventar |
| P1 | Biblioteca de motivos por familia + conservación | Es el motor real del censo |
| P1 | Calibrar promotor con TSS/RNA-seq de K-12 | Evita umbrales mágicos |
| P2 | 3D solo para top N, con scores | Evita falsa precisión |
| P2 | Rutas de rescate como módulos | No forzar explicaciones |
| P3 | RNA-seq y lógica | Solo si hay datos y perturbaciones |

---

## 6. Criterios de aceptación

- [ ] Cada etiqueta (`ON`, `OFF`, `represor`, `activador`, `MerR`) tiene un
      score y una calibración, no solo una regla.
- [ ] Cada descarte tiene `reason`.
- [ ] El backtest K-12 reporta precisión y recall por clase.
- [ ] Los controles negativos están incluidos (motivo barajado, PWM de otra
      familia, promotor no regulado, ligando *decoy*).
- [ ] El 3D se reporta como score, no como veredicto.
- [ ] Las rutas de rescate no se usan para forzar una conclusión.
- [ ] La coherencia metabólica es una capa, no un juez.
- [ ] Ningún `confirmed` sin experimento.

---

## 7. Frase final

**El diseño tiene la arquitectura correcta: filtrar barato, reservar el 3D para
pocos, y tener rutas de rescate. Lo que hay que evitar es convertir heurísticas
útiles (bits, spacer, choque, Fisher) en verdades binarias. El pipeline debe
producir hipótesis jerarquizadas con incertidumbre, no etiquetas definitivas.**
