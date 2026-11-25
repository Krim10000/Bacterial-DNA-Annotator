# K-12 como laboratorio de aprendizaje: patrón general y subcasos

> **Pregunta:** ¿usamos K-12 para aprender el patrón general y los subcasos, y
> luego aplicamos eso a los 400 genomas?
>
> **Opinión corta:** sí, pero con una condición. K-12 es el mejor lugar para
> **aprender la gramática** del sistema, no para **definir la verdad** de todas
> las bacterias. K-12 da el patrón general y muchos subcasos; los 400 genomas
> dicen cuáles son universales, cuáles son de linaje y cuáles son excepciones.
>
> Documento relacionado: [`CENSUS_REPRESSORES.md`](CENSUS_REPRESSORES.md).

---

## 1. Por qué K-12 es el lugar correcto para aprender

K-12 no es una bacteria cualquiera. Es el genoma bacteriano mejor anotado del
mundo para regulación:

| Dato | Valor |
|---|---:|
| TF curados | 243 |
| TF con efector anotado | 93 (38 %) |
| Represores decisivos | 36 (27 inductor / 9 correpresor) |
| Interacciones TF–promotor con secuencia | 3.444 (216 TF) |
| PWM de TF | 82 |
| Operones | 2.601 |
| Genoma | 4,6 Mb, ~4.400 genes |

Eso permite hacer algo que en los 400 genomas no se puede: **ver la respuesta
correcta**. K-12 es el único lugar donde podemos comparar el censo contra una
verdad curada, medir precisión y aprender qué evidencia basta.

Pero K-12 es también una bacteria concreta: gamma-proteobacteria, GC ~50,8 %,
genoma mediano, con una historia de investigación centrada en *lac*, *trp*,
*ara* y *lexA*. No es una muestra aleatoria de la diversidad bacteriana.

**Conclusión:** K-12 es una lente, no la verdad. Sirve para aprender la
gramática; los 400 genomas sirven para saber si la gramática es universal.

---

## 2. Qué patrón general se puede aprender de K-12

El patrón general no es una secuencia universal. Es la **lógica del sistema**:

1. **Los sitios no están al azar.** Se concentran en regiones reguladoras,
   típicamente solapando el promotor o el TSS.
2. **La geometría manda.** Un represor oligomérico reconoce repeticiones
   simétricas o semisimétricas; el espaciador y la orientación importan.
3. **El efector cambia el equilibrio DBD–ADN.** No cambia "la función" en
   abstracto: estabiliza una conformación que une o que no une operador.
4. **El regulón es coherente.** Un represor de metabolismo de azúcares tiende a
   regular transporte/catabolismo de azúcares; un sensor de metal regula
   homeostasis de metal.
5. **La conservación es señal.** Los sitios funcionales se conservan en
   ortólogos; los *near-sites* no.
6. **La arquitectura del operón importa.** Divergón, distancia al TSS, número de
   genes y presencia de atenuadores cambian la interpretación.

Estos seis puntos son el "patrón general". Se pueden aprender de K-12 y se
pueden testear en los 400.

---

## 3. Qué subcasos hay que aprender de K-12

Un TF puede tener **varias etiquetas a la vez**. No se fuerza una sola clase.

| Subcaso | Ejemplo K-12 | Qué enseña | Implicación para el censo |
|---|---|---|---|
| Represor simple | TetR (inhibidor de tetO) | Un operador principal, un efector | Caso base |
| Múltiples operadores | LacI (O1/O2/O3), TrpR, PurR | Varios sitios, cooperatividad, *looping* | No truncar a top-1 |
| Inductor | LacI, TetR | El efector estabiliza el estado libre | `direction=inducer` |
| Correpresor | PurR, TrpR | El efector estabiliza el estado unido | `direction=corepressor` |
| Sensor de metal | Fur, Zur, CueR | El efector no es un metabolito | GBA puede ser ciego |
| Sensor de xenobiótico | MarR, MprA | El efector no está en la ruta | No forzar química del operón |
| Arquitectónico | H-NS, Fis, IHF, HU | No hay operador discreto | Excluir del censo de operadores |
| Regulador dual | AraC, XylR | Activa o reprime según contexto | Permitir ambos signos |
| Autoproteólisis | LexA | La inactivación no es un ligando | No buscar efector |
| Antirepresor | KilAC, AntA/AntB | No une ADN como represor | Excluir |
| Ambiguo | MetJ, IclR, ExuR, ZntR, NikR, AcrR, AllR | El efector aparece en ambas conformaciones | No forzar dirección |
| Huérfano | TF sin efector | No hay etiqueta | `status=unknown` |

**Lección:** el "subcaso" no es un error del modelo. Es una **región distinta del
espacio de parámetros**. El censo debe reconocerla y bajar la confianza, no
inventar una respuesta.

---

## 4. Qué NO puede enseñar K-12

| Límite | Por qué | Consecuencia |
|---|---|---|
| Diversidad filogenética | Es una gamma-proteobacteria | Un patrón puede ser de linaje, no universal |
| GC y composición | ~50,8 % GC | Motivos y espaciadores pueden cambiar en genomas AT/GC extremos |
| Sigma factors | Repertorio de *E. coli* | Otros filos usan otros sigma y promotores |
| Familias ausentes | Algunas familias no están en K-12 | No hay prior para ellas |
| n pequeño de dirección | 36 decisivos, 9 correpresores | No alcanza para entrenar un clasificador general |
| Un TF por operón | Regulones superpuestos | El censo debe permitir varios TF por operón |
| K-12 no es representativo | Es el genoma más estudiado | Sesgo de "lo que se publicó", no de "lo que existe" |

**Regla:** K-12 se usa para **inducir reglas y calibrar**, no para reportar la
accuracy final. La accuracy final se mide en *held-out* por familia y por genoma.

---

## 5. Modelo mental: gramática + parámetros + excepciones

En vez de "un patrón general y casos raros", conviene pensar en tres capas:

```text
Capa 1 — Gramática (universal)
  el represor une un sitio en una región reguladora;
  el efector modula la afinidad por ADN;
  el regulón tiene coherencia funcional.

Capa 2 — Parámetros por familia (aprendidos de K-12 + 400)
  motivo, longitud, simetría, espaciador, orientación;
  número típico de operadores;
  clase química del efector.

Capa 3 — Excepciones por linaje o por TF (subcasos)
  H-NS arquitectónico;
  LexA autoproteolítico;
  MarR xenobiótico;
  TF dual;
  familia huérfana.
```

El censo aplica la gramática a todos, los parámetros de familia a cada familia,
y las excepciones solo cuando la evidencia las activa. Eso evita dos errores:

- forzar un modelo único para todo (falso negativo en subcasos);
- llenar el modelo de excepciones *ad hoc* (sobreajuste a K-12).

---

## 6. Protocolo concreto

### Paso 1 — Atlas K-12

Construir `data/regulation/k12_atlas.tsv` con una fila por TF:

```text
tf, family, class, operator_count, operator_arrangement,
promoter_overlap, effector_class, direction, mechanism,
evidence, confidence, source
```

### Paso 2 — Etiquetar subcasos

Aplicar las etiquetas de la tabla §3. Un TF puede tener varias. Marcar
`ambiguous` y `orphan` explícitamente; no inventar.

### Paso 3 — Inducir reglas por familia

Para cada familia con ≥3 TF en K-12:

- motivo consenso y PWM;
- número y disposición típica de operadores;
- solapamiento con promotor;
- clase de efector y dirección;
- evidencia mínima para `high`.

### Paso 4 — Testear en los 400

- Aplicar las reglas a todos los genomas.
- Medir por clase, no en agregado.
- Si una regla falla sistemáticamente en un filo, se marca como **subcaso de
  linaje**, no como error.
- Si falla en todas partes, la regla estaba mal.

### Paso 5 — Volver a K-12

Cada subcaso nuevo se documenta con su ejemplo K-12 y su frecuencia en los 400.
El atlas crece, pero la gramática no se rompe.

---

## 7. Criterios de aceptación

- [ ] Los 243 TF de K-12 están en el atlas, con clase y evidencia.
- [ ] Los 36 decisivos tienen `direction` y `mechanism` documentados.
- [ ] Los 7 ambiguos están marcados como ambiguos, no forzados.
- [ ] Las reglas por familia tienen ≥3 casos o se declaran `insufficient`.
- [ ] El test en los 400 es por clase, no solo global.
- [ ] Los subcasos de linaje se reportan como tales.
- [ ] Ningún `confirmed` sin experimento.

---

## 8. Opinión final

**Sí a usar K-12 para aprender.** Es el único lugar donde tenemos la respuesta
curada y suficiente densidad de datos para inducir reglas. Pero K-12 no es el
patrón universal: es el **laboratorio de calibración** y el **catálogo de
subcasos**.

La frase que resume la estrategia:

> **K-12 enseña la gramática y los subcasos; los 400 genomas dicen qué parte de
> esa gramática es universal y qué parte es solo de *E. coli*.**
