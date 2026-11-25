# ¿Vale la pena usar rutas metabólicas como insumo de clasificación?

> **Respuesta corta:** **no como clasificador principal ni como feature fuerte.**
> Como prior débil de desempate aporta muy poco y puede empeorar. Sirve más
> para interpretar y priorizar que para decidir el sigma.

---

## 1. Prueba hecha

Sobre K-12 gold, TSS con gen y producto anotado:

| Dato | Valor |
|---|---:|
| TSS gold | 3.314 |
| TSS con producto | 1.425 (43,0 %) |
| Sigma evaluables | 6 |

Se añadió un prior binario por palabras clave de función:

- σ28: flagelo/quimiotaxis;
- σ32: heat shock/chaperonas;
- σ38: estacionaria/osmoprotección;
- σ54: nitrógeno/formato/transporte;
- σ24: envelope/ECF;
- σ70: ribosoma/traducción.

Score final: `z + w × prior`.

| w | Accuracy | Δ vs sin prior | Cambios de top-1 |
|---:|---:|---:|---:|
| 0,00 | 0,5144 | — | 0 |
| 0,25 | **0,5179** | **+0,35 pp** | 17 |
| 0,50 | 0,5172 | +0,28 pp | 47 |
| 1,00 | 0,5144 | 0,00 pp | 107 |
| 1,50 | 0,5109 | −0,35 pp | 150 |
| 2,00 | 0,5039 | −1,05 pp | 194 |

En el subconjunto de **margen bajo** (top-1 − top-2 < 1, n=1.004):

| w | Accuracy | Δ |
|---:|---:|---:|
| 0,0 | 0,4602 | — |
| 0,5 | 0,4641 | +0,40 pp |
| 1,0 | 0,4602 | 0,00 pp |

**Lectura:** el prior apenas mueve la aguja; con peso alto empeora. No justifica
una iteración completa del clasificador.

---

## 2. Por qué aporta tan poco

1. **Cobertura baja:** solo 43 % de los TSS tienen gen/producto.
2. **Anotación derivada de homología:** la función del gen no es independiente
   del pipeline; puede haber circularidad.
3. **Los regulones sigma no son estrictamente metabólicos:** transporte,
   reguladores, genes hipotéticos, estrés, envelope.
4. **Un sigma controla varios módulos:** σ38 no es “solo estacionaria”; σ54 no
   es “solo nitrógeno”.
5. **El cuello está en otra capa:** localización 49,5 % y motivo σ54/σ38
   28–49 %. La ruta no arregla eso.
6. **Rutas ubicuas:** ribosoma, traducción y metabolismo central aparecen en
   muchas clases; no discriminan.

---

## 3. Dónde sí sirve

| Uso | ¿Vale? | Por qué |
|---|---|---|
| Decidir el top-1 | **No** | Aporta +0,35 pp y cambia pocos casos |
| Desempatar top-1/top-2 con margen bajo | **Débil** | +0,40 pp; solo como evidencia |
| Subir confianza de un candidato ya fuerte | **Sí** | Contexto coherente refuerza la hipótesis |
| Interpretar un regulón | **Sí** | Explica qué módulo controla |
| Detectar un `novel_candidate` | **Sí, con cautela** | Coherencia funcional + conservación |
| Filtrar duro | **No** | Pierde genes no metabólicos y regulones mixtos |

---

## 4. Decisión

No se incorpora como feature del clasificador en esta versión.

Se mantiene:

- `pathway_support` / `pathway_hits` como **evidencia** en
  `sigma_promoters.tsv`;
- el uso en la vista Sigma y en el regulón para interpretación;
- la posibilidad de subir confianza cuando hay contexto coherente **y** la
  señal PWM ya es fuerte.

No se hace:

- sumar un bonus de ruta al z-score;
- usar la ruta para cambiar el sigma asignado;
- declarar falso positivo por falta de ruta.

---

## 5. Cuándo revisar

Vale la pena otra iteración si se consigue:

1. un mapa gen→ruta/EC curado para K-12 y Subtilis (no solo palabras clave);
2. cobertura >70 % de los TSS con gen;
3. un gold con negativos;
4. validación *held-out* por especie;
5. expresión o conservación que haga la ruta más informativa.

Hasta entonces, la ruta es **contexto**, no clasificador.

---

## 6. Frase final

**La ruta metabólica no es un mal insumo, pero en este pipeline no es el cuello
ni la señal que desempata: aporta +0,35 pp y puede empeorar con peso alto. Se
usa para interpretar y reforzar, no para decidir el sigma.**
