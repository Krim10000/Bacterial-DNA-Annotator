# Buscador de complejos proteicos (subunidades)

Responde una pregunta que el código genético **no** responde: dado un genoma
anotado, ¿qué proteínas son candidatas a formar parte del mismo complejo, con
qué se combinan y con qué evidencia?

No predice la interfaz ni la estequiometría. Ordena hipótesis y las agrupa en
candidatos a complejo; la prueba sigue siendo experimental (co-IP, MS nativo,
crosslinking, cryo-EM) o estructural (PDB/ComplexPortal/AlphaFold-Multimer).

## La biología, en dos líneas

El código genético traduce codones a aminoácidos: da **una cadena**. La
estructura cuaternaria no está escrita de forma explícita en el gen; está
implícita en la superficie de la proteína.

Y en bacterias el ensamblaje **no ocurre en la transcripción**: ocurre en la
traducción (co-traduccional) o después. Lo que sí está acoplado en bacterias es
transcripción→traducción (no hay núcleo), y lo que el **operón** garantiza es
que las subunidades se produzcan juntas en tiempo y espacio. Por eso el contexto
genómico es una señal útil, pero no una prueba de interacción física.

Requiere `numpy` para aplicar; `scikit-learn` (extra `[full]`) solo para
ajustar y validar.

## Uso

```bash
# 1) gold de K-12 + caché de STRING (descarga ~10 MB, una vez)
python scripts/complexes/fetch_gold_k12.py

# 2) ajustar con el gold y aplicar al genoma
python scripts/complexes/build_candidates.py --fit --apply

# 3) validar (5-fold agrupada por complejo, ~2-3 min)
python scripts/complexes/validate_k12.py --save-oof
```

Salidas en `data/complexes/` (todo `data/` está en `.gitignore`, es regenerable):

| fichero | qué es |
|---|---|
| `gold_k12.tsv` | complejos curados de ComplexPortal mapeados a b-numbers |
| `model_k12.json` | pesos, media, escala y procedencia del modelo |
| `model_k12.genomic.json` | respaldo sin STRING (solo señales genómicas) |
| `pairs_ecoli_k12.tsv` | pares puntuados con la evidencia declarada |
| `clusters_ecoli_k12.tsv` | candidatos a complejo con confianza |
| `clusters_ecoli_k12_evidencia.tsv` | clusters etiquetados documentado/subcomplejo/parcial/nuevo |
| `evidence_summary.json` | conteo por evidencia × confianza |
| `VALIDATION.md` | informe de validación (se regenera) |

## Señales por par

| señal | fuente | qué aporta |
|---|---|---|
| `same_operon` | `find_operons` | co-transcripción, estequiometría |
| `adjacent`, `gap_bp`, `dist_bp` | genoma | vecindad física |
| `synteny` | `k12_synteny.tsv` | adyacencia conservada en 400 genomas |
| `same_ec` | anotación | misma reacción |
| `string_exp` | STRING `experiments` | interacción experimental |
| `string_nbr` | STRING `neighborhood` | vecindario conservado entre especies |
| `string_fus` | STRING `fusion` | fusión de dominios (Rosetta) |
| `string_cooc` | STRING `cooccurrence` | perfil filogenético |
| `string_coex` | STRING `coexpression` | coexpresión |
| `string_db` | STRING `database` | **excluida por defecto** |

### Por qué se excluye `database`

El canal `database` de STRING agrega complejos curados de bases de datos, entre
ellas ComplexPortal. Usarlo para validar contra ComplexPortal sería circular.
El modelo por defecto lo excluye (`--include-string-db` permite medir el techo:
R@1 sube de 0.917 a 0.921, poco — la señal no depende de esa vía).

## ¿Documentado o sospecha?

`annotate_evidence.py` compara cada cluster con ComplexPortal y lo etiqueta:

| etiqueta | criterio |
|---|---|
| `documentado` | reproduce el complejo curado (Jaccard ≥0.8) |
| `subcomplejo` | todos sus genes están en un complejo curado, pero es un trozo |
| `parcial` | ≥50 % de sus genes en un mismo complejo curado |
| `nuevo` | no se parece a nada del catálogo → sospecha |

Y separa la evidencia del catálogo: `experimental` (`ECO:0000353`) vs `inferida`
(ortología, parálogos, conocimiento previo). De los 273 complejos usables de
K-12, **133 son experimentales y 140 inferidos**: "curado" no es "medido".

En la salida de K-12: 86 documentados, 33 subcomplejos, 84 parciales y 757
nuevos (335 de ellos con confianza alta). "Nuevo" significa *no catalogado en
ComplexPortal*, no falso: el catálogo cubre una fracción de los complejos reales.

Es **anotación, no validación**: el catálogo se usó para entrenar. Para validar,
`validate_k12.py`.

## Validación en *E. coli* K-12

Gold: **ComplexPortal** (taxid 83333), 324 complejos → **273 usables** (≥2 genes
distintos y todas las subunidades mapeadas a b-number vía el `db_xref` UniProt
del GenBank). De esos 273, 133 son evidencia experimental y 140 son inferencia
curada: la validación mezcla ambos y así se reporta. Universo: **731 genes**, **266.815 pares**, **1.773 positivos**
(prevalencia 0,66 %). Los negativos son pares del universo que no comparten
complejo (mundo cerrado). CV 5-fold agrupada por complejo; métricas fuera de
fold.

| features | PR-AUC | R@1 | R@3 | R@5 | MRR |
|---|---|---|---|---|---|
| genómica | 0.426 | 0.765 | 0.824 | 0.837 | 0.797 |
| + STRING sin `database` | **0.700** | **0.917** | **0.973** | **0.980** | **0.945** |
| + STRING con `database` | 0.736 | 0.921 | 0.971 | 0.980 | 0.947 |

`R@k`: para cada proteína, fracción en que alguno de sus **k mejores partners**
pertenece a su complejo curado. Es la pregunta "¿con qué se combina?".

El top-1 cambia según la evidencia del gold: **94,0 %** de acierto en genes de
complejos experimentales (301 genes) vs **90,2 %** en complejos inferidos (449).
El modelo rinde más donde hubo experimento; no es un detalle, porque 140 de los
273 complejos del gold son inferencia curatorial.

Pares individuales (modelo sin db): precisión 34 % / recall 95 % con score ≥0.5.
La precisión por par es baja **a propósito**: STRING conecta proteínas de la
misma ruta que no forman complejo. La especificidad no viene del par, viene del
agrupamiento.

Clusters (top-2 mutuo sobre scores OOF): **209** candidatos, pureza **79,9 %**
(todos sus pares son del mismo complejo curado), cobertura Jaccard ≥0.5
**83,2 %**, recuperación exacta **42,9 %**, Jaccard medio **0,725**.

| confianza | n | pureza |
|---|---|---|
| alta | 171 | 90,6 % |
| media | 24 | 50,0 % |
| baja | 14 | 0,0 % |

Usar `alta`. `media` con cautela. `baja` es ruido y se conserva solo para
auditoría (quitarlas bajaría la cobertura de 83 % a 72 %).

### Casos documentados

| complejo | genes | Jaccard | ¿recuperado? |
|---|---|---|---|
| GroEL-GroES | 2 | 1.00 | sí |
| AcrAB-TolC | 3 | 1.00 | sí |
| ATP sintasa (`atp`) | 8 | 0.875 | sí |
| ARN polimerasa σ70 | 5 | 0.60 | sí |
| BAM | 5 | 0.60 | sí |
| ADN polimerasa III, núcleo | 3 | 0.43 | no |
| divisoma | 7 | 0.29 | no |
| 30S / 50S ribosomal | 21 / 34 | 0.10 / 0.11 | no |

## Límites (léelos antes de creerle)

- **Mundo cerrado.** Solo se evalúan genes que ya están en algún complejo
  curado. Un gen fuera de ComplexPortal no se puntúa en la validación.
- **Complejos grandes se fragmentan.** El top-2 mutuo no puede mantener unido un
  anillo de 30 subunidades: el ribosoma sale en trozos de 2-6. Un umbral más
  laxo no arregla esto, lo empeora (encadena rutas enteras: 14 clusters, uno
  gigante).
- **Mismo operón ≠ mismo complejo.** `sucABCD` codifica 2-oxoglutarato
  deshidrogenasa (`sucAB`) **y** succinil-CoA sintetasa (`sucCD`); el clúster
  predicho los mezcla.
- **Complejo ≠ holoenzima.** El núcleo de la ADN pol III (`dnaE`, `dnaQ`,
  `holE`) se mezcla con el holoenzima (`dnaN`, `dnaX`, ...). Los dos existen; la
  frontera es discutible.
- **`score` no es una probabilidad calibrada.** El modelo entrena con clases
  balanceadas y la prevalencia real es ~0,7 %.
- **La CV agrupa por complejo, pero un gen puede estar en dos complejos** y
  aparecer en train y test. Mide recuperación, no generalización a genes nuevos.
- **Sin STRING no hay milagro.** Si el genoma no tiene red de STRING, el script
  cambia solo a `model_k12.genomic.json` (R@1 0,76 en K-12). Recupera complejos
  codificados en operones, no los dispersos, y la confianza "alta" pasa a
  significar "operón", no "interacción".

## Aplicar a otro genoma

```bash
python scripts/complexes/build_candidates.py --apply-only \
  --genome GCF_000006985.1 \
  --genes data/regulation/GCF_000006985.1/genes.tsv \
  --operons data/operon_finder/GCF_000006985.1_operons.tsv \
  --string data/complexes/raw/string_<taxid>_full.txt.gz \
  --string-species <taxid> \
  --model data/complexes/model_k12.json
```

El modelo está ajustado en K-12: en otro genoma es un prior, no una
calibración. `--string` es opcional; sin él, `string_*` queda en 0.

## Catálogo de herramientas externas

[SUGERENCIAS.md](SUGERENCIAS.md) ordena los siguientes pasos por ROI: transferencia
por ortología, consistencia de masa, perfil filogenético, Foldseek y coevolución,
cada uno con nivel de evidencia, score y criterio de aceptación.

[TOOLS.md](TOOLS.md) evalúa 16 herramientas (ComplexPortal, EcoCyc, PDB/PDBe,
IntAct, BioGRID, Foldseek-Multimer, AlphaFold-Multimer, Boltz-2, Operon-mapper,
UniProt, …) con nivel de evidencia A/B/C/D, score/umbral, licencia, coste y
veredicto. Incluye el contrato de puntuación y la trampa de circularidad.

## Tests

```bash
python -m pytest tests/test_complexes_k12.py -q
```

Cubren el parser de GenBank (la continuación de `/product` que empieza con
dígito no debe contaminar la location), el mapeo de ComplexPortal (ignora
ligandos ChEBI), las features y el filtro top-k mutuo.
