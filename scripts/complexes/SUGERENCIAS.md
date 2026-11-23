# Siguientes pasos para el buscador de complejos

> **Integración 2026-11-23:** el pipeline ya llama a `build_candidates.py`
> mediante `webapp/complex_model.py`; el modelo K-12 se aplica como prior y sus
> clusters salen como `complex_candidate` amarillos en la GUI.
>
> Estado de esta lista:
> 1. transferencia por ortología — **implementada** (`transfer_orthology.py`);
> 2. masa/estequiometría — **parcial** (masa observada y `expected_n`, sin masa
>    documentada porque el gold local no la trae);
> 3. perfil filogenético — **implementado** (`phylo_jacc`);
> 4. Foldseek-Multimer — **interfaz implementada**
>    (`foldseek_multimer.py`), falta la DB de ensamblajes PDB;
> 5. coevolución — **pendiente** (falta MSA/CPU).

Documento de trabajo: qué agregar, en qué orden, con qué evidencia y con qué
criterio se acepta. Nada de esto está implementado todavía salvo lo indicado.

## Punto de partida

- Pipeline: `pairs_*` (score de par) → `clusters_*` (top-2 mutuo) →
  `annotate_evidence.py` (documentado/subcomplejo/parcial/nuevo).
- Validación K-12: R@1 0.917 sin el canal `database` de STRING; cobertura de
  clusters 83.2 %; pureza 79.9 %.
- Catálogo: 273 complejos usables, **133 experimentales** (`ECO:0000353`) y
  **140 inferidos**.
- Límites conocidos: sin STRING el modelo baja a R@1 0.765; los complejos
  grandes se fragmentan; "mismo operón" no es "mismo complejo".

## Criterio de prioridad

Ordené por **(valor para el producto) / (coste)**, con un requisito: que la
evidencia quede etiquetada como A/B/C/D y con score, igual que el resto del
módulo. Nada de "herramienta que dice sí/no" sin umbral.

| # | idea | nivel | score / umbral | coste | insumos ya disponibles |
|---|---|---|---|---|---|
| 1 | Transferencia por ortología | B | fracción de subunidades ≥0.8 × identidad | bajo | `data/terminators/k12_vs_all.tsv`, `.cache/synteny/hits.tsv` |
| 2 | Consistencia de masa | A/B | error de masa ±15 % | bajo | `gold_k12.tsv` (`properties`), GenBank |
| 3 | Perfil filogenético | D | Jaccard de presencia ≥0.5 | bajo | `k12_vs_all.tsv` |
| 4 | *Subunit propensity* por proteína | D | score combinado 0-1 | bajo | derivado de `pairs_*` |
| 5 | Foldseek-Multimer vs PDB | B/C | TM ≥0.7 fuerte, ≥0.5 mismo fold | medio | binario + DB en `data/db/` |
| 6 | Coevolución de interfaz | C | score de acoplamiento | alto | requiere MSA por par |
| 7 | Exclusividad de parálogos | D | presencia mutuamente excluyente | bajo | `k12_vs_all.tsv` |

---

## 1. Transferencia del catálogo por ortología (nivel B)

**Qué.** Proyectar los 273 complejos curados de K-12 a los otros 400 genomas
usando ortología. Cada genoma tendría su propio catálogo "documentado por
homología", en vez de solo "sospecha".

**Por qué.** Hoy el catálogo documentado es solo K-12. Es la mayor ganancia de
cobertura y encaja con el objetivo del proyecto (400 genomas).

**Insumos.** Ya existen:

- `.cache/synteny/hits.tsv` — 400 genomas → K-12, 222.770 hits,
  columnas `qseqid, sseqid, pident, evalue, bitscore`;
- `data/terminators/k12_vs_all.tsv` — K-12 → 400 genomas, 335.109 hits.

**Método.**

1. Ortología por *reciprocal best hit* (RBH): A→B y B→A.
2. Filtros mínimos: identidad ≥30 %, `evalue ≤1e-5`, y cobertura de alineamiento
   ≥70 % (hoy `hits.tsv` no trae longitud de alineamiento: hay que re-exportar
   con `--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore`).
3. Para cada genoma y cada complejo: subunidades con ortólogo / total.
4. Emitir solo si ≥2 subunidades mapeadas.

**Salida.** `data/complexes/transfer_<genoma>.tsv`:

```
complex_id  complex_name  n_subunits  n_mapped  frac_mapped
identity_media  genes_k12  genes_genoma  evidence_class  eco
```

**Score.** `score = frac_mapped × identidad_media`, con etiqueta
`documentado_inferido` (`ECO:0005544`). Umbral de salida: `frac_mapped ≥0.5`;
"alta" si `≥0.8` y `identidad_media ≥0.5`.

**Validación.** Dos niveles:

- contra K-12 dejando fuera el genoma evaluado (no hay fuga);
- si ComplexPortal tiene la especie, comparar con su `complextab/<taxid>.tsv`.

**Criterio de aceptación.** ≥60 % de los genomas con ≥50 complejos a
`frac_mapped ≥0.8`, y precisión ≥0.8 contra el catálogo de la especie cuando
exista.

**Riesgo.** El best-hit solo (sin RBH) infla la cobertura: un homólogo lejano
puede colar. Por eso el RBH y el umbral de identidad no son opcionales.

---

## 2. Consistencia de masa y estequiometría (nivel A/B)

**Qué.** Para un cluster que coincide con un complejo curado, comparar la masa
sumada de sus subunidades con la masa documentada.

**Por qué.** Convierte "Jaccard 0.875" en algo accionable: "falta una
subunidad". Ejemplo real: ATP sintasa con 7 de 8 genes.

**Insumos.**

- `gold_k12.tsv`, columna `properties`: trae `Molecular weight = 206 kDa`
  (a veces varias especies: "A secondary 400-500 kDa species…").
- Longitudes de proteína: `translation` del GenBank o `(end-start+1)/3`.

**Método.**

1. Parsear `Molecular weight\s*=\s*([\d.,]+)\s*kDa` (tomar la primera especie).
2. Masa de cada subunidad ≈ `n_aa × 110 Da`.
3. `error = |masa_cluster − masa_doc| / masa_doc`.
4. Si `error > 0.15` y `n_genes < n_subunits`, listar las subunidades curadas
   ausentes y sugerirlas.

**Salida.** Columnas nuevas en `clusters_*_evidencia.tsv`:
`masa_kda`, `masa_doc_kda`, `error_masa`, `sugiere_faltan`.

**Score.** `confianza_masa = 1 - min(1, error/0.5)`. No reemplaza el score del
modelo: lo acompaña.

**Criterio de aceptación.** En los clusters que reproducen un complejo curado
(Jaccard ≥0.8), el error de masa debe ser <15 % en ≥90 % de los casos.

**Riesgo.** La masa documentada puede incluir ligandos o especies alternativas;
por eso solo se aplica a coincidencias documentadas, no a candidatos nuevos.

---

## 3. Perfil filogenético con los 400 genomas (nivel D)

**Qué.** Genes que aparecen y desaparecen juntos tienden a formar complejo
(*phylogenetic profiling*, Pellegrini 1999).

**Por qué.** Es la única señal fuerte que **no depende de STRING**, justo la
debilidad del modelo actual en genomas sin red.

**Insumos.** `data/terminators/k12_vs_all.tsv` (presencia de cada gen de K-12 en
cada genoma).

**Medición de factibilidad (ya hecha).**

| métrica | pares dorados | pares al azar |
|---|---|---|
| Jaccard medio de presencia | 0.257 | 0.105 |
| Jaccard ≥0.5 | 18.9 % | 4.2 % |
| Jaccard ≥0.9 | 2.7 % | 1.5 % |

Enriquecimiento 4.5× en J≥0.5: señal real pero moderada. Cobertura mediana por
gen: 18 de 400 genomas (son especies diversas, no cepas).

**Método.**

1. Matriz 400 × 3761 de presencia (`k12_vs_all.tsv`).
2. `phylo_jacc(a,b) = |P(a)∩P(b)| / |P(a)∪P(b)|`.
3. Agregar como feature `phylo_jacc` al modelo y reajustar.
4. Reportar el lift en el modelo **genómico** (sin STRING), que es donde importa.

**Salida.** `data/complexes/phylo_profiles.tsv` y la feature en
`model_k12.genomic.json`.

**Criterio de aceptación.** Subir el R@1 genómico de 0.765 a ≥0.80 sin usar
STRING, y no degradar el modelo completo.

**Riesgo.** Muchos genes tienen perfil escaso; hay que tratar `0/0` como "sin
dato", no como Jaccard 0.

---

## 4. *Subunit propensity* por proteína (nivel D, para la UI)

**Qué.** Un score por proteína que responda "¿es subunidad de algo?" sin exigir
un partner concreto.

**Por qué.** Hoy la respuesta es a nivel de par; la pregunta original es a nivel
de proteína. Sirve para ordenar la tabla y para priorizar qué mandar a Boltz.

**Features.** Mejor score de par, número de partners con score ≥0.9, contexto
de operón (¿comparte operón con ≥2 genes?), longitud, TM/signal peptide,
dominios de interacción.

**Salida.** `data/complexes/protein_propensity.tsv` con
`locus_tag, propensity, n_partners, mejor_partner, evidencia`.

**Criterio de aceptación.** AUC ≥0.8 para separar "gen en algún complejo curado"
vs "gen fuera", en validación por gen (no por par).

**Riesgo.** Es derivado del mismo modelo: puede ser circular. Etiquetarlo como
prior de UI, no como evidencia.

---

## 5. Foldseek-Multimer contra PDB (nivel B/C)

**Qué.** Buscar el candidato en ensamblajes del PDB y devolver TM-score.

**Por qué.** Pasa de "hipótesis" a "existe una estructura resuelta parecida".
No requiere GPU.

**Insumos.** Binario y DB en `data/db/foldseek_swissprot*` (hoy es DB de
monómeros: hace falta una DB de **complejos**, p. ej. PDB biological
assemblies).

**Método.** `foldseek easy-complexsearch <modelo.pdb> <db_complejos> <out> --monomer-include-mode 0`.

**Score.** TM ≥0.7 evidencia fuerte; 0.5-0.7 mismo fold; <0.5 descartar.

**Riesgo.** Licencia GPL-3.0: no redistribuir en un SaaS propietario sin revisar.

---

## 6. Coevolución de interfaz (nivel C, apuesta)

**Qué.** Para el top de pares, calcular acoplamiento de residuos (EVcouplings,
GREMLIN, DCA) y usarlo como evidencia de interfaz.

**Por qué.** Es la señal que mejor separa "misma ruta" de "interfaz real", que
es justo donde el par falla (precisión 46.7 % a score ≥0.9).

**Coste.** Alto: MSA por par, alineamientos profundos, GPU/CPU intensiva. Solo
para los ~20 pares que ya pasaron el filtro genómico.

**Criterio de aceptación.** Mejorar la precisión de pares a score alto sin bajar
el recall, medido en el mismo universo de 266.815 pares.

---

## 7. Evidencia negativa: exclusividad de parálogos (nivel D)

**Qué.** Pares que nunca coocurren en los 400 genomas (patrón mutuamente
excluyente) no deberían unirse en un cluster.

**Por qué.** Evita fusionar parálogos que ocupan el mismo sitio en complejos
distintos. Es el contrapeso del perfil filogenético.

**Método.** `P(a)∩P(b)=∅` y ambos con perfil ≥10 genomas → penalización.
Score: `-log(1 + n_genomas_sin_coocurrencia)` acotado.

---

## Lo que no vale la pena (y por qué)

- **AlphaFold DB**: modelos monómero; no responde "con quién".
- **KEGG MODULE / MetaCyc**: módulos funcionales, no complejos físicos.
- **CORUM**: mamíferos.
- **MCL / ClusterONE**: algoritmos, no evidencia; ya probamos MCL (Jaccard 0.43
  vs 0.72 del top-2 mutuo).
- **Operon-mapper**: mejora marginal de una feature que STRING ya aproxima;
  solo si se quiere independencia de STRING.

## Orden recomendado

1. **Transferencia por ortología** — pan-genoma, nivel B, coste bajo.
2. **Masa/estequiometría** — mejor retorno por hora; convierte Jaccard en
   "falta `atpI`".
3. **Perfil filogenético** — barato e independiente de STRING.
4. **Foldseek-Multimer** — evidencia estructural con score.
5. **Coevolución** — solo si 1-4 ya están y se quiere subir la precisión fina.

## Reglas que no se negocian

- Nada documentado se usa a la vez como feature y como validación: circularidad.
- Toda evidencia lleva nivel (A/B/C/D) y score con umbral.
- Ortología ≠ prueba: nivel B, siempre.
- Los candidatos nuevos se etiquetan `nuevo`, no `falso`.
