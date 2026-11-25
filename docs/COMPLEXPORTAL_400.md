# ComplexPortal en las 400 especies: estado real y plan

Auditoría: **2026-11-24 21:13 -03** (2026-11-25 00:13 UTC).
Alcance: los 400 genomas GCF del panel más `ecoli_k12` como referencia.
No se descargó ningún TSV de especie; solo se leyó el listado FTP de
`complextab/`.

## Respuesta corta

| pregunta | respuesta |
|---|---|
| genomas inventariados | **401** (400 GCF + `ecoli_k12`) |
| con especie local | **401/401** |
| con taxid local | **401/401** (400/400 del panel) |
| con catálogo ComplexPortal local | **1/401**: K-12 (`83333`) |
| panel con catálogo local | **0/400** |
| catálogos de organismo publicados por ComplexPortal | **28** |
| panel con entrada en ComplexPortal | **0/400** |
| especies del panel con un catálogo curado disponible hoy | **0** |

No falta una descarga: **no hay catálogo ComplexPortal para ninguna de las 400
especies del panel**. ComplexPortal publica 28 catálogos de organismo, casi
todos de organismos modelo; ninguno coincide con los taxids del panel, ni por
taxid de ensamblado ni por `species_taxid`. El único catálogo local del proyecto
es el de *E. coli* K-12, que ya está en `data/complexes/gold_k12.tsv`.

## Evidencia local

Salida real de
[`inventory_complexportal_400.py`](../scripts/complexes/inventory_complexportal_400.py)
sobre este checkout:

```text
genomas: 401 (panel GCF: 400, referencias: 1)
con especie: 401
con taxid: 401 (panel: 400)
con catálogo local: 1 (panel: 0)
estados: ready=1, raw_only=0, catalog_missing=400, taxid_missing=0, taxonomy_missing=0
```

Los detalles por genoma están en
[`complexportal_inventory.tsv`](../reports/complexportal_inventory.tsv) y
[`complexportal_inventory.md`](../reports/complexportal_inventory.md). Columnas:
`genome, species, taxid, local_catalog, status, reason`.

Hechos comprobados:

- 400 JSON en `data/taxonomy/full/`, con 400 especies y 400 taxids distintos.
  Estados de taxonomía: 382 `assigned`, 18 `provisional`; todos tienen taxid.
- El taxid usado es el **del ensamblado** (`ncbi_taxid`), no el
  `species_taxid`. El README de taxonomía ya documenta que NCBI no siempre
  propaga la especie; usar el del ensamblado evita atribuir un catálogo de otra
  cepa.
- `ecoli_k12` no es un GCF del panel y no tiene JSON en `data/taxonomy/`. Su
  especie se lee del GBFF (`data/test_ecoli.gbff`, línea `ORGANISM`) y su taxid
  `83333` se lee de la cabecera del TSV crudo
  `data/complexes/raw/complexportal_83333.tsv` (columna
  `Taxonomy identifier`). Si ese TSV falta, el inventario deja el taxid en
  `not_available`; no lo rellena a mano.
- K-12: `data/complexes/raw/complexportal_83333.tsv` tiene 324 complejos;
  `data/complexes/gold_k12.tsv` tiene **273** usables (43 homómeros descartados,
  7 con subunidades sin mapear, según
  [`manifest.json`](../data/complexes/manifest.json)).
- Los 400 `data/complexes/transfer_<genoma>.tsv` y los 400
  `data/complexes/model_<genoma>.complexes.tsv` **no son catálogos
  ComplexPortal**. Son transferencia por ortología y clusters del modelo K-12;
  el inventario no los cuenta como catálogo local.

## Auditoría remota: ComplexPortal no cubre el panel

Fuente:
`https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/`
(listado leído una vez; no se descargaron los TSV).

| comprobación | resultado |
|---|---:|
| ficheros `complextab/<taxid>.tsv` publicados | **28** |
| taxids del panel (400) presentes | **0** |
| taxids del panel por `species_taxid` presentes | **0** |
| K-12 `83333` presente | **sí** |

Los 28 taxids del listado son:

```text
562, 3702, 6239, 6523, 7227, 7787, 7788, 7955, 8355, 8732, 9031, 9606,
9615, 9823, 9913, 9940, 9986, 10090, 10116, 83333, 208964, 243277,
284812, 559292, 694009, 1235996, 1263720, 2697049
```

Es una foto de esa fecha. ComplexPortal puede añadir organismos; el inventario
local no consulta la red, así que hay que repetir esta comprobación antes de
planificar una descarga. Un `catalog_missing` en el TSV significa **sin
catálogo local**, no “ComplexPortal no lo tiene”: para eso está esta sección.

## Plan por taxid

El registro por taxid es
[`complexportal_inventory.tsv`](../reports/complexportal_inventory.tsv). No se
duplican aquí los 400 taxids porque todos caen en la misma clase y su acción es
idéntica:

| clase | taxids | n | acción |
|---|---|---|---:|---|
| K-12 | `83333` | 1 | usar `gold_k12.tsv` como anotación/entrenamiento; no como validación independiente |
| panel con entrada remota | ninguno | 0 | descargar `complextab/<taxid>.tsv`, mapear y filtrar incompletos |
| panel sin entrada remota | los 400 del TSV | 400 | **no descargar**; usar transferencia K-12 u otra fuente, marcada como inferida |

Pasos concretos:

1. Re-ejecutar el inventario para confirmar genomas, taxids y catálogos locales.
2. Revisar el listado remoto. Si aparece un taxid del panel:
   - descargar **solo ese** TSV;
   - mapear UniProt→gen con `functional.tsv` (ver abajo);
   - decidir `gold` (complejo completo) o parcial (`--allow-incomplete`);
   - guardarlo como `gold_<taxid>.tsv`, no como validación si se usó para
     ajustar el modelo.
3. Para los 400 sin entrada: no hay nada que descargar de ComplexPortal. La
   alternativa ya existente es:
   - `transfer_<genoma>.tsv` (ortología K-12, nivel B);
   - `model_<genoma>.complexes.tsv` (modelo genómico K-12);
   - evidencia externa (PDB/IntAct/BioGRID, literatura, masa, sintenía).
4. Si se amplía la cobertura, priorizar especies con más impacto en las brechas
   metabólicas y con mapeo UniProt de buena calidad; no las 400 a la vez.

## Mapeo UniProt→locus

### Qué hay en los datos

- `data/regulation/<genoma>/genes.tsv`:
  - `gene_id`: ID de proteína, p. ej. `NZ_AP019822.1__NZ_AP019822.1_1`;
  - `locus_tag`: a menudo **vacío** en los 400 GCF;
  - `acc`: accesión UniProt/Swiss-Prot cuando hay hit (p. ej. `A7FPR6`);
  - `locus_tag` y `gene_id` son claves estables dentro del genoma; en el panel
    el `gene_id` es la clave que existe.
- `data/annotations_ab_initio/<genoma>/functional.tsv`:
  - sin cabecera en este proyecto;
  - columna 1 = `qseqid` = el mismo `gene_id`;
  - columna 2 = `sp|ACC|...` o `tr|ACC|...`.
  - Existe para **400/400** genomas del panel.
- K-12 no tiene `functional.tsv`, pero su `genes.tsv` sí trae `acc`; el
  inventario lo marca como `genes.tsv acc`.

### Cómo mapear

1. Parsear la columna de miembros de ComplexPortal
   (`Identifiers (and stoichiometry) of molecules in complex`). El lector
   existente usa `ACC(n)` y descarta líneas que no encajan
   (`complexes_lib._PART`).
2. Construir `ACC → gene_id`:
   - desde `functional.tsv` (columna 2 → columna 1), o
   - desde `genes.tsv` (`acc` → `gene_id`/`locus_tag`).
3. Unir con `genes.tsv` por `gene_id`; usar `locus_tag` si no está vacío y, si
   está vacío, `gene_id` como identificador de locus.
4. Descartar el complejo si falta alguna subunidad (`drop_incomplete=True`) para
   el gold; conservarlo como parcial solo con `--allow-incomplete`.
5. Reportar `total`, `kept`, `unmapped` y `homomers`; no convertir un parcial en
   gold completo.

### Límite actual de `fetch_complexportal.py`

`uniprot_map()` devuelve `ACC → qseqid`, y en el panel `qseqid` es el `gene_id`,
no el `locus_tag`. El script llama `locus_tags` a esa columna. Como `gene_id` es
una clave estable y `gene_names()` también busca `locus_tag or gene_id`, el
resultado sirve para K-12 y para el panel, pero el nombre `locus_tag` es
engañoso cuando `locus_tag` está vacío. Antes de publicar un `gold_<taxid>` de
otra especie conviene:

- hacer el join explícito `gene_id → locus_tag` y guardar ambos;
- no asumir que el identificador de ComplexPortal es un b-number;
- probar con un genoma y verificar 3–5 complejos conocidos.

## Circularidad

- El gold K-12 se usó para construir el modelo K-12 (`gold_complexes=273`,
  `n_pos=1773` en `model_k12.json`) y también para validarlo. Evaluar con ese
  mismo gold es circular; no es una validación independiente.
- La transferencia K-12 → 400 y el modelo genómico K-12 son **anotaciones
  inferidas** (`documentado_inferido`, `ECO:0005544`), no evidencia
  experimental de la especie.
- Añadir catálogos ComplexPortal de otras especies no elimina la circularidad
  por sí solo:
  - si se usan para ajustar el modelo, dejan de ser test;
  - si se usan para test, hay que excluirlos del entrenamiento y separar por
    especie o por complejo, no por pares de genes;
  - un catálogo parcial no puede validar un complejo completo.
- El canal `database` de STRING (que incluye complejos curados) ya se excluye
  por defecto al ajustar; mantenerlo así si se compara contra ComplexPortal.
- El inventario en sí no introduce circularidad: solo registra qué catálogo
  local existe. El riesgo aparece al usarlo como entrenamiento y como test a la
  vez.

## Coste

| tarea | coste | base |
|---|---|---|
| ejecutar el inventario | < 1 s, 0 red | 401 `genes.tsv` + 400 JSON |
| comprobar el listado remoto | 1 request (~76 KB) | listado FTP |
| descargar 1 catálogo | ~0,5 MB, < 1 min | K-12: 556 KB |
| mapear 1 genoma | CPU segundos | 400 genomas ya tienen `functional.tsv` |
| descargar los 400 | **no aplica** | 0/400 publicados; serían 404 |
| curaduría manual | el coste dominante; no medido | requiere revisar composición y evidencia |

Hipotético: si ComplexPortal publicara los 400, serían ~200 MB de TSV
(400 × 0,5 MB) y minutos de CPU para mapear, pero semanas/meses de revisión
manual. El cuello no es la descarga. Con la cobertura actual el coste de
“ComplexPortal para las 400” es cero porque no hay catálogos que descargar; el
coste real es la transferencia/ortología y su validación independiente.

## Reproducir

```bash
python scripts/complexes/inventory_complexportal_400.py
python -m pytest tests/test_complexportal_inventory.py -q
```

El script no descarga nada y no toca `app/` ni el navegador. Los tests usan
árboles sintéticos y comprueban, entre otras cosas, que un genoma sin taxid
queda `not_available` y que K-12 no inventa `83333` si falta el TSV crudo.

## Límites

- La foto remota es de una fecha; el inventario local no consulta la red.
- `taxid` es el del ensamblado. Un catálogo de especie podría no coincidir con
  la cepa; no se fuerza la equivalencia.
- `catalog_missing` no distingue “no existe en ComplexPortal” de “existe pero no
  está descargado”; para eso está la sección remota.
- Un `gold_<taxid>` con complejos incompletos no es un catálogo completo.
- No se ha medido la calidad de los catálogos porque, fuera de K-12, no hay
  ninguno.
