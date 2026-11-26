# GO MF/CC por accesión UniProt

## Objetivo

Anotar **Molecular Function (MF)** y **Cellular Component (CC)** en las 442
especies sin depender de que exista un GAF de esa especie.

## Principio

La anotación se transfiere por **accesión UniProt**, no por ortología:

```text
gen → proteína → accesión UniProt → GO MF/CC
```

Si una proteína de una especie nueva tiene una accesión con GO, se usa esa
anotación. No importa que el GAF sea de otra cepa o especie.

## Insumos

```text
data/regulation/<genoma>/genes.tsv
  columna acc (UniProt)
data/labels/external/goa/by_code/*.gaf.gz
data/labels/external/go/go-basic.obo
```

## Salidas

```text
data/labels/go_mf_cc/<genoma>/gene_go.tsv
data/labels/go_mf_cc/<genoma>/go_summary.json
```

`gene_go.tsv`:

```text
genome
gene_id
gene_name
acc
go_id
go_name
aspect          F | C
namespace       molecular_function | cellular_component
evidence        IEA | IDA | EXP | ...
evidence_class  electronic | curated | experimental
source          GAF file
taxon
assigned_by
```

## Cobertura actual

```text
442 genomas
1.522.305 genes
275.094 genes con al menos un MF/CC
1.533.774 filas GO
```

K-12:

```text
2.788 genes anotados
9.807 términos MF
6.240 términos CC
```

B. subtilis:

```text
2.904 genes anotados
7.161 términos MF
3.831 términos CC
```

## Pipeline

```text
1. Construir índice SQLite accesión → GO (una vez).
2. Para cada genoma, leer genes y accesiones.
3. Cruzar accesión → GO.
4. Escribir gene_go.tsv y go_summary.json.
5. El navegador dibuja una pista GO MF/CC.
```

Comandos:

```bash
python scripts/labels/14_build_go_mf_cc.py
python scripts/labels/15_add_bsu_go.py
```

Módulo:

```text
modules/go_mf_cc/module.py
order = 31
```

## Especie nueva

Si el usuario integra un genoma:

1. `genes.tsv` debe traer `acc` de DIAMOND/Swiss-Prot.
2. Se ejecuta el módulo `go_mf_cc`.
3. Se escribe en `data/labels/go_mf_cc/<genoma>/`.
4. El navegador muestra la pista.
5. Si no hay accesión o GO: `no_data`.

No se descarga nada automáticamente.

## GUI

En el navegador:

```text
casilla "GO MF/CC"
```

Dibuja una línea por gen anotado:

- azul: MF y CC;
- violeta: solo MF/CC parcial;
- tooltip: número de términos y GAF fuente.

En la pestaña **Etiquetas**:

- genes anotados;
- términos MF/CC;
- fuente por GAF;
- tabla y descarga.

## Límites

- Es una transferencia por accesión, no una validación experimental.
- `IEA` es evidencia electrónica.
- No es ortología: si no hay accesión, no hay anotación.
- No es BP todavía.
- La fuente siempre se muestra; no se oculta.

## Cita

GO: `go-basic.obo` releases/2026-07-26, GO Consortium, CC BY 4.0.
GOA/UniProt: UniProt Consortium, CC BY 4.0.
Ver `docs/LABELS_CITATIONS.md`.
