# Cobertura GOA descargada

## Qué se descargó

```text
193 GAF de GO (UniProt/mod)
349 MB
182 taxid únicos
taxdump NCBI para mapear nombres
```

Archivos:

```text
data/labels/external/goa/by_code/
data/labels/external/goa/goa_map.tsv
data/labels/external/goa/goa_fetch_summary.json
```

## Mapeo al catálogo

```text
442 genomas
  9 con GOA por especie
433 no_data
```

Motivo:

- GO GAF cubre 193 especies, mayoría eucariotas;
- nuestros 442 son bacterias diversas;
- el taxid de cepa no coincide con la cepa tipo;
- el mapeo por nombre de especie recupera solo 9.

K-12 y B. subtilis están cubiertos. El resto queda `no_data`, no se inventa.

## Con y sin evidencia electrónica

GAF no tiene un campo "microscopía electrónica". Tiene **códigos de evidencia**:

```text
experimental: EXP, IDA, IPI, IMP, IGI, IEP
curado:       TAS, NAS, IC
electrónico:  IEA
```

`IDA`/`EXP` pueden incluir microscopía, pero GAF no distingue la técnica.

Vistas posibles:

```text
allow_iea=False  → excluye IEA (curado/experimental)
allow_iea=True   → incluye IEA
experimental     → solo EXP/IDA/IPI/IMP/IGI/IEP
```

La opción `allow_iea` ya existe en `modules/labels`.

## Límite

No se puede decir "GOA con microscopía electrónica" desde GAF: solo
"experimental" vs "IEA". Si se necesita la técnica, hay que ir al artículo
PMID de cada anotación.
