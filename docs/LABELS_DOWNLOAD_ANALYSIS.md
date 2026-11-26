# Descarga y análisis de etiquetas externas

Fecha: 2026-11-26.

## Descargado

| Recurso | Tamaño | Cobertura |
|---|---:|---|
| `go-basic.obo` | 31 MB | global |
| GOA `ECOLI-uniprot` | 931 KB | K-12 |
| GOA `ECOLI-mod` | 875 KB | K-12 |
| GOA `BACSU-uniprot` | 433 KB | B. subtilis |
| KEGG `gene→KO/pathway` | 41 MB | 339/442 genomas |

KEGG:

```text
442 genomas
339 con código KEGG
338 códigos únicos
103 sin match
1 error: gjp pathway (HTTP 400)
```

## Cobertura de etiquetas

### K-12

```text
genes:              4.534
GAF IDs:            4.095
genes con etiqueta:   570
operones:           2.544
operones con label:   440
KEGG operones:        935
KEGG pathways:        140
```

### B. subtilis

```text
genes:              4.134
GAF IDs:            2.348
genes con etiqueta:   372
operones:           1.334
operones con label:   244
KEGG operones:        471
KEGG pathways:        129
```

## Perfil K-12 por sigma

Fuente de dianas: `NetworkSigmaTU.tsv`.
Etiquetas: léxico + GO + KEGG.

| Sigma | Dianas | Señal de proceso | Señal KEGG |
|---|---:|---|---|
| sigma24 | 37 | envelope 8 | LPS biosynthesis, CAMP resistance, glycerophospholipid |
| sigma28 | 14 | flagellum 6 | chemotaxis, flagellar assembly, two-component |
| sigma32 | 46 | heat_shock 12 | integrated stress response, RNA degradation, aminoacyl-tRNA |
| sigma38 | 116 | stationary 13, envelope 15 | glycolysis, carbon metabolism, biofilm |
| sigma54 | 30 | nitrogen 8 | two-component, quorum sensing, amino acid metabolism |

## Perfil B. subtilis (descriptivo)

TSS gold mapeados a operones; **no es validación independiente**.

| Sigma | TSS | Señal |
|---|---:|---|
| sigma28 | 51 | flagellum 9, sporulation 10 |
| sigma32 | 10 | heat_shock 1, envelope 2 |
| sigma38 | 125 | sporulation 39, nitrogen 22, envelope 18 |

El mapeo TSS→operón y la mezcla BSGatlas/DBTBS limitan la interpretación.

## Hallazgos

1. KEGG aporta señal clara para sigma24, sigma28, sigma32, sigma38 y sigma54 en K-12.
2. GO aporta cobertura, pero con GAF sin IEA el recall de procesos es bajo.
3. El léxico sigue siendo mejor para heat shock/flagelo; KEGG para metabolismo.
4. Las tres fuentes son complementarias.
5. B. subtilis no tiene un sigmulón independiente local; su perfil es descriptivo.

## Límites

- K-12 usa el mismo universo de operones para dianas y etiquetas.
- GO GAF sin IEA reduce cobertura.
- KEGG no cubre 103 genomas.
- B. subtilis gold es BSGatlas/Nicolas.
- No hay negativos.

## Siguiente

1. Combinar léxico + GO + KEGG en `labels`.
2. Propagar GO.
3. Añadir sigmulón B. subtilis.
4. Evaluar prior contextual sobre sigma.
5. Aplicar a 339 genomas KEGG y 103 con fallback.

## GO con y sin IEA

Comparación de GAF:

| Genoma | Genes sin IEA | Genes con IEA | Nuevos genes |
|---|---:|---:|---:|
| K-12 | 4.095 | 4.482 | 163 |
| B. subtilis | 2.348 | 3.439 | 212 |

Etiquetas que aparecen solo al incluir IEA:

| Genoma | nitrogen | envelope | stationary | heat_shock | flagellum | sporulation |
|---|---:|---:|---:|---:|---:|---:|
| K-12 | 35 | 74 | 44 | 9 | 3 | 0 |
| B. subtilis | 102 | 55 | 2 | 8 | 5 | 54 |

IEA añade principalmente:

- metabolismo de nitrógeno;
- envoltura;
- stationary;
- esporulación en B. subtilis.

Por eso `labels` tiene `allow_iea` como opción explícita y no lo activa por defecto.
