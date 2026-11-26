# Piloto de etiquetas K-12 / B. subtilis

## Qué se construyó

`webapp/labels.py` + `scripts/labels/01_build_labels.py`.

Léxico controlado de procesos:

```text
heat_shock
flagellum
nitrogen
stationary
envelope
sporulation
```

No usa gold de sigma. Etiqueta genes/operones por producto y nombre, excluyendo
el propio factor sigma y anti-sigma.

## Cobertura

### K-12

```text
operones: 2.544
heat_shock:  58 (2,3 %)
flagellum:   22 (0,9 %)
nitrogen:    58 (2,3 %)
stationary:  22 (0,9 %)
envelope:   178 (7,0 %)
sporulation:250 (9,8 %)
```

### B. subtilis

```text
operones: 1.334
heat_shock:  13 (1,0 %)
flagellum:   17 (1,3 %)
nitrogen:    43 (3,2 %)
stationary:  11 (0,8 %)
envelope:    24 (1,8 %)
sporulation:295 (22,1 %)
```

## Validación K-12 contra sigmulón curado

Fuente: `data/regulondb/NetworkSigmaTU.tsv`.

| Sigma | Etiqueta | Diana | Mapeados | Con etiqueta | Recall | Fondo | OR | FDR |
|---|---|---|---:|---:|---:|---:|---:|---:|
| sigma32 | heat_shock | 94 | 46 | 9 | 19,6 % | 2,3 % | 12,2 | 2e-6 |
| sigma28 | flagellum | 29 | 14 | 5 | 35,7 % | 0,9 % | 82,1 | <1e-6 |
| sigma54 | nitrogen | 54 | 30 | 4 | 13,3 % | 2,3 % | 7,0 | 0,005 |
| sigma38 | stationary | 201 | 116 | 10 | 8,6 % | 0,9 % | 19,0 | <1e-6 |
| sigma24 | envelope | 82 | 37 | 5 | 13,5 % | 7,0 % | 2,1 | 0,11 |

Lectura:

- `heat_shock`, `flagellum`, `nitrogen` y `stationary` están enriquecidos
  frente al fondo.
- `envelope` no alcanza significación con este léxico: es demasiado amplio.
- El recall es bajo: el léxico es incompleto; no sustituye GO/KEGG.
- La validación es de **proceso**, no de identidad de sigma.

## B. subtilis

No hay un sigmulón independiente local equivalente a `NetworkSigmaTU`.
Los resultados son **descriptivos**, no validación.

## Límites

- Léxico, no ontología completa.
- Cobertura funcional dependiente de anotación.
- `envelope` es demasiado general.
- No hay GO/KEGG cargados.
- No hay negativos.
- La validación K-12 usa el mismo universo de operones, no un held-out.

## Siguiente paso

1. Sustituir el léxico por GO/KEGG/EC.
2. Añadir términos padre y propagación.
3. Calcular FDR por operón.
4. Validar con regulones independientes.
5. Aplicar a 400 con cobertura y abstención.
