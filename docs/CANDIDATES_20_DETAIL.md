# Detalle de 20 candidatos

Archivo: `data/metabolism/CANDIDATES_20_DETAIL.tsv`

## KARI: 10 candidatos

- Brecha: `(2R,3R)-2,3-dihidroxi-3-metilpentanoato`.
- Dominio: `KARI_N`.
- Homólogo: `P06168`.
- Producto anotado: `A8ZTR0` KARI.
- Identidad: 55–60%.
- Cobertura: ~1.0.
- Consistencia: alta.

**Lectura:** candidato fuerte a cierre de KARI. No confirmado.

## L-cysteinylglycine: 10 candidatos

- Brecha: `L-cysteinylglycine`.
- Dominio: `G_glu_transpept`.
- Homólogo: `P15557` Acylase ACY 1.
- Identidad: 31–35%.
- Cobertura: ~1.0.
- Producto anotado: `A6T9C8` oxamate amidohidrolasa.

**Lectura:** hay conflicto entre el dominio/homólogo (acilasa/GGT-like) y la anotación del producto. No se debe cerrar como GGT sin revisión.

## Prioridad

1. KARI: revisar contexto y reacción.
2. L-cys: revisar si la acilasa puede hidrolizar L-cysteinylglycine.
3. No marcar `confirmed`.
