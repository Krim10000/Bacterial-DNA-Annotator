# Validación de KARI

## Candidatos

- 10 genomas.
- Brecha: `(2R,3R)-2,3-dihidroxi-3-metilpentanoato`.
- Rhea master: **13493**.
- EC: **1.1.1.86**.
- Dirección candidata: `RL`.
- Dominio: `KARI_N`.
- Homólogo: `P06168` (ILV5_YEAST, KARI).
- Query anotado: `A8ZTR0` ILVC_DESOH, KARI.

## Evidencia

| genoma | longitud | identidad | cobertura | contexto ilv |
|---|---:|---:|---:|---|
| GCF_000021905.1 | 351 | 57,2% | 0,99 | sí |
| GCF_001017655.1 | 366 | 55,4% | 0,94 | **no** |
| GCF_000455605.1 | 348 | 58,3% | 0,99 | sí |
| GCF_009688965.1 | 351 | 55,3% | 0,99 | sí |
| GCF_048565275.1 | 351 | 55,3% | 0,99 | sí |
| GCF_003966975.1 | 355 | 55,3% | 0,98 | sí |
| GCF_000199675.1 | 351 | 55,3% | 0,99 | sí |
| GCF_018292205.1 | 348 | 60,5% | 1,00 | sí |
| GCF_003815875.1 | 349 | 56,8% | 0,99 | sí |
| GCF_023921165.1 | 348 | 57,9% | 1,00 | sí |

**9/10** están en un clúster génico con `ilvD`, `ilvB`, `ilvH`, `ilvA`, `leu` o `cimA`.

## Interpretación

- La anotación del query, el dominio, el homólogo y el contexto génico coinciden.
- Es el **candidato más fuerte** de los 508.
- El caso aislado `GCF_001017655.1` sigue siendo candidato, pero con menos contexto.
- No hay validación experimental.
- No se debe marcar `confirmed`.

## Veredicto

`candidate_strong` para los 9 con contexto ilv.  
`candidate` para el aislado.

Ninguno es `confirmed`.
