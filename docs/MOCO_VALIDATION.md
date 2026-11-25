# Validación MoCF / MoaE

## Candidatos

Total: **18**, no 13.

| grupo | candidatos | brecha | dominio |
|---|---:|---|---|
| MoaE | 5 | molybdopterin | MoaE |
| MoCF | 13 | adenylyl-molybdopterin | MoCF_biosynth |

## Contexto génico

- **13/18** están en un clúster con genes `moa`, `mob`, `moe`, `moc`, `mod` o `mog`.
- **5/18** están aislados.

### MoaE

- 4/5 con contexto mo.
- query mayormente `Y903_SYNY3` no caracterizado.
- vecinos `MoaA`, `MoaC`, `MoaE`.
- gap `molybdopterin`.

### MoCF

- 9/13 con contexto mo.
- query anotado como `MoaA2`, `MoaB`, `MoaC` o `MoaE`.
- gap `adenylyl-molybdopterin`.
- dominio `MoCF_biosynth` es amplio.

## Interpretación

- MoaE: candidato razonable en los 4 con contexto.
- MoCF: hay señal, pero el dominio es amplio y algunas anotaciones no coinciden con la brecha exacta.
- Los 5 aislados necesitan más revisión.
- No hay confirmación experimental.

## Veredicto

| subgrupo | veredicto |
|---|---|
| MoaE + contexto mo | `candidate` |
| MoaE aislado | `candidate_weak` |
| MoCF + contexto mo | `candidate` |
| MoCF aislado | `candidate_weak` |

Ninguno `confirmed`.
