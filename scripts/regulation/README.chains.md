# `scripts/regulation` — constructor de cadenas Rhea

Módulo del agente `chain-builder`. Contrato:
[`docs/AGENT_CHAIN_BUILDER.md`](../../docs/AGENT_CHAIN_BUILDER.md).
Salidas en [`data/chains/`](../../data/chains/README.md).

## Scripts

| Script | Qué hace |
|---|---|
| `build_chains.py` | lee Rhea/KEGG/gaps, construye el grafo, busca rutas y escribe `data/chains/` |
| `validate_chains.py` | valida TSV/JSON contra `rhea_reactions.pkl` y escribe `validation.json` |
| `run.sh` | ejecuta ambos en orden |

## Uso

```bash
python scripts/regulation/build_chains.py            # ~30 s
python scripts/regulation/validate_chains.py         # exit 0 si todo cuadra
bash scripts/regulation/run.sh                       # ambos
```

Opciones:

- `--root PATH` — raíz del proyecto (por defecto `.`).
- `--out-dir PATH` — directorio de salida (por defecto `data/chains`).
- `--include-unmapped-rhea` — usa las 17.177 reacciones maestras no-transporte en
  vez de las 4.215 mapeadas a KEGG (más cobertura, más ruido).
- `--strict-cofactors` — exige que todo co-sustrato esté en la cadena o en el
  pool; produce menos rutas pero más estrictas.

## Decisiones de diseño

1. **Solo reacciones maestras `UN` y no transporte.** Los IDs `RHEA:*` que
   aparecen en los TSV son los mismos que usan `refined_gaps`/`real_gaps`.
2. **Grafo colapsado metabolito↔metabolito.** Cada reacción produce aristas en
   ambos sentidos; es una simplificación declarada (no hay estequiometría).
3. **Moneda y pool.** Ver `CURRENCY` y `COFACTOR_POOL` en `build_chains.py`.
   Cambiarlos cambia las rutas; por eso están en el código y en `summary.json`.
4. **Sin aleatoriedad.** Todo se ordena antes de buscar y la salida es
   determinista (probado con dos `PYTHONHASHSEED` distintos).
5. **No escribe fuera de `data/chains/`.**

## Tests rápidos

No hay suite formal; el validador es el test de aceptación:

- falla si un `RHEA:*` no existe en el pickle, si `length` no cuadra, si la
  cadena de metabolitos no conecta con la reacción o si los totales del JSON no
  coinciden con los TSV.

Para probar el camino completo en un directorio temporal:

```bash
python scripts/regulation/build_chains.py --out-dir /tmp/chains_test
python scripts/regulation/validate_chains.py --data-dir /tmp/chains_test
```
