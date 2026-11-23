# Hipótesis precalculadas en GUI y pipeline

## Dataset

`UNCHAR_AFDB_PRESENT_filters.tsv`:

- 1.192 accesiones con AFDB.
- 195 con F2.
- 176 con F2+F3+F4.
- 3.374 clusters propagados.
- 0 errores.

## GUI

Nuevo modo:

**🔬 Hipótesis estructurales**

Muestra:

- tabla de 1.192;
- filtros por confianza/EC/cofactor;
- detalle F1–F4;
- revisión;
- descarga TSV/JSON/PDB.

No modifica la anotación base.

## Pipeline

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --precomputed-hypotheses data/metabolism/targets/UNCHAR_AFDB_PRESENT_filters.tsv
```

Se cargan solo las filas del genoma procesado y se añaden a
`*.structural_hypotheses.tsv/.json`.

## Regla

- `supported_hypothesis` sigue siendo hipótesis.
- No hay `confirmed`.
- No cierra brechas.
- F5 propio no se inventa.

## Hallazgos

Casos con tres fuentes:

- P9WN14: glucósido hidrolasa probable (`1gro`, Ca/Mg).
- P0DUF3: metiltransferasa probable (`3exr`, `1aqj`, MTA/SAH).
- P76585: posible enzima (`1nww`, `4pmh`).
- P37443: posible enzima Zn/Mn.
- P71353: posible ATPasa/ligasa (`2om2`, ATP/VGL).
- P9WIY6: posible nucleótido-enzima (`1d4e`, CMP/K).
