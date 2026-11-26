# Inventario de sigma por especie

## Regla

No existe un panel global de seis sigma. Cada especie declara sus propios
factores.

El inventario **no asigna** un subtipo funcional exacto a partir del nombre.
Solo conserva:

```text
native_name
protein_id / locus_tag
product
family: sigma70_like | sigma54_like | ecf_like | unknown
evidence
confidence
```

`family` es una pista amplia de dominio/anotación, no una equivalencia con
`sigma32`, `sigma38`, etc.

## API

```python
from webapp.sigma_inventory import load_inventory, summary

factores = load_inventory("data", "subtilis")
s = summary(factores)
s["n_factors"], s["families"], s["names"]
```

CLI:

```bash
python scripts/sigma/00_sigma_inventory.py --genome subtilis \
  --out data/sigma/inventory/subtilis.tsv
```

## Ejemplos

```text
E. coli K-12: 7 sigma candidatos
B. subtilis: 19 sigma candidatos
```

No se fuerza que ambas especies tengan los mismos.

## Límites

- Se basa en producto/nombre anotado; no sustituye HMMER/Pfam.
- `anti-sigma`, reguladores sigma-dependientes y subunidades del core se excluyen.
- `sigKC`, `sigO`, etc. se conservan con su nombre nativo.
- `family` no es una llamada funcional.
- La detección por dominios y la ortología quedan como siguiente capa.
