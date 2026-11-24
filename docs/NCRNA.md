# ncRNA con Rfam/Infernal (CPU)

El pipeline puede buscar ncRNA con `cmsearch` contra un subconjunto Rfam local.
No usa red ni GPU.

## Entrada

```text
data/rfam_models/bac_ncRNA.cm
```

Panel curado con:

- RNase P bacteriano;
- SRP bacteriano;
- tmRNA;
- 6S;
- CsrB/CsrC.

No es Rfam completo. Se genera desde un `Rfam.cm` local:

```bash
bash scripts/rrna/00_build_ncrna_panel.sh Rfam.cm
```

## Salida

```text
<genoma>.ncrna.tsv
```

Y las features `ncRNA` se añaden al GFF3.

Si faltan `cmsearch` o el panel:

```text
not_available
```

No se inventan ncRNA.

## Uso

```bash
python -m webapp.pipeline genoma.fna --outdir results/
# desactivar:
python -m webapp.pipeline genoma.fna --outdir results/ --no-ncrna
```

## Límites

- `cmsearch` detecta familia/modelo, no función ni expresión.
- Un `no_hit` no prueba ausencia: el panel es parcial.
- No sustituye Rfam completo, Infernal completo ni anotación experimental.
- No cierra brechas.
