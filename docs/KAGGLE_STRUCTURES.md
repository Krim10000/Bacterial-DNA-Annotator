# Estructuras ESMFold de Kaggle (`data/Kaggler`)

## Qué contiene

ZIP con:

```text
<lote>/metrics.tsv
<lote>/pdb/<uniprot>.pdb
```

`metrics.tsv` trae:

```text
uniprot status length plddt ptm seconds error
```

Son estructuras **predichas por ESMFold**, no experimentales.

## Ingesta

```bash
python scripts/structure/ingest_kaggle_zips.py --dry-run
python scripts/structure/ingest_kaggle_zips.py
```

- escanea todos los ZIP;
- cataloga por accesión;
- no duplica una accesión;
- si dos lotes tienen la misma accesión, conserva la de mayor pLDDT;
- copia el PDB a `data/metabolism/targets/esmfold_kaggle/`;
- catálogo: `data/metabolism/targets/esmfold_catalog.tsv`.

Es reejecutable: si añades `set_006.zip`, se ingiere sin repetir lo anterior.

## Uso en la GUI

La ficha del gen busca, en orden:

1. `pdb` declarado;
2. AFDB en `afdb25k/`;
3. ESMFold en `esmfold_kaggle/`.

Si hay PDB, muestra descarga y visor 3D. Si no hay PAE, lo dice.

## Límites

- ESMFold no da PAE.
- pLDDT/pTM bajos no significan que la proteína no exista.
- Una estructura predicha no confirma función, interfaz ni ensamblaje.
- No es estructura cuaternaria.
