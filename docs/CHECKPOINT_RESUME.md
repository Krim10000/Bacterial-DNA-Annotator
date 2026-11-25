# Checkpoint y reanudación

## Lo que no se puede

Un kernel de Kaggle/Colab **no puede**:

- descargar automáticamente a tu disco local;
- reiniciarse solo tras morir;
- conservar `/kaggle/working` después de cerrar la sesión.

## Lo que sí se puede

- guardar checkpoint cada N proteínas;
- copiarlo a almacenamiento persistente;
- reanudar sin repetir PDB ya hechos;
- borrar del lote lo completado.

## Opción recomendada: Colab + Google Drive

```python
from google.colab import drive
drive.mount('/content/drive')
CKPT = Path('/content/drive/MyDrive/dna_poc_esmfold')
CKPT.mkdir(parents=True, exist_ok=True)
```

Checkpoint:

```python
import shutil, datetime

def checkpoint():
    zip_base = str(OUT)
    shutil.make_archive(zip_base, "zip", root_dir=OUT.parent, base_dir=OUT.name)
    dst = CKPT / f"esmfold_{datetime.datetime.now():%Y%m%d_%H%M}.zip"
    shutil.copy(zip_base + ".zip", dst)
    print("checkpoint:", dst)
```

Llamarlo cada 20 proteínas desde `run_rows`.

## Kaggle

`/kaggle/working` no sobrevive.

Opciones:

1. Mini-lotes de 10–20 proteínas.
2. Descargar ZIP manualmente.
3. Subir el ZIP como dataset en la siguiente sesión.
4. Descomprimir en `OUT`; `run_rows` salta PDB existentes.

## Reanudar

```python
# Subir el checkpoint como dataset y descomprimir
!unzip -q /kaggle/input/checkpoint/esmfold_*.zip -d /kaggle/working/
```

El runner debe saltar:

```python
if target.exists() and target.stat().st_size > 1000:
    continue
```

## Borrar completados del TSV

```python
pending = [r for r in rows if not (OUT / "pdb" / f"{r['uniprot']}.pdb").exists()]
```

## Regla

- No lanzar 257 en una sola celda.
- Mini-lotes de 10–20.
- Checkpoint después de cada mini-lote.
- Descargar/copiar antes de cerrar.
