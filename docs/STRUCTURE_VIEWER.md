# Visor 3D de estructuras

La GUI tiene una pestaña **🧊 Estructura 3D** y un visor embebido en las fichas
de estructura/ESMFold.

## Cómo se usa

1. Pestaña **Estructura 3D**:
   - subir un PDB/mmCIF, o
   - indicar una ruta local.
2. En cierre por estructura/ESMFold:
   - se dibuja el PDB generado o seleccionado.
3. En cinco filtros/evidencia:
   - se dibuja el PDB evaluado.

Dependencias: `py3Dmol` y `stmol` (declaradas en `requirements-annotator.txt`).

## Qué muestra

- opción **Colorear por B-factor/RMSX**; útil para modelos por slice de RMSX.

- proteína en cartoon, color por cadena/conformación;
- heteroátomos en stick;
- descarga del PDB.

## Qué NO significa

Una estructura 3D:

- no confirma función;
- no confirma sustrato;
- no confirma interfaz ni estequiometría;
- no cierra una brecha.

Una estructura experimental, AFDB o ESMFold se muestran igual; la diferencia
está en la procedencia, no en la imagen.

## Límites

- no calcula pLDDT/ipTM por sí mismo;
- no valida la geometría;
- no muestra density maps;
- no sustituye PyMOL/ChimeraX;
- no prueba interacción.
