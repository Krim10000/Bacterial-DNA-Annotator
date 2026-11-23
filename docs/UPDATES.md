# Actualizaciones

La aplicación no consulta internet al arrancar. Para saber si hay una versión
nueva hay tres vías, de la más cómoda a la más manual.

## 1. Suscribirse a las releases de GitHub

1. Abre https://github.com/Krim10000/Bacterial-DNA-Annotator
2. Arriba a la derecha, pulsa **Watch** → **Custom**.
3. Marca **Releases**.
4. Guarda.

GitHub enviará un correo cuando se publique una release nueva.

También puedes usar el feed RSS:

```text
https://github.com/Krim10000/Bacterial-DNA-Annotator/releases.atom
```

## 2. Mirar la versión y la release

- El **pie de la GUI** muestra la versión activa: `Anotador v1.1.3`.
- El **badge Release** del README muestra la última versión publicada.
- La página de releases:
  https://github.com/Krim10000/Bacterial-DNA-Annotator/releases

Si el badge o la release muestran una versión mayor que el pie de la GUI, hay
una actualización.

## 3. Actualizar una instalación desde el código

```bash
cd /ruta/al/repositorio
git pull origin v1.0
bash install.sh
```

`install.sh` reutiliza el entorno virtual, actualiza dependencias y vuelve a
crear el lanzador. No borres `data/` ni `out_fl/` si quieres conservar los
resultados.

Para una instalación desde ZIP:

1. Descarga el ZIP de la última release.
2. Descomprímelo en otra carpeta.
3. Copia `data/` y `out_fl/` si los necesitas.
4. Ejecuta `bash install.sh`.

## 4. Zenodo

Cada versión archivada tiene DOI propio. El **DOI de concepto** siempre apunta a
la última versión:

- Concepto: https://doi.org/10.5281/zenodo.24469731
- Versión 1.1.3: https://doi.org/10.5281/zenodo.24469732

## Lo que no hace la app

- No busca actualizaciones automáticamente.
- No envía datos del genoma a GitHub ni a Zenodo.
- No actualiza archivos por sí sola.

Una comprobación manual desde la GUI sería posible, pero enviaría una petición a
`api.github.com`; por ahora se prefiere que sea una decisión explícita del
usuario mediante el botón de GitHub o el feed RSS.
