# Empaquetado Bioconda + Galaxy

Estado: **preparado localmente, sin publicar y sin usar red**. Esta carpeta y
este documento dejan la receta y el wrapper listos para revisión, pero no se
ha enviado nada a Bioconda ni al ToolShed de Galaxy.

## Archivos

| Archivo | Para qué |
| --- | --- |
| `packaging/bioconda/meta.yaml` | Receta Bioconda del paquete Python. |
| `packaging/bioconda/build.sh` | Instalación con `pip` y ajuste del paquete `ghost`. |
| `packaging/galaxy/bacterial_dna_annotator.xml` | Wrapper mínimo de Galaxy. |
| `tests/test_packaging.py` | Validación offline de YAML/XML y coherencia. |

## Fuente de verdad

- `pyproject.toml`: nombre `dna-poc-annotator`, versión `1.1.3`, licencia MIT,
  dependencias y entry point `dna-annotate = webapp.pipeline:main`.
- `requirements-annotator.txt`: se revisó para detectar dependencias de
  ejecución que no están en `pyproject.toml`.
- `install.sh`: confirma que el entorno soportado es Python >= 3.10 y que los
  binarios externos (tRNAscan-SE, barrnap, DIAMOND, HMMER, etc.) son opcionales.
- `LICENSE`: licencia MIT del código.

La receta usa las dependencias reales de `pyproject.toml` y añade `requests`,
que `webapp/esmfold_tools.py` y `webapp/structure_evidence.py` importan sin
condición pero que falta en `pyproject.toml`. Las dependencias opcionales
(`scikit-learn`, `orcapy`, `ViennaRNA`, etc.) no se incluyen en el paquete
base; se documentan como extras.

## Bioconda

### Fuente local, sin SHA256 inventado

`meta.yaml` usa `source.path: ../..` para construir desde este checkout. Es
deliberado: no hay red y no se inventa un `sha256`. Para una release real hay
que crear el tarball, calcular su hash con `sha256sum` y sustituir `path` por
`url` + `sha256`. La receta incluye un comentario con ese paso, pero no un
hash falso.

### Construcción (requiere conda-build y red, no ejecutada aquí)

En un equipo con `conda-build` y los canales configurados:

```bash
conda build packaging/bioconda --channel conda-forge --channel bioconda
```

No se ejecutó en este entorno porque no hay `conda`, `mamba` ni
`conda-build`, y la tarea prohíbe usar red. Por tanto la receta está
**validada estructuralmente, no construida**.

### Verificación offline ejecutada

```bash
.venv/bin/python -m pytest tests/test_packaging.py -q

.venv/bin/python -c "import pathlib, yaml; d=yaml.safe_load(pathlib.Path('packaging/bioconda/meta.yaml').read_text()); print(d['package'])"

.venv/bin/python -c "import xml.etree.ElementTree as ET; r=ET.parse('packaging/galaxy/bacterial_dna_annotator.xml').getroot(); print(r.tag, r.get('id'), r.get('version'))"
```

La suite comprueba, entre otras cosas: YAML y XML parseables, versión y
nombre contra `pyproject.toml`, entry point, dependencias declaradas, licencia,
ausencia de `sha256`, sintaxis de `build.sh` y que las opciones usadas por el
wrapper existen en `python -m webapp.pipeline --help`.

### Prueba local del wheel (sin red)

```bash
.venv/bin/python -m pip wheel . --no-deps --no-build-isolation -w dist/packaging

ROOT=$PWD
SITE=$(mktemp -d)
"$ROOT/.venv/bin/python" -m pip install --no-deps --target "$SITE" dist/packaging/*.whl
cp -r "$ROOT/scripts/ghost/ghost" "$SITE/ghost"
(cd /tmp && PYTHONPATH="$SITE" "$ROOT/.venv/bin/python" -c "import webapp.pipeline")
rm -rf "$SITE"
```

El `cp` reproduce lo que hace `build.sh`. Es necesario porque
`pyproject.toml` empaqueta `webapp*` y `modules*`, pero no
`scripts/ghost`; `webapp/annotator.py` importa `ghost` al cargarse.

## Galaxy

El wrapper declara el paquete conda `dna-poc-annotator` y ejecuta:

```bash
python -m webapp.pipeline input.fna --outdir out \
  --threads "${GALAXY_SLOTS:-1}" \
  --no-trna --no-rrna \
  --no-functional --no-membrane --no-oric --no-metabolic
```

Salidas: `<genoma>.gff3`, `<genoma>.candidates.tsv` y
`<genoma>.summary.json`. El wrapper es mínimo a propósito.

### Limitaciones del wrapper

- Desactiva tRNA/rRNA porque `tRNAscan-SE` y `barrnap` no están en el paquete
  Python y el pipeline falla si no existen.
- Desactiva DIAMOND/Swiss-Prot, membrana, ORCA y la red metabólica porque
  necesitan el árbol `data/` y binarios externos que no se empaquetan.
- Por lo anterior **no genera `gaps.tsv`** ni tablas de evidencia estructural.
- No se ejecutó con `planemo` ni Galaxy: no hay `planemo`/Galaxy instalados y
  la tarea prohíbe red. La validación es de XML y de opciones del CLI.

Para una integración real: instalar el paquete conda en el entorno de la
herramienta, copiar el XML al directorio de herramientas de Galaxy y recargar
la configuración. El ToolShed no se ha tocado.

## Publicación (pendiente, no hecha)

Bioconda:

1. Crear un tarball de release inmutable.
2. Calcular su `sha256` real con `sha256sum` y ponerlo en `source`.
3. Añadir `extra.recipe-maintainers` con cuentas que acepten mantenerlo.
4. Pasar `conda build` y el lint de `bioconda-recipes` en un entorno con red.
5. Abrir el PR a `bioconda-recipes`.

Galaxy:

1. Instalar `planemo` y validar con `planemo lint` y `planemo test` usando
   datos de prueba reales.
2. Publicar en el ToolShed solo después de que el paquete conda exista.

No se ha publicado ni se ha usado red en ningún paso.

## Limitaciones conocidas

- La receta no se construyó: no hay `conda-build` en el entorno.
- `source.path` no es una fuente inmutable; es solo para prueba local.
- No hay `sha256` porque no se puede calcular sin un tarball de release.
- No se verificó la disponibilidad de `py3dmol` ni `stmol` en los canales
  conda; hay que comprobarlo antes de enviar la receta.
- `webapp/annotator.py` depende de `ghost`, que `pyproject.toml` no empaqueta.
  `build.sh` lo copia a site-packages; lo ideal es corregir el empaquetado
  aguas arriba.
- `requests` falta en `pyproject.toml`; la receta lo añade.
- Los datos (`data/`), modelos, binarios externos y resultados no se incluyen.
  La licencia MIT cubre el código; revisar `LICENSE-SCOPE.md`, `DATA_NOTICE.md`
  y `THIRD_PARTY_LICENSES.md` para datos y componentes de terceros.
- Las pruebas son estáticas (parseo y coherencia) más `--help`; no sustituyen
  una construcción conda ni una ejecución Galaxy de extremo a extremo.
