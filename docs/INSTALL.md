# Instalación y uso

🌐 **Español** · [English](INSTALL.en.md)

## Forma rápida (Ubuntu)

```bash
git clone <URL-del-repositorio> DNA_TOKEN_BAC
cd DNA_TOKEN_BAC
bash install.sh --with-data
```

Esto crea `.venv`, instala dependencias, descarga Swiss-Prot y Rhea (~120 MB) y
añade **«Anotador de ADN»** al menú de aplicaciones. Después abre la app desde
el menú o con `./lanzar.sh`.

Opciones:

```bash
bash install.sh                 # entorno + lanzador
bash install.sh --with-data     # + Swiss-Prot y Rhea
bash install.sh --with-tools    # + tRNAscan-SE, barrnap, HMMER, DIAMOND (sudo)
bash install.sh --full          # + dependencias completas del pipeline
bash install.sh --dry-run       # muestra lo que haría
bash install.sh --uninstall     # quita el lanzador
```

Requisitos: Python 3.10+, 8 GB de RAM (16 recomendado), 10 GB libres. **No hace
falta GPU** para anotar, buscar brechas ni usar la GUI. Ver
[`GPU_STATUS.md`](GPU_STATUS.md).

## Instalación manual

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-annotator.txt   # GUI + anotación
# o todo el pipeline:
.venv/bin/pip install -r requirements.txt
```

Herramientas del sistema (opcionales, la app detecta las que falten):

```bash
sudo apt install trnascan-se infernal hmmer diamond-aligner barrnap ncbi-blast+
```

## CLI

```bash
.venv/bin/python -m webapp.pipeline genoma.fna --outdir results/ --threads 8
# o, si instalaste el paquete:
dna-annotate genoma.fna --outdir results/ --threads 8
```

Salidas: GFF3, proteínas, candidatos, brechas, informe Markdown y
`summary.json`.

## GUI

```bash
.venv/bin/python -m streamlit run webapp/app.py
```

Abre <http://localhost:8501>. El selector de idioma está en la barra lateral:
español, inglés y chino simplificado. También `?lang=en` o `?lang=zh`.

Sin datos del proyecto puedes probar **🧬 Navegador de genoma → `demo`** (o el
botón **Explorar la demo del navegador**): 150 kb anotados desde cero, incluidos
en el repositorio.

## Datos

```bash
python -m webapp.setup_data --data-root data            # Swiss-Prot + Rhea
python -m webapp.setup_data --data-root data --pfam     # + Pfam (~600 MB)
python -m webapp.fetch_kegg --data-root data            # catálogo KEGG
```

Sin datos, la GUI avisa y muestra el comando. Los modos que exploran los 400
genomas necesitan que esos datos estén generados; ver
[`../data/README.md`](../data/README.md).

## Instalación editable (desarrollo)

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest tests -q
```

La instalación editable mantiene acceso a `scripts/ghost/` (motor de ORFs).
Un wheel aislado necesitaría empaquetar también ese motor; por ahora el paquete
incluye `webapp` y sus datos de interfaz (locales, componente, muestras).

## Problemas frecuentes

| síntoma | solución |
|---|---|
| `streamlit: command not found` | usa `.venv/bin/python -m streamlit run …` |
| Faltan Swiss-Prot/Rhea | `python -m webapp.setup_data --data-root data` |
| DIAMOND/tRNAscan aparecen en ❌ | instala las herramientas del sistema o usa `--with-tools` |
| No abre el navegador | entra a `http://localhost:8501` a mano |
| El menú no muestra la app | cierra sesión y vuelve a entrar, o ejecuta `update-desktop-database ~/.local/share/applications` |
| Error de GPU/ESMFold | es opcional; la anotación funciona en CPU |
