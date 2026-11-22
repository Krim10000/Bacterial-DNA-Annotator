🌐 [Español](INSTALL.md) · **English**

# Installation and usage

## Quick start (Ubuntu)

```bash
git clone <URL-del-repositorio> DNA_TOKEN_BAC
cd DNA_TOKEN_BAC
bash install.sh --with-data
```

This creates `.venv`, installs dependencies, downloads Swiss-Prot and Rhea
(~120 MB) and adds **«Anotador de ADN»** ("DNA annotator") to the applications
menu. Then open the app from the menu or with `./lanzar.sh`.

Options:

```bash
bash install.sh                 # entorno + lanzador
bash install.sh --with-data     # + Swiss-Prot y Rhea
bash install.sh --with-tools    # + tRNAscan-SE, barrnap, HMMER, DIAMOND (sudo)
bash install.sh --full          # + dependencias completas del pipeline
bash install.sh --dry-run       # muestra lo que haría
bash install.sh --uninstall     # quita el lanzador
```

Requirements: Python 3.10+, 8 GB of RAM (16 recommended), 10 GB free. **No GPU
is needed** to annotate, search for gaps or use the GUI. See
[`GPU_STATUS.en.md`](GPU_STATUS.en.md).

## Manual installation

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-annotator.txt   # GUI + anotación
# o todo el pipeline:
.venv/bin/pip install -r requirements.txt
```

System tools (optional, the app detects whichever are missing):

```bash
sudo apt install trnascan-se infernal hmmer diamond-aligner barrnap ncbi-blast+
```

## CLI

```bash
.venv/bin/python -m webapp.pipeline genoma.fna --outdir results/ --threads 8
# o, si instalaste el paquete:
dna-annotate genoma.fna --outdir results/ --threads 8
```

Outputs: GFF3, proteins, candidates, gaps, a Markdown report and
`summary.json`.

## GUI

```bash
.venv/bin/python -m streamlit run webapp/app.py
```

Open <http://localhost:8501>. The language selector is in the sidebar: Spanish,
English and Simplified Chinese. You can also use `?lang=en` or `?lang=zh`.

Without project data you can try **🧬 Genome browser → `demo`** (or the
**🧬 Explore the browser demo** button): 150 kb annotated from scratch, bundled
in the repository.

## Data

```bash
python -m webapp.setup_data --data-root data            # Swiss-Prot + Rhea
python -m webapp.setup_data --data-root data --pfam     # + Pfam (~600 MB)
python -m webapp.fetch_kegg --data-root data            # catálogo KEGG
```

Without data, the GUI warns you and shows the command. The modes that explore
the 400 genomes need that data to be generated; see
[`../data/README.md`](../data/README.md).

## Editable installation (development)

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest tests -q
```

The editable install keeps access to `scripts/ghost/` (ORF engine). An
isolated wheel would need to package that engine as well; for now the package
includes `webapp` and its interface data (locales, component, samples).

## Common problems

| symptom | fix |
|---|---|
| `streamlit: command not found` | use `.venv/bin/python -m streamlit run …` |
| Swiss-Prot/Rhea missing | `python -m webapp.setup_data --data-root data` |
| DIAMOND/tRNAscan show up as ❌ | install the system tools or use `--with-tools` |
| The browser does not open | go to `http://localhost:8501` manually |
| The menu does not show the app | log out and back in, or run `update-desktop-database ~/.local/share/applications` |
| GPU/ESMFold error | it is optional; annotation works on CPU |
