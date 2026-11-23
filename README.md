# 🧬 Bacterial DNA annotator

🌐 [Español](README.es.md) · **English** · [中文](README.zh.md)

[![CI](https://github.com/Krim10000/Bacterial-DNA-Annotator/actions/workflows/ci.yml/badge.svg)](https://github.com/Krim10000/Bacterial-DNA-Annotator/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/Krim10000/Bacterial-DNA-Annotator)](https://github.com/Krim10000/Bacterial-DNA-Annotator/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24469732.svg)](https://doi.org/10.5281/zenodo.24469732)
[![Version](https://img.shields.io/badge/version-1.1.3-blue)](CHANGELOG.md)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20|%203.11%20|%203.12-blue)](pyproject.toml)
[![Languages](https://img.shields.io/badge/languages-ES%20|%20EN%20|%20ZH-purple)](#-languages)

Annotates bacterial genomes **from scratch**, finds genes nobody had annotated
before, rebuilds metabolic pathways, looks for the gaps that keep an organism
from living, and proposes operons and promoters. Everything runs on your own
computer and is driven from a web interface.

Built so a high-school student can use it without fear and a researcher with
their own data can use it too: every result states **where it comes from**, with
what **confidence**, and what it does **not** prove.

> **You do not need a GPU.** Annotation, gaps, operons and the browser all run on
> CPU. The GPU is only used for local ESMFold and for training the DNA model,
> which are optional. See [`docs/GPU_STATUS.en.md`](docs/GPU_STATUS.en.md).

## 🧭 Scope at a glance

The goal is to go from a DNA sequence to a species, its proteins, metabolic
routes, operons and regulatory loops — always with evidence and limits.

| layer | status |
|---|---|
| Species / taxonomy | ✅ implemented |
| Genes, proteins, unannotated ORFs | ✅ implemented |
| Functions / EC / UniProt | ✅ implemented for homologs |
| Metabolic routes and gaps | ✅ implemented |
| Operons and promoters | ✅ implemented |
| Regulatory loops | 🟡 partial: repressor candidates and operators; direction/magnitude of effectors is not |
| Unknown-protein reduction (400 genomes) | 🟡 partial: cheap filters and structural gating; full structure set pending |
| DNA language model | ⛔ experimental, not evaluated |

Full status, GPU limits and Kaggle policy:
[`docs/ROADMAP.en.md`](docs/ROADMAP.en.md).

![Annotation results](webapp/screenshot_v2_results.png)

## ⚡ Installing on Ubuntu (3 steps)

```bash
git clone https://github.com/Krim10000/Bacterial-DNA-Annotator.git DNA_TOKEN_BAC
cd DNA_TOKEN_BAC
bash install.sh --with-data
```

`install.sh` creates a virtual environment, installs the dependencies, downloads
Swiss-Prot and Rhea (~120 MB) and adds **«Anotador de ADN» to the applications
menu**. After that:

- open **Anotador de ADN** from Applications, or
- run `./lanzar.sh` (double-clicking works too).

Installer options:

| command | what it does |
|---|---|
| `bash install.sh` | environment + launcher (without downloading data) |
| `bash install.sh --with-data` | + Swiss-Prot and Rhea |
| `bash install.sh --with-tools` | + tRNAscan-SE, barrnap, HMMER, DIAMOND (asks for sudo) |
| `bash install.sh --full` | + every dependency of the pipeline |
| `bash install.sh --dry-run` | shows what it would do, without touching anything |
| `bash install.sh --uninstall` | removes the launcher from the menu |

Windows or macOS? Use Docker: [`docs/INSTALL_DOCKER.md`](docs/INSTALL_DOCKER.md).

Updates and how to subscribe to releases: [`docs/UPDATES.md`](docs/UPDATES.md).

On GNOME, the first time you run a `.desktop` file you may have to mark it as
trusted (right-click → *Allow Launching*). The installer already leaves it in
`~/.local/share/applications`.

### Manual alternative

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-annotator.txt
.venv/bin/python -m streamlit run webapp/app.py
```

Open <http://localhost:8501>.

### Included demo (without downloading data)

In **🧬 Genome browser** choose the **`demo`** genome, or press
**🧬 Explore the browser demo** on the start screen. It is 150 kb of
*Chlorobaculum tepidum* annotated from scratch: genes predicted with pyrodigal
and operons grouped by proximity and strand. It is there to try out the browser
and the cards you get when you click; it is **not** experimental evidence and is
not used in the project metrics. The other modes need the project data.

## ⏱️ Requirements and runtimes

| | minimum | recommended |
|---|---|---|
| System | Ubuntu 22.04/24.04, macOS or WSL | Ubuntu 24.04 |
| Python | 3.10 | 3.11 or 3.12 |
| RAM | 8 GB | 16 GB |
| Disk | 10 GB | 25 GB with all databases |
| CPU | 4 cores | 8 cores |
| GPU | **not needed** | optional, for ESMFold |
| Network | only for the first download | — |

Approximate times on an 8-core laptop without a GPU:

| step | time |
|---|---|
| Annotating a ~5 Mb genome | 1–3 min |
| Searching for unannotated ORFs | ~5 s per genome |
| Metabolic network and gaps | ~2 min per genome |
| Rule-based operons | ~8 s per genome |
| Termination hairpins | 1–10 min per genome |
| Comparing two genomes | 1–5 min |
| ESMFold (Neurosnap API) | 1–2 min per protein (~0.04 credits) |
| Training the DNA model | hours on GPU (not needed for annotation) |

The app shows this same table under **⏱️ What do I need and how long does it
take?**, along with a glossary for anyone who has never heard of DIAMOND, HMMER
or an operon.

## 🌐 Languages

The interface is available in **Spanish, English and Simplified Chinese**. The
selector sits at the top of the sidebar and your choice is remembered on the
machine. You can also link directly:

```
http://localhost:8501/?lang=en
http://localhost:8501/?lang=zh
```

**Data is not translated**: gene products, pathways and pipeline labels are
shown exactly as the source produces them. To add a language, copy
`webapp/locales/en.json`, translate the values and add the code in
`webapp/i18n.py`; `tests/test_i18n.py` checks that no key is missing.

## 🖥️ The interface

Six modes in the sidebar:

| mode | what it is for |
|---|---|
| 🧬 **Annotate FASTA** | you upload one or more genomes and download GFF3, proteins, candidates and deserts |
| 🧪 **Metabolic pathways** | you explore the already computed analysis of 400 genomes (gaps, closures, operons, mobility, membrane, taxonomy…) |
| 🧬 **Genome browser** | you see a window of the genome with tracks and **click on a gene** to read its card |
| 📥 **Reports and downloads** | gathers the available TSV/JSON/Markdown files, with size, source and license, and builds a ZIP of results |
| 🧩 **Modules** | each module declares its options, status, time and self-test; the interface builds itself |
| 🙏 **Credits and tools** | who made each piece, with link, citation and license |

### DNA is clickable

In the browser, every gene, operon, promoter, hairpin and metabolic suspicion has
a card: clicking it shows the product, coordinates, strand, EC, confidence,
−10/−35 boxes, ΔG, neighbouring genes and the interpretation of that clue. The
detail does not replace the table: it complements it.

You can also search by coordinate (`1000000-1050000`), gene (`thrA`), product
(`protease`), EC or reaction (`RHEA:12345`), and share a link to a region with
`?genoma=...&contig=...&pos=...`.

### Downloads for research

Every table has its TSV/CSV download button; reports download as Markdown and
structured results as JSON. The GFF3 files can be opened in Artemis, IGV or NCBI.
Under **Gaps** you will find the downloads for review, evidence and mapping
candidates; under **Completeness**, the essential-pathway report.

Files carry **evidence, confidence and status** columns. There are no
"confirmed" results without an EC, context or experiment.

## 🔬 What it computes

```
FASTA
  ├─ CDS (pyrodigal), tRNA (tRNAscan-SE), rRNA (barrnap)
  ├─ Unannotated ORFs (ghost engine) + optional ESM-2
  ├─ Annotation deserts (≥3 kb without signals)
  ├─ Metabolic network (DIAMOND + Swiss-Prot + Rhea)
  │    └─ gaps, closure by domain/homology, transport
  ├─ Operons: deterministic rules, promoters, hairpins
  ├─ Mobile/defence, membrane, taxonomy, chains, allostery
  └─ Report + GFF3 + downloadable TSV/JSON
```

Every stage is idempotent and can be run on its own:

```bash
bash run_pipeline.sh doctor        # check the environment
bash run_pipeline.sh components    # reference truth and components
bash run_pipeline.sh rules         # rule-based operons (8 s)
bash run_pipeline.sh hairpins      # termination hairpins
bash run_pipeline.sh closure       # gap closure
```

With no arguments it does **not** run everything: you have to ask for a stage.

## 🧩 Modules

The **Modules** view discovers the folders in `modules/` and builds the selector,
options, status, estimated time and self-test button by itself. Adding a module
means dropping `modules/<id>/module.py` with its `MODULE`; no `app.py` edit is
needed. Every module can be tested without the GUI:

```bash
python -m modules --list
python -m modules --describe annotation
python -m modules --selftest annotation
```

Full contract in [`docs/MODULES.en.md`](docs/MODULES.en.md).

## 📦 Data: what the repository includes and what it does not

This repository **does not redistribute** RegulonDB, KEGG, DoriC, STRING, DBTBS
or BATTER: they have academic or subscription licenses. They are downloaded on
your machine with the scripts shown and stay out of git. The databases that can
be redistributed (Rfam, Pfam, UniProt, PDB, AlphaFold DB, Rhea, ChEBI) are not
versioned either, because of their size: they are downloaded once.

- Full terms: [`DATA_NOTICE.md`](DATA_NOTICE.md).
- How to regenerate each database: [`docs/INSTALL.en.md`](docs/INSTALL.en.md) and
  [`docs/EXTERNAL_SOURCES.md`](docs/EXTERNAL_SOURCES.md).
- Demo data: `webapp/samples/demo.fna` and `demo1b.fna`.

The app works without project data **to annotate a FASTA**. The modes that
explore the 400 genomes need that data to be generated or downloaded; if it is
missing, the app says so and shows you the command.

## 🖥️ GPU

| stage | GPU | status in v1.0 |
|---|---|---|
| Annotation, gaps, operons, browser | no | ✅ included |
| Structural filters (Foldseek, P2Rank, AlphaFill) | no | ✅ included |
| ESMFold through Neurosnap | no (remote) | ☁️ optional, consumes credit |
| Local ESMFold | ≥16 GB VRAM | ⛔ not included |
| DNA model training/evaluation | A10G | ⚠️ experimental, not evaluated |
| MM/GBSA, ABFE, molecular dynamics | GPU | ⛔ outside v1.0 |

Details, evidence and costs in [`docs/GPU_STATUS.en.md`](docs/GPU_STATUS.en.md).

## 🙏 Credits

**Author:** Camilo Rojas Bravo — Biotechnology Engineer, graduated from the University of Santiago, Chile (USACH); currently without an institutional affiliation — Camilo.Rojas@usach.cl.

**Development assistance:** DeepSeek Harness (DSH) with the DeepSeek model,
**November 2026**, on code review, interface translation, documentation, GUI
testing and packaging. Scientific decisions, data and validation are the
author's responsibility.

This project builds on third-party programs and databases: pyrodigal, DIAMOND,
HMMER, Infernal, tRNAscan-SE, barrnap, ViennaRNA, Streamlit, Rhea, RegulonDB,
UniProt, Pfam, Rfam, Foldseek, P2Rank, AlphaFill, M-CSA, ESM-2 and many more.
Each one has its authorship, its citation and its license on the **Credits**
page and in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

## 📄 Citation and license

- Cite with [`CITATION.cff`](CITATION.cff) or with
  [`AUTHORS.md`](AUTHORS.md#cómo-citar). DOI for this version:
  [10.5281/zenodo.24469732](https://doi.org/10.5281/zenodo.24469732).
- Original code: **MIT** ([`LICENSE`](LICENSE)).
- Dependencies: they keep their own license. **pyrodigal, orcapy, DIAMOND,
  tRNAscan-SE and barrnap are GPL-3.0**; any distribution that bundles them
  (Docker image, `.deb`, service) must comply with the GPL-3.0.
- Data: [`DATA_NOTICE.md`](DATA_NOTICE.md).

## 📁 Structure

```
modules/           self-discovered modules (options, status, test, run)
webapp/            Streamlit interface, i18n, browser, credits
  locales/         es.json, en.json, zh.json
  components/      clickable genome component (no build)
  samples/         two demo FASTA files
scripts/           340 scripts by stage (see scripts/README.md)
docs/              installation, pipeline, GPU, decisions, sources
tools/             development utilities (extract texts, index)
tests/             i18n, import and GUI tests
install.sh         Ubuntu installer + launcher
lanzar.sh          double-click startup
run_pipeline.sh    stage orchestrator
```

## ⚠️ Known limits

- Operon prediction from sequence tops out at **F1 ≈ 0.78–0.81**; two independent
  methods converge there.
- **Hairpins are a clue, not a border**: 70 of the 375 curated K-12 terminators
  lie inside their operon and ~10 % are Rho-dependent.
- The **3′** side of operons is poorly covered: TTS cover 16 % of the boundaries
  compared with 61 % for TSS.
- Curated components exist only for *E. coli* K-12; for everything else they are
  predicted from the sequence.
- The **DNA model is not evaluated** in this version (no GPU). It is not used for
  annotation.
- A gap **does not prove** that an enzyme is missing: it may be an unannotated
  gene, transport, an alternative route or a false positive. The app says so on
  every card.

## 🇪🇸 Resumen en español

Anotador de genomas bacterianos con GUI en Streamlit (español, inglés y chino
simplificado), navegador de genoma clicable, análisis de brechas metabólicas y
descargas en GFF3/TSV/JSON. Funciona en CPU; la GPU no es necesaria para anotar.
No redistribuye datos de RegulonDB, KEGG, DoriC, STRING, DBTBS ni BATTER: se
descargan localmente bajo sus condiciones. El código original es MIT; las
herramientas de terceros conservan su licencia (varias son GPL-3.0). Autor:
Camilo Rojas Bravo (Ingeniero en Biotecnología, titulado de la USACH;
sin afiliación institucional actual; Camilo.Rojas@usach.cl).
Asistencia de desarrollo:
DeepSeek Harness con DeepSeek, noviembre de 2026.

## 🔭 Roadmap

- **Protein query**: protein FASTA → homologs, domains and function.
- **Unknown-protein reduction**: finish the GPU-free filters and broad
  categories for the 400 genomes.
- **GPU/Kaggle**: first-pass ESMFold/co-folding on Kaggle; optional, never
  required to annotate a genome.
- **Packaging**: macOS installer, WSL2 guide and PyPI.

Details and honest limits: [`docs/ROADMAP.en.md`](docs/ROADMAP.en.md).

## ⭐ Did it help?

If this project saved you time, **give it a star on GitHub**:
[github.com/Krim10000/Bacterial-DNA-Annotator](https://github.com/Krim10000/Bacterial-DNA-Annotator).
It costs nothing and helps other people find it.

---

Made in Chile, with real data and real biology questions. If you find a bug, open
an *issue* with the genome, the stage and the exact command.
