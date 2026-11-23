---
title: "Bacterial DNA Annotator: a modular, evidence-aware pipeline for genome annotation and metabolic gap exploration"
tags:
  - bioinformatics
  - bacterial genomics
  - genome annotation
  - metabolic gap filling
  - operon prediction
  - reproducible research
authors:
  - name: Camilo Rojas Bravo
    affiliation: 1
affiliations:
  - name: "Sin afiliación institucional actual; ex estudiante de la Universidad de Santiago de Chile (USACH), Santiago, Chile"
    index: 1
date: 23 November 2026
bibliography: paper.bib
---

# Summary

**Bacterial DNA Annotator** is an open-source pipeline and web application for
annotating a bacterial genome from scratch and exploring what the organism can
and cannot do metabolically. It combines gene prediction, non-coding RNA
detection, unannotated ORF discovery, metabolic network reconstruction, gap
analysis, operon and promoter prediction, and a clickable genome browser in a
single reproducible workflow. The application runs on CPU and is available in
Spanish, English, and Simplified Chinese.

Unlike a genome viewer tied to a fixed dataset, the project is designed around
the method: the 400 genomes used during development are the validation base,
while the target is a new genome. Every result declares its source, evidence,
confidence, and status. The tool never promotes a prediction to `confirmed`
without EC, contextual, or experimental support.

# Statement of need

Bacterial genome annotation is often split between command-line tools that
require bioinformatics expertise and web viewers that hide how results were
produced. This creates two problems: students cannot see the reasoning behind a
call, and researchers cannot easily audit or reuse the evidence. Metabolic gap
lists are especially prone to over-interpretation: an absent enzyme may be an
unannotated gene, a transport event, an alternative route, or a false positive.

Bacterial DNA Annotator addresses this gap with three design choices:

1. **Evidence-aware outputs.** Every gap, operon, and structural hypothesis
   carries its evidence, confidence, and status, with explicit limitations.
2. **A modular architecture.** Each domain (annotation, metabolism, regulation,
   structure, taxonomy, membrane, mobile elements, targets, validation, ghost
   genes) declares its options, requirements, estimated time, self-test, and
   execution. The GUI discovers modules automatically and the same modules can
   be tested from the command line.
3. **A low barrier to entry.** A bundled 150 kb demo, plain-language glossary,
   per-step time estimates, and an Ubuntu installer with a desktop launcher let
   a student start without a GPU or a cluster.

The software is intended for teaching, method development, and exploratory
research. It does not replace curated databases or experimental validation; it
makes the evidence and its limits visible.

# Architecture and reproducibility

The pipeline is organized as Python modules with deterministic outputs and
idempotent stages. Inputs can be a FASTA file, a genome identifier from the
project, or a precomputed table. Outputs are written as GFF3, TSV, JSON, and
Markdown, and can be downloaded from the GUI or the command line. The public
repository contains source code, tests, documentation, a demo genome, and the
module contract, but not third-party data with academic or subscription
licenses.

The GUI is built with Streamlit; the genome browser uses a small, build-free
component that renders SVG tracks and returns the clicked feature. The code has
unit and integration tests, a path audit for the script tree, and continuous
integration covering Python 3.10–3.12 and a Playwright smoke test of the GUI.

# Acknowledgements

The author thanks the developers of pyrodigal, DIAMOND, HMMER, Infernal,
tRNAscan-SE, barrnap, ViennaRNA, Streamlit, Rhea, RegulonDB, UniProt, Pfam,
Rfam, Foldseek, P2Rank, AlphaFill, M-CSA, and ESM-2, whose work makes this
project possible. The core methods rely on gene prediction [@hyatt2010prodigal],
protein search [@buchfink2015diamond], profile HMMs [@eddy2011hmmer],
reaction knowledge [@bansal2022rhea], curated regulation
[@salgado2024regulondb], RNA folding [@lorenz2011viennarna], and protein
language models [@lin2023esm]. Development assistance was provided by DeepSeek Harness (DSH)
with the DeepSeek model in November 2026 for code review, interface
translation, documentation, GUI testing, and packaging. Scientific decisions,
data, and validation remain the author's responsibility.

# References
