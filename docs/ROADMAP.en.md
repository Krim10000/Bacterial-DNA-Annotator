🌐 [Español](ROADMAP.md) · **English**

# Vision, real status and next steps

## The idea

From a DNA sequence, the user should be able to answer, with evidence and
explicit limits:

```text
ADN
 → especie y taxonomía
 → genes y proteínas
 → función probable (UniProt/EC/Pfam)
 → rutas metabólicas y brechas
 → operones y promotores
 → reguladores y bucles de regulación
 → hipótesis estructurales
```

It is not an “annotator that knows everything”. It is a system that **reduces
uncertainty layer by layer** and says what is measured, what is inference and
what cannot be concluded.

## What is implemented

| capability | status | evidence |
|---|---|---|
| Species/taxonomy | ✅ | ANI, GTDB/NCBI, confidence and conflicts |
| CDS, tRNA, rRNA | ✅ | pyrodigal, tRNAscan-SE, barrnap |
| Proteins and unannotated ORFs | ✅ | `proteins.faa`, ghost engine |
| Function by homology | ✅ | DIAMOND + Swiss-Prot, EC, UniProt |
| Metabolic pathways and gaps | ✅ | Rhea, transport, closures, evidence |
| Operons and promoters | ✅ | rules, TSS/TTS, PWMs, RegulonDB as validation |
| Clickable browser | ✅ | card per gene/operon/gap |
| Modules and ES/EN/ZH GUI | ✅ | 10 modules, self-test, summary |
| Docker and DOI | ✅ | Compose, Zenodo `10.5281/zenodo.24469732` |

## What is partial

### Repressors and regulatory loops

Implemented:

- genome → repressor candidate pipeline;
- family by domain, homolog by DIAMOND and inherited effector;
- PWM sites and relative score;
- 92 % recall in *E. coli* and 98.7 % of sites versus random.

Not implemented:

- effector direction (inducer/corepressor);
- magnitude of the effect (ΔΔG, absolute affinity);
- validation that the loop closes in a real context.

Reason: it requires a GPU (MM/GBSA, ABFE, co-folding) and, above all,
experimental data (EMSA/ITC). The reopening criterion is in
[`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md): first MM/GBSA of LacI,
PurR and TrpR; if it does not rank them, no more is spent.

### Reduction of unknown proteins across the 400 genomes

Implemented:

- deduplication, clustering and manifests;
- cheap filters before structure;
- F1 local Foldseek;
- F2 EnzyMM/M-CSA;
- F3 P2Rank and F4 AlphaFill;
- propagation to clusters and prioritization;
- **gating by pLDDT/pTM and early-exit** in `structural_filters.py`.

Pending:

- splitting the cheap and expensive queues into independent commands;
- completing the 2,570 structures ≤600 aa and the 169 >600 aa;
- assigning broad categories (`enzyme_hypothesis`, `transporter_hypothesis`,
  `regulator_hypothesis`, `unknown`) automatically;
- manual review and final propagation.

The rule is: **no filter closes a function**. It only reduces candidates.

### GPU: Kaggle as a first filter

Kaggle is used as a **first GPU filter**, not as a dependency of the app:

- it allows testing ESMFold/co-folding and large models without a local GPU;
- sequences are uploaded manually and only if the user accepts;
- results are imported as tables/structures;
- it does not run at startup or in the GUI;
- the local pipeline keeps working on CPU.

This is deliberate: Kaggle is not reproducible as a permanent service and must
not be a requirement for annotating a genome.

## Upcoming

| priority | piece | cost | depends on |
|---|---|---|---|
| 1 | `protein_query`: protein → homologs/function | 1–2 days | Swiss-Prot + DIAMOND |
| 2 | Full GPU-free reduction of unknowns | 3–5 days | current tables |
| 3 | Automatic broad categories | 2–3 days | filters + context |
| 4 | GPU reopening of effector direction | budget | Kaggle/Modal + bench |
| 5 | macOS installer / WSL2 guide | 1–2 days | test on a clean machine |
| 6 | PyPI / TestPyPI | 1 day | packaging |
| 7 | JOSS / bio.tools | 1–3 days | DOI already available |
| 8 | Reproducible benchmark vs. Bakta/Prokka/DFAST | 2–4 days | 5–6 genomes |
| 9 | Export SBML draft (CarveMe/COBRApy) | 1 week | Rhea mapping |

## Benchmark

The comparison with 12+ projects (Prokka, Bakta, DRAM, CarveMe, gapseq,
Operon-mapper, BV-BRC, etc.) is in [`BENCHMARK.md`](BENCHMARK.md).
Bottom line: base annotation is table stakes; our advantage is local
integration with evidence and the ES/EN/ZH GUI; the main gap is curated
functional databases, packaging and a public benchmark.

## What we do not promise

- A gap does not prove that an enzyme is missing.
- An HTH domain does not confirm a regulon.
- ESMFold does not confirm function.
- AlphaFill does not confirm substrate.
- P2Rank does not confirm active site.
- No result is `confirmed` without an experiment, EC or strong context.
- The custom DNA model is experimental and has not been evaluated.

## Related documents

- [`UNKNOWN_REDUCTION_PLAN.md`](UNKNOWN_REDUCTION_PLAN.md)
- [`FILTER_STRATEGY.md`](FILTER_STRATEGY.md)
- [`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md)
- [`PROTEINAS_REGULADORAS.md`](PROTEINAS_REGULADORAS.md)
- [`PIPELINE_GENOMA.md`](PIPELINE_GENOMA.md)
- [`GPU_STATUS.md`](GPU_STATUS.md)
- [`BENCHMARK.md`](BENCHMARK.md)
- [`INTEGRATIONS.md`](INTEGRATIONS.md)
- [`ORIC.md`](ORIC.md)
- [`TRANSITION.md`](TRANSITION.md)
- [`PIPELINE_PARITY.md`](PIPELINE_PARITY.md)
- [`COMPLEX_GAPS.md`](COMPLEX_GAPS.md)
- [`COMPLEX_CPU.md`](COMPLEX_CPU.md)
- [`COMPLEX_TRANSFER.md`](COMPLEX_TRANSFER.md)
