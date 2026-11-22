# Integración de hipótesis estructurales

## Principio

La estructura **no sustituye** la anotación. Se añade una capa separada:

```text
anotación base
+ hipótesis estructural revisable
+ decisión de revisión
+ evidencia descargable
```

## Pipeline

### Manifest de estructuras

`manifest.tsv`:

```text
genome  protein_id  uniprot  pdb
GCF_...  NZ_...      P77161   /ruta/2055.pdb
```

### Comando

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --structure-manifest manifest.tsv \
  --structural-filters
```

Salidas:

- `GENOMA.structural_hypotheses.tsv`
- `GENOMA.structural_hypotheses.json`
- `GENOMA.structural_filters.tsv`
- `GENOMA.report.md`
- `GENOMA.summary.json`

### Campos

`genome, protein_id, uniprot, pdb, hypothesis, ec, cofactor, pocket,
template_mcsa, foldseek_top, f1_foldseek, f2_enzymm, f3_pocket,
f4_alphafill, f5_active_site, n_sources, confidence, review_decision,
limitation, evidence_json`

## Confianza

| combinación | etiqueta |
|---|---|
| F2 + F3 + F4 | `supported_hypothesis` |
| F2 + F5 | `high_priority` |
| F2 solo | `orientative` |
| F1 + (F3 o F4) | `structural_only` |
| sin fuentes | `uncharacterized` |

## GUI

Nueva pestaña:

**🔬 Hipótesis estructural**

Muestra:

- tabla de proteínas;
- hipótesis;
- EC;
- cofactor;
- bolsillo;
- plantilla M-CSA;
- homólogo Foldseek;
- número de fuentes;
- confianza;
- decisión.

Al seleccionar una proteína:

- F1 Foldseek;
- F2 EnzyMM/M-CSA;
- F3 P2Rank;
- F4 AlphaFill;
- F5 sitio activo;
- fuente y límite;
- evidencia JSON;
- PDB descargable;
- decisión de revisión.

## Decisiones

```text
pending
reviewed_orientative
reviewed_supported
conflict
rejected
needs_validation
```

Se guardan en la sesión y se exportan como:

```text
GENOMA_review_decisions.tsv
```

No se modifica la anotación base.

## Regla

- `supported_hypothesis` no es `structural_supported`.
- `high_priority` no es cierre.
- `confirmed` solo con experimento.
- `no_template`/`no_pocket` no eliminan.

## Coste

- Foldseek, EnzyMM, P2Rank, AlphaFill, UniProt: $0.
- CPU local.
- Sin GPU.
- Sin créditos Neurosnap.
