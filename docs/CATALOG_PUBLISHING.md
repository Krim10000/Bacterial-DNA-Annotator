# Publicar un catálogo de hipótesis priorizadas

## Principio

Publicar:

- una tabla maestra por cluster;
- evidencia por proteína;
- modelos PDB;
- código y versiones;
- limitaciones explícitas.

No publicar como “funciones confirmadas”.

## Dónde

| componente | destino |
|---|---|
| Tablas TSV/CSV | **Zenodo** (DOI) |
| Código | GitHub + release |
| PDB/mmCIF | **ModelArchive** o Zenodo |
| Preprint | bioRxiv |
| Familias validadas | UniProt/InterPro/Pfam, tras revisión |
| Versión congelada | Zenodo + GitHub release |

## Formato principal

`hypotheses.tsv`:

```text
catalog_version
cluster_id
uniprot
genome
protein_id
sequence_length
sequence_sha256
plddt
ptm
pfam_domains
tm_helices
signal_peptide
f1_foldseek
f2_enzymm
f2_template
f2_ec
f3_pocket
f4_alphafill
f5_active_site
category
confidence
review_status
evidence_sources
limitations
pdb_path
model_license
```

Una fila por **representante de cluster**.

## Tablas auxiliares

```text
clusters.tsv        cluster -> proteínas/genomas
propagation.tsv     cluster -> proteínas
evidence.jsonl      evidencia completa por proteína
review.tsv          decisiones humanas
data_dictionary.md  vocabularios controlados
```

## Vocabularios

`category`:

```text
enzyme_hypothesis
transporter_hypothesis
regulator_hypothesis
membrane_protein
disordered
unknown
```

`confidence`:

```text
orientative
supported_hypothesis
high_priority
conflict
uncharacterized
```

`review_status`:

```text
pending
reviewed
rejected
needs_validation
```

## Evidencia JSON

```json
{
  "cluster_id": "...",
  "uniprot": "...",
  "plddt": 0.74,
  "ptm": 0.68,
  "f1": {"top": "AF-...", "bits": 1200},
  "f2": {"template": "1abc", "ec": "1.1.1.1"},
  "f3": {"pockets": 2},
  "f4": {"ligands": ["NAD"]},
  "f5": {"active_positions": "170;171"},
  "sources": ["Foldseek", "EnzyMM", "P2Rank"],
  "limitations": "hipótesis computacional"
}
```

## Flujo de publicación

1. Congelar versión.
2. Calcular hashes.
3. Generar tabla maestra.
4. Escribir `data_dictionary.md`.
5. Empaquetar:
   - README;
   - LICENSE;
   - tablas;
   - evidencia;
   - modelos;
   - código/commit.
6. Subir a Zenodo → DOI.
7. Crear release en GitHub.
8. Depositar PDB en ModelArchive.
9. Publicar preprint.
10. Enviar familias validadas a UniProt/InterPro/Pfam.

## Qué no hacer

- No subir 2.570 PDB a GitHub.
- No llamar `confirmed`.
- No publicar sin licencia.
- No omitir limitaciones.
- No enviar familias a Pfam sin validación.

## Primer catálogo

Publicar el piloto de **1.192** ya:

- 176 hipótesis triples;
- 3.374 clusters propagados;
- 508 candidatos;
- PDB/evidencia.

Después actualizar con 2.570 y 169 largas.

## Licencias

- Tablas: CC BY 4.0.
- Código: MIT.
- PDB/modelos: revisar términos de ESMFold/AlphaFold/PDB.
- AlphaFill/M-CSA: citar y respetar términos.
