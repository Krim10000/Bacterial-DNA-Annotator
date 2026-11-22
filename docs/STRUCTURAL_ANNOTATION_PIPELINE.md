# Pipeline de anotación estructural automática

## Objetivo

Reducir candidatos antes de docking/validación, añadiendo evidencia estructural con procedencia.

## Etapas

```text
candidato
  → PDB (AFDB o ESMFold)
  → filtro pLDDT
  → Foldseek PDB100/AFDB
  → EnzyMM/PyJess (plantillas M-CSA)
  → bolsillo (fpocket/P2Rank) [opcional]
  → mapeo de sitio activo local
  → fusión de evidencia
  → revisión
  → docking solo finalistas
```

## 1. Entrada

- PDB del ORF exacto (ESMFold) o AFDB del homólogo.
- Metadatos: genoma, gap, UniProt, pLDDT.

## 2. Filtro de calidad

- pLDDT ≥70 para considerar.
- pLDDT ≥80 para evidencia fuerte.
- Registra el valor; no lo oculta.

## 3. Homología estructural

- Foldseek contra PDB100/AFDB.
- Guarda top hits, identidad, bits, cobertura.
- Busca PDB holo con ligando/cofactor.

## 4. EnzyMM/PyJess

- Escanea el PDB contra plantillas M-CSA.
- Salida: template, M-CSA ID, UniProt, EC, RMSD, residuos.
- `no_template` es un resultado, no un fallo.
- Fuente: `EnzyMM/PyJess + M-CSA`.

## 5. Bolsillo (opcional)

- fpocket o P2Rank.
- Solo apoyo: debe co-localizar con residuos catalíticos.
- No asigna función.

## 6. Sitio activo local

- UniProt features / mapa estructural existente.
- Conservación y distancia.
- No confundir binding site con active site.

## 7. Fusión de evidencia

| señal | peso |
|---|---|
| EnzyMM hit con ≥2 residuos | alto |
| residuo catalítico local conservado | alto |
| Foldseek holo con ligando | medio |
| bolsillo co-localizado | apoyo |
| cofactor compatible | apoyo |
| contexto genómico | alto |

## Estados

- `no_template`: sin plantilla; no descarta.
- `template_match`: hit geométrico.
- `template_plus_active_site`: hit + residuo local.
- `multi_source_support`: ≥2 fuentes independientes.
- `structural_supported`: requiere ≥2 residuos catalíticos + bolsillo/cofactor + revisión.
- `confirmed`: solo experimento.

## GUI

En el cierre ESMFold:

1. PDB generado.
2. Botón **Escanear con EnzyMM/PyJess**.
3. Tabla: template, M-CSA, EC, RMSD, residuos.
4. Aviso: no confirma función.
5. Botón separado **Cruzar Foldseek + AlphaFill**.
6. Resultados exportables.

## Pipeline

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --esmfold "neuros nap" \
  --esmfold-gaps 3 \
  --structure-evidence \
  --structural-annotation
```

Salidas:

- `GENOMA.enzymm_scan.tsv`
- `GENOMA.structural_evidence.tsv`
- `GENOMA.esmfold_closure.tsv`
- `GENOMA.report.md`
- `GENOMA.summary.json`

## Scripts

```bash
python scripts/structure/11_enzymm_scan.py \
  --pdb candidato.pdb --out resultados.tsv --jobs 4
```

## Coste/tiempo

- EnzyMM: $0, local, CPU; ~30–40 s/PDB; plantillas 16 MB.
- Foldseek: $0, local; segundos/estructura.
- fpocket/P2Rank: $0, local; segundos.
- Docking: $0, CPU; 3–8 min.
- No requiere GPU.

## Límites

- EnzyMM no cubre todas las enzimas.
- Un hit no prueba reacción.
- AFDB es modelo.
- El bolsillo no es sitio activo.
- La decisión final es manual.

## Validación

- `1LYZ`: 5 matches M-CSA 203 (EC 3.2.1.17), RMSD 0.45–0.65 Å.
- `2055`: `no_template`.
- `P77161` AFDB: `no_template`.
- `no_template` no se interpreta como ausencia de función.
