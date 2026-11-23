# Transferencia de complejos K-12 y evidencia CPU

## Qué se implementó

`scripts/complexes/transfer_orthology.py` proyecta los **273 complejos
curados de K-12** a los 400 genomas mediante *reciprocal best hit* (RBH).

```text
gold_k12.tsv
  + k12_vs_all.tsv      (K-12 → genoma)
  + hits.tsv            (genoma → K-12)
  → RBH por genoma
  → frac_mapped, identidad, masa, perfil filogenético
```

Filtros:

- identidad ≥30 %;
- evalue ≤1e-5;
- RBH por bitscore;
- ≥2 subunidades mapeadas para emitir el complejo.

Salida por genoma:

```text
data/complexes/transfer_<genoma>.tsv
```

Resumen:

```text
reports/complex_transfer_summary.tsv/json
```

## Resultado en los 400

| métrica | valor |
|---|---:|
| genomas con transferencia | 400 |
| complejos transferidos (≥2 subunidades) | **19.085** |
| con `frac_mapped ≥0.8` | **8.408** |
| genomas con al menos un complejo | 400 |

Ejemplo, `GCF_000006985.1_ASM698v1`:

- 59 complejos transferidos;
- 50S: 13/34 subunidades, `frac=0.382`;
- identidad media 54,1 %;
- masa observada 235,9 kDa;
- `phylo_jacc=0.193`.

## Columnas nuevas

`transfer_<genoma>.tsv`:

```text
complex_id  complex_name  n_subunits  n_mapped  frac_mapped
identity_mean  score  genes_k12  genes_genoma
observed_n  expected_n  stoich_status  mass_kda  phylo_jacc
evidence_class  eco  limitation
```

- `stoich_status`: `complete` si `observed_n ≥ expected_n`; `incomplete` si hay
  `expected_n`; `unknown` si no.
- `mass_kda`: suma de `n_aa × 110 Da` de las subunidades observadas.
- `phylo_jacc`: Jaccard medio de presencia de las subunidades en los 400.

## Consistencia de masa

El `gold_k12.tsv` local no trae masa documentada: `properties` contiene
`HeteroN-mer`, no `Molecular weight = ... kDa`. Por eso:

- se calcula **masa observada**;
- se compara **número de subunidades** (`expected_n` vs `observed_n`);
- **no** se calcula error de masa contra una masa documentada que no existe.

Si se importa ComplexPortal con `Molecular weight`, el campo `mass_error` se
puede añadir sin cambiar el esquema.

## Perfil filogenético

`phylo_jacc` mide si las subunidades aparecen juntas en los 400 genomas:

- Jaccard de presencia por par;
- promedio por complejo;
- no se usa como prueba de interacción; es una señal de asociación.

## Límites

- RBH sin cobertura de alineamiento: los TSV actuales no traen `length/qlen`.
- K-12 es un **prior**, no una calibración para otro genoma.
- Los complejos grandes se fragmentan.
- `frac_mapped` bajo no significa que el complejo no exista: puede faltar
  anotación/ortólogo.
- La transferencia es `documentado_inferido` (`ECO:0005544`), no experimental.
- No cierra brechas.

## Transferencia contra brechas

`tools/audit_transfer_400.py` cruza las EC/Rhea de los complejos transferidos
con los metabolitos de las brechas de los 400:

| métrica | valor |
|---|---:|
| brechas totales | 8.655 |
| brechas con transferencia | **22** |
| brechas con transferencia fuerte (`frac≥0.8`) | **14** |
| brechas sin asociación | 8.633 |

Ejemplos fuertes:

- `CPX-1092` succinyl-CoA synthetase (varios genomas);
- `CPX-3446` TrpAB tryptophan synthase;
- `CPX-4384` putrescine ABC transporter.

Estas 22 se anotan como `complex_suspected`/`amber`, no como cerradas. La
mayoría de las brechas sigue sin una EC/Rhea que las conecte a un complejo.

Detalle: `reports/transfer_gap_details.tsv`.

## Sospecha fuerte en pipeline/GUI

La transferencia con `frac_mapped ≥0.8` se integra como
`complex_status="complex_suspected"` y `complex_flag="yellow"` cuando el genoma
está en los 400. No se queda solo en el navegador: `webapp/complex_transfer.py`
la carga en la tabla de complejos y en `protein_paths.tsv`.

No se convierte en `confirmed` ni cierra la brecha. La nota dice
"transferencia por ortología, no interacción demostrada".

## Visualización

El navegador carga `data/complexes/transfer_<genoma>.tsv` si no hay una tabla
`complexes` propia:

- azul: `complex_inferred_complete`;
- ámbar: `complex_inferred_incomplete`.

La GUI muestra los candidatos del modelo K-12 y los de firma; la transferencia
es una capa de 400 que se puede añadir al navegador.

## Pendiente

- Foldseek-Multimer: interfaz lista (`scripts/complexes/foldseek_multimer.py`),
  pero falta una base de ensamblajes PDB; sin DB devuelve `not_configured`.
- Coevolución de interfaz: requiere MSA por par y cómputo intensivo; no se
  ejecutó.
- Importar masas documentadas de ComplexPortal. El gold local no trae
  `Molecular weight`, así que hoy se reporta masa observada y estequiometría.
