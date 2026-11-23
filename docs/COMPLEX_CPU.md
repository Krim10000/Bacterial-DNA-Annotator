# Complejos con CPU: evaluación restringida por vecindad génica

**Estado 2026-11-23:** fase 1 implementada en `webapp/complexes.py` (18 firmas,
`complexes.tsv`, columnas en proteínas y brechas, tabla GUI), auditoría de los
400 en `tools/audit_complexes_400.py` y pista `complexes` en el navegador.
No hay todavía template PDB/Foldseek para interfaz cuaternaria.

## Principio

La estructura cuaternaria no se puede predecir de forma fiable con CPU desde
una lista global de proteínas. Pero sí se puede **evaluar qué complejos son
plausibles** si se cumplen tres condiciones:

1. las subunidades están codificadas en genes **físicamente cercanos**;
2. la arquitectura de dominios/familia es compatible con un complejo conocido;
3. el conjunto de subunidades está completo o se declara incompleto.

La proximidad génica es **necesaria pero no suficiente**. Un operón no es un
complejo, y un complejo puede tener genes dispersos. Por eso el resultado se
etiqueta como `complex_candidate`, no como `confirmed`.

## Flujo CPU posterior a la primera pasada

```text
proteínas + función + coordenadas
  → clústeres de vecindad (contig/hebra/distancia/operón)
  → firma de complejo conocida o arquitectura de dominios
  → subunidades esperadas vs presentes
  → estado de completitud
  → plantilla estructural CPU (PDB/AlphaFold DB/Foldseek)
  → tarjeta y pista en el navegador
```

### 1. Clústeres de vecindad

Para un genoma nuevo:

- mismo contig;
- misma hebra;
- distancia entre genes ≤ umbral (por defecto 100–300 pb para operón;
  hasta 10 kb solo para vecindad laxa);
- sin tRNA/rRNA intercalado;
- sin terminador conocido entre ellos, si existe;
- si hay tabla de operones, usarla como agrupación primaria.

No se permite emparejar genes de contigs distintos ni buscar homólogos globales
para “armar” un complejo.

### 2. Firmas de complejo

Diccionario curado, no inventado. Ejemplos:

| complejo | subunidades típicas | evidencia |
|---|---|---|
| ATP sintasa | `atpA..atpI` | genes contiguos, PFAM ATPase |
| NADH deshidrogenasa | `nuoA..nuoN` | operón `nuo` |
| succinato deshidrogenasa | `sdhA..sdhD` | operón `sdh` |
| citocromo oxidasa | `cyoA..cyoD`, `ccoN...` | operón |
| piruvato deshidrogenasa | `aceE`, `aceF`, `lpdA` | operón |
| nitrogenasa | `nifH`, `nifD`, `nifK` | operón |
| sistema de secreción | `secA/Y/E/G`, `tatA/B/C` | operón |
| ribosoma | `rps*`, `rpl*` | vecindad y familias |
| transportador ABC | `ABC_tran` + dominios TM | arquitectura |
| pili/fimbrias | `fim`, `pil` | operón |

La firma debe declarar:

- `complex_id`;
- `subunits_expected`;
- `essential`;
- `source` (Complex Portal, EcoCyc/MetaCyc, KEGG module, PDB, literatura);
- `license`.

Si no hay firma conocida, no se inventa un nombre: se emite
`complex_candidate` con `architecture_only`.

### 3. Completitud

| estado | significado |
|---|---|
| `complex_complete` | todas las subunidades esenciales están en el clúster |
| `complex_incomplete` | falta al menos una subunidad esencial |
| `complex_candidate` | hay arquitectura/vecindad, pero no firma completa |
| `complex_ambiguous` | varias interpretaciones posibles |
| `monomer_sufficient` | no se espera complejo |
| `not_evaluated` | sin datos suficientes |

Nunca se cierra una brecha metabólica con esto; solo se cambia la explicación y
la confianza.

## Estructura cuaternaria con CPU

CPU puede hacer **evaluación por plantilla**, no predicción de novo:

1. buscar en PDB un complejo homólogo (Foldseek/MMseqs) para las subunidades;
2. verificar que las cadenas del template correspondan a las subunidades;
3. transferir interfaz, estequiometría y cofactores si la identidad es
   suficiente;
4. marcar `template_complex` con identidad, cobertura y PDB;
5. si no hay template, `structural_complex_unknown`.

No es AlphaFold-Multimer. Para predecir una interfaz nueva hace falta GPU
(AlphaFold-Multimer, Boltz, etc.) y validación.

## Amarillo = sospecha, no cierre

Cada fila lleva `complex_flag`:

- `green`: complejo completo;
- `yellow`: complejo incompleto o candidato, **se sospecha** que falta una
  subunidad;
- `orange`: ambiguo;
- `grey`: no evaluado.

Una brecha amarilla no se cierra: se explica como “posible complejo
incompleto”. El texto `complex_note` lo dice y lista las subunidades faltantes.
No se cambia el estado metabólico a `confirmed`.

## Visualización en el navegador

### Pista `complexes`

- span desde el inicio de la primera subunidad hasta el fin de la última;
- color por estado:
  - verde: `complex_complete`;
  - ámbar: `complex_incomplete`;
  - gris: `complex_candidate`;
  - rojo: `complex_ambiguous`.
- clic → ficha:
  - nombre/firma;
  - subunidades y genes;
  - esperadas, presentes y ausentes;
  - evidencia (`operon`, `domain`, `synteny`, `template`);
  - confianza y límite;
  - enlace a cada proteína.

### Panel de completitud

Además de “% del genoma anotado”, mostrar:

```text
Proteínas: 1.301
  A función nombrada: 1.100
  B no caracterizada: 201
  C sin hit: 708
Complejos:
  completos: 12
  incompletos: 8
  candidatos: 15
Brechas con complejo no evaluado: 42
```

Así el usuario ve la diferencia entre “no hay gen” y “no se evaluó el complejo”.

### Ficha de brecha

Añadir:

- `complex_status`;
- `complex_id` si existe;
- enlace “ver complejo en el navegador”;
- texto: “esta reacción puede requerir un complejo; no se evaluó” cuando
  corresponda.

## Esquema de salida

`complexes.tsv`:

```text
complex_id  name  contig  start  end  subunits  expected  missing
status  evidence  template_pdb  template_identity  confidence  limitation
```

`protein_paths.tsv`:

- `complex_id` y `complex_role` (`subunit`, `candidate`, `none`);
- no cambia A/B/C; el complejo es una capa adicional.

## Límites

- Operón ≠ complejo.
- Vecindad ≠ interacción física.
- Un complejo puede tener genes dispersos; la restricción lo omitirá.
- Un operón puede codificar funciones no relacionadas.
- El template PDB puede ser de otra especie; la transferencia debe declarar
  identidad/cobertura.
- Sin GPU no se predice una interfaz nueva.
- No se debe usar pLDDT de monómero como evidencia de interfaz.

## Fases

| fase | contenido | esfuerzo |
|---|---|---|
| 1 | clústeres de vecindad + 10–15 firmas curadas + `complexes.tsv` | 2–3 días |
| 2 | pista y ficha en el navegador + panel de completitud | 1–2 días |
| 3 | template PDB/Foldseek e interfaz | 1 semana |
| 4 | AlphaFold-Multimer/Boltz opcional | GPU |

## Decisión

La primera fase debe ser **CPU y conservadora**: solo complejos con firma
conocida y genes cercanos. Es mejor un `complex_incomplete` honesto que un
complejo inventado.
