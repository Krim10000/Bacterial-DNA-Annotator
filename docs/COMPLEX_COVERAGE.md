# Cobertura de nombres “subunit” fuera del catálogo K-12

Pregunta: de los 105 genes K-12 cuyo producto dice *subunit/component* pero que
no están en el gold curado, ¿cuántos cubre el pipeline?

## Medición

Script: `tools/audit_subunit_coverage.py`
Resultado: `reports/subunit_coverage.json`.

| métrica | valor |
|---|---:|
| nombres fuera del catálogo | 105 |
| cubiertos por el modelo K-12 (candidato) | **88 (83,8 %)** |
| sin cobertura | **17 (16,2 %)** |
| evidencia del modelo | 81 `nuevo`, 7 `parcial` |

“Cubierto” significa que el gen aparece en un cluster candidato del modelo. **No**
significa que sea una subunidad confirmada ni que el complejo exista.

## Los 17 sin cobertura

Son casos que el modelo actual no agrupa:

- PTS: `yadI`, `nagE`, `glvB`, `glvC`;
- translocón Sec: `yfgM`;
- RNase P: `rnpB` (RNA catalítico);
- citocromo bd accesorio: `cydH`;
- hidrogenasa 4: `hyfA`;
- septo: `ymgF`;
- fimbria: `elfA`;
- conector de dos componentes: `safA`;
- flippase: `arnE`;
- maduración ribosomal: `era`, `lepA`;
- transportador ABC: `phnE1`.

## Cómo cubrirlos

### 1. Modelo como candidato, no como catálogo

Los 88 ya están cubiertos como `nuevo`/`parcial`. Se etiquetan amarillo. No hace
falta forzar el catálogo, pero sigue siendo una hipótesis.

### 2. Firmas curadas para familias recurrentes

Añadir firmas para:

- PTS (`ptsI`, `ptsH`, `crr`, `nagE`, `glvB/C`);
- citocromo bd y accesorios (`cydA/B`, `cydH`);
- translocón Sec (`secY/E/G`, `secA`, `yfgM`);
- RNase P (`rnpA`, `rnpB`);
- hidrogenasas;
- transportadores ABC por dominio, no por nombre.

### 3. Dominios/HMM

Un transportador ABC o un PTS EII se reconoce por arquitectura de dominios aunque
el nombre no diga “subunit”. Es la vía más general para C y B.

### 4. Catálogo parcial

`fetch_complexportal.py --allow-incomplete` conserva complejos con subunidades no
mapeadas. Se pueden usar como **sospecha parcial**, no como gold completo.

### 5. Otras fuentes

- EcoCyc/MetaCyc (licencia local, K-12);
- ComplexPortal por taxid;
- PDB/IntAct/BioGRID para interacciones;
- UniProt “Subunit” con NER;
- masa/estequiometría y Foldseek-Multimer.

## A/B/C no es lo mismo que complejo

- **A:** función nombrada por homología.
- **B:** homólogo no caracterizado/DUF.
- **C:** sin hit Swiss-Prot.
- **`subunit_hint`:** el nombre del producto menciona *subunit/component*.

Son ejes independientes. Un gen puede ser A y no ser subunidad, o C y ser una
subunidad real. Añadir `subunit_hint` **no cambió** A/B/C.

Los 88 cubiertos por el modelo son sobre todo A; los 17 sin cobertura pueden ser
A, B o C. Para cubrirlos hay que usar dominio/contexto/estructura, no renombrar
A/B/C.
