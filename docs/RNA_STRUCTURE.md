# Estructura secundaria de ARN (tRNA y horquillas)

`webapp/rna_structure.py` predice **estructura secundaria** en formato
dot-bracket para una secuencia de ARN corta. Es una **predicción
computacional**, no una estructura 3D y no una validación experimental.

- No genera coordenadas atómicas.
- No descarga PDB ni consulta estructuras experimentales.
- No dibuja nada en la GUI (la API queda lista para conectarla).
- No inventa pares: si no hay algoritmo disponible, la fuente es
  `not_available` y el dot-bracket queda vacío.

## Dependencia y estado

La fuente preferida es **ViennaRNA** (`RNAfold`). En este entorno ya está
instalada en `.venv`:

```bash
.venv/bin/python -c "import RNA; print(RNA.__version__)"
# 2.7.2
```

`requirements.txt` declara `ViennaRNA>=2.6` y `pyproject.toml` la incluye en el
extra `full`. Si faltara, se puede instalar en el entorno virtual (rueda CPU):

```bash
.venv/bin/python -m pip install "ViennaRNA>=2.6"
```

No hace falta GPU ni red para plegar. Sin ViennaRNA, la API usa el respaldo
Nussinov (si se permite) o devuelve `not_available`.

## API

```python
from webapp.rna_structure import predict_structure

resultado = predict_structure(secuencia)
resultado.dot_bracket        # "((((....))))"
resultado.energy_kcal_mol    # MFE en kcal/mol (solo vienna)
resultado.score              # número de pares de bases
resultado.source             # "vienna" | "nussinov" | "not_available"
resultado.pairs              # ((0, 9), (1, 8), ...) índices 0-based
resultado.reason             # por qué se usó ese backend, o ""
resultado.as_dict()          # dict JSON-serializable para la GUI
```

Firma:

```python
predict_structure(
    sequence: str,
    *,
    backend: str = "auto",        # "auto" | "vienna" | "nussinov"
    allow_fallback: bool = True,
) -> RnaStructure
```

Comportamiento:

| `backend` | ViennaRNA presente | ViennaRNA ausente o falla |
| --- | --- | --- |
| `auto` | `vienna` | `nussinov` si `allow_fallback`; si no, `not_available` |
| `vienna` | `vienna` | `nussinov` si `allow_fallback`; si no, `not_available` |
| `nussinov` | `nussinov` | `nussinov` |

La secuencia se normaliza: se quitan espacios/saltos de línea, se pasa a
mayúsculas y `T` se convierte en `U`. Solo se aceptan `A`, `C`, `G`, `U`.
Una base ambigua (`N`, `R`, `Y`, ...) se rechaza con `ValueError`: no se puede
emparejar sin inventar.

### Campos de `RnaStructure`

| Campo | Contenido |
| --- | --- |
| `sequence` | ARN normalizado (mayúsculas, T→U) |
| `dot_bracket` | `.` sin par, `(` y `)` para cada par; vacío en `not_available` |
| `energy_kcal_mol` | MFE de RNAfold en kcal/mol; `None` en Nussinov/`not_available` |
| `score` | número de pares de bases; `None` en `not_available` |
| `source` | `vienna`, `nussinov` o `not_available` |
| `pairs` | pares `(i, j)` 0-based con `i < j` sobre `sequence` |
| `reason` | explicación del backend o del fallo; `""` si todo fue normal |

## Ejemplos (ViennaRNA 2.7.2, 37 °C)

### Horquilla corta

```python
>>> r = predict_structure("GGGGAAAACCCC")
>>> r.dot_bracket
'((((....))))'
>>> round(r.energy_kcal_mol, 2)
-5.4
>>> r.score, r.source
(4, 'vienna')
```

### tRNA-Phe de *S. cerevisiae* (76 nt)

```python
>>> trna = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
>>> r = predict_structure(trna)
>>> r.dot_bracket
'(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'
>>> round(r.energy_kcal_mol, 2)
-22.4
>>> r.score, r.source
(21, 'vienna')
```

La estructura tiene tres bucles de horquilla internos (compatibles con D,
anticodón y T) y un tallo aceptor. Es una predicción MFE: no garantiza el
cloverleaf biológico ni incluye modificaciones de bases.

### Secuencia sin pares

```python
>>> predict_structure("AAAAAAAA").dot_bracket
'........'
```

### Respaldo forzado (sin energía)

```python
>>> r = predict_structure("GGGGAAAACCCC", backend="nussinov")
>>> r.dot_bracket, r.energy_kcal_mol, r.score, r.source
('((((....))))', None, 4, 'nussinov')
```

### Sin herramienta

Si ViennaRNA no está instalado y se desactiva el respaldo:

```python
>>> r = predict_structure(trna, backend="vienna", allow_fallback=False)
>>> r.source, r.dot_bracket, r.pairs
('not_available', '', ())
```

No se rellena con una estructura "plausible" inventada.

## Algoritmos

### `vienna` (preferido)

- `RNA.fold(secuencia)` de ViennaRNA: mínima energía libre (MFE).
- Parámetros por defecto de la biblioteca: 37 °C, modelo Turner 2004.
- Energía en kcal/mol (`energy_kcal_mol`).
- Sin restricciones de pares, sin pseudoknots y sin modelar ligandos,
  Mg²⁺, modificaciones ni plegamiento co-transcripcional.
- El dot-bracket devuelto se revalida antes de confiar en él (longitud,
  balance de paréntesis y pares canónicos).

### `nussinov` (respaldo)

Programación dinámica clásica de máximo número de pares:

- Pares permitidos: `AU`, `UA`, `GC`, `CG`, `GU`, `UG` (incluye wobble).
- Sin pseudoknots.
- Tallo-bucle mínimo de 3 nt.
- `dp[i][j]` = máximo de pares en `sequence[i:j+1]`; se prueban dejar `i`
  sin parear, dejar `j` sin parear, parear `(i, j)` y todas las
  bifurcaciones `k`.
- Desempate determinista: `i` sin parear > `j` sin parear > par `(i, j)` >
  primera bifurcación.
- `score` es el **número de pares**, no una energía. No es comparable con
  `energy_kcal_mol`.
- Coste O(n³); por eso solo se aplica hasta
  `FALLBACK_MAX_LENGTH = 200` nt. Por encima devuelve `not_available`.

El respaldo no es un modelo termodinámico: en tRNAs puede devolver más pares
que ViennaRNA y una topología distinta. Sirve para comprobar que una
horquilla puede formarse cuando falta la herramienta, no para reportar ΔG.

## Determinismo

- No hay aleatoriedad ni semilla.
- La misma secuencia, versión de ViennaRNA y parámetros producen el mismo
  dot-bracket.
- Nussinov tiene un orden de desempate fijo (documentado arriba).
- Una versión distinta de ViennaRNA puede cambiar la MFE y, en casos
  ambiguos, el dot-bracket: la versión usada debe registrarse junto al
  resultado si se compara entre ejecuciones.

## Límites (leer antes de usar)

1. **No es estructura 3D.** El dot-bracket describe pares de bases, no
   coordenadas, dominios, interacciones terciarias ni un modelo molecular.
2. **No valida nada experimentalmente.** Un resultado plausible no confirma
   la estructura in vivo; no sustituye SHAPE, DMS, cristalografía, cryo-EM
   ni NMR.
3. **No descarga PDB** ni usa plantillas estructurales.
4. **Sin pseudoknots** en ninguna de las dos fuentes.
5. **Sin bases modificadas ni ambiguas.** Se rechazan `N` y códigos IUPAC en
   lugar de adivinar.
6. **Sin control de temperatura, iones o restricciones.** ViennaRNA usa sus
   valores por defecto; no hay parámetros expuestos para Mg²⁺, ligandos,
   temperatura o pares forzados.
7. **Sin plegamiento co-transcripcional.** Es una MFE de la secuencia
   completa, no la estructura cinética de un transcrito naciente.
8. **El respaldo Nussinov no da energía** y puede diferir del cloverleaf de
   un tRNA. Para tRNA usar `vienna`.
9. **Límites de longitud:** 1000 nt para la API (`MAX_SEQUENCE_LENGTH`) y
   200 nt para el respaldo (`FALLBACK_MAX_LENGTH`). tRNA y horquillas están
   muy por debajo.
10. **No integra la GUI todavía.** No hay pestaña, ficha ni dibujo; la API y
    `as_dict()` están listas para que el root las conecte.
11. **No clasifica terminadores.** Una horquilla predicha no decide fronteras
    de operón ni demuestra terminación; la clase/confianza de una horquilla
    sigue siendo una pista separada.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_rna_structure.py -q
```

Cubren:

- horquilla conocida y tRNA-Phe con ViennaRNA;
- determinismo;
- respaldo Nussinov y ausencia de energía;
- simulación de ViennaRNA ausente (`not_available`, sin pares inventados);
- fallo de `RNA.fold`;
- secuencia vacía, bases ambiguas, longitudes y parámetros inválidos;
- validación de dot-bracket y pseudoknots.
