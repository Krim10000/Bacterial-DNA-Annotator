# De MF/CC a Biological Process

## Regla

**MF y CC no implican BP automáticamente.**

Una proteína puede tener una función molecular y una localización, pero eso no
prueba en qué proceso participa, ni en qué condiciones, ni en qué organismo.

## Cuándo sí

Se puede asignar BP si hay una de estas evidencias:

1. Anotación GO BP directa con evidencia experimental.
2. Anotación GO BP curada.
3. Ortólogo curado con BP y reciprocidad suficiente.
4. Ruta completa: enzima + reacción + ruta + contexto.
5. Operón/regulón conservado con evidencia.

## Cuándo no

No asignar BP fuerte por:

- una enzima aislada;
- una reacción reversible compartida;
- una localización CC;
- un dominio Pfam;
- una homología remota;
- una ruta incompleta;
- un contexto de operón solo.

## Evidencia

| Fuente | Código | Fuerza |
|---|---|---|
| GO BP experimental | EXP/IDA/IMP | alta |
| GO BP curado | TAS/NAS | media-alta |
| ortología | ISS | media |
| ruta/EC | IEA/ISS | media-baja |
| contexto | IEA/ISS | baja |
| sin datos | — | no asignar |

## Ejemplo: ciclo de Krebs

```text
enzima + reacción TCA + ruta completa
→ GO:0006099 tricarboxylic acid cycle
evidence = inferred_from_pathway
confidence = medium
```

Pero:

- fumarasa, aconitasa, succinato deshidrogenasa pueden tener otras rutas;
- una sola enzima no demuestra que el ciclo funcione;
- una ruta incompleta no justifica el BP completo.

## Ejemplo: rRNA

```text
rRNA + GO:0005840 ribosome (CC)
+ GO:0003735 structural constituent of ribosome (MF)
+ GO:0006412 translation (BP)
```

El BP se asigna por tipo de feature curado, no por inferir de CC/MF.

## Decisión operativa

```text
BP directo       → usar
BP por ruta      → usar con confidence medium y source=ruta
BP por contexto  → solo hipótesis, no anotación firme
sin evidencia    → no_data
```

## Cita

Ver `docs/LABELS_CITATIONS.md` y la documentación de GO/ECO.
