# Sets automáticos de ~20 minutos

## Criterio

- Excluidas las **587** ya resueltas.
- Restantes: **1.983**.
- Ordenadas de **más corta a más larga**.
- Empaquetadas por tiempo estimado, objetivo ~19 min.
- 49 sets.

## Modelo de tiempo

Calibrado con las 587 métricas reales:

```text
segundos ≈ 10,2 − 0,130·L + 0,000625·L²
```

## Tamaño

| sets | proteínas | tiempo |
|---|---:|---:|
| 49 | 1.983 | ~15,4 h total |

Sets iniciales:

- set_001: 290 proteínas cortas;
- set_002: 282;
- set_003: 183;
- luego sets más pequeños, hasta las más largas.

## Uso

1. Subir el ZIP como dataset.
2. Abrir el notebook.
3. Ejecutar `run_next()`.
4. Descargar el ZIP del set.
5. Repetir.

## Regla

- Un set por ejecución.
- No cargar el modelo dos veces.
- Descargar antes de seguir.
- No mezclar sets.

## Por qué así

- Las cortas se procesan rápido.
- Los sets grandes iniciales siguen dando ~20 min.
- Los sets finales tienen pocas proteínas largas.
- Se evita un bloque gigante de 1.983.
