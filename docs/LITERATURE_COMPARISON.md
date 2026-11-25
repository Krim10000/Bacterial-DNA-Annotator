# Comparación con literatura de DUFs y proteoma oscuro

## Referencias recientes

- Interacciones de proteínas con DUF a escala, AlphaFold3 + STRING:
  https://www.biorxiv.org/content/biorxiv/early/2026/08/20/2026.08.19.745780.full.pdf
- ECLIPSE / DPPS: priorización de proteoma oscuro:
  https://academic.oup.com/bioinformatics/article/42/8/btag613/8761985
- Revisión de predicción funcional en bacterias/fagos:
  https://journals.asm.org/doi/pdf/10.1128/mmbr.00022-25
- PARSE, hipótesis funcionales para folds nuevos:
  https://www.pnas.org/doi/pdf/10.1073/pnas.2513219122

## Comparación

| aspecto | nuestro plan | literatura reciente |
|---|---|---|
| escala taxonómica | 400 genomas | hasta 12.535 organismos |
| proteínas | 54.031 uncharacterized | pan-proteomas completos |
| estructura | ESMFold monómero + AFDB | AlphaFold2/3, MSA, complejos |
| familia | clustering + Foldseek | grafos all-vs-all, SSN, CLANS |
| dominio | Pfam/DUF | CATH/TED, ECOD, InterPro |
| contexto | operón/vecindad | GCsnap, STRING, co-expresión |
| evolución | plan: orthogroups/árbol | filogenómica y ganancia/pérdida |
| interacción | pares curados | STRING + AF3 a gran escala |
| benchmark | pendiente | controles, precisión/recall |
| validación | pendiente | experimentos en familias top |

## Nivel actual

### Sí está al nivel

- enfoque por familias, no proteínas sueltas;
- integración secuencia + estructura + contexto;
- uso de Foldseek/M-CSA/P2Rank/AlphaFill;
- trazabilidad y límites;
- priorización de candidatos;
- catálogo de hipótesis.

### Le falta para un paper funcional fuerte

1. **AlphaFold2/3/ColabFold** para familias top.
2. **MSA** y modelado de complejos.
3. **Grafo all-vs-all** de familias y SSN.
4. **CATH/TED/ECOD**.
5. **STRING/co-expresión** y GCsnap.
6. **Benchmark** con controles conocidos.
7. **Validación experimental**.
8. **Escala taxonómica** mayor que 400.

## Conclusión

- Como **recurso computacional y catálogo de hipótesis**: sólido.
- Como **paper de caracterización funcional**: todavía no.
- Conceptualmente alineado con 2026, pero por detrás en modelado, interacciones y validación.

## Qué haría para subir nivel

```text
1. Terminar estructuras.
2. Benchmark contra controles.
3. Familias all-vs-all + CATH.
4. Evolución/contexto.
5. AF2/AF3 en top.
6. Colaboración experimental.
7. Publicar recurso + preprint.
```
