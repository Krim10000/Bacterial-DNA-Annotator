# Índice de literatura sigma

**Estado:** 2026-11-26. Revisión en curso.
Cada documento verifica DOI/PMID/URL; las afirmaciones sin verificación se marcan.

| Documento | Estado | Líneas | Referencias | Pregunta |
|---|---|---:|---:|---|
| `SIGMA_BIOLOGY.md` | terminado | 960 | 52 + 11 Pfam | Familias, dominios, ortología, conteos |
| `SIGMA_PREDICTION.md` | terminado | 910 | 37 | Predictores de promotor/sigma |
| `SIGMA_DATABASES.md` | terminado | 696 | 15 bases | TSS/sigma curados, licencias |
| `GENERALIZATION.md` | terminado | 1129 | 57 | LOSO/LOCO, transferencia, negativos |
| `ALT_SIGMA.md` | terminado | 1016 | 48 | Arquitecturas y motivos alternativos |
| `SIGMA_CONTEXT.md` | terminado | 2150 | 74 | Contexto/regulón |
| `SIGMA_STRUCTURE.md` | en curso | — | — | Estructura, dominios, bolsillos |

## Síntesis

### Lo que la literatura sostiene

1. No hay un panel universal de seis sigma.
2. Los nombres nativos no son ortología.
3. Los conteos varían: 1 en *M. genitalium*, 7 en *E. coli*, 17–18 en *B. subtilis*, 24 en *P. aeruginosa*, 65–66 en *S. coelicolor*.
4. PF04997 no es sigma; es Rpb1.
5. σ54 exige geometría −24/−12 y activador.
6. σ38 se solapa con σ70.
7. ECF tiene subgrupos y una base puede cambiar el regulón.
8. Solo RegulonDB (K-12) y DBTBS (B. subtilis) tienen TSS+sigma; BSGatlas es derivado.
9. El “gold DBTBS” local es BSGatlas: 1.777/2.226 Nicolas-only; 430 estrictos.
10. Las 400 no tienen gold sigma-TSS; no se puede reportar accuracy.
11. El top-1 premia el prior; macro-F1/cobertura son necesarios.
12. La generalización requiere LOSO/LOCO y negativos.
13. El contexto prioriza, no prueba. Predice mejor el módulo/condición que la identidad del sigma.
14. GBA está confundido por multifuncionalidad y densidad de anotación.
15. La cobertura de contexto es limitada: ~37–50 % de genes en el benchmark de *M. genitalium*; operones con sensibilidad 30–50 % en *E. coli*.
16. RegulonDB y EcoCyc no son independientes: comparten curación y códigos ECO.
17. La circularidad tiene al menos 10 formas; hace falta provenance holdout.
18. En nuestros datos, K-12 no gana (escala 0); B. subtilis sube +13,70 pp, pero gold y tablas comparten BSGatlas/DBTBS.

### Lo que no se puede afirmar

- Que un número de sigma coincida con la literatura.
- Que un nombre local equivalga a una familia global.
- Que una PWM única modele todos los promotores.
- Que una asociación contextual sea función.
- Que una estructura/bolsillo decida el promotor.
- Que haya precisión en las 400 sin gold.

## Implicaciones para el producto

El producto no es “contar sigma como la literatura”. Es:

```text
inventario por especie
+ orthogroup
+ estado conocido/desconocido
+ evidencia de promotor
+ contexto/regulón
+ hipótesis funcional
+ confianza
+ abstención
```

## Bibliografía adicional

Basta de búsqueda amplia por ahora. Falta una revisión en curso:

```text
SIGMA_STRUCTURE.md
```

Después solo búsquedas dirigidas a huecos concretos:

- ECF por subgrupo;
- σ38 vs σ70;
- σ54 activadores;
- gold TSS multi-especie;
- validación LOSO/LOCO.

## Cómo citar

Cada documento contiene DOI/PMID verificados. Para afirmaciones del producto,
citar el documento de literatura y la sección; no citar de memoria.
