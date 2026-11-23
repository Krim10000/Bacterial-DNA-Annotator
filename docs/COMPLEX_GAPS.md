# Brechas metabólicas y estructura cuaternaria

## Respuesta corta

**Sí, puede explicar una parte de las brechas.** El pipeline actual evalúa
proteínas individuales: DIAMOND asigna cada proteína a UniProt/Rhea y construye
la red con esas reacciones. **No exige que estén todas las subunidades de un
complejo obligado, no modela estequiometría y no ensambla el complejo.**

Pero no es la única causa y probablemente no es la mayoría. Hay que separar:

```text
capacidad enzimática presente
  ├─ proteína individual suficiente
  ├─ complejo obligado incompleto (subunidad ausente/no anotada)
  ├─ complejo presente pero mal anotado (EC/Rhea en otra subunidad)
  ├─ homólogo remoto no detectado por Swiss-Prot
  ├─ enzima presente pero reacción no mapeada a Rhea
  ├─ ruta alternativa / transporte / compartimento
  └─ brecha real o definición incompleta de la ruta
```

## Qué hace hoy el pipeline

`build_network()`:

1. DIAMOND de cada proteína contra Swiss-Prot;
2. accesión UniProt → Rhea;
3. reacciones con al menos un productor;
4. brecha = reacción sin productor.

Limitaciones:

- no hay tabla de complejos;
- no hay subunidades obligadas ni estequiometría;
- un complejo puede quedar como brecha aunque todas sus subunidades estén
  presentes si ninguna mapea la reacción;
- una subunidad puede mapear la reacción aunque falten las demás;
- los filtros estructurales usan monómeros (AlphaFold DB/ESMFold);
- no hay AlphaFold-Multimer ni co-folding de complejos.

Para no sobreinterpretar, cada brecha ahora incluye:

| campo | valor actual |
|---|---|
| `complex_status` | `not_evaluated` |
| `quaternary_status` | `unknown` |
| `complex_evidence` | vacío |
| `complex_note` | explica que no se ensamblan complejos |

## Ejemplos donde la estructura cuaternaria importa

- ATP sintasa (varias subunidades, EC 7.1.2.2 / 3.6.3.14).
- Citocromo c oxidasa y NADH deshidrogenasa.
- Succinato deshidrogenasa.
- Nitrogenasa.
- Piruvato deshidrogenasa y 2-oxoglutarato deshidrogenasa.
- Ácido graso sintasa.
- Complejos de transporte de membrana.

En estos casos, la actividad catalítica puede residir en una subunidad, pero la
reacción fisiológica requiere el complejo. La anotación individual puede ser
correcta y la red seguir mostrando una brecha.

## Auditoría de los 400

`tools/audit_complexes_400.py` sobre 400 genomas:

- 2.218 complejos por vecindad;
- 639 completos y 1.579 incompletos;
- 8.655 brechas;
- **0 brechas asociadas por EC/Rhea** en esta primera pasada.

El 0 no prueba que los complejos no importen: la asociación exige que una EC del
complejo esté mapeada a Rhea y que su producto sea exactamente el metabolito de
la brecha. Muchas subunidades no tienen EC/Rhea asignada. Por eso las brechas
quedan `not_evaluated` y no se cierran.

## Primer resultado real

En `GCF_000006985.1_ASM698v1` (158 brechas, 6 complejos por vecindad):

- 0/158 brechas pudieron asociarse a un complejo por EC/Rhea;
- los complejos detectados (ATP sintasa, NADH deshidrogenasa, nitrogenasa,
  citocromo bc1, Fts) no explican las brechas de este genoma;
- las 158 quedan `not_evaluated`, no `closed`.

Interpretación: en este genoma la estructura cuaternaria **no** es la
explicación principal de las brechas. No se fuerza una asociación.

## Cómo diagnosticarlo sin GPU

1. **Mapa reacción → complejo** con Complex Portal, EcoCyc/MetaCyc (si hay
   licencia), KEGG modules o literatura.
2. **Subunidades esperadas** por complejo y su EC/UniProt.
3. **Comprobar presencia** de cada subunidad en el proteoma.
4. **Usar contexto**: operón, sintenía y conservación entre los 400.
5. **Etiquetar** la brecha:
   - `complex_complete`: todas las subunidades esenciales presentes;
   - `complex_incomplete`: falta al menos una;
   - `complex_ambiguous`: no se puede decidir;
   - `monomer_sufficient`: no requiere complejo.
6. **No cerrar** la brecha con esto: solo cambiar la explicación y la confianza.

Diseño detallado de evaluación CPU y visualización:
[`COMPLEX_CPU.md`](COMPLEX_CPU.md).

## Qué implementar

### Fase CPU, 2–4 días

- tabla `data/metabolism/complexes.tsv` con `complex_id`, `subunit`,
  `uniprot`, `ec`, `essential`, `source`, `license`;
- función `annotate_complexes(gaps, proteins, ...)`;
- campos `complex_status`, `complex_evidence`, `quaternary_status`;
- salida en `gaps.tsv` y en el informe;
- tests con 3–5 complejos conocidos.

### Fase estructura, GPU opcional

- AlphaFold-Multimer/Boltz para complejos candidatos;
- interfaz con PDB/AlphaFold DB;
- no usar pLDDT de monómero como evidencia de complejo.

## Lo que no hay que hacer

- No decir “la enzima está ausente” si el complejo está incompleto.
- No decir “el complejo está presente” solo porque una subunidad aparece.
- No cerrar una brecha por una predicción de estructura.
- No usar un monómero para concluir una interfaz.

## Relación con los tres caminos

Una subunidad de complejo suele caer en:

- **A** si tiene función nombrada;
- **B** si es uncharacterized/DUF;
- **C** si no tiene hit Swiss-Prot.

El pipeline ya clasifica A/B/C; lo que falta es la capa de complejo. Se puede
añadir `complex_subunit` como subcategoría de A/B/C y no como un cuarto camino
independiente.

## Próximo paso recomendado

Antes de tocar la estructura cuaternaria, hacer un **diagnóstico de 20–30
brechas** en 3 genomas:

1. ¿Cuántas tienen una subunidad conocida de complejo?
2. ¿Cuántas tienen todas las subunidades pero ninguna mapea Rhea?
3. ¿Cuántas no tienen ningún homólogo?
4. ¿Cuántas son transporte/ruta alternativa?

Eso evita implementar un modelo de complejos para un problema que quizá sea
mayoritariamente de homología o de base de datos.
