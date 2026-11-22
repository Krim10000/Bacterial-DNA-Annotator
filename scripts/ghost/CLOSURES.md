# Cierre de brechas refinadas

## Método

1. Brecha refinada: productor **enzimático y no-transporte**, metabolito **consumido**
   en la red del genoma, no ión/ROS/currency, y conservada como `candidate_loss`.
2. Evidencia del candidato: dominio Pfam → EC → reacción Rhea.
3. Cierre si la reacción de la enzima está entre las productoras del metabolito,
   en **cualquier genoma miembro** del cluster.

Resultado: **8.655 brechas refinadas** en 309 genomas → **1 cierre** con evidencia
de dominio.

## Cierre encontrado

| campo | valor |
|---|---|
| cluster | `c050517` (3 genomas) |
| dominio | Glyoxalase (Pfam PF00903), EC **4.4.1.5** |
| genoma con la brecha | `GCF_000199675.1_ASM19967v1` |
| locus candidato | `NC_014960.1:824153-824347(+)` |
| brecha | **(R)-S-lactoilglutationa** (CHEBI:57474) |
| reacción | RHEA:19069 |
| descripción original | Uncharacterized protein YyaH (*Bacillus subtilis*) |

La glioxalasa cataliza la interconversión metilglioxal + glutatión ⇄
(R)-S-lactoilglutationa. El genoma no puede producir el compuesto y el cluster
conservado aporta el dominio enzimático.

**Matiz:** Rhea publica RHEA:19069 como maestro reversible y el sentido del
registro es la hidrólisis; el flujo de la traza lo cuenta como productor. Hay que
verificar la dirección curada antes de dar el cierre por definitivo.

## Casi-cierres (evidencia real, reacción distinta)

| candidato | dominio | brecha del genoma | por qué no cierra |
|---|---|---|---|
| c201454 | MoaC (Mo-cofactor) | molybdate | la brecha es transporte de un ión, no la reacción de MoaC |
| c201455 | MobB | molybdate | ídem |
| c037043 | RibD_C (riboflavina) | adenylyl-molibdopterina | la brecha es Mo-cofactor, no riboflavina |
| c050517 | Glyoxalase | molybdopterin | brecha no relacionada |
| c112479 | Glyoxalase | molybdopterin | brecha no relacionada |

## Lectura

El refinamiento funcionó: de 24.224 brechas crudas a 8.655 biosintéticas.
Pero el empalme dominio→reacción es estricto y solo un cluster cumple las tres
condiciones (genoma con la brecha + dominio de la enzima + reacción productora).

Los demás son **genes reales que no son la enzima de esa brecha**, exactamente
la distinción que planteábamos.

## Etiqueta del cierre

- `label`: `gap_closure_candidate`
- `label_source`: Pfam-A → EC → Rhea + conservación
- `label_evidence`: `pfam=Glyoxalase; EC=4.4.1.5; RHEA:19069; genoma=GCF_000199675.1; ng=4; direccion=UN`
- `label_confidence`: `medium` (reacción reversible; verificar dirección)
- `label_interpretation`: proteína conservada con el dominio de la enzima que produce el
  metabolito faltante; pendiente de validar dirección y especificidad.
