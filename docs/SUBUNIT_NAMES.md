# ¿El nombre “subunit” indica que es una subunidad?

No. Es una pista débil. Medimos el problema en K-12 contra el catálogo curado
ComplexPortal, no contra verdad experimental.

## Medición K-12

| métrica | valor |
|---|---:|
| genes | 4.523 |
| subunidades en complejos curados | 731 |
| genes con `subunit/component` en el producto | **513** |
| de esos, dentro de un complejo curado | 408 |
| de esos, **fuera** del catálogo | **105 (20,5 %)** |
| precisión del nombre vs catálogo | **79,5 %** |
| recall del nombre vs catálogo | **55,8 %** |
| subunidades curadas sin `subunit/component` | **44,2 %** |

Script: `tools/audit_subunit_names.py`; resultado:
`reports/subunit_name_audit.json`.

## Cómo interpretar el 20,5 %

**No** significa que 20,5 % de las proteínas llamadas “subunit” no sean
subunidades. Significa que **no están en el catálogo curado usado como
referencia**. Muchas son subunidades reales de complejos que el catálogo no
incluye o no mapea:

- `dnaX` — DNA polymerase III subunit tau;
- `holA` — DNA polymerase III subunit delta;
- `b4506` — 50S ribosomal subunit protein L36B;
- `yadH` — ABC transporter membrane subunit;
- `pxpB`, `pxpC` — 5-oxoprolinase components.

Por eso el 20,5 % es una **tasa de no-coincidencia con el catálogo**, no una
tasa de falsos “subunit”.

Para medir la tasa real de nombres falsos haría falta verdad experimental
(ComplexPortal ampliado, co-IP, MS nativo, PDB) y una definición de “ser
subunidad”. Hoy no la tenemos.

## El caso “chain” es peor

Si se busca `chain` sin restringir, aparecen falsos positivos claros:

```text
branched chain amino acid transporter
long-chain acyl-CoA thioesterase
peptide chain release factor
```

Con `chain` genérico: 18 nombres, solo 6 dentro del catálogo y 12 fuera
(66,7 % fuera). Por eso el pipeline usa solo `subunit/component` como pista,
nunca `chain` genérico.

## Lo que hace el pipeline ahora

Añade dos columnas a `protein_paths.tsv`:

```text
subunit_hint = named_subunit | none
subunit_evidence = product_name_only
```

No se usa para asignar complejo ni para cerrar brechas. Es una pista para el
usuario y una feature débil para el modelo, no una etiqueta.

## Cómo mejorar la pipeline

1. **No usar el nombre como etiqueta dura.** Peso bajo, siempre con
   `product_name_only`.
2. **Combinar con vecindad/operón/sintenia.** El nombre solo no basta; la
   vecindad aporta el contexto de “se producen juntos”.
3. **Catálogo más amplio.** ComplexPortal por taxid, no solo K-12.
4. **Transferencia por ortología.** Ya implementada; nivel B, no prueba.
5. **Dominios/HMM.** Un transportador ABC tiene dominios y arquitectura
   reconocibles aunque su nombre no diga “subunit”.
6. **Masa/estequiometría.** Detecta “falta una subunidad”.
7. **Modelo K-12.** Ya integrado; calibrar en K-12 y no asumir calibración en
   otras especies.
8. **Estructura.** Foldseek-Multimer/PDB; hoy falta DB de ensamblajes.
9. **Validación experimental.** La única forma de medir el falso “subunit” real
   es contrastar con complejos medidos, no con nombres.

## Regla

```text
nombre “subunit” ≠ subunidad
nombre sin “subunit” ≠ no subunidad
```

La pertenencia a complejo se decide por evidencia combinada y se etiqueta como
sospecha (`complex_suspected`/`complex_candidate`), no como `confirmed`.
