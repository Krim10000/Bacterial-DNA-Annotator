# Límites de CPU y ruta GPU

## Regla de interpretación

Ninguna capa cierra una brecha metabólica por sí sola. Las etiquetas dicen qué
se sabe y qué falta:

```text
A/B/C          función putativa por homología
domain         hipótesis de familia/dominio
structural     modelo disponible o evaluado
complex        sospecha de ensamblaje
gap            brecha; no se cierra sin evidencia
```

## Qué hace CPU

| capa | CPU | evidencia | límite |
|---|---|---|---|
| A/B/C | DIAMOND/Swiss-Prot | homología | no función confirmada |
| B estructural | busca AFDB local, manifiesto | modelo disponible | no pliega ni filtra |
| C dominio | hmmscan panel Pfam | familia/dominio | panel curado, no Pfam completo |
| complejos | firmas + K-12 + transferencia | vecindad/ortología | no estructura cuaternaria |
| membrana | hmmscan + hidropatía/señal/lipobox | heurística | no TMHMM/SignalP/LipoP |
| red | Rhea/EC | reacción putativa | no cierra brechas |
| oriC | ORCA/GC skew/dnaA/DoriC | estado | DoriC opcional |

## Qué no puede CPU

- plegar miles de proteínas;
- predecir interfaz/estequiometría;
- calcular ipTM/pDockQ;
- resolver coevolución de interfaz;
- demostrar función por estructura;
- cerrar una brecha sin reacción/EC/experimento.

## Ruta GPU

GPU acelera estructura, no verdad experimental:

1. **ESMFold local** para el top de B/C.
2. **AlphaFold-Multimer/ColabFold** para pares candidatos de complejo.
3. **Boltz-2** para co-folding/ipTM/pDockQ.
4. **Foldseek-Multimer** para buscar el ensamblaje en PDB.
5. **Coevolución** con MSA profunda (CPU/GPU).

Toda salida GPU se etiqueta:

```text
structural_evidence = predicted_model
confidence = pLDDT/ipTM/pDockQ
closure = no
```

### Cómo seguir en GPU

```bash
# ESMFold local para brechas priorizadas
python -m webapp.pipeline genoma.fna --outdir results/ \
  --esmfold local --esmfold-gaps 3 --structure-evidence --structural-filters

# AlphaFold-Multimer/Boltz: solo el top de pares candidatos
# (scripts/modelado; requiere GPU y MSA)
```

No se debe correr sobre 647.588 C ni sobre todos los pares. Priorizar:

- B/C con dominio/contexto;
- complejos con `frac_mapped ≥0.8`;
- brechas con soporte en 400;
- pares con evidencia de operón/sintenia.

## Cómo leer la GUI

- `structural_status=afdb_model_available`: hay modelo, no está evaluado.
- `structural_status=no_afdb_model`: falta modelo; no es “sin estructura”.
- `domain`: hipótesis de dominio; no función.
- `complex_suspected`: sospecha; no cierre.
- `complex_incomplete`: falta subunidad; no cierre.
