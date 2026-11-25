# Qué hacer cuando estén las 2.570 estructuras

## Orden

```text
QC → filtros masivos → clasificación → complejos → docking selectivo → validación
```

## 1. QC y congelar tabla

- pLDDT/pTM.
- Excluir <0,60 / pTM <0,40.
- Recortar colas de baja confianza.
- Tabla maestra por cluster.
- Versión congelada.

## 2. Filtros estructurales masivos

- F1 Foldseek local.
- F2 EnzyMM/M-CSA.
- F3 P2Rank.
- F4 AlphaFill.
- F5 solo con mapa propio.
- pLDDT alto para interpretar.

## 3. Clasificación estructural

- Foldseek clustering.
- Familias CATH/ECOD.
- Dominios y límites.
- Comparación con PDB.
- Detección de fold conocido/nuevo.

## 4. Clase funcional amplia

- enzima;
- transportador;
- regulador;
- membrana;
- desordenada;
- desconocida.

Con Pfam/TMHMM/contexto, no solo estructura.

## 5. Complejos

Solo pares curados:

- P44104/P44105;
- P77795/YdcV-U-S;
- Q58207/LolD;
- P37494/ComA;
- P55680/Y4WA.

AlphaFold-Multimer/ColabFold:

- ipTM;
- PAE;
- estequiometría;
- no confirmar interacción in vivo.

## 6. Ligandos y docking

- AlphaFill/P2Rank.
- Docking solo si hay sustrato/cofactor conocido.
- No docking ciego masivo.
- Vina/Gnina para top.

## 7. Validación experimental

- clonar/expresar top;
- ensayo enzimático;
- complementación;
- termoforesis;
- SEC-MALS para complejos;
- cristalografía/criptoEM solo si hay caso fuerte.

## 8. Publicación

- tabla maestra;
- evidencia JSON;
- PDB;
- limitaciones;
- DOI.

## Límite

Las 169 >600 aa siguen aparte.
