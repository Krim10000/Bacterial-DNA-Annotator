# Comparar dos genomas

GUI: pestaña **🔀 Comparar genomas**.
Módulo: `webapp/compare_ui.py`.

## Qué hace

- lee dos GFF3;
- agrupa CDS multiparte;
- empareja por `locus_tag`, `ID` o coordenadas exactas;
- muestra compartidos y exclusivos;
- descarga TSV.

## Qué no hace

- no lee ADN/proteína;
- no calcula identidad;
- no infiere ortología;
- no compara función;
- no usa red.

Una coincidencia por ID/coordenadas es solo de anotación: puede ser un contador
local, un fragmento o una diferencia real. No es homología.
