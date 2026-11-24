# Sprint UX: navegador, ficha de gen y taxonomía

## Hecho

- Especie activa visible **encima del menú**, con separación conceptual entre
  genoma y catálogo.
- Glosario y límites en la barra lateral.
- Navegador:
  - tabla “Elementos en pantalla” con contador y explicación;
  - clic/“Ver en el dibujo” con zoom y resaltado;
  - “Todos los elementos en pantalla” + exportación;
  - vista circular del contig con la ventana lineal marcada.
- Ficha de gen:
  - ADN, traducción y PDB descargables;
  - visor 3D si hay PDB;
  - A/B/C, confianza, EC y complejo explicados.
- Hipótesis estructurales:
  - sin `st.json` sobre texto plano;
  - F1–F5 explicados in situ;
  - tabla priorizada, ficha y visor 3D.
- Árbol de vida con los 400 y posición del genoma nuevo.
- El demo ya no bloquea el cambio de genoma.
- Nota de GPU sin referencia al modelo de ADN.
- Reportes de auditoría añadidos al centro de descargas.

## Límites

- Algunos textos nuevos están en español; la traducción completa de módulos
  sigue pendiente.
- La longitud del contig en la vista circular se estima desde la anotación.
- El árbol taxonómico es textual, no Plotly/SVG.
- La ficha de gen solo descarga ADN si el FASTA local está disponible.
- El catálogo de reportes no incluye todos los JSON.

## Pendiente

1. Agrupar el menú en Anotar / Explorar / Catálogo / Reportes.
2. Glosario contextual dentro de cada módulo.
3. Reporte descargable por módulo (no solo centro global).
4. Tree of life gráfico, no solo texto.
5. Completar i18n de visor, comparador y reportes.
