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

## Integrado después

1. Menú agrupado en Anotar / Explorar / Catálogo / Reportes / Créditos.
2. `help_ui` contextual dentro de cada módulo.
3. Reportes descargables por módulo desde `downloads_ui`.
4. Árbol de vida SVG además del textual.
5. i18n de visor 3D, comparador, fichas, árbol y reportes.

## Límites que siguen

- El árbol SVG tiene zoom, reset y pantalla completa; sigue siendo un SVG generado localmente.
- La tabla de reportes sigue siendo un catálogo global, aunque haya descargas
  por módulo.
- El contenido de los reportes de auditoría está en español.
- La vista circular estima la longitud desde la anotación.
