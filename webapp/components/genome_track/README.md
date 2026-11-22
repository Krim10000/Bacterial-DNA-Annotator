# genome_track

Componente Streamlit mínimo (HTML + JS, sin build) que dibuja el SVG del
navegador de genoma y devuelve el elemento pulsado.

- Entrada (`args`): `svg` (cadena SVG), `show_hint` (bool).
- Salida: `{kind, id, ts}` donde `kind` es `gene`, `operon`, `promoter`,
  `hairpin`, `tss`, `tts` o `suspect`.
- Protocolo: `streamlit:componentReady`, `streamlit:render`,
  `streamlit:setComponentValue`, `streamlit:setFrameHeight`.

No usa `streamlit-component-lib` ni npm: se mantiene como un solo archivo para
que el repositorio no necesite un paso de compilación.
