# ARN y oriC en el navegador

## tRNA, rRNA y ncRNA

Al seleccionar un gen de tipo ARN:

- se muestra la secuencia en orientación 5′→3′;
- se descarga como FASTA de ARN (`.rna.fna`);
- **no se traduce** a proteína;
- no se muestra un visor 3D vacío ni se inventa una estructura.

No hay modelos 3D de ARN locales en el proyecto. La ficha lo declara.

## Horquillas

Al seleccionar una horquilla:

- coordenadas de la región;
- hebra;
- clase/confianza;
- ΔG, tallo, bucle, poli-U y genes vecinos;
- secuencia de ARN (T→U), si el genoma está disponible;
- aviso de que el plegamiento no se dibuja.

No se dibuja un tallo-bucle inventado: no hay estructura 3D de ARN local.

## oriC

`load_oric` conserva:

```text
pos       punto ORCA
start/end región DoriC si existe
source    ORCA/DoriC
confidence
span_source
```

La pista dibuja una **banda** para la región y una línea para el punto ORCA.
La ficha muestra:

- inicio, fin y longitud;
- punto ORCA;
- confianza;
- fuente;
- secuencia de la región, si existe;
- aviso de que oriC no tiene hebra codificante: la dirección no define un
  transcrito.

Si solo hay un punto, se declara `ORCA (punto)` y no se inventa un intervalo.
