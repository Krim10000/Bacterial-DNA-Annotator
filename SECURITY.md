# Seguridad

## Qué cubre este proyecto

La aplicación se ejecuta **en el equipo del usuario**. Los FASTA que se suben no
se envían a ningún servidor, salvo que el usuario active explícitamente una
llamada a una API externa (por ejemplo, ESMFold por Neurosnap). En ese caso la
secuencia sale del equipo y el usuario debe autorizarlo.

## Riesgos conocidos

- **Subida de archivos**: la GUI acepta FASTA de hasta 200 MB. Un archivo
  malformado puede consumir memoria; el límite de Streamlit y el tamaño del
  genoma acotan el riesgo.
- **Comandos externos**: el pipeline invoca binarios (`diamond`, `hmmscan`,
  `tRNAscan-SE`, …). No ejecutes el proyecto como root.
- **APIs externas**: las claves (`ESM_API_KEY`) se leen de `.env` o del entorno.
  No las pegues en un *issue* ni en un chat. El archivo `.env` está ignorado por
  git.
- **Datos**: no se redistribuyen bases con licencia académica; ver
  `DATA_NOTICE.md`.

## Cómo reportar una vulnerabilidad

No abras un *issue* público. Escribe a **Camilo.Rojas@usach.cl** con:

- versión (`cat VERSION` o `pip show dna-poc-annotator`),
- sistema operativo,
- pasos para reproducir,
- impacto que crees que tiene.

Se responde en un plazo razonable y se coordina la divulgación.

## Versiones soportadas

La rama `v1.0` recibe correcciones de seguridad. Las versiones anteriores no.
