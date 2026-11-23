# Presentar el proyecto

Material breve para distintas audiencias. No prometas el modelo de ADN como
evaluado ni presentes una brecha como prueba de ausencia.

## Una frase

> Anotador open-source de genomas bacterianos que, además de encontrar genes,
> explica qué rutas metabólicas le faltan a la bacteria, con evidencia,
> confianza y una interfaz clara en tres idiomas.

## Pitch de 30 segundos

> Cuando secuenciamos una bacteria nueva, saber qué genes tiene es solo el
> principio: la pregunta real es qué puede hacer y qué le falta. Bacterial DNA
> Annotator anota el genoma desde cero, reconstruye su red metabólica, busca las
> brechas y propone operones y promotores. Todo corre en CPU, con una demo
> incluida y una GUI que un estudiante puede usar sin saber qué es DIAMOND.
> Cada resultado declara su evidencia y sus límites; no convierte una sospecha
> en un hecho.

## Charla de 5 minutos

1. **Problema (1 min):** anotar no es entender; las listas de brechas se
   sobreinterpretan.
2. **Método (1.5 min):** anotación → red metabólica → brechas → regulación;
   evidencia/confianza/estado.
3. **Producto (1 min):** demo, GUI ES/EN/ZH, navegador clicable, descargas.
4. **Validación (1 min):** 400 genomas como base, benchmarks, tests, CI.
5. **Límites y siguiente paso (0.5 min):** modelo experimental, datos no
   redistribuibles, TTS poco cubiertos.

## Para audiencia de secundaria

- Empieza con la analogía: el ADN es el manual; los genes, las instrucciones;
  una brecha, una receta a la que le falta un ingrediente.
- Muestra la demo y haz clic en un gen.
- No empieces por DIAMOND ni por F1: empieza por la pregunta biológica.
- Cierra con «la herramienta no adivina: dice qué sabe y qué no».

## Para audiencia de investigación

- Enfatiza la trazabilidad: `evidence`, `confidence`, `status`.
- Muestra el contrato de módulos y la CLI reproducible.
- Aclara que los 400 genomas son validación, no el producto.
- Presenta los límites: F1 de operones ≈ 0,78–0,81, cobertura de TTS, modelo
  no evaluado.
- Ofrece el DOI de Zenodo cuando esté disponible.

## Texto para LinkedIn

> Publiqué la v1.0 de Bacterial DNA Annotator, un anotador de genomas
> bacterianos que no se queda en la lista de genes: reconstruye la red
> metabólica, busca brechas y propone regulación, mostrando siempre la evidencia
> y los límites. Corre en CPU, incluye una demo y una GUI en español, inglés y
> chino. Es open source (MIT) y está pensado para estudiantes e investigadores.
> Feedback bienvenido: https://github.com/Krim10000/Bacterial-DNA-Annotator

## Texto para Reddit / Mastodon

> [OC] Bacterial DNA Annotator v1.0: a CPU-friendly, open-source pipeline and
> GUI for bacterial genome annotation, metabolic gaps, and operon prediction.
> Includes a 150 kb demo, evidence/confidence labels, and ES/EN/ZH interface.
> I’d love feedback from people who annotate genomes or teach bioinformatics.
> Code: https://github.com/Krim10000/Bacterial-DNA-Annotator

## Dónde enviarlo

| canal | qué enviar | requisito |
|---|---|---|
| GitHub | topics, release, repo fijado | ninguno |
| Zenodo | DOI de la release | conectar cuenta |
| bio.tools | ficha de herramienta | cuenta |
| JOSS | `paper.md` + revisión | software con tests/docs/licencia |
| BOSC/ISMB | póster o charla | call abierto |
| SoIBio / USACH | charla o póster | inscripción |
| r/bioinformatics, Biostars | post con demo | enlace y contexto |
| Streamlit forum | caso de uso | enlace |

## Checklist antes de presentar

- [ ] Release `v1.0.0` publicada y CI verde.
- [ ] Topics y descripción del repo.
- [ ] DOI de Zenodo (si está listo).
- [ ] Captura o GIF de 20–30 s de la demo.
- [ ] Enlace a `docs/INSTALL.md` y a la demo.
- [ ] Frase explícita sobre datos no redistribuibles y modelo experimental.
- [ ] Créditos y licencias visibles.
