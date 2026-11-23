# Benchmark de proyectos similares

**Fecha:** 2026-11-23. **Método:** metadatos de GitHub (estrellas, última
actividad, licencia) y lectura de README/documentación. Los servicios web sin
API se marcan como *no verificado*. Las estrellas y fechas cambian.

## Ficha del proyecto evaluado

| dato | valor |
|---|---|
| Repositorio | `Krim10000/Bacterial-DNA-Annotator` |
| Versión | 1.1.3 |
| DOI | 10.5281/zenodo.24469732 |
| Estrellas / forks | 1 / 0 (al 2026-11-23) |
| Descargas del DOI | 0 (al 2026-11-23) |
| Licencia | MIT (GitHub la detectaba como NOASSERTION; se corrigió el texto) |
| Interfaz | ES / EN / ZH |
| Módulos | 10 |
| SBML / FBA | no |

**Consecuencia honesta:** no hay adopción que citar. Cualquier comparación debe
basarse en funciones y benchmarks, no en popularidad.

## Anotación, función y estructura

| proyecto | qué hace / qué no | entrada | GUI/CLI | GPU | licencia | ★ / última act. | diferenciador |
|---|---|---|---|---|---|---|---|
| **Prokka** | anotación procariota rápida; no metabolismo/operones/GUI | FASTA | CLI | no | GPL-3.0 | 1005 / 2026-01 | estándar de facto en anotación base |
| **Bakta** | anotación con UniRef, sORF, AMR, ncRNA, oriC, dbxref; no metabolismo/operones | FASTA/GFF/proteína | CLI + web | no | GPL-3.0 | 696 / 2026-10 | DB curada y submission-ready |
| **PGAP** | pipeline NCBI de referencia; no metabolismo/operones/GUI local | FASTA | CLI/CWL + web | no | NOASSERTION | 392 / 2026-08 | calidad para submission |
| **DFAST** | anotación + envío DDBJ; no metabolismo/operones | FASTA/GFF | CLI + web | no | sin licencia declarada | 87 / 2026-11 | submission |
| **Prodigal / pyrodigal** | gene calling; no anotación funcional | FASTA | CLI/librería | no | GPL-3.0 | 556 / 2024; 191 / 2026 | dependencia nuestra, no competidor |
| **eggNOG-mapper** | función por ortología; no genoma→GFF ni metabolismo | proteína/CDS | CLI + web | no | AGPL-3.0 | 810 / 2026-08 | mucho más profundo en función |
| **InterProScan** | dominios/familias; no genoma completo | proteína | CLI + web | no | Apache-2.0 | 397 / 2026-11 | estándar de dominios |
| **Foldseek** | búsqueda estructural; no anotación | PDB/mmCIF | CLI + web | opcional | GPL-3.0 | 1332 / 2026-11 | lo usamos como filtro |
| **ESM / ESMFold** | LM y plegado; no anotación de genoma | secuencia | Python/CLI | sí | MIT | 4175 (repo archivado) | opcional vía API; AlphaFold DB primero |
| **AlphaFold DB** | estructuras precomputadas; no anotación | UniProt | web/API | no | no verificado | servicio | recurso estructural masivo y gratuito |

## Metabolismo, regulación y plataformas

| proyecto | qué hace / qué no | entrada | GUI/CLI | GPU | licencia | ★ / última act. | diferenciador |
|---|---|---|---|---|---|---|---|
| **CarveMe** | GEM draft + gap-filling; no anotación/operones/GUI | proteína | CLI | no | NOASSERTION | 224 / 2025-09 | modelo **simulable** (SBML/FBA) |
| **gapseq** | rutas + redes + gap-filling; no GUI/operones | FASTA/GFF | CLI (R) | no | GPL-3.0 | 231 / 2026-10 | más profundo en rutas |
| **DRAM** | anotación metabólica de MAGs; no operones/GUI | FASTA/GFF | CLI | no | GPL-3.0 | 346 / 2026-11 | el más cercano en anotación + metabolismo |
| **antiSMASH** | clusters de metabolitos secundarios; no metabolismo central | FASTA/GFF | CLI + web | no | AGPL-3.0 | 272 / 2026-11 | BGCs, un gap real nuestro |
| **Operon-mapper** | operones con confianza; no anotación completa | FASTA/GFF | web | no | no verificado | servicio | referencia específica de operones |
| **RegulonDB** | red curada de *E. coli* K-12 | consultas | web | no | datos académicos | servicio | ground truth de K-12 |
| **iModulonDB** | iModulones de transcriptoma; no anotación | compendios | web | no | no verificado | servicio | regulación a nivel de módulo |
| **GhostKOALA / BlastKOALA** | KO + KEGG Mapper; no estructura | proteína | web | no | KEGG propietario | servicio | KO oficial; Kofam local sería mejor |
| **BV-BRC** | plataforma integrada: anotación, comparativa, metagenómica | FASTA/reads | web + CLI + API | no | gratis con registro | servicio | competidor “todo en uno” en nube |
| **KBase** | narratives + apps + FBA | FASTA/reads | web | no | gratis con registro | servicio | potencia de cómputo |
| **Galaxy** | workflows | cualquiera | web + CLI | no | NOASSERTION | 1901 / 2026-11 | distribución |
| **MicroScope** | anotación experta + curación | FASTA/GBK | web | no | acceso por proyecto | servicio | curación manual |

## Hallazgos

1. **La anotación base es tabla stakes.** Prokka/Bakta usan la misma familia
   (pyrodigal, tRNAscan-SE, Infernal/Rfam) y tienen más adopción. No competimos
   ahí; integramos capas posteriores.
2. **La ventaja real es la integración local/CPU con evidencia por hallazgo.**
   Ningún proyecto revisado ofrece en una instalación local: anotación +
   metabolismo/brechas + operones/promotores + regulación candidata + estructura
   + taxonomía + GUI ES/EN/ZH.
3. **Estamos por detrás en base funcional curada y versionada.** Bakta usa
   UniRef y dbxref; nosotros Swiss-Prot/Rhea. La sensibilidad en homólogos
   remotos y la reproducibilidad de la versión de DB son una brecha.
4. **Estamos por detrás en empaquetado y benchmark público.** Prokka/Bakta/DRAM
   están en Bioconda/Galaxy; nosotros tenemos install.sh + Docker, pero no una
   tabla de métricas publicada.
5. **No somos comparables a CarveMe/gapseq/ModelSEED en FBA.** Ellos producen
   SBML y simulan crecimiento; nosotros explicamos brechas. Presentarnos como
   alternativa a CarveMe sería deshonesto.
6. **Los operones necesitan un benchmark replicable.** El F1 0,78–0,81 no es
   comparable con Operon-mapper sin los mismos criterios y genes.
7. **ES/EN/ZH es un diferenciador real y barato.** Ningún proyecto del set
   anuncia interfaz en español o chino.
8. **Higiene:** EGAPx es eucariota y RAST clásico está caído; no sirven como
   comparadores bacterianos. Nuestro repo tiene 1★ y 0 descargas: hay que ganar
   credibilidad con benchmarks, no con afirmaciones.

## Recomendaciones priorizadas

| # | acción | impacto | esfuerzo |
|---|---|---|---|
| R1 | Benchmark reproducible contra Bakta, Prokka y DFAST en 5–6 genomas: CDS, tRNA/rRNA, % hypothetical, tiempo y RAM | alto | bajo-medio |
| R2 | AlphaFold DB antes de ESMFold; ESMFold solo como fallback | alto | bajo |
| R3 | Exportar borrador SBML y mapear Rhea → ModelSEED/BiGG | alto | medio |
| R4 | Receta Bioconda + wrapper Galaxy, además de Docker | alto | medio |
| R5 | Validar operones contra RegulonDB con P/R/F1 por clase e intervalos | medio | medio |

**No hacer:** reescribir el pipeline, competir con PGAP en submission, prometer
FBA o presentar EGAPx/RAST como comparadores bacterianos.

Plan de integración concreta: [`INTEGRATIONS.md`](INTEGRATIONS.md).

## Fuentes

- https://github.com/Krim10000/Bacterial-DNA-Annotator
- https://github.com/tseemann/prokka
- https://github.com/oschwengers/bakta
- https://github.com/ncbi/pgap
- https://github.com/nigyta/dfast_core
- https://github.com/eggnogdb/eggnog-mapper
- https://github.com/ebi-pf-team/interproscan
- https://github.com/steineggerlab/foldseek
- https://github.com/cdanielmachado/carveme
- https://github.com/jotech/gapseq
- https://github.com/BortonWrightonLabs/DRAM
- https://github.com/antismash/antismash
- https://biocomputo.ibt.unam.mx/operon_mapper/
- https://regulondb.ccg.unam.mx/
- https://www.bv-brc.org/
- https://www.kbase.us/
- https://galaxyproject.org/
- https://alphafold.ebi.ac.uk/
