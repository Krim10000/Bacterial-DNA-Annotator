# Créditos del pipeline de sigma

> Para la sección **🙏 Créditos y herramientas**. Lista qué debe acreditarse,
> qué ya existe en `webapp/credits.py` y qué hay que añadir o actualizar.
>
> Regla: no acreditar herramientas que el módulo **no** usa. El pipeline de
> sigma no usa MEME/FIMO, AlphaFold, Rosetta, Vina ni DM.

---

## 1. Créditos directos del módulo

| Herramienta / dato | Por qué | Licencia | Acción |
|---|---|---|---|
| **NumPy** | Puntaje PWM, arrays y fondo nulo | BSD-3-Clause | Ya existe; ampliar `use` |
| **Biopython** | Lectura de GBFF/GenBank cuando falta el FASTA | Biopython License (tipo BSD) | **Añadir** |
| **DIAMOND** | Agrupación de ortólogos sigma (all-vs-all) | GPL-3.0 | Ya existe; ampliar `use` |
| **HMMER** | Dominios Pfam del inventario sigma | BSD-3-Clause | Ya existe; ampliar `use` |
| **Pfam** | Dominios `Sigma70_r2/r4`, `Sigma54`, etc. | CC0 / InterPro | Ya existe; ampliar `use` |
| **RegulonDB** | PWMs de sigma, promotores y TSS curados de K-12 | Académica / no comercial | Ya existe; ampliar `use` |
| **NCBI RefSeq/GenBank** | Genomas y anotaciones; fallback GBFF | Políticas NCBI | **Añadir** |
| **UniProt / Swiss-Prot** | Anotación funcional usada por el contexto del operón | CC BY 4.0 | Ya existe; ampliar `use` |
| **Streamlit** | Interfaz del módulo | Apache-2.0 | Ya existe; ampliar `use` |
| **pandas** | Tabla de evidencia en la GUI | BSD-3-Clause | Ya existe; ampliar `use` |
| **Python** | Lenguaje del módulo | PSF | Ya existe en créditos generales |

La agrupación usa **DIAMOND**; no se usa MMseqs2 ni `diamond cluster` como
dependencia separada. El código ejecuta `diamond makedb` + `diamond blastp`.

---

## 2. Créditos heredados del pipeline

El módulo lee datos ya calculados por el proyecto. Estos créditos ya deberían
estar en la sección, pero conviene que el `use` mencione sigma:

| Pieza | Uso en sigma |
|---|---|
| **pyrodigal** | `proteins.faa` de los genomas |
| **HMMER + Pfam** | Inventario de dominios sigma |
| **DIAMOND + Swiss-Prot** | Anotación funcional y homólogos |
| **RegulonDB** | PWMs de promotores sigma |
| **NCBI** | Genomas, GFF/GBFF, proteínas |
| **UniProt** | Producto/función de genes del operón |

---

## 3. Entradas listas para `webapp/credits.py`

> No duplicar las entradas existentes; actualizar `use` o añadir solo las que
> falten. Los textos están en español, como el resto de la sección.

### 3.1 Nuevas

```python
dict(
    group="libs", name="Biopython",
    what="Lee FASTA, GenBank y estructuras. Es la pieza que permite recuperar "
         "la secuencia desde un GBFF cuando falta el FASTA.",
    author="Biopython Contributors (Peter J. A. Cock y colaboradores)",
    url="https://biopython.org/",
    cite="Cock P.J.A. et al., Biopython: freely available Python tools for "
         "computational molecular biology and bioinformatics. Bioinformatics "
         "2009;25:1422-1423. doi:10.1093/bioinformatics/btp163",
    use="Lectura de GBFF/GenBank y FASTA en el módulo sigma; fallback de "
        "secuencia cuando `genes_meta.json` no trae `seq_source`.",
    license="Biopython License (tipo BSD)",
),

dict(
    group="data", name="NCBI RefSeq / GenBank",
    what="Genomas, anotaciones y proteínas de referencia de las 400 especies.",
    author="National Center for Biotechnology Information (NCBI)",
    url="https://www.ncbi.nlm.nih.gov/",
    cite="(datos; citar el assembly/accession usado y las políticas NCBI)",
    use="Genomas y anotaciones de los 400; fallback GBFF para recuperar "
        "secuencia. Los accessions se conservan en `data/annotations/`.",
    license="Datos públicos; respetar las políticas de uso de NCBI",
),
```

### 3.2 Actualizar `use` de las existentes

**DIAMOND**

```python
use="Anotación funcional contra Swiss-Prot, búsqueda de homólogos y "
    "agrupación all-vs-all de factores sigma por ortólogo.",
```

**HMMER**

```python
use="Detección de dominios Pfam en cierres de brechas y en el inventario "
    "de factores sigma (Sigma70_r2/r4, Sigma54, etc.).",
```

**Pfam**

```python
use="Dominios proteicos para cierres de brechas y para clasificar factores "
    "sigma por familia.",
```

**RegulonDB**

```python
use="Verdad de referencia de operones, promotores, TSS y terminadores de K-12; "
    "PWMs de σ70/σ54/σ32/σ38/σ28/σ24 usadas por el módulo sigma.",
```

**NumPy**

```python
use="Puntaje de PWMs, arrays de secuencia y fondo nulo del módulo sigma.",
```

**Streamlit / pandas**

```python
use="Interfaz y tabla de evidencia del módulo sigma; descarga de TSV/JSON.",
```

---

## 4. Referencias científicas de sigma

Si la sección tiene un bloque de referencias, añadir estas. Si no, van en
[`SIGMA_MODULE.md`](SIGMA_MODULE.md).

| Tema | Referencia |
|---|---|
| Sigma factors, revisión general | Gruber T.M., Gross C.A. Annu Rev Microbiol. 2003;57:441-466. doi:10.1146/annurev.micro.57.030502.090913 |
| Sigma factors, estructura y genómica | Feklístov A., Sharon B.D., Darst S.A., Gross C.A. Annu Rev Microbiol. 2014;68:357-376. doi:10.1146/annurev-micro-092412-155737 |
| Sigma/anti-sigma | Paget M.S. Biomolecules. 2015;5:1245-1265. doi:10.3390/biom5031245 |
| ECF sigma factors | Staroń A. et al. Mol Microbiol. 2009;74:557-581. doi:10.1111/j.1365-2958.2009.06870.x |
| σ54 y bEBP | Bush M., Dixon R. Microbiol Mol Biol Rev. 2012;76:497-529. doi:10.1128/MMBR.00006-12 |
| RegulonDB | Salgado H. et al. Nucleic Acids Res. 2024;52:D255-D264. doi:10.1093/nar/gkad1072 |
| DIAMOND | Buchfink B., Xie C., Huson D.H. Nat Methods. 2015;12:59-60. doi:10.1038/nmeth.3176 |
| DIAMOND a escala | Buchfink B., Reuter K., Drost H.G. Nat Methods. 2021;18:366-368. doi:10.1038/s41592-021-01101-x |
| HMMER | Eddy S.R. PLoS Comput Biol. 2011;7:e1002195. doi:10.1371/journal.pcbi.1002195 |
| Pfam | Mistry J. et al. Nucleic Acids Res. 2021;49:D412-D419. doi:10.1093/nar/gkaa913 |
| Biopython | Cock P.J.A. et al. Bioinformatics. 2009;25:1422-1423. doi:10.1093/bioinformatics/btp163 |
| NumPy | Harris C.R. et al. Nature. 2020;585:357-362. doi:10.1038/s41586-020-2649-2 |

---

## 5. Lo que NO debe acreditarse para este módulo

- **MEME Suite / FIMO**: no se usan; el escaneo es propio con NumPy.
- **AlphaFold / RoseTTAFold / Boltz**: no se usan en sigma.
- **AutoDock Vina / fpocket / Rosetta / GROMACS**: no se usan en sigma.
- **MMseqs2**: no se usa; la agrupación es DIAMOND.
- **KEGG**: no se usa directamente; el contexto es por producto/EC del gen.
- **ChEBI**: no se usa en sigma (sí en efectores/alostería).

---

## 6. Frase para la sección

> **Anotador de sigma:** PWMs y promotores curados de **RegulonDB**; dominios
> **Pfam** detectados con **HMMER**; agrupación de ortólogos con **DIAMOND**;
> secuencias de **NCBI**; anotación funcional de **UniProt/Swiss-Prot**; cálculo
> numérico con **NumPy**; lectura de GBFF con **Biopython**; interfaz con
> **Streamlit** y **pandas**.

---

## 7. Estado de integración (2026-11-25)

Integrado en el repositorio:

- `modules/sigma/module.py` — autodescubierto por `python -m modules --list`;
- `scripts/sigma/10_annotate_sigma.py` y `20_group_sigma.py`;
- `webapp/module_ui.py` — preview `outputs["preview"]`;
- `webapp/app.py` — atajo `🧪 Anotador de sigmas`;
- `run_pipeline.sh` — etapa `sigma`;
- `docs/MODULES.md` — fila `sigma`.

Créditos actualizados en `webapp/credits.py`:

- añadidos **Biopython** y **NCBI RefSeq/GenBank**;
- ampliado el uso de **DIAMOND**, **HMMER**, **Pfam**, **RegulonDB**, **NumPy**, **Streamlit** y **pandas**;
- añadidos `SHORT` para Biopython, NCBI y sigma;
- actualizada la nota de licencias.

Verificación:

```bash
python -m modules --list | grep sigma
python -m modules --selftest sigma
python scripts/sigma/10_annotate_sigma.py --selftest
python -m pytest tests/test_sigma_annotator.py tests/test_modules.py -q
```
