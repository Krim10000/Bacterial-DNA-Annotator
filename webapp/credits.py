"""Créditos de las herramientas, datos y librerías que usa el proyecto.

Cada entrada declara: qué hace, quién lo hizo, dónde está, su referencia y en
qué parte de este proyecto se usa. Los datos viven aquí y no en la interfaz,
para poder revisarlos sin tocar la presentación.

Regla de honestidad: si una referencia no se pudo verificar, se escribe
"(referencia por confirmar)" en vez de inventarla.
"""
from __future__ import annotations

from i18n import tr
from repo import REPO_URL

#: Autoría del proyecto. Se muestra en la GUI, en README y en CITATION.cff.
AUTHOR = {
    "name": "Camilo Rojas Bravo",
    "email": "Camilo.Rojas@usach.cl",
    "affiliation": "Ingeniero en Biotecnología, titulado de la Universidad de Santiago de Chile (USACH); sin afiliación institucional actual",
    "role": "Diseño, implementación y validación científica",
}

#: Ayuda de desarrollo con IA. Se declara de forma explícita y verificable.
AI_HELP = {
    "tool": "DeepSeek Harness (DSH) con el modelo DeepSeek",
    "period": "noviembre de 2026",
    "scope": "revisión de código, traducción de la interfaz, documentación, "
             "pruebas de la GUI y empaquetado para publicación",
    "limit": "Las decisiones científicas, los datos y la validación de los "
             "resultados son responsabilidad del autor.",
}

GROUPS = [
    ("tools", "🧰 Herramientas de análisis",
     "Programas que calculan algo sobre el ADN o las proteínas."),
    ("data", "🗃️ Bases de datos",
     "Colecciones de conocimiento biológico que consultamos."),
    ("libs", "🧱 Librerías y plataformas",
     "Piezas de software general que sostienen la aplicación."),
    ("models", "🤖 Modelos de lenguaje y estructura",
     "Modelos entrenados con millones de secuencias."),
]

CREDITS: list[dict[str, str]] = [

    {
        "group": "data",
        "name": "DBTBS",
        "what": "Regulación curada de Bacillus subtilis: operones, promotores y "
                "sitios de unión de factores de transcripción",
        "author": "Sierro N., Makita Y., de Hoon M., Nakai K. — Human Genome "
                  "Center, Universidad de Tokio",
        "url": "https://dbtbs.hgc.jp/",
        "cite": "Sierro N. et al., DBTBS: a database of transcriptional "
                "regulation in Bacillus subtilis containing upstream "
                "intergenic conservation information. Nucleic Acids Research "
                "2008;36:D93-D96. Ver también Ishii T. et al., 2001",
        "use": "Segunda especie curada (B. subtilis) para comparar el buscador "
               "de operones fuera de E. coli. Se obtiene de sus tablas públicas "
               "porque el archivo plano del formulario devuelve 404.",
        "license": "Acceso libre con registro (afiliación y propósito de uso). "
                   "No se redistribuye: lo descarga fetch_dbtbs_tables.py.",
    },
    # ---------------------------------------------------------- herramientas
    dict(
        group="tools", name="pyrodigal",
        what="Encuentra los genes que producen proteína (CDS) en un genoma.",
        author="Martin Larralde; método original de Prodigal (Oak Ridge National Laboratory)",
        url="https://github.com/althonos/pyrodigal",
        cite="Prodigal: Hyatt D. et al., BMC Bioinformatics 2010;11:119",
        use="Anotación de CDS en los genomas subidos y en los 400 del proyecto.",
        license="GPL-3.0",
    ),
    dict(
        group="tools", name="DIAMOND",
        what="Compara proteínas muy rápido contra una base de proteínas conocidas.",
        author="Benjamin Buchfink, Chao Xie, Daniel H. Huson",
        url="https://github.com/bbuchfink/diamond",
        cite="Buchfink B., Xie C., Huson D.H., Nature Methods 2015;12:59-60",
        use="Anotación funcional contra Swiss-Prot, búsqueda de homólogos de "
            "enzimas y agrupación all-vs-all de factores sigma por ortólogo.",
        license="GPL-3.0",
    ),
    dict(
        group="tools", name="HMMER",
        what="Busca familias de proteínas comparando su perfil, no la secuencia exacta.",
        author="Sean R. Eddy (HHMI / Harvard)",
        url="http://hmmer.org/",
        cite="Eddy S.R., PLoS Computational Biology 2011;7:e1002195",
        use="Detección de dominios Pfam en cierres de brechas y en el inventario "
            "de factores sigma (Sigma70_r2/r4, Sigma54, etc.).",
        license="BSD-3-Clause",
    ),
    dict(
        group="tools", name="Infernal",
        what="Busca ARN por su estructura, no solo por su secuencia.",
        author="Eric P. Nawrocki, Sean R. Eddy",
        url="http://eddylab.org/infernal/",
        cite="Nawrocki E.P., Eddy S.R., Bioinformatics 2013;29:2933-2935",
        use="Detección de ARN ribosómico contra Rfam.",
        license="BSD-3-Clause",
    ),
    dict(
        group="tools", name="tRNAscan-SE",
        what="Localiza los genes de ARN de transferencia y su anticodón.",
        author="Patricia P. Chan, Todd M. Lowe (UC Santa Cruz)",
        url="http://lowelab.ucsc.edu/tRNAscan-SE/",
        cite="Chan P.P. et al., Nucleic Acids Research 2021;49:9077-9096",
        use="Anotación de tRNA en genomas subidos y en los 400.",
        license="GPL-3.0",
    ),
    dict(
        group="tools", name="barrnap",
        what="Encuentra los genes de ARN ribosómico de forma rápida.",
        author="Torsten Seemann",
        url="https://github.com/tseemann/barrnap",
        cite="(sin publicación; citar el repositorio)",
        use="Prefiltro de rRNA antes de la búsqueda con Infernal.",
        license="GPL-3.0",
    ),
    dict(
        group="tools", name="ViennaRNA",
        what="Calcula la estructura secundaria del ARN y su energía libre (ΔG).",
        author="Ronny Lorenz e Ivo Hofacker (Universidad de Viena)",
        url="https://www.tbi.univie.ac.at/RNA/",
        cite="Lorenz R. et al., Algorithms for Molecular Biology 2011;6:26",
        use="Energía de las horquillas de terminación.",
        license="Licencia propia de ViennaRNA (uso libre con condiciones: ver el sitio)",
    ),
    dict(
        group="tools", name="ORCA (orcapy)",
        what="Predice el origen de replicación de cromosomas bacterianos circulares.",
        author="Zoya van Meel, Jasmijn A. Baaijens",
        url="https://github.com/ZoyavanMeel/ORCA",
        cite="van Meel Z., Baaijens J.A., ORCA: predicting replication origins in circular prokaryotic chromosomes. bioRxiv (preprint, 2024), doi:10.1101/2024.03.28.587133",
        use="Predicción de oriC en los 400 genomas.",
        license="Consultar en el repositorio",
    ),
    dict(
        group="tools", name="fastANI",
        what="Mide el parecido global entre dos genomas (ANI).",
        author="Chirag Jain y colaboradores",
        url="https://github.com/ParBLiSS/FastANI",
        cite="Jain C. et al., Nature Communications 2018;9:5114",
        use="Taxonomía del piloto: comparación contra genomas de referencia.",
        license="Consultar en el repositorio",
    ),

    # ---------------------------------------------------------------- datos
    dict(
        group="data", name="Rhea",
        what="Base de datos de reacciones químicas enzimáticas, con sus direcciones.",
        author="SIB Swiss Institute of Bioinformatics y EMBL-EBI",
        url="https://www.rhea-db.org/",
        cite="Bansal P. et al., Nucleic Acids Research 2022;50:D693-D700",
        use="Reacciones, ecuaciones y dirección de cada paso metabólico.",
        license="Abierta (consultar condiciones en el sitio)",
    ),
    dict(
        group="data", name="ChEBI",
        what="Identificadores únicos para moléculas químicas.",
        author="EMBL-EBI",
        url="https://www.ebi.ac.uk/chebi/",
        cite="Hastings J. et al., Nucleic Acids Research 2016;44:D1214-D1219",
        use="Nombres e identificadores de las moléculas de cada brecha.",
        license="Abierta (CC BY 4.0)",
    ),
    dict(
        group="data", name="RegulonDB",
        what="Regulación transcripcional de E. coli K-12: operones, promotores y terminadores curados.",
        author="Centro de Ciencias Genómicas, UNAM (Heladia Salgado y colaboradores)",
        url="https://regulondb.ccg.unam.mx/",
        cite="Salgado H. et al., Nucleic Acids Research 2024;52:D255-D264",
        use="Verdad de referencia de operones, promotores, TSS y terminadores de "
            "K-12; PWMs de σ70/σ54/σ32/σ38/σ28/σ24 usadas por el módulo sigma.",
        license="Académica / no comercial; los datos NO se redistribuyen en el repositorio",
    ),
    dict(
        group="data", name="NCBI RefSeq / GenBank",
        what="Genomas, anotaciones y proteínas de referencia de las 400 especies.",
        author="National Center for Biotechnology Information (NCBI)",
        url="https://www.ncbi.nlm.nih.gov/",
        cite="(datos; citar el assembly/accession usado y las políticas NCBI)",
        use="Genomas y anotaciones de los 400; fallback GBFF para recuperar "
            "secuencia y evidencia del módulo sigma.",
        license="Datos públicos; respetar las políticas de uso de NCBI",
    ),
    dict(
        group="data", name="STRING",
        what="Puntuaciones de asociación funcional entre proteínas.",
        author="Lars J. Jensen y colaboradores (Universidad de Copenhague)",
        url="https://string-db.org/",
        cite="Szklarczyk D. et al., Nucleic Acids Research 2023;51:D638-D646",
        use="Evidencia de acoplamiento funcional en las reglas de operones.",
        license="Abierta para uso académico (consultar condiciones)",
    ),
    dict(
        group="data", name="DoriC",
        what="Base de datos de orígenes de replicación en procariotas.",
        author="Feng Gao y colaboradores (Universidad de Tianjin)",
        url="https://tubic.org/doric/",
        cite="Dong M.-J., Luo H., Gao F., Nucleic Acids Research 2022;50:D117-D120",
        use="Verdad de referencia del origen de replicación (oriC).",
        license="Académica; no se redistribuye en el repositorio",
    ),
    dict(
        group="data", name="Gene Ontology (GO)",
        what="Ontología jerárquica de funciones, procesos y componentes celulares.",
        author="Gene Ontology Consortium",
        url="https://geneontology.org/",
        cite="Ashburner M. et al., Nature Genetics 2000;25:25-29; "
             "The Gene Ontology Consortium, Nucleic Acids Research "
             "2021;49:D325-D334. Ontología local: go-basic.obo, "
             "data-version releases/2026-07-26.",
        use="Etiquetas de proceso y propagación jerárquica para genes/operones "
            "y como prior contextual de sigma.",
        license="CC BY 4.0",
    ),
    dict(
        group="data", name="KEGG",
        what="Enciclopedia de rutas metabólicas con nombre y mapa.",
        author="Minoru Kanehisa y Susumu Goto (Universidad de Kioto)",
        url="https://www.kegg.jp/",
        cite="Kanehisa M., Goto S., Nucleic Acids Research 2000;28:27-30",
        use="Nombres de rutas y mapas de reacciones.",
        license="Académica / no comercial",
    ),
    dict(
        group="data", name="UniProt / Swiss-Prot",
        what="Proteínas con función revisada por expertos.",
        author="UniProt Consortium",
        url="https://www.uniprot.org/",
        cite="The UniProt Consortium, Nucleic Acids Research 2023;51:D523-D531",
        use="Anotación funcional de proteínas y homologías.",
        license="CC BY 4.0",
    ),
    dict(
        group="data", name="Pfam",
        what="Familias de proteínas representadas como perfiles.",
        author="Pfam Consortium (EMBL-EBI / Sanger)",
        url="https://www.ebi.ac.uk/interpro/",
        cite="Mistry J. et al., Nucleic Acids Research 2021;49:D412-D419",
        use="Dominios proteicos para cierres de brechas y para clasificar "
            "factores sigma por familia.",
        license="Abierta (consultar InterPro)",
    ),
    dict(
        group="data", name="Rfam",
        what="Familias de ARN con su estructura consenso.",
        author="Rfam Consortium (EMBL-EBI)",
        url="https://rfam.org/",
        cite="Kalvari I. et al., Nucleic Acids Research 2021;49:D192-D200",
        use="Detección de rRNA y otros ARN no codificantes.",
        license="Abierta (consultar el sitio)",
    ),
    dict(
        group="data", name="eQuilibrator",
        what="Estima la energía libre (ΔG′) de reacciones bioquímicas.",
        author="Elad Noor, Ron Milo y colaboradores (Instituto Weizmann)",
        url="https://equilibrator.weizmann.ac.il/",
        cite="Beber M.E. et al., Nucleic Acids Research 2022;50:D603-D609",
        use="Dirección y termodinámica de los cierres de brechas.",
        license="Código abierto (consultar el proyecto)",
    ),
    dict(
        group="data", name="TECRDB",
        what="Datos termodinámicos medidos de reacciones enzimáticas.",
        author="Robert N. Goldberg, Yadu B. Tewari, Talapady N. Bhat (NIST)",
        url="https://trc.nist.gov/TECRDB/",
        cite="Goldberg R.N., Tewari Y.B., Bhat T.N., Bioinformatics 2004;20:2874-2877",
        use="Validación de los valores de energía libre (ΔG′).",
        license="Académica (NIST)",
    ),

    # ------------------------------------------------------------ librerías
    dict(
        group="libs", name="Streamlit",
        what="Convierte código de Python en una página web interactiva.",
        author="Streamlit Inc. (hoy parte de Snowflake)",
        url="https://streamlit.io/",
        cite="(software; sin publicación)",
        use="Toda la interfaz web, incluido el módulo de sigmas y promotores.",
        license="Apache-2.0",
    ),
    dict(
        group="libs", name="pandas",
        what="Tablas de datos para Python.",
        author="pandas Development Team",
        url="https://pandas.pydata.org/",
        cite="(software; sin publicación)",
        use="Lectura y filtrado de todas las tablas; tabla de evidencia del módulo sigma.",
        license="BSD-3-Clause",
    ),
    dict(
        group="libs", name="NumPy",
        what="Cálculo numérico y matrices.",
        author="NumPy Developers",
        url="https://numpy.org/",
        cite="(software; sin publicación)",
        use="Puntaje de PWMs, arrays de secuencia y fondo nulo del módulo sigma.",
        license="BSD-3-Clause",
    ),
    dict(
        group="libs", name="Biopython",
        what="Lee FASTA, GenBank y estructuras; permite recuperar secuencia desde GBFF.",
        author="Biopython Contributors (Peter J. A. Cock y colaboradores)",
        url="https://biopython.org/",
        cite="Cock P.J.A. et al., Bioinformatics 2009;25:1422-1423",
        use="Lectura de GBFF/GenBank y FASTA en el módulo sigma; fallback de "
            "secuencia cuando genes_meta.json no trae seq_source.",
        license="Biopython License (tipo BSD)",
    ),
    dict(
        group="libs", name="scikit-learn",
        what="Aprendizaje automático clásico: clasificación, validación, métricas.",
        author="scikit-learn Developers",
        url="https://scikit-learn.org/",
        cite="Pedregosa F. et al., Journal of Machine Learning Research 2011;12:2825-2830",
        use="Clasificadores de evaluación y comparación de métodos.",
        license="BSD-3-Clause",
    ),
    dict(
        group="libs", name="NetworkX",
        what="Grafos y redes.",
        author="NetworkX Developers",
        url="https://networkx.org/",
        cite="(software; sin publicación)",
        use="Redes metabólicas y de regulación.",
        license="BSD-3-Clause",
    ),
    dict(
        group="libs", name="SciPy",
        what="Funciones científicas: estadística, optimización, álgebra.",
        author="SciPy Developers",
        url="https://scipy.org/",
        cite="(software; sin publicación)",
        use="Apoyo estadístico y numérico.",
        license="BSD-3-Clause",
    ),

    # -------------------------------------------------------------- modelos
    dict(
        group="models", name="ESM-2 / ESMFold",
        what="Modelo de lenguaje entrenado con millones de proteínas: dice si una secuencia «parece» una proteína real.",
        author="Meta AI (Fundamental AI Research)",
        url="https://github.com/facebookresearch/esm",
        cite="Lin Z. et al., Science 2023;379:1123-1130",
        use="Puntuación de ORFs candidatos en genomas subidos.",
        license="Código MIT; pesos ESM-2 MIT (verificar antes de redistribuir); pesos ESMFold con condiciones no comerciales",
    ),
    dict(
        group="tools", name="Foldseek",
        what="Busca proteínas por parecido estructural 3D contra PDB/AFDB.",
        author="van Kempen M. et al.; Steinegger Lab",
        url="https://github.com/steineggerlab/foldseek",
        cite="van Kempen M. et al., Nature Biotechnology 2024;42:243-246",
        use="Filtro F1: homología estructural de candidatos.",
        license="GPL-3.0",
    ),
    dict(
        group="tools", name="EnzyMM / PyJess",
        what="Busca arreglos geométricos de residuos catalíticos con plantillas M-CSA.",
        author="Raymund Hackett; Zeller Group (LUMC); PyJess de Martin Larralde",
        url="https://github.com/RayHackett/enzymm",
        cite="Hackett R. et al., Protein Science 2026, doi:10.1002/pro.70744",
        use="Filtro F2: anotación estructural automática contra M-CSA.",
        license="MIT",
    ),
    dict(
        group="tools", name="P2Rank",
        what="Predice cavidades y posibles sitios de unión a ligando en una estructura.",
        author="Radoslav Krivák, David Hoksza",
        url="https://github.com/rdk/p2rank",
        cite="Krivák R., Hoksza D., BMC Bioinformatics 2018;19:314",
        use="Filtro F3: bolsillo y co-localización con residuos catalíticos.",
        license="MIT",
    ),
    dict(
        group="data", name="AlphaFill",
        what="Transplanta ligandos y cofactores desde estructuras PDB homólogas a modelos AFDB.",
        author="Hekkelman M. et al. (NKI/EMBL-EBI)",
        url="https://alphafill.eu/",
        cite="Hekkelman M. et al., Nature Methods 2023;20:205-213",
        use="Filtro F4: ligandos/cofactores compatibles.",
        license="Abierta para uso académico (consultar el sitio)",
    ),
    dict(
        group="data", name="M-CSA",
        what="Atlas curado de mecanismos y residuos catalíticos de enzimas.",
        author="Thornton Group, EMBL-EBI",
        url="https://www.ebi.ac.uk/thornton-srv/m-csa/",
        cite="Ribeiro A.J.M. et al., Nucleic Acids Research 2018;46:D618-D623",
        use="Plantillas de sitios catalíticos para EnzyMM.",
        license="Abierta para uso académico",
    ),
    dict(
        group="data", name="AlphaFold DB",
        what="Estructuras predichas para proteínas UniProt.",
        author="DeepMind y EMBL-EBI",
        url="https://alphafold.ebi.ac.uk/",
        cite="Varadi M. et al., Nucleic Acids Research 2024;52:D368-D375",
        use="Estructuras de homólogos para precribado estructural.",
        license="CC BY 4.0 (consultar términos)",
    ),
    dict(
        group="data", name="PDB",
        what="Estructuras experimentales de proteínas y complejos.",
        author="wwPDB Consortium",
        url="https://www.rcsb.org/",
        cite="wwPDB Consortium, Nucleic Acids Research 2019;47:D464-D474",
        use="Homólogos holo con ligandos/cofactores para Foldseek.",
        license="Abierta (CC0)",
    ),
]

LICENSES_NOTE = [
    ("Uso libre", "Streamlit, pandas, NumPy, scikit-learn, NetworkX, SciPy, HMMER, "
                  "Infernal, pyrodigal, DIAMOND, tRNAscan-SE, barrnap, Biopython."),
    ("Uso abierto con condiciones", "ViennaRNA, Rhea, ChEBI, STRING, UniProt, Pfam, "
                                    "Rfam, eQuilibrator, Gene Ontology (CC BY 4.0)."),
    ("Uso académico o no comercial", "RegulonDB, KEGG, DoriC, TECRDB y los pesos de "
                                     "ESMFold."),
    ("No se redistribuye en el repositorio", "RegulonDB, Pfam, STRING, DoriC y "
                                             "DBTBS: se "
                                             "descargan o se generan localmente y se "
                                             "guardan fuera del control de versiones."),]


def render() -> None:
    """Página de créditos. Importa Streamlit aquí para que el módulo de datos
    pueda leerse (o probarse) sin depender de la interfaz."""
    import streamlit as st

    st.title(tr("🙏 Créditos y herramientas"))
    st.caption(
        "Este proyecto no sería posible sin el trabajo de otras personas: "
        "programas, bases de datos y modelos que llevan años construyéndose. "
        "Aquí está cada pieza, quién la hizo y en qué parte de este proyecto se usa. "
        "Los nombres propios y las citas se conservan en su idioma original."
    )

    st.subheader(tr("👤 Autoría"))
    st.markdown(
        f"**{AUTHOR['name']}** — {tr(AUTHOR['role'])}  \n"
        f"🏛️ {tr(AUTHOR['affiliation'])}  \n"
        f"✉️ [{AUTHOR['email']}](mailto:{AUTHOR['email']})"
    )
    st.markdown(
        f"**{tr('Asistencia de desarrollo:')}** {AI_HELP['tool']}, "
        f"{tr(AI_HELP['period'])}.  \n"
        f"{tr('Se usó para')} {tr(AI_HELP['scope'])}.  \n"
        f"{tr(AI_HELP['limit'])}"
    )
    st.caption(
        tr("Si usas este proyecto en un trabajo académico, cita el repositorio "
           "y las herramientas de las que dependen tus resultados. "
           "El archivo CITATION.cff tiene la referencia lista para usar.")
    )

    for key, title, subtitle in GROUPS:
        items = [c for c in CREDITS if c["group"] == key]
        if not items:
            continue
        st.subheader(tr(title))
        st.caption(tr(subtitle))
        for c in items:
            with st.container(border=True):
                st.markdown(f"**{c['name']}** — {tr(c['what'])}")
                st.markdown(f"👤 {c['author']}")
                st.markdown(f"🔗 [{c['url']}]({c['url']})")
                st.markdown(f"📄 {c['cite']}")
                st.markdown(f"🧭 **{tr('En este proyecto:')}** {tr(c['use'])}")
                st.caption(f"⚖️ {tr('Licencia:')} {tr(c['license'])}")

    st.subheader("⚖️ Nota sobre licencias")
    st.caption(
        "No todo lo que usa el proyecto se puede redistribuir. Los datos de uso "
        "académico se descargan en tu equipo y quedan fuera del control de "
        "versiones."
    )
    for title, detail in LICENSES_NOTE:
        st.markdown(f"- **{tr(title)}:** {tr(detail)}")

    st.info(
        "Si usas este proyecto en un trabajo académico, cita también las "
        "herramientas y bases de datos de las que dependen tus resultados, no "
        "solo este repositorio."
    )
    st.markdown(
        f"⭐ **{tr('¿Te sirvió?')}** {tr('Dale una estrella en GitHub:')} "
        f"[{REPO_URL}]({REPO_URL})"
    )


# --------------------------------------------------------------------------
# Líneas cortas de acreditación para usar bajo cada resultado.
# Formato: "🧰 Calculado con **herramienta** (autor, año) — enlace"
# --------------------------------------------------------------------------
SHORT: dict[str, str] = {
    "pyrodigal": "🧰 *Genes calculados con **pyrodigal** (Larralde; método Prodigal, Hyatt et al., 2010) — [github.com/althonos/pyrodigal](https://github.com/althonos/pyrodigal)*",
    "diamond": "🧰 *Homologías calculadas con **DIAMOND** (Buchfink et al., 2015) — [github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)*",
    "hmmer": "🧰 *Dominios calculados con **HMMER** (Eddy, 2011) sobre **Pfam** (Mistry et al., 2021) — [hmmer.org](http://hmmer.org/)*",
    "infernal": "🧰 *ARN buscado con **Infernal** (Nawrocki & Eddy, 2013) contra **Rfam** (Kalvari et al., 2021) — [eddylab.org/infernal](http://eddylab.org/infernal/)*",
    "trnascan": "🧰 *tRNA anotados con **tRNAscan-SE** (Chan et al., 2021) — [lowelab.ucsc.edu/tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)*",
    "barrnap": "🧰 *rRNA prefiltrado con **barrnap** (Seemann) — [github.com/tseemann/barrnap](https://github.com/tseemann/barrnap)*",
    "vienna": "🧰 *Energía calculada con **ViennaRNA** (Lorenz et al., 2011) — [tbi.univie.ac.at/RNA](https://www.tbi.univie.ac.at/RNA/)*",
    "fastani": "🧰 *Parecido entre genomas calculado con **fastANI** (Jain et al., 2018) — [github.com/ParBLiSS/FastANI](https://github.com/ParBLiSS/FastANI)*",
    "orca": "🧰 *Origen predicho con **ORCA** (van Meel & Baaijens, 2024) — [github.com/ZoyavanMeel/ORCA](https://github.com/ZoyavanMeel/ORCA)*",
    "doric": "🗃️ *Verdad de referencia: **DoriC** (Dong et al., 2022) — [tubic.org/doric](https://tubic.org/doric/)*",
    "regulondb": "🗃️ *Datos curados de **RegulonDB** (Salgado et al., 2024) — [regulondb.ccg.unam.mx](https://regulondb.ccg.unam.mx/)*",
    "rhea": "🗃️ *Reacciones de **Rhea** (Bansal et al., 2022) — [rhea-db.org](https://www.rhea-db.org/)*",
    "chebi": "🗃️ *Moléculas identificadas con **ChEBI** (Hastings et al., 2016) — [ebi.ac.uk/chebi](https://www.ebi.ac.uk/chebi/)*",
    "string": "🗃️ *Acoplamiento funcional de **STRING** (Szklarczyk et al., 2023) — [string-db.org](https://string-db.org/)*",
    "kegg": "🗃️ *Nombres de rutas de **KEGG** (Kanehisa & Goto, 2000) — [kegg.jp](https://www.kegg.jp/)*",
    "foldseek": "🧰 *Homología estructural con **Foldseek** (van Kempen et al., 2024) — [github.com/steineggerlab/foldseek](https://github.com/steineggerlab/foldseek)*",
    "enzymm": "🧰 *Plantillas catalíticas con **EnzyMM/PyJess** y **M-CSA** (Hackett et al., 2026) — [enzymm.readthedocs.io](https://enzymm.readthedocs.io)*",
    "p2rank": "🧰 *Bolsillos predichos con **P2Rank** (Krivák & Hoksza, 2018) — [github.com/rdk/p2rank](https://github.com/rdk/p2rank)*",
    "alphafill": "🗃️ *Ligandos/cofactores de **AlphaFill** (Hekkelman et al., 2023) — [alphafill.eu](https://alphafill.eu)*",
    "mcsa": "🗃️ *Plantillas de **M-CSA** (Ribeiro et al., 2018) — [ebi.ac.uk/thornton-srv/m-csa](https://www.ebi.ac.uk/thornton-srv/m-csa/)*",
    "afdb": "🗃️ *Estructuras de **AlphaFold DB** (Varadi et al., 2024) — [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/)*",
    "biopython": "🧱 *Lectura de GBFF/GenBank con **Biopython** (Cock et al., 2009) — [biopython.org](https://biopython.org/)*",
    "ncbi": "🗃️ *Genomas y anotaciones de **NCBI RefSeq/GenBank** — [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/)*",
    "sigma": "🧪 *Factores sigma y promotores: **RegulonDB** (PWMs), **Pfam/HMMER**, **DIAMOND** y **NumPy***",
    "uniprot": "🗃️ *Proteínas de **UniProt/Swiss-Prot** (UniProt Consortium, 2023) — [uniprot.org](https://www.uniprot.org/)*",
    "equilibrator": "🗃️ *Energía estimada con **eQuilibrator** (Beber et al., 2022) y **TECRDB** (Goldberg et al., 2004)*",
    "esm": "🤖 *Puntuación con **ESM-2** (Lin et al., 2023) — [github.com/facebookresearch/esm](https://github.com/facebookresearch/esm)*",
    "streamlit": "🧱 *Interfaz construida con **Streamlit** — [streamlit.io](https://streamlit.io/)*",
}


def line(*names: str) -> str:
    """Une varias acreditaciones, omitiendo las que no existan."""
    return " · ".join(SHORT[n] for n in names if n in SHORT)
