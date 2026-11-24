"""Requisitos, tiempos estimados y estado de la GPU.

Todo lo que un usuario nuevo necesita saber **antes** de pulsar Anotar:
qué hardware hace falta, cuánto tarda cada paso, qué necesita GPU y qué
significa cada herramienta. Está separado de ``app.py`` para que los textos se
puedan revisar y traducir sin tocar la interfaz.

Los tiempos son rangos observados en un portátil de 8 núcleos sin GPU y en la
documentación del proyecto; se muestran siempre como aproximados.
"""
from __future__ import annotations

import shutil
import subprocess

from i18n import tr

#: (paso, qué hace en lenguaje llano, tiempo aproximado, ¿GPU?, requisito)
STEPS = [
    ("1. Anotar un genoma", "Encuentra genes, tRNA y rRNA en un FASTA.",
     "1–3 min por genoma de ~5 Mb", "No",
     "pyrodigal, tRNAscan-SE y barrnap (o solo pyrodigal)"),
    ("2. Buscar ORFs no anotados", "Mira en los huecos entre genes conocidos.",
     "~5 s por genoma; ~6 h los 400", "No", "Python puro"),
    ("3. Comparar genomas subidos", "Agrupa proteínas parecidas al 95 %.",
     "1–5 min por pareja", "No", "DIAMOND"),
    ("4. Red metabólica y brechas", "Conecta enzimas con reacciones y busca huecos.",
     "~2 min por genoma", "No", "DIAMOND + Swiss-Prot + Rhea (~120 MB)"),
    ("5. Operones por reglas", "Decide fronteras de operones sin modelo.",
     "~8 s por genoma", "No", "GFF + TSS/TTS curados"),
    ("6. Horquillas de terminación", "Pliega el ARN y mide su estabilidad.",
     "~1–10 min por genoma", "No", "ViennaRNA"),
    ("7. Estructura (Foldseek/AlphaFill)", "Busca parecidos 3D y ligandos.",
     "~1–5 min por proteína", "No", "Foldseek local + APIs"),
    ("8. ESMFold (opcional)", "Predice la estructura 3D del ORF exacto.",
     "~1–2 min por proteína", "Sí (o API)",
     "GPU ≥16 GB, o Neurosnap (~0,04 créditos/proteína)"),
]

SYSTEM = [
    ("Sistema operativo", "Ubuntu 22.04 o 24.04 (también funciona en macOS y WSL)."),
    ("Python", "3.10, 3.11 o 3.12."),
    ("Memoria RAM", "8 GB mínimo; 16 GB recomendado para el pipeline completo."),
    ("Disco", "10 GB libres para la app; ~25 GB si descargas todas las bases."),
    ("CPU", "4 núcleos; 8 recomendado."),
    ("GPU", "No hace falta para anotar, buscar brechas ni usar la GUI. "
            "Solo ESMFold local (opcional) la necesita."),
    ("Navegador", "Firefox o Chrome actualizados."),
    ("Red", "Solo para descargar bases de datos y herramientas la primera vez."),
]

GLOSSARY = [
    ("FASTA", "Archivo de texto con la secuencia de ADN: líneas que empiezan por "
              "«>» con el nombre y, debajo, las letras A, T, G, C."),
    ("CDS", "Tramo de ADN que se traduce a una proteína. Es lo que busca "
            "pyrodigal."),
    ("DIAMOND", "Buscador ultrarrápido de parecidos entre proteínas. Compara tu "
                "proteína con millones de proteínas ya conocidas: si se parece "
                "mucho a una enzima conocida, probablemente hace lo mismo. Es "
                "como buscar una frase en una biblioteca por similitud, no por "
                "coincidencia exacta."),
    ("Swiss-Prot", "Base de proteínas revisadas a mano; es la referencia contra "
                   "la que compara DIAMOND."),
    ("HMMER / Pfam", "Busca familias de proteínas por su patrón, no por la "
                     "secuencia exacta. Sirve para detectar dominios."),
    ("Rhea", "Base de reacciones químicas: qué entra, qué sale y qué enzima "
             "cataliza cada paso."),
    ("KEGG", "Mapa de rutas metabólicas con nombres y números EC."),
    ("RegulonDB", "Regulación curada de E. coli K-12: operones, promotores y "
                  "factores de transcripción."),
    ("Brecha", "Paso químico que la bacteria necesita pero para el que no "
               "encontramos la enzima. No prueba que falte: puede estar en un "
               "gen no anotado, entrar del exterior o usar otra ruta."),
    ("Operón", "Grupo de genes contiguos que se copian juntos a un mismo ARN."),
    ("ESM-2", "Modelo de lenguaje de proteínas: puntúa si una secuencia parece "
              "una proteína real. Corre en CPU."),
    ("ESMFold", "Predice la estructura 3D de una proteína a partir de su "
                "secuencia. Necesita GPU o un servicio externo."),
]



STATUS_GLOSSARY = [
    ("A / B / C", "A: homología a función nombrada; B: homólogo no caracterizado/DUF; "
     "C: sin hit Swiss-Prot. No son función confirmada ni ausencia de función."),
    ("strong / moderate / weak / none", "Fuerza de la evidencia de homología. "
     "strong no significa función confirmada; none no significa que no exista."),
    ("complex_candidate", "Candidato de complejo por modelo/vecindad. No demuestra "
     "interacción física ni estequiometría."),
    ("complex_suspected", "Sospecha fuerte por ortología recíproca; nivel B. "
     "No es confirmación experimental."),
    ("complex_incomplete", "Se detectó parte de un complejo; falta al menos una "
     "subunidad. No significa que el complejo no exista."),
    ("not_evaluated / not_available", "No se evaluó o no hay dato. No es un "
     "resultado negativo."),
    ("gap / brecha", "Paso metabólico sin enzima detectada. No prueba ausencia: "
     "puede ser transporte, gen no anotado, ruta alternativa o falso positivo."),
    ("supported_hypothesis", "Varias fuentes computacionales coinciden. Sigue "
     "siendo hipótesis; no cierra una brecha."),
    ("afdb_model_available / filter_f1", "Hay modelo o filtro local, pero no "
     "necesariamente se evaluó. Un campo vacío no es un resultado negativo."),
]


def render_glossary() -> None:
    """Glosario corto y reutilizable, visible desde la barra lateral."""
    import streamlit as st
    for term, meaning in GLOSSARY + STATUS_GLOSSARY:
        st.markdown(f"**{term}** — {tr(meaning)}")

def gpu_status() -> dict:
    """Estado de GPU sin importar torch (evita 2 s de arranque)."""
    out = {"nvidia_smi": bool(shutil.which("nvidia-smi")), "gpu": None,
           "note": ""}
    if out["nvidia_smi"]:
        try:
            r = subprocess.run(
                ["nvidia-smi", "--query-gpu=name,memory.total",
                 "--format=csv,noheader"],
                capture_output=True, text=True, timeout=5, check=False,
            )
            line = (r.stdout or "").strip().splitlines()
            if line:
                out["gpu"] = line[0]
        except (OSError, subprocess.SubprocessError):
            pass
    if not out["gpu"]:
        out["note"] = tr("No se detectó una GPU NVIDIA utilizable. La anotación, "
                         "las brechas y la GUI funcionan igual en CPU; ESMFold "
                         "local y el entrenamiento quedan fuera de esta máquina.")
    return out


def render_requirements(expanded: bool = False) -> None:
    """Panel de requisitos, tiempos y glosario para usuarios nuevos."""
    import pandas as pd
    import streamlit as st

    with st.expander(tr("⏱️ ¿Qué necesito y cuánto tarda?"), expanded=expanded):
        st.markdown(
            tr("**No necesitas GPU para usar el anotador.** La GPU solo se usa "
               "para ESMFold local, que es opcional. Los tiempos son aproximados "
               "y dependen del tamaño del genoma y de tu computador.")
        )
        st.markdown("**" + tr("Requisitos del sistema") + "**")
        st.dataframe(
            pd.DataFrame(SYSTEM, columns=[tr("Componente"), tr("Mínimo / recomendado")]),
            hide_index=True, width="stretch",
        )
        st.markdown("**" + tr("Tiempo estimado por paso") + "**")
        df = pd.DataFrame(
            [(tr(a), tr(b), tr(c), tr(d), tr(e)) for a, b, c, d, e in STEPS],
            columns=[tr("Paso"), tr("Qué hace"), tr("Tiempo aproximado"),
                     tr("¿GPU?"), tr("Requiere")],
        )
        st.dataframe(df, hide_index=True, width="stretch")
        st.caption(
            tr("Los 400 genomas del proyecto ya están calculados: explorarlos en "
               "«Rutas metabólicas» es instantáneo. Lo que tarda es analizar un "
               "genoma nuevo o regenerar los datos.")
        )

        g = gpu_status()
        if g["gpu"]:
            st.success(tr("GPU detectada: {gpu}. ESMFold local puede funcionar si "
                          "los pesos están instalados.").format(gpu=g["gpu"]))
        else:
            st.info("🖥️ " + g["note"])
        st.caption(
            tr("La GPU es opcional; el núcleo de anotación funciona en CPU.")
        )

        with st.expander(tr("📖 Glosario para principiantes")):
            for term, meaning in GLOSSARY:
                st.markdown(f"**{term}** — {tr(meaning)}")
