"""Ayuda contextual in situ para la interfaz.

Cada tema vive en :data:`HELP` como un bloque por idioma (español, inglés y
chino). ``render_help(topic)`` dibuja un expander «❓ ¿Qué significa esto?» con
la explicación completa: qué es el concepto, cómo se calcula, cómo se lee la
etiqueta y qué **no** demuestra. La idea es que la persona pueda entender la
pantalla sin abrir la documentación.

Las traducciones EN/ZH son equivalentes del mismo contenido, no un resumen
aparte. El español es la versión de referencia.

Ejemplo de uso::

    from help_ui import render_help

    render_help("transporte", key_prefix="gaps_", expanded=False)

``key_prefix`` solo es necesario cuando el mismo tema se dibuja más de una vez
en la misma página: se usa para que las claves de Streamlit no colisionen.
"""
from __future__ import annotations

from typing import Any

import streamlit as st

#: Etiqueta única del expander que pide la interfaz.
EXPANDER_LABEL = "❓ ¿Qué significa esto?"

#: Orden de las pestañas de idioma. El español siempre es la primera.
LANGUAGES: tuple[tuple[str, str], ...] = (
    ("es", "Español"),
    ("en", "English"),
    ("zh", "中文"),
)

#: Temas de ayuda disponibles. Es el contrato que consumen la GUI y los tests.
HELP: dict[str, dict[str, dict[str, Any]]] = {
    "transporte": {
        "es": {
            "title": "Transporte: qué es un transportador y qué significa high/medium/low",
            "summary": (
                "Un **transportador** es una proteína de membrana que mueve un "
                "metabolito (sustrato) de un lado a otro de la membrana, sola o "
                "acoplada a iones o a ATP. Aquí se infiere por homología con "
                "transportadores curados y por la reacción Rhea exacta del "
                "metabolito; no se mide el flujo ni la actividad en laboratorio."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Qué se exige:** que el metabolito aparezca en los dos lados "
                "de una reacción Rhea de transporte y que exista un homólogo "
                "Swiss-Prot bacteriano ligado a esa misma reacción.",
                "**high:** homólogo bacteriano con pident ≥ 40 %, qcov ≥ 80 % y "
                "e ≤ 1e-30. Es la señal más fuerte, pero sigue siendo una "
                "inferencia computacional.",
                "**medium:** cumple los umbrales mínimos (pident ≥ 35 %, "
                "qcov ≥ 70 %, e ≤ 1e-5) pero no alcanza los de high.",
                "**low:** hay un hit, pero es eucariota o no alcanza los "
                "umbrales (t1_review); no habilita importación.",
                "**Sin reacción exacta:** se marca no_transport o "
                "transport_reaction_only; no se inventa un transportador.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "Una etiqueta high **no es una función confirmada**: no prueba "
                "transporte, dirección ni acoplamiento iónico.",
                "Un alto parecido no descarta otro sustrato: la especificidad "
                "puede divergir entre homólogos.",
                "La capa de membrana por HMM usa además e-values y motivos; no "
                "es TMHMM.",
            ],
        },
        "en": {
            "title": "Transport: what a transporter is and what high/medium/low mean",
            "summary": (
                "A **transporter** is a membrane protein that moves a metabolite "
                "(substrate) across the membrane, either on its own or coupled to "
                "ions or ATP. Here it is inferred from homology to curated "
                "transporters and from the exact Rhea transport reaction of the "
                "metabolite; flux and activity are not measured."
            ),
            "details_heading": "How to read it",
            "details": [
                "**What is required:** the metabolite on both sides of an exact "
                "Rhea transport reaction, plus a bacterial Swiss-Prot homolog "
                "linked to that same reaction.",
                "**high:** bacterial homolog with pident ≥ 40%, qcov ≥ 80% and "
                "e ≤ 1e-30. This is the strongest signal, but still a "
                "computational inference.",
                "**medium:** meets the minimum thresholds (pident ≥ 35%, "
                "qcov ≥ 70%, e ≤ 1e-5) but not the high ones.",
                "**low:** a hit exists, but it is eukaryotic or below the "
                "thresholds (t1_review); it does not enable import.",
                "**No exact reaction:** the gap is marked no_transport or "
                "transport_reaction_only; no transporter is invented.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "A high label **is not a confirmed function**: it does not prove "
                "transport, direction, or ion coupling.",
                "High similarity does not rule out another substrate: "
                "specificity can diverge between homologs.",
                "The membrane HMM layer also uses e-values and motifs; it is "
                "not TMHMM.",
            ],
        },
        "zh": {
            "title": "转运：什么是转运蛋白，high/medium/low 是什么意思",
            "summary": (
                "**转运蛋白**是位于膜上的蛋白质，可将代谢物（底物）从膜的一侧"
                "运到另一侧，可单独完成，也可与离子或 ATP 偶联。这里依据与已整理"
                "转运蛋白的同源性以及代谢物对应的精确 Rhea 转运反应来推断，"
                "并未测定运输通量或实验活性。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**判定条件：** 代谢物必须在同一条 Rhea 转运反应的两侧出现，"
                "并且存在与该反应关联的细菌 Swiss-Prot 同源蛋白。",
                "**high：** 细菌同源蛋白，pident ≥ 40%，qcov ≥ 80%，"
                "e ≤ 1e-30。这是最强信号，但仍属计算推断。",
                "**medium：** 达到最低阈值（pident ≥ 35%，qcov ≥ 70%，"
                "e ≤ 1e-5），但未达到 high 标准。",
                "**low：** 存在命中，但为真核来源或未达阈值（t1_review），"
                "不能据此判断可输入。",
                "**没有精确反应：** 标记为 no_transport 或 "
                "transport_reaction_only；系统不会凭空生成转运蛋白。",
            ],
            "limits_heading": "局限",
            "limits": [
                "high 标签**不是已确认的功能**：不能证明转运、方向或离子偶联。",
                "高相似度不能排除其他底物：同源蛋白的特异性可能已分化。",
                "膜蛋白 HMM 层还使用 e 值和基序，并非 TMHMM。",
            ],
        },
    },
    "membrana": {
        "es": {
            "title": "Membrana: señal péptido, hélices TM y lipobox",
            "summary": (
                "La capa de membrana combina un panel Pfam con tres heurísticas "
                "de secuencia. Sirve para priorizar proteínas de membrana, "
                "exportadas o lipoproteínas; **no sustituye TMHMM, SignalP ni "
                "LipoP**."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Señal péptido:** región N-terminal hidrofóbica con posibles "
                "residuos pequeños en −1/−3 (sitio de corte A-X-A). Sugiere "
                "exportación por la vía Sec; es una heurística débil.",
                "**TM (transmembrana):** segmentos de unos 18–22 aa con "
                "hidropatía alta (Kyte-Doolittle). Varias hélices TM proponen "
                "una proteína integral de membrana.",
                "**Lipobox:** motivo [LVI][ASTM][AG]C en los primeros 40 aa, con "
                "una región hidrofóbica aguas arriba. Sugiere lipoproteína y "
                "anclaje a membrana; no prueba la modificación lipídica.",
                "**Tipos y confianza:** transporter, secretion, lipoprotein, tm "
                "y signal_peptide; high/medium/low combinan el hit Pfam, el "
                "e-value y las señales de secuencia.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "**No es TMHMM:** no calcula topología ni el número exacto de "
                "hélices con un modelo entrenado.",
                "**No es SignalP ni LipoP:** puede producir falsos positivos y "
                "falsos negativos.",
                "Una señal débil no demuestra que la proteína no sea de "
                "membrana.",
            ],
        },
        "en": {
            "title": "Membrane: signal peptide, TM helices and lipobox",
            "summary": (
                "The membrane layer combines a Pfam panel with three sequence "
                "heuristics. It helps prioritize membrane, exported, or "
                "lipoprotein candidates; it **does not replace TMHMM, SignalP, "
                "or LipoP**."
            ),
            "details_heading": "How to read it",
            "details": [
                "**Signal peptide:** hydrophobic N-terminal region with small "
                "residues possibly at −1/−3 (A-X-A cleavage site). It suggests "
                "Sec-pathway export; this is weak heuristic evidence.",
                "**TM (transmembrane):** segments of about 18–22 aa with high "
                "hydropathy (Kyte-Doolittle). Several TM helices suggest an "
                "integral membrane protein.",
                "**Lipobox:** [LVI][ASTM][AG]C motif within the first 40 aa, "
                "with a hydrophobic region upstream. It suggests a lipoprotein "
                "and membrane anchoring; it does not prove lipid modification.",
                "**Types and confidence:** transporter, secretion, lipoprotein, "
                "tm, and signal_peptide; high/medium/low combine the Pfam hit, "
                "the e-value, and the sequence signals.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "**It is not TMHMM:** it does not compute topology or the exact "
                "number of helices with a trained model.",
                "**It is not SignalP or LipoP:** it can produce false positives "
                "and false negatives.",
                "A weak signal does not prove that the protein is not a membrane "
                "protein.",
            ],
        },
        "zh": {
            "title": "膜蛋白：信号肽、TM 螺旋与 lipobox",
            "summary": (
                "膜蛋白层将 Pfam 面板与三种序列启发式方法结合，用于优先筛选"
                "膜蛋白、分泌蛋白或脂蛋白候选；**它不能替代 TMHMM、SignalP "
                "或 LipoP**。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**信号肽：** N 端疏水区，−1/−3 位常为小残基（A-X-A 切割位点），"
                "提示经 Sec 途径输出；这是较弱的启发式证据。",
                "**TM（跨膜段）：** 约 18–22 个氨基酸的高疏水性片段"
                "（Kyte-Doolittle）。多个 TM 螺旋提示整合膜蛋白。",
                "**Lipobox：** N 端前 40 个氨基酸内的 [LVI][ASTM][AG]C 基序，"
                "其上游有疏水区；提示脂蛋白和膜锚定，但不能证明脂质修饰。",
                "**类型与置信度：** transporter、secretion、lipoprotein、tm、"
                "signal_peptide；high/medium/low 综合 Pfam 命中、e 值和序列信号。",
            ],
            "limits_heading": "局限",
            "limits": [
                "**并非 TMHMM：** 未用训练模型计算拓扑和精确螺旋数。",
                "**并非 SignalP 或 LipoP：** 可能产生假阳性和假阴性。",
                "信号弱并不证明该蛋白不是膜蛋白。",
            ],
        },
    },
    "operones": {
        "es": {
            "title": "Operones: distancia, reglas y lo que no demuestran",
            "summary": (
                "Un **operón** es un grupo de genes que se transcriben juntos "
                "desde un promotor. Aquí hay dos fuentes: operones curados de "
                "RegulonDB (K-12) y operones putativos construidos por reglas de "
                "distancia en genomas nuevos. Ninguna prueba cotranscripción."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Regla de distancia:** mismo contig, misma hebra y hueco "
                "≤ 300 pb entre CDS consecutivos. Es una condición necesaria, no "
                "una prueba de cotranscripción.",
                "**Reglas R1–R5:** atenuador, terminador fuerte, promotor σ70, "
                "distancia y ausencia de señal deciden unir o cortar; quedan como "
                "hipótesis con evidencia rule_based.",
                "**TSS:** sitio de inicio de transcripción. El promotor putativo "
                "es una ventana fija de 50 pb upstream con consenso −35 TTGACA / "
                "−10 TATAAT (espaciador 15–19); no es un TSS experimental.",
                "**TTS/terminador:** sitio de terminación. En el catálogo curado "
                "se muestra si el operón tiene terminador documentado (cobertura "
                "baja); en las reglas, un terminador fuerte puede cortar el grupo.",
                "**Catálogo curado:** cada operón puede traer promotor, TSS, "
                "cajas −10/−35, sigma, TF, terminador y PMIDs; la GUI muestra la "
                "completitud para no asumir piezas ausentes.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "Distancia y motivo no demuestran cotranscripción ni expresión.",
                "Un operón puede quedar partido o unido de más: la regla no "
                "conoce los TSS/TTS reales del genoma nuevo.",
                "Los componentes curados son de E. coli K-12 y no se transfieren "
                "automáticamente a otra especie.",
            ],
        },
        "en": {
            "title": "Operons: distance, rules, and what they do not prove",
            "summary": (
                "An **operon** is a group of genes transcribed together from one "
                "promoter. There are two sources here: curated RegulonDB operons "
                "(K-12) and putative operons built with distance rules on new "
                "genomes. Neither proves cotranscription."
            ),
            "details_heading": "How to read it",
            "details": [
                "**Distance rule:** same contig, same strand, and a gap ≤ 300 bp "
                "between consecutive CDSs. It is necessary, but not proof of "
                "cotranscription.",
                "**Rules R1–R5:** attenuator, strong terminator, σ70 promoter, "
                "distance, and absence of signal decide whether to join or split; "
                "they remain hypotheses with rule_based evidence.",
                "**TSS:** transcription start site. The putative promoter is a "
                "fixed 50 bp upstream window matching the −35 TTGACA / −10 TATAAT "
                "consensus (15–19 bp spacer); it is not an experimental TSS.",
                "**TTS/terminator:** transcription termination site. The curated "
                "catalog shows whether a documented terminator exists (low "
                "coverage); in the rules, a strong terminator can split the group.",
                "**Curated catalog:** each operon may include promoter, TSS, "
                "−10/−35 boxes, sigma factor, TF, terminator, and PMIDs; the GUI "
                "shows completeness so missing pieces are not assumed.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "Distance and motifs do not prove cotranscription or expression.",
                "An operon can be split or merged incorrectly: the rule does not "
                "know the real TSS/TTS of the new genome.",
                "Curated components come from E. coli K-12 and are not "
                "automatically transferred to another species.",
            ],
        },
        "zh": {
            "title": "操纵子：距离、规则以及它们不能证明什么",
            "summary": (
                "**操纵子**是从同一启动子一起转录的一组基因。本界面有两个来源："
                "RegulonDB（K-12）整理的操纵子，以及在新基因组上按距离规则构建的"
                "推定操纵子。两者都不能证明共转录。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**距离规则：** 同一 contig、同一链，相邻 CDS 间隔 ≤ 300 bp。"
                "这是必要条件，不是共转录的证据。",
                "**R1–R5 规则：** 衰减子、强终止子、σ70 启动子、距离和无信号"
                "决定合并或拆分；结果仍是带 rule_based 证据的假设。",
                "**TSS：** 转录起始位点。推定启动子是首个基因上游 50 bp 的固定"
                "窗口，匹配 −35 TTGACA / −10 TATAAT 一致性序列（间隔 15–19 bp），"
                "并非实验 TSS。",
                "**TTS/终止子：** 转录终止位点。整理目录显示是否有已记录的"
                "终止子（覆盖率低）；在规则中，强终止子可拆分基因簇。",
                "**整理目录：** 每个操纵子可含启动子、TSS、−10/−35 框、sigma "
                "因子、TF、终止子和 PMID；界面显示完整度，避免假定缺失部件存在。",
            ],
            "limits_heading": "局限",
            "limits": [
                "距离和基序不能证明共转录或表达。",
                "操纵子可能被错误拆分或合并：规则不知道新基因组的真实 TSS/TTS。",
                "整理组件来自大肠杆菌 K-12，不会自动转移到其他物种。",
            ],
        },
    },
    "brechas": {
        "es": {
            "title": "Brechas metabólicas: qué son y por qué no se cierran",
            "summary": (
                "Una **brecha** es una reacción que la ruta necesita pero para la "
                "que no se encontró una enzima candidata en la anotación (o no se "
                "pudo asignar a un gen). No es necesariamente un error: puede "
                "faltar anotación, ser una reacción espontánea, un transportador, "
                "un complejo o una ruta que el organismo no tiene."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Cómo se detecta:** se compara el conjunto de reacciones "
                "esperadas de Rhea/KEGG con las que los genes anotados pueden "
                "catalizar. Si no hay candidato, queda como brecha.",
                "**Por qué no se cierra:** falta homología, el gen está mal "
                "anotado o truncado, la enzima requiere un complejo o cofactor, la "
                "reacción es no enzimática, o el sustrato/producto está mal "
                "asignado.",
                "**Estados:** validated (soporte alto en el panel de 400 genomas), "
                "consensus (soporte medio), novel (sin soporte) y anomalous "
                "(entrada/salida, compuesto exógeno o química espontánea).",
                "**Dirección:** queda direction_unspecified hasta tener ΔG o "
                "fisiología; Rhea no fija qué lado es sustrato y cuál producto.",
                "**Transporte aparte:** no_transport, transport_reaction_only "
                "(T2), t1_review y t1_bacterial (T1). Una brecha puede cerrarse "
                "por importación, no por una enzima interna.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "Una brecha no prueba que falte la función: la ausencia de hit no "
                "es evidencia negativa.",
                "No se cierra por parecido estructural ni por vecindad de "
                "complejo; eso solo genera hipótesis.",
                "Los estados resumen evidencia computacional, no validación "
                "experimental.",
            ],
        },
        "en": {
            "title": "Metabolic gaps: what they are and why they stay open",
            "summary": (
                "A **gap** is a reaction the pathway needs but for which no "
                "candidate enzyme was found in the annotation (or it could not be "
                "assigned to a gene). It is not necessarily an error: annotation "
                "may be missing, the reaction may be spontaneous, or it may be "
                "handled by a transporter, a complex, or a pathway the organism "
                "lacks."
            ),
            "details_heading": "How to read it",
            "details": [
                "**How it is detected:** the expected Rhea/KEGG reaction set is "
                "compared with those the annotated genes can catalyze. With no "
                "candidate, the reaction stays as a gap.",
                "**Why it is not closed:** homology is missing, the gene is "
                "misannotated or truncated, the enzyme needs a complex or "
                "cofactor, the reaction is non-enzymatic, or the substrate/product "
                "is misassigned.",
                "**States:** validated (high support in the 400-genome panel), "
                "consensus (medium support), novel (no support), and anomalous "
                "(entry/exit, exogenous compound, or spontaneous chemistry).",
                "**Direction:** it stays direction_unspecified until ΔG or "
                "physiology is available; Rhea does not fix which side is the "
                "substrate.",
                "**Transport is separate:** no_transport, "
                "transport_reaction_only (T2), t1_review, and t1_bacterial (T1). A "
                "gap may close through import rather than an internal enzyme.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "A gap does not prove the function is missing: no hit is not "
                "negative evidence.",
                "Structural similarity or complex neighborhood does not close a "
                "gap; it only generates hypotheses.",
                "The states summarize computational evidence, not experimental "
                "validation.",
            ],
        },
        "zh": {
            "title": "代谢缺口：是什么，为什么没有闭合",
            "summary": (
                "**缺口**是通路需要、但在注释中找不到候选酶的某个反应"
                "（或无法分配到基因）。它不一定是错误：可能是注释缺失、"
                "自发反应、转运蛋白、复合体，或该生物根本没有这条通路。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**检测方式：** 将 Rhea/KEGG 期望的反应集合与已注释基因可能"
                "催化的反应比较；没有候选即记为缺口。",
                "**未闭合原因：** 缺少同源序列、基因注释错误或截短、酶需要"
                "复合体或辅因子、反应非酶促，或底物/产物指派错误。",
                "**状态：** validated（400 基因组面板中支持度高）、consensus"
                "（中等支持）、novel（无支持）、anomalous（进出点、外源化合物"
                "或自发化学）。",
                "**方向：** 在获得 ΔG 或生理证据前保持 direction_unspecified；"
                "Rhea 不规定哪一侧是底物。",
                "**转运单独评估：** no_transport、transport_reaction_only"
                "（T2）、t1_review、t1_bacterial（T1）。缺口可通过输入闭合，"
                "而不一定靠内部酶。",
            ],
            "limits_heading": "局限",
            "limits": [
                "缺口不能证明功能缺失：没有命中不是阴性证据。",
                "结构相似或复合体邻接不能闭合缺口，只产生假设。",
                "这些状态只概括计算证据，不是实验验证。",
            ],
        },
    },
    "complejos": {
        "es": {
            "title": "Complejos: vecindad, modelo y transferencia",
            "summary": (
                "Un **complejo** es un conjunto de subunidades que funcionan "
                "juntas. Aquí se sospecha por tres vías: genes vecinos en el "
                "mismo contig/hebra, un modelo entrenado con complejos de K-12 y "
                "transferencia por ortología desde K-12. Ninguna confirma "
                "interacción física ni estructura cuaternaria."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Vecindad:** agrupa CDS cercanos en una ventana de distancia. "
                "La proximidad es necesaria para muchos operones, pero no "
                "demuestra que las proteínas se unan.",
                "**Modelo:** aplica un modelo de K-12 a genes/operones del genoma "
                "nuevo y propone complex_candidate. Es una sospecha estadística, "
                "no un ensayo.",
                "**Transferencia:** si hay ortólogos recíprocos de las subunidades, "
                "el complejo se infiere por ortología (nivel B, no experimental).",
                "**Colores:** verde = firma curada completa; azul = transferido "
                "inferido; amarillo = sospecha/incompleto; naranja = ambiguo "
                "(varios complejos compatibles); gris = no evaluado.",
                "**Amarillo no confirma:** indica que falta una subunidad, que el "
                "candidato no está validado o que la evidencia es débil. No cierra "
                "la brecha.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "No se predicen interfaces, estequiometría ni estructura "
                "cuaternaria de novo.",
                "La vecindad puede compartir operón sin interacción física, y "
                "subunidades esenciales pueden estar lejos.",
                "Un complejo incompleto puede reflejar anotación faltante, no "
                "ausencia real.",
            ],
        },
        "en": {
            "title": "Complexes: neighborhood, model, and transfer",
            "summary": (
                "A **complex** is a set of subunits that work together. Here it is "
                "suspected through three routes: neighboring genes on the same "
                "contig/strand, a model trained on K-12 complexes, and transfer by "
                "orthology from K-12. None confirms physical interaction or "
                "quaternary structure."
            ),
            "details_heading": "How to read it",
            "details": [
                "**Neighborhood:** groups nearby CDSs within a distance window. "
                "Proximity is necessary for many operons, but it does not prove "
                "that the proteins bind each other.",
                "**Model:** applies a K-12 model to genes/operons of the new genome "
                "and proposes a complex_candidate. It is a statistical suspicion, "
                "not an assay.",
                "**Transfer:** when reciprocal orthologs of the subunits exist, "
                "the complex is inferred by orthology (level B, not experimental).",
                "**Colors:** green = complete curated signature; blue = inferred "
                "transfer; yellow = suspected/incomplete; orange = ambiguous "
                "(several compatible complexes); gray = not evaluated.",
                "**Yellow does not confirm:** it means a subunit is missing, the "
                "candidate is not validated, or the evidence is weak. It does not "
                "close the gap.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "Interfaces, stoichiometry, and quaternary structure are not "
                "predicted de novo.",
                "Neighborhood can mean a shared operon without physical "
                "interaction, and essential subunits can be far apart.",
                "An incomplete complex may reflect missing annotation rather than "
                "true absence.",
            ],
        },
        "zh": {
            "title": "复合体：邻接、模型与转移",
            "summary": (
                "**复合体**是共同发挥功能的一组亚基。这里通过三条途径推测："
                "同一 contig/链上的邻近基因、用 K-12 复合体训练的模型，以及从 "
                "K-12 的直系同源转移。它们都不能确认物理互作或四级结构。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**邻接：** 将距离窗口内的相邻 CDS 分组。邻近是许多操纵子的"
                "必要条件，但不能证明蛋白彼此结合。",
                "**模型：** 把 K-12 模型应用于新基因组的基因/操纵子，给出 "
                "complex_candidate。这是统计猜测，不是实验。",
                "**转移：** 若各亚基存在双向最佳直系同源，则按同源关系推断"
                "复合体（B 级，非实验）。",
                "**颜色：** 绿色 = 整理特征完整；蓝色 = 推断转移；黄色 = "
                "可疑/不完整；橙色 = 模糊（多个兼容复合体）；灰色 = 未评估。",
                "**黄色不代表确认：** 表示缺少亚基、候选未验证或证据弱，"
                "不能闭合缺口。",
            ],
            "limits_heading": "局限",
            "limits": [
                "不预测界面、化学计量比或从头预测四级结构。",
                "邻接可能只是同一操纵子而无物理互作；必需亚基也可能相距很远。",
                "不完整复合体可能只是注释缺失，而非真的不存在。",
            ],
        },
    },
    "regulacion": {
        "es": {
            "title": "Regulación: reglas frente a evidencia experimental",
            "summary": (
                "La capa de regulación propone operones y promotores putativos con "
                "reglas de secuencia y distancia. Son hipótesis (rule_based), no "
                "mediciones de expresión ni sitios de unión confirmados. La "
                "evidencia experimental solo aparece cuando el catálogo curado la "
                "aporta, por ejemplo RegulonDB en K-12."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Promotor putativo:** ventana de 50 pb upstream del primer gen "
                "con el mejor par −35 TTGACA / −10 TATAAT, espaciador 15–19 y "
                "puntuación de coincidencia. No modela sigma alternativo, "
                "UP-element, −10 extendida ni promotores internos.",
                "**Operón putativo:** mismo contig/hebra y hueco ≤ 300 pb; puede "
                "unir genes que no se cotranscriben o partir un operón real.",
                "**Confianza:** solo low o medium; refleja cuántas señales "
                "coinciden, no una probabilidad calibrada.",
                "**Evidencia experimental:** TSS, sitios de unión de TF, "
                "promotores y terminadores documentados con PMID se muestran "
                "aparte y no se infieren para genomas nuevos.",
                "**Predicho no es observado:** un motivo puede ser falso positivo "
                "y un gen puede regularse sin un motivo evidente.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "El cálculo no usa red ni bases externas; no consulta "
                "RegulonDB/DBTBS para genomas nuevos.",
                "No mide expresión, ocupancia ni efecto regulatorio.",
                "Un promotor putativo no convierte un gen en regulado ni cierra "
                "una brecha.",
            ],
        },
        "en": {
            "title": "Regulation: rules versus experimental evidence",
            "summary": (
                "The regulation layer proposes putative operons and promoters with "
                "sequence and distance rules. They are hypotheses (rule_based), "
                "not expression measurements or confirmed binding sites. "
                "Experimental evidence appears only when the curated catalog "
                "provides it, for example RegulonDB for K-12."
            ),
            "details_heading": "How to read it",
            "details": [
                "**Putative promoter:** 50 bp upstream window of the first gene "
                "with the best −35 TTGACA / −10 TATAAT pair, a 15–19 bp spacer, "
                "and a match score. It does not model alternative sigma factors, "
                "UP elements, extended −10, or internal promoters.",
                "**Putative operon:** same contig/strand and gap ≤ 300 bp; it can "
                "join genes that are not cotranscribed or split a real operon.",
                "**Confidence:** only low or medium; it reflects how many signals "
                "agree, not a calibrated probability.",
                "**Experimental evidence:** TSSs, TF binding sites, promoters, and "
                "terminators documented with PMIDs are shown separately and are "
                "not inferred for new genomes.",
                "**Predicted is not observed:** a motif can be a false positive, "
                "and a gene can be regulated without an obvious motif.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "The calculation uses no network or external databases; it does "
                "not query RegulonDB/DBTBS for new genomes.",
                "It does not measure expression, occupancy, or regulatory effect.",
                "A putative promoter does not make a gene regulated or close a "
                "gap.",
            ],
        },
        "zh": {
            "title": "调控：规则与实验证据",
            "summary": (
                "调控层用序列和距离规则提出推定操纵子和启动子。它们是 "
                "rule_based 假设，不是表达测量或已确认的结合位点。只有整理目录"
                "提供时（如 K-12 的 RegulonDB）才出现实验证据。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**推定启动子：** 首基因上游 50 bp 窗口，寻找最佳 "
                "−35 TTGACA / −10 TATAAT 对，间隔 15–19 bp 并打分。不建模"
                "替代 sigma、UP 元件、延伸 −10 或内部启动子。",
                "**推定操纵子：** 同一 contig/链且间隔 ≤ 300 bp；可能合并"
                "并不共转录的基因，或拆开真实操纵子。",
                "**置信度：** 只有 low 或 medium，反映信号重合程度，"
                "不是校准概率。",
                "**实验证据：** 带 PMID 的 TSS、TF 结合位点、启动子和终止子"
                "单独显示，不会为新基因组推断。",
                "**预测调控 ≠ 观测调控：** 基序可能假阳性，基因也可能无明显"
                "基序仍受调控。",
            ],
            "limits_heading": "局限",
            "limits": [
                "计算不使用网络或外部数据库；对新基因组不查询 "
                "RegulonDB/DBTBS。",
                "不测量表达、占据率或调控效应。",
                "推定启动子不会把基因变成受调控基因，也不能闭合缺口。",
            ],
        },
    },
    "hipotesis": {
        "es": {
            "title": "Hipótesis estructurales: F1–F5 y sus etiquetas",
            "summary": (
                "La pestaña de hipótesis combina cinco filtros computacionales. "
                "Cada uno aporta una señal distinta; juntos priorizan qué revisar, "
                "pero ninguno confirma función, sustrato, reacción ni cierre de "
                "brecha."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**F1 Foldseek:** encuentra proteínas con pliegue similar en "
                "PDB/AlphaFold. La forma puede conservarse aunque la reacción "
                "cambie.",
                "**F2 EnzyMM/M-CSA:** plantilla enzimática con EC; sugiere familia "
                "y tipo de reacción.",
                "**F3 P2Rank:** predice un bolsillo en la superficie; la geometría "
                "no identifica el sustrato.",
                "**F4 AlphaFill:** transfiere ligandos/cofactores desde estructuras "
                "homólogas; apoya qué químico podría unirse.",
                "**F5 Sitio activo:** posiciones curadas transferidas desde la "
                "plantilla; permite revisar residuos concretos, no mide actividad.",
                "**supported_hypothesis:** F2 + F3 + F4, tres fuentes coinciden y "
                "no hay conflicto. **conflict:** la evidencia contradice la "
                "anotación base o dos fuentes discrepan; hay que revisar antes de "
                "usar la hipótesis.",
                "**Otras etiquetas:** high_priority (F2 + F5), orientative (solo "
                "F2), structural_only (F1 + F3 o F4) y uncharacterized (sin señal).",
            ],
            "limits_heading": "Límites",
            "limits": [
                "Ninguna etiqueta es un ensayo enzimático: no confirma EC, "
                "sustrato ni función in vivo.",
                "Tres señales pueden coincidir por sesgo estructural y no por la "
                "misma reacción.",
                "Un conflicto no se resuelve solo ni cambia la anotación base.",
            ],
        },
        "en": {
            "title": "Structural hypotheses: F1–F5 and their labels",
            "summary": (
                "The hypotheses tab combines five computational filters. Each "
                "contributes a different signal; together they prioritize what to "
                "review, but none confirms function, substrate, reaction, or gap "
                "closure."
            ),
            "details_heading": "How to read it",
            "details": [
                "**F1 Foldseek:** finds proteins with a similar fold in "
                "PDB/AlphaFold. Shape can be conserved even when the reaction "
                "changes.",
                "**F2 EnzyMM/M-CSA:** enzymatic template with an EC number; it "
                "suggests a family and reaction type.",
                "**F3 P2Rank:** predicts a surface pocket; geometry does not "
                "identify the substrate.",
                "**F4 AlphaFill:** transfers ligands/cofactors from homologous "
                "structures; it supports which chemical could bind.",
                "**F5 Active site:** curated positions transferred from the "
                "template; it allows residue-level review but does not measure "
                "activity.",
                "**supported_hypothesis:** F2 + F3 + F4, three sources agree and "
                "there is no conflict. **conflict:** the evidence contradicts the "
                "base annotation or two sources disagree; review before using the "
                "hypothesis.",
                "**Other labels:** high_priority (F2 + F5), orientative (F2 only), "
                "structural_only (F1 + F3 or F4), and uncharacterized (no signal).",
            ],
            "limits_heading": "Limits",
            "limits": [
                "No label is an enzyme assay: it does not confirm EC, substrate, "
                "or in vivo function.",
                "Three signals can agree because of structural bias rather than "
                "the same reaction.",
                "A conflict is not resolved automatically and does not change the "
                "base annotation.",
            ],
        },
        "zh": {
            "title": "结构假设：F1–F5 与标签",
            "summary": (
                "假设页面结合五个计算过滤器。每个过滤器提供不同证据；它们共同"
                "用于排序复核优先级，但任何一个都不能确认功能、底物、反应或"
                "缺口闭合。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**F1 Foldseek：** 在 PDB/AlphaFold 中寻找折叠相似蛋白。"
                "形状可保守而反应已改变。",
                "**F2 EnzyMM/M-CSA：** 带 EC 的酶模板，提示酶家族和反应类型。",
                "**F3 P2Rank：** 预测表面口袋；几何形状不能鉴定底物。",
                "**F4 AlphaFill：** 从同源结构转移配体/辅因子，支持可能结合的"
                "化学物。",
                "**F5 活性位点：** 从模板转移的整理位点，可用具体残基复核，"
                "但不测活性。",
                "**supported_hypothesis：** F2 + F3 + F4，三来源一致且无冲突。"
                "**conflict：** 证据与基础注释矛盾或来源不一致，使用前必须复核。",
                "**其他标签：** high_priority（F2 + F5）、orientative（仅 F2）、"
                "structural_only（F1 + F3 或 F4）、uncharacterized（无信号）。",
            ],
            "limits_heading": "局限",
            "limits": [
                "任何标签都不是酶学实验：不确认 EC、底物或体内功能。",
                "三个信号可能因结构偏倚而一致，并非同一反应。",
                "冲突不会自动解决，也不会改变基础注释。",
            ],
        },
    },
    "taxonomia": {
        "es": {
            "title": "Taxonomía: precalculada, sin ANI en vivo",
            "summary": (
                "La taxonomía se lee de un JSON local ya precalculado (normalmente "
                "GTDB r232; si no, NCBI) y se normaliza a una tabla de rangos. Esta "
                "capa no consulta NCBI/GTDB, no ejecuta fastANI/GTDB-Tk y **no "
                "valida la especie en vivo**."
            ),
            "details_heading": "Cómo interpretarlo",
            "details": [
                "**Qué contiene:** rangos domain → phylum → class → order → "
                "family → genus → species → strain, con confidence, method "
                "(gtdb/ncbi) y source.",
                "**Cómo se elige:** prefiere GTDB r232 si el JSON lo trae; si no, "
                "NCBI. No fusiona ambos linajes: cada fila conserva su método y las "
                "discrepancias quedan en conflicts.",
                "**Confianza:** se copia del JSON; si falta, unknown. No se "
                "recalcula ni se convierte en una probabilidad.",
                "**Rangos ausentes:** no se escriben y no son evidencia negativa; "
                "una especie solo aparece si el JSON la declara.",
                "**ANI:** si el panel incluye grupos ANI, son datos precalculados "
                "(por ejemplo ani_groups); no se ejecuta un ANI nuevo al abrir la "
                "vista.",
            ],
            "limits_heading": "Límites",
            "limits": [
                "No valida la especie por ANI, marcadores ni red en tiempo real.",
                "Un rango ausente no significa que el organismo no tenga ese "
                "rango.",
                "Los conflictos NCBI/GTDB no se resuelven automáticamente: hay que "
                "revisarlos.",
            ],
        },
        "en": {
            "title": "Taxonomy: precomputed, with no live ANI",
            "summary": (
                "Taxonomy is read from a precomputed local JSON (usually GTDB "
                "r232; otherwise NCBI) and normalized into a rank table. This layer "
                "does not query NCBI/GTDB, does not run fastANI/GTDB-Tk, and **does "
                "not validate the species live**."
            ),
            "details_heading": "How to read it",
            "details": [
                "**What it contains:** ranks domain → phylum → class → order → "
                "family → genus → species → strain, with confidence, method "
                "(gtdb/ncbi), and source.",
                "**How it is chosen:** it prefers GTDB r232 when the JSON provides "
                "it; otherwise NCBI. It does not merge the two lineages: each row "
                "keeps its method, and discrepancies remain in conflicts.",
                "**Confidence:** copied from the JSON; if missing, unknown. It is "
                "not recalculated or turned into a probability.",
                "**Missing ranks:** they are not written and are not negative "
                "evidence; a species appears only if the JSON declares it.",
                "**ANI:** if the panel includes ANI groups, they are precomputed "
                "data (for example ani_groups); no new ANI runs when the view "
                "opens.",
            ],
            "limits_heading": "Limits",
            "limits": [
                "It does not validate the species by ANI, markers, or network in "
                "real time.",
                "A missing rank does not mean the organism lacks that rank.",
                "NCBI/GTDB conflicts are not resolved automatically: they must be "
                "reviewed.",
            ],
        },
        "zh": {
            "title": "分类学：预计算，无实时 ANI",
            "summary": (
                "分类学从本地预计算 JSON 读取（通常 GTDB r232，否则 NCBI），"
                "并规范化为等级表。本层不查询 NCBI/GTDB，不运行 fastANI/"
                "GTDB-Tk，也**不实时验证物种**。"
            ),
            "details_heading": "如何解读",
            "details": [
                "**内容：** domain → phylum → class → order → family → genus → "
                "species → strain 等级，附 confidence、method（gtdb/ncbi）"
                "和 source。",
                "**选择方式：** JSON 有 GTDB r232 时优先使用，否则用 NCBI。"
                "不合并两个谱系：每行保留其方法，差异记录在 conflicts。",
                "**置信度：** 直接从 JSON 复制；缺失时为 unknown。不重新计算，"
                "也不当作概率。",
                "**缺失等级：** 不写入，且不是阴性证据；只有 JSON 声明的物种"
                "才会出现。",
                "**ANI：** 若面板包含 ANI 分组，那是预计算数据（如 ani_groups）；"
                "打开视图时不会重新运行 ANI。",
            ],
            "limits_heading": "局限",
            "limits": [
                "不通过 ANI、标记基因或网络实时验证物种。",
                "缺失等级不代表该生物没有该等级。",
                "NCBI/GTDB 冲突不会自动解决，需要人工复核。",
            ],
        },
    },
}

#: Orden estable de temas (también usado por la GUI y los tests).
TOPICS: tuple[str, ...] = tuple(HELP)

_LANG_CODES = tuple(code for code, _ in LANGUAGES)


def get_help(topic: str, lang: str = "es") -> dict[str, Any]:
    """Devuelve el bloque de ayuda de ``topic`` en ``lang`` (sin tocar Streamlit)."""
    if topic not in HELP:
        raise KeyError(f"tema de ayuda desconocido: {topic!r}; "
                       f"disponibles: {', '.join(TOPICS)}")
    if lang not in _LANG_CODES:
        raise KeyError(f"idioma no soportado: {lang!r}; "
                       f"disponibles: {', '.join(_LANG_CODES)}")
    return HELP[topic][lang]


def _render_block(block: dict[str, Any]) -> None:
    """Dibuja un bloque de idioma ya validado."""
    st.markdown(f"**{block['title']}**")
    st.markdown(block["summary"])
    st.markdown(f"**{block['details_heading']}**")
    st.markdown("\n".join(f"- {item}" for item in block["details"]))
    st.markdown(f"**{block['limits_heading']}**")
    st.markdown("\n".join(f"- {item}" for item in block["limits"]))


def render_help(topic: str, *, key_prefix: str = "",
                expanded: bool = False) -> None:
    """Dibuja la ayuda in situ de ``topic`` en un expander.

    Parámetros
    ----------
    topic:
        Una clave de :data:`HELP` (``transporte``, ``membrana``, ``operones``,
        ``brechas``, ``complejos``, ``regulacion``, ``hipotesis``,
        ``taxonomia``).
    key_prefix:
        Prefijo opcional para la clave de Streamlit del expander. Úsalo cuando
        el mismo tema se dibuje más de una vez en una página.
    expanded:
        Si el expander empieza abierto.

    Un tema desconocido no rompe la página: muestra un aviso con los temas
    disponibles y vuelve sin dibujar el expander.
    """
    entry = HELP.get(topic)
    if entry is None:
        st.warning(
            f"No hay ayuda para «{topic}». Temas disponibles: "
            f"{', '.join(TOPICS)}."
        )
        return

    key = f"{key_prefix}help_{topic}" if key_prefix else None
    with st.expander(EXPANDER_LABEL, expanded=bool(expanded), key=key):
        tabs = st.tabs([label for _, label in LANGUAGES])
        for tab, (code, _) in zip(tabs, LANGUAGES):
            with tab:
                _render_block(entry[code])
