"""Textos explicativos de la interfaz, en un solo sitio.

Separar los textos del código permite revisarlos y traducirlos sin tocar la
lógica. El tono está pensado para un estudiante de secundaria: frases cortas,
sin jerga, con una analogía cuando ayuda, y definiendo el término técnico en la
misma frase.

Convención de claves: ``<sección>_<uso>`` (``intro``, ``read``, ``m1``…).

Traducción: ``T`` guarda el texto **en español** y lo resuelve en el idioma
activo al leerlo (``T["clave"]``). Las traducciones viven en
``webapp/locales/<código>.json`` y usan el texto español como clave, así que
este archivo sigue siendo la fuente de verdad en español.
"""
from __future__ import annotations

from i18n import tr


class _Textos(dict):
    """Diccionario que traduce al vuelo; el español es la fuente."""

    def __getitem__(self, key: str) -> str:
        return tr(super().__getitem__(key))

    def get(self, key: str, default=None):
        return self[key] if key in self else default


T: dict[str, str] = _Textos()

# ---------------------------------------------------------------- general
T["app_subtitle"] = (
    "Esta aplicación lee el ADN de una bacteria y te dice qué genes tiene, "
    "qué reacciones químicas puede hacer y qué le falta para vivir. "
    "Todo se ejecuta en tu computador: nada se envía a internet."
)
T["mode_label"] = "¿Qué quieres hacer?"
T["mode_help"] = (
    "**Anotar FASTA** analiza un genoma que subas tú. "
    "**Rutas metabólicas** muestra el análisis ya calculado de los 400 genomas "
    "del proyecto, para explorarlos sin esperar. "
    "**Reportes y descargas** reúne los resultados que puedes reutilizar. "
    "**Módulos** muestra las opciones y autopruebas de cada módulo."
)
T["mode_annotate"] = "🧬 Anotar FASTA"
T["mode_pathways"] = "🧪 Rutas metabólicas"
T["mode_browser"] = "🧬 Navegador de genoma"
T["mode_reports"] = "📥 Reportes y descargas"
T["mode_modules"] = "🧩 Módulos"
T["mode_credits"] = "🙏 Créditos y herramientas"
T["sidebar_tools"] = (
    "Herramientas detectadas en tu equipo. Si una aparece con ❌, la casilla que "
    "la usa queda desactivada."
)
T["sidebar_options"] = (
    "Marca lo que quieres calcular. Todo lo marcado se ejecuta solo al pulsar "
    "**Anotar**."
)
T["sidebar_metabolism"] = (
    "Opcional: reconstruir la red de reacciones químicas y buscar sus huecos."
)
T["sidebar_params"] = "Ajustes finos. Los valores por defecto funcionan bien."
T["upload_label"] = "Arrastra uno o varios archivos FASTA (.fa, .fna, .fasta, .gz)"
T["upload_help"] = (
    "Un FASTA es un archivo de texto con la secuencia de ADN. Puedes subir "
    "varios genomas a la vez: se procesan uno por uno."
)

# ------------------------------------------------------------- anotación
T["opt_cds"] = "**CDS (genes que producen proteína)** — busca los tramos de ADN que empiezan por un codón de inicio, terminan en un codón de parada y se traducen a una proteína. Sin esto no hay nada más que analizar."
T["opt_trna"] = "**tRNA** — los ARN que llevan los aminoácidos a la fábrica de proteínas. Saber cuántos hay y de qué tipo ayuda a detectar genomas incompletos."
T["opt_rrna"] = "**rRNA** — el ARN que forma los ribosomas, las máquinas que construyen proteínas. En genomas grandes puede tardar un minuto."
T["opt_ghost"] = "**ORFs no anotados** — busca posibles genes en los huecos entre genes conocidos. La idea: si una bacteria necesita una enzima y no la encontramos, puede estar escondida en un tramo que nadie anotó."
T["opt_esm"] = "**ESM-2** — un modelo de lenguaje de proteínas. Puntúa si una secuencia parece una proteína real o una casualidad estadística. Corre en CPU, sin tarjeta gráfica."
T["opt_deserts"] = "**Desiertos** — tramos de 3.000 letras o más sin genes conocidos. No son basura: son regiones que todavía no sabemos leer."
T["opt_compare"] = "**Comparar genomas subidos** — si subes varios, agrupa los genes parecidos para saber cuáles se conservan entre ellos. Conservarse suele significar que importan."
T["opt_metabolic"] = "**Red metabólica y brechas** — anota las proteínas contra una base de enzimas conocidas (Swiss-Prot) y las conecta con las reacciones de Rhea. Requiere descargar Swiss-Prot una vez."
T["param_min_region"] = "Longitud mínima de un hueco entre genes para analizarlo. Valores bajos generan mucho ruido."
T["param_min_aa"] = "Longitud mínima de una posible proteína. Por debajo de 60 aminoácidos es difícil distinguir un gen real de una casualidad."
T["param_esm_top"] = "A cuántos candidatos se les pasa el modelo ESM. Puntuarlos todos tarda mucho; los mejores suelen estar entre los primeros."
T["bell"] = "Avisa con un sonido cuando termine, para no tener que mirar la pantalla."
T["demo"] = "Carga dos genomas de ejemplo para probar la aplicación sin subir nada."
T["run_button"] = "Empezar el análisis"
T["run_help"] = "El análisis es automático: no hay que decidir nada mientras corre."

# --------------------------------------------------------------- resultado
T["res_contigs"] = "**Contigs** — trozos de ADN ensamblados. Una bacteria con el genoma completo suele tener 1 o 2; muchos contigs indican un ensamblaje fragmentado."
T["res_size"] = "**Tamaño** — cuántas letras (pares de bases) tiene el genoma. Una bacteria típica ronda los 2 a 6 millones."
T["res_gc"] = "**GC** — porcentaje de las letras G y C. Varía mucho entre especies (25% a 70%) y sirve para saber si el genoma parece completo y de qué tipo de bacteria se trata."
T["res_cds"] = "**CDS** — número de genes que producen proteína. Un genoma de 4 millones de letras suele tener entre 3.000 y 4.000."
T["res_rna"] = "**tRNA / rRNA** — cuántos ARN de transferencia y ribosómicos se encontraron. Muy pocos tRNA puede indicar un genoma incompleto."
T["res_deserts"] = "**Desiertos** — regiones largas sin genes conocidos. No son errores: son zonas que aún no sabemos interpretar."
T["res_coverage"] = (
    "Porcentaje del genoma que quedó cubierto por alguna anotación. "
    "Un valor alto significa que entendimos la mayor parte; un valor bajo "
    "indica muchos tramos desconocidos, no que el análisis fallara."
)
T["res_ghost"] = (
    "Posibles genes en zonas sin anotar. Están ordenados por **rank**, una "
    "puntuación de 0 a 1: cuanto más alto, más señales apoyan que sea un gen real."
)
T["res_deserts_tab"] = (
    "Regiones de 3.000 letras o más sin genes conocidos. En un genoma subido no "
    "tenemos datos de origen de replicación ni de elementos móviles, así que la "
    "etiqueta lo indica en cada fila."
)
T["res_gaps"] = (
    "Huecos de la red metabólica: moléculas que la bacteria debería producir y "
    "para las que no encontramos la enzima. Cada uno se investiga por capas."
)
T["res_downloads"] = (
    "Descarga los resultados para usarlos en otro programa: el GFF3 (mapa de "
    "genes), las proteínas, los candidatos y los desiertos."
)
T["res_features"] = "Todo lo anotado, una fila por elemento."
T["res_next"] = (
    "Siguiente paso: revisa «Genes fantasma» y «Desiertos». Si algo te interesa, "
    "descarga el informe para conservar la evidencia."
)

# ---------------------------------------------------------------- catálogo
T["catalog_intro"] = (
    "Un **catálogo** es una lista ordenada. Aquí agrupamos los huecos "
    "metabólicos por tipo de ruta. Un **hueco** (o brecha) es un paso químico "
    "que la bacteria necesita pero para el que no encontramos la enzima: como "
    "una receta a la que le falta un ingrediente y no sabes cuál es."
)
T["catalog_read"] = (
    "Cada tarjeta es un paso de la ruta: 🟢 el gen está · 🔴 falta · ⚪ no hay "
    "datos. Varios 🔴 seguidos suelen significar que la bacteria usa otra ruta "
    "o que consigue la molécula del exterior."
)
T["catalog_scope"] = (
    "Mostramos solo un trozo de la ruta (hasta 8 pasos alrededor del hueco). "
    "La ruta completa puede ser mucho más larga y se vería como una maraña."
)

# -------------------------------------------------------------------- KEGG
T["kegg_intro"] = (
    "**KEGG** es una enciclopedia de rutas metabólicas con nombre propio, "
    "organizadas como mapas. Aquí vemos cada reacción del mapa y si la bacteria "
    "tiene la enzima que la realiza."
)
T["kegg_read"] = (
    "🟢 presente · 🔴 falta · ⚪ sin datos. Si falta casi todo el mapa, "
    "probablemente esa ruta no exista en esta bacteria; si faltan uno o dos "
    "pasos sueltos, pueden ser genes que no detectamos."
)
T["kegg_pages"] = (
    "Los mapas tienen decenas de reacciones. Los mostramos por bloques para que "
    "se puedan leer."
)

# ------------------------------------------------------------------ brechas
T["gaps_intro"] = (
    "Una **brecha** es una molécula que la bacteria debería poder fabricar "
    "(porque otras rutas la necesitan) pero para la que no aparece la enzima. "
    "El trabajo es explicar cada brecha: ¿falta el gen de verdad, está escondido, "
    "o la bacteria consigue la molécula del exterior? Cada brecha declara sus "
    "fuentes para que puedas juzgar."
)
T["gaps_m_total"] = "**Brechas** — cuántos huecos quedan tras aplicar los filtros. Más brechas no es peor: una bacteria que vive en un entorno rico puede importar muchas moléculas en vez de fabricarlas."
T["gaps_m_support"] = "**Soporte ≥10** — brechas cuya enzima existe en 10 o más de los 400 genomas de referencia. Un soporte alto indica que la ruta es común y que el gen probablemente esté, aunque no lo detectemos."
T["gaps_m_novel"] = "**Novel** — brechas que ningún genoma de referencia resuelve. Son las más interesantes y también las más inciertas: pueden ser una ruta rara o un error de anotación."
T["gaps_m_transport"] = "**Con transporte** — brechas que se explican porque la bacteria importa la molécula en vez de fabricarla. Deja de ser un problema."
T["gaps_m_closure"] = "**Con cierre** — brechas donde encontramos un trozo de proteína (un dominio) con la forma de la enzima perdida. Es una pista fuerte, no una prueba."
T["gaps_e1"] = "**Trazado y soporte** — de dónde sale la molécula según la base de datos Rhea y en cuántos genomas de referencia aparece una enzima que la produzca."
T["gaps_e2"] = "**Transporte** — si la bacteria no fabrica la molécula, ¿la trae de fuera? Buscamos una proteína parecida a transportadores conocidos y la reacción de transporte que le corresponde."
T["gaps_e3"] = "**Dirección y termodinámica** — una reacción química puede ir en los dos sentidos. Comprobamos si la dirección que la célula necesita es energéticamente favorable o está prohibida."
T["gaps_e4"] = "**Revisión manual** — lo que una persona revisó y decidió, con su motivo escrito. Si no hay revisión, lo dice."
T["gaps_e5"] = "**Cierre candidato** — un trozo de proteína del genoma con la forma de la enzima que falta. Se llama cierre porque «cierra» el hueco, pero es una hipótesis hasta confirmarla."
T["gaps_e6"] = "**HMM** — búsqueda de familias de proteínas por su perfil. Un resultado aquí refuerza que la enzima existe aunque la anotación normal no la haya visto."
T["gaps_support_help"] = "Número de los 400 genomas de referencia que tienen alguna enzima capaz de producir esta molécula. Cuanto más alto, más común es la ruta y menos probable que el gen esté realmente ausente."
T["gaps_review"] = (
    "Capa de revisión: una clasificación de cada brecha con la acción propuesta. "
    "No modifica la tabla base, solo añade una opinión revisable."
)
T["gaps_queue"] = (
    "Cola priorizada: qué brechas conviene revisar primero, ordenadas por "
    "impacto y evidencia."
)

# ------------------------------------------------------------------ cierres
T["closures_intro"] = (
    "Un **cierre** es una propuesta de solución para una brecha: una proteína "
    "del genoma que podría ser la enzima que falta. Cada cierre dice qué "
    "dominio tiene, en qué gen está y qué tan sólida es la evidencia."
)
T["closures_m_total"] = "**Cierres** — cuántas brechas tienen al menos una proteína candidata a ser la enzima perdida."
T["closures_m_forward"] = "**forward** — la reacción va en el sentido que la bacteria necesita y la energía lo permite. Es el caso más sólido."
T["closures_m_reversible"] = "**reversible** — la reacción puede ir en los dos sentidos; sirve, pero no garantiza el sentido que hace falta."
T["closures_m_reverse"] = "**reverse** — la reacción, tal como la proponemos, iría en sentido contrario al necesario. Hay que revisarla."
T["closures_m_unknown"] = "**unknown** — no hay dato energético. El cierre queda como hipótesis sin comprobar."
T["closures_limit"] = (
    "Un cierre `candidate` no es `confirmed`: significa «encontramos una proteína "
    "que encaja», no «demostramos que hace esa reacción». La energía se estima en "
    "condiciones estándar (pH 7, fuerza iónica 0,1 M, 25 °C)."
)
T["closures_evidence"] = (
    "Evidencia del cierre elegido: dominio encontrado, gen candidato, soporte en "
    "los 400 genomas, dirección de la reacción y resultado de la revisión manual."
)

# ----------------------------------------------------------------- operones
T["operons_intro"] = (
    "Un **operón** es un grupo de genes que se leen juntos, como una sola "
    "instrucción. Comparten un interruptor (el **promotor**) y se copian en un "
    "mismo mensaje. Esto permite a la bacteria encender varios genes a la vez. "
    "Aquí ves cada pieza que conocemos de cada operón de *E. coli* K-12."
)
T["operons_m_total"] = "**Operones** — cuántos grupos de genes se leen juntos según la base de datos RegulonDB."
T["operons_m_multi"] = "**Multigénicos** — operones con dos o más genes. Son los interesantes: los de un solo gen no explican nada en grupo."
T["operons_m_term"] = "**Con terminador** — operones donde se conoce la señal que marca el final del mensaje. Solo una parte está documentada."
T["operons_m_complete"] = "**Completos (7/7)** — operones de los que se conocen todas las piezas: promotor, inicio de transcripción, cajas de unión, factor sigma, reguladores, terminador y evidencia publicada."
T["operons_pieces"] = (
    "Cada columna es una pieza: ✓ está documentada · — no se conoce. "
    "Que falte una pieza no significa que no exista, solo que nadie la ha descrito todavía."
)
T["operons_tu"] = (
    "Una **TU** (unidad de transcripción) es un mensaje copiado desde un "
    "promotor. Un mismo operón puede tener varios promotores internos, y por eso "
    "se copia en trozos distintos según la situación."
)
T["operons_why"] = (
    "Esto importa porque al buscar operones en un genoma nuevo suponemos «un "
    "promotor por grupo», y casi la mitad de los operones conocidos tienen más "
    "de uno."
)

# ------------------------------------------------------------------- reglas
T["rules_intro"] = (
    "Estas no son predicciones estadísticas: son **reglas escritas** que "
    "deciden, entre cada dos genes, si siguen en el mismo mensaje o empiezan uno "
    "nuevo. Cada decisión muestra su motivo, así que puedes discutirla."
)
T["rules_r1"] = "**R1 · péptido líder** — si el gen anterior es muy corto (un péptido líder), la horquilla que viene detrás no separa genes: regula la producción (atenuación). No se corta."
T["rules_r1b"] = "**R1b · atenuador funcional** — si hay una horquilla fuerte pero los dos genes trabajan juntos (según STRING), se interpreta como pausa interna, no como final del grupo."
T["rules_r2"] = "**R2 · terminador** — horquilla fuerte con cola de uracilos pegada al final del gen anterior: el mensaje termina aquí."
T["rules_r3"] = "**R3 · promotor** — aparece una señal de inicio de transcripción justo antes del gen siguiente: empieza un mensaje nuevo."
T["rules_r4"] = "**R4 · distancia** — los genes están muy separados. Es la regla más simple y la que más se equivoca."
T["rules_r5"] = "**R5 · sin señal** — no hay evidencia fuerte; se unen por defecto."
T["rules_m_total"] = "**Operones** — grupos resultantes de aplicar las reglas."
T["rules_m_borders"] = "**Fronteras decididas** — pares de genes consecutivos que se evaluaron."
T["rules_m_cuts"] = "**Cortes** — cuántas de esas fronteras se decidieron como inicio de un operón nuevo."

# ---------------------------------------------------------------- horquillas
T["hairpins_intro"] = (
    "Una **horquilla** es un pliegue del ARN: la molécula se dobla sobre sí "
    "misma y forma un lazo, como un cordón de zapato que se enrosca. Ese lazo "
    "puede frenar a la máquina que copia el ADN y hacer que se detenga. "
    "Es una **pista**, no una prueba: una de cada cinco horquillas reales corta "
    "en medio de un operón (atenuación), y algunas terminaciones no usan "
    "horquilla en absoluto."
)
T["hairpins_m_total"] = "**Regiones con horquilla** — huecos entre genes donde se encontró una estructura de lazo."
T["hairpins_m_a"] = "**Clase A** — horquilla fuerte: tallo largo, rico en G y C, con cola de uracilos. Son las más parecidas a un terminador real."
T["hairpins_m_b"] = "**Clase B** — horquilla intermedia: se encuentra tanto en terminadores reales como en regiones sin función, así que no distingue."
T["hairpins_m_c"] = "**Clase C** — estructura débil. Casi siempre es ruido."
T["hairpins_class_help"] = "La clase resume la fuerza de la estructura. A es la más fiable; C la menos."
T["hairpins_dg"] = (
    "**ΔG** (energía libre de plegamiento, en kcal/mol) — cuánta energía se "
    "libera al formarse el lazo. Un valor negativo grande indica una estructura "
    "estable; cuanto más negativo, más difícil es deshacerla."
)

# ------------------------------------------------------------------ móviles
T["mobile_intro"] = (
    "Aquí se anotan los elementos que se mueven o defienden a la bacteria: "
    "**elementos móviles** (trozos de ADN que se copian y se pegan en otro "
    "sitio), **profagos** (virus bacterianos dormidos dentro del genoma), "
    "**sistemas de defensa**, **CRISPR** (el sistema inmune de la bacteria), "
    "**AMR** (genes de resistencia a antibióticos) y **virulencia** (capacidad "
    "de causar daño)."
)
T["mobile_limit"] = (
    "Estos datos vienen precalculados para los genomas del proyecto. La fuente "
    "no incluye el contig ni la hebra, así que las coordenadas son relativas al "
    "genoma completo. Para un genoma subido harían falta herramientas "
    "específicas que no están integradas."
)
T["mobile_types"] = (
    "Filtra por tipo si solo te interesa una categoría. Los totales de cada "
    "tarjeta corresponden a la selección."
)

T["mobile_typed_intro"] = (
    "Versión más detallada de la pestaña anterior, probada en 10 genomas: cada "
    "elemento lleva su tipo concreto, la coordenada dentro de su contig y un "
    "nivel de confianza."
)
T["mobile_typed_limit"] = (
    "Es un piloto: las llamadas se basan en familias de proteínas y reglas "
    "sencillas, no en las herramientas especializadas de referencia. La "
    "confianza baja significa «puede ser, revísalo»."
)
T["membrane_intro"] = (
    "La **membrana** es la piel de la bacteria: una capa de grasa que separa el "
    "interior del exterior. Aquí se anotan las proteínas que viven en ella o la "
    "atraviesan: partes hidrofóbicas (TM), señales de envío (péptido señal), "
    "lipoproteínas, transportadores y sistemas de secreción."
)
T["membrane_limit"] = (
    "Piloto: las hélices transmembrana y el péptido señal se estiman con reglas "
    "sencillas de hidrofobicidad, no con los predictores especializados. La "
    "mayoría de las filas quedan en confianza baja a propósito."
)

# ------------------------------------------------------------------ taxonomía
T["taxonomy_intro"] = (
    "**Taxonomía** es ordenar la vida por parentesco. Para saber qué especie es "
    "un genoma comparamos su ADN con genomas de referencia y calculamos el "
    "**ANI** (identidad nucleotídica media): el porcentaje de coincidencia entre "
    "dos genomas. Por encima del 95% suele ser la misma especie."
)
T["taxonomy_m_ncbi"] = "**Especie NCBI** — el nombre que da la base de datos oficial de secuencias."
T["taxonomy_m_gtdb"] = "**Especie GTDB** — el nombre según GTDB, una clasificación alternativa muy usada en microbiología."
T["taxonomy_m_ani"] = "**ANI** — parecido del ADN con la referencia más cercana. «sin ANI» significa que no había ninguna referencia comparable, no que el parecido sea 0%."
T["taxonomy_m_conf"] = "**Confianza** — cuánto coinciden las dos clasificaciones. Si discrepan, el nombre es provisional."
T["taxonomy_limit"] = (
    "Límites: el ANI se calcula con fastANI, sin el flujo completo de GTDB-Tk; "
    "en varios genomas no hay una referencia independiente; y el 95% no es una "
    "frontera universal. «Provisional» no es una verdad definitiva."
)
T["taxonomy_lineage"] = "Lista de grupos desde lo más general (dominio) hasta lo más concreto (especie)."

# ------------------------------------------------------------- explicaciones
T["family_intro"] = (
    "No todas las brechas son «falta una enzima». Algunas se explican porque la "
    "molécula entra del exterior, otras porque el gen existe pero está mapeado "
    "mal, y otras porque la ruta usa un atajo. Aquí se agrupan por familia de "
    "ruta y se propone una explicación y un siguiente paso."
)
T["family_limit"] = (
    "Las explicaciones son hipótesis de trabajo: se apoyan en la evidencia "
    "disponible y están pensadas para que una persona las revise."
)

# -------------------------------------------------------------------- mapeo
T["mapping_intro"] = (
    "A veces el gen está, pero la reacción no aparece en la red. Eso se llama "
    "**brecha de mapeo**: la enzima existe y hace el trabajo, pero nuestro "
    "sistema no la conectó con la reacción correcta. Aquí se busca el homólogo "
    "más parecido en genomas de referencia para decidir si el gen existe de "
    "verdad."
)
T["mapping_m_total"] = "**Brechas** — huecos de esta familia que se están investigando."
T["mapping_m_strong"] = "**homolog_strong** — hay una proteína muy parecida a la enzima esperada: probablemente el gen está y solo falta conectarlo."
T["mapping_m_weak"] = "**homolog_weak** — parecido parcial: puede ser la enzima o una proteína emparentada con otra función."
T["mapping_m_low"] = "**homolog_low** — parecido apenas suficiente para sugerir algo. Requiere revisión."
T["mapping_m_none"] = "**no_homolog** — no aparece ninguna proteína parecida. La enzima probablemente no está."
T["mapping_corrections"] = (
    "Correcciones propuestas y su validación: si al aplicar el mapeo la brecha "
    "desaparece de la red, la corrección era buena."
)

# --------------------------------------------------------------- completitud
T["completeness_intro"] = (
    "**Completitud** es qué fracción de las reacciones de una ruta está "
    "presente en el genoma. Si una ruta tiene 10 pasos y encontramos 7, su "
    "completitud es 0,7. No es lo mismo que esencialidad: una ruta incompleta "
    "puede funcionar si la célula importa el producto."
)
T["completeness_m_genomes"] = "**Genomas** — cuántos genomas se comparan en esta ruta."
T["completeness_m_cov"] = "**Cobertura media** — fracción media de reacciones presentes. 1,0 sería la ruta completa en todos los genomas."
T["completeness_m_gaps"] = "**Gaps** — reacciones de la ruta que faltan en al menos un genoma."
T["completeness_m_ess"] = "**Gaps esenciales** — reacciones que, si faltan, probablemente rompen la ruta entera."

# ------------------------------------------------------------------ cadenas
T["chains_intro"] = (
    "Una **cadena** es un camino de reacciones que conecta una molécula que "
    "sobra con otra que falta. Si la bacteria no puede hacer un paso, quizá "
    "pueda llegar al mismo resultado por otro camino: eso es una ruta "
    "alternativa."
)
T["chains_m_total"] = "**Cadenas** — caminos encontrados para esta familia de rutas."
T["chains_m_alt"] = "**Alternativas** — cadenas que evitan el paso que falta."
T["chains_m_noco"] = "**Sin co-sustratos** — cadenas que no necesitan moléculas auxiliares, así que son más fáciles de ejecutar para la célula."
T["chains_len"] = "**Longitud máxima** — cuántos pasos puede tener el camino. Más pasos significa más suposiciones encadenadas."
T["chains_limit"] = (
    "Las cadenas son posibilidades topológicas, no flujos reales: no dicen a qué "
    "velocidad ocurrirían ni si la célula las usa. `reversible=unknown` significa "
    "que no hay dato energético."
)

# --------------------------------------------------------------- transporte
T["transport_intro"] = (
    "Un **transportador** es una puerta en la membrana: una proteína que deja "
    "entrar o salir una molécula concreta. Si una brecha se explica porque la "
    "molécula entra de fuera, aquí está la puerta."
)
T["transport_m_t1"] = "**T1** — transportadores con reacción de transporte exacta en Rhea y sustrato conocido. Es el nivel más fiable."
T["transport_m_genomes"] = "**Genomas** — cuántos genomas distintos aportan al menos un transportador válido."
T["transport_m_high"] = "**Alta confianza** — coincidencia fuerte entre el transportador encontrado y el esperado."
T["transport_limit"] = (
    "T2 y T1_review no son importación demostrada: son parecidos que aún no "
    "prueban que la molécula entre. Solo T1 se cuenta como evidencia."
)

# ------------------------------------------------------------------- origen
T["origin_intro"] = (
    "El **origen de replicación** (oriC) es el punto del ADN donde empieza la "
    "copia del genoma. Localizarlo importa porque los genes cerca del origen se "
    "copian más veces y suelen estar más activos. Se busca con varios métodos "
    "independientes y se declara cuáles coinciden."
)
T["origin_m_state"] = "**Estado** — `validated` (dos métodos coinciden) · `consensus` (coinciden sin referencia externa) · `anomalous` (los métodos discrepan) · `unresolved` (sin candidato)."
T["origin_m_pos"] = "**Posición oriC** — número de letra del genoma donde se sitúa el origen."
T["origin_m_conf"] = "**Confianza** — de 0 a 1. Resume cuánta evidencia apoya esa posición."
T["origin_states"] = (
    "Estados: `validated` dos métodos coinciden · `consensus` coinciden sin "
    "verdad externa · `lineage_specific` el linaje no tiene el marcador habitual "
    "· `anomalous` discrepan · `unresolved` sin candidato."
)
T["origin_triage"] = (
    "Triaje de las anomalías: se añadió una comprobación extra (sesgo GC y "
    "cajas DnaA). Esa comprobación acierta el 57% de las veces en genomas ya "
    "validados, así que **no reclasifica sola**: solo deja la evidencia a la vista."
)

# ------------------------------------------------------------- rutas curadas
T["curated_intro"] = (
    "Rutas de referencia paso a paso, con las enzimas conocidas y el estado de "
    "cada una en el genoma. Sirven para entender una ruta concreta sin perderse "
    "entre cientos de reacciones."
)
T["curated_read"] = (
    "🟢 la enzima está · 🔴 falta · ⚪ sin datos. Si falta una enzima pero el "
    "producto aparece igual, probablemente la bacteria lo consigue de otra forma."
)

# ------------------------------------------------------------------- varios
T["download_tsv"] = "Descarga la tabla en formato TSV (se abre con Excel o cualquier editor de texto)."
T["filter_hint"] = "Los filtros no cambian los datos: solo eligen qué filas ves."

# ------------------------------------------------- ayudas de widgets repetidos
T["help_genoma"] = (
    "Elige de qué genoma quieres ver los datos. Los del proyecto son 400 "
    "bacterias de referencia, cada una con su análisis ya calculado."
)
T["help_estado"] = (
    "Filtra por nivel de confianza. `validated` = dos métodos independientes "
    "coinciden; `consensus` = un método con respaldo; `candidate` = una sola "
    "evidencia; `anomalous` = los métodos se contradicen."
)
T["help_tipo"] = (
    "Filtra por categoría del elemento. Deja todas marcadas para ver el "
    "conjunto completo."
)
T["help_clase"] = (
    "Fuerza de la horquilla: A es una estructura clara (se detecta bien), "
    "B es intermedia (no se distingue del azar) y C es débil."
)
T["help_ruta"] = (
    "Ruta metabólica concreta. Cada paso indica si la enzima está presente "
    "en este genoma, falta o no hay datos."
)
T["help_bloque"] = (
    "Las listas largas se muestran por páginas. Cambia el número para ver "
    "el siguiente grupo."
)
T["help_operon"] = (
    "Operón curado de *E. coli* K-12. Al elegirlo se muestran sus unidades "
    "de transcripción con el promotor y los factores que lo regulan."
)
T["help_brecha"] = (
    "Brecha concreta. Abajo se despliegan las evidencias que la explican: "
    "qué se intentó y por qué se llegó a esa conclusión."
)
