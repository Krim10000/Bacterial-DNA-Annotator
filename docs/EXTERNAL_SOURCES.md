# Fuentes externas: DBTBS y BATTER-TPE

Fecha: 2026-11-21. Dos encargos independientes: una segunda especie curada y un
segundo detector de terminadores. Ninguno de los dos salió como estaba previsto,
y en ambos casos la conclusión se apoya en lo que respondió cada servidor.

---

## Resumen

| encargo | resultado | qué quedó |
|---|---|---|
| **DBTBS** (*B. subtilis*) | **bloqueado por diseño**: descarga detrás de un formulario de registro + `robots.txt` prohíbe `/cgi-bin/` | integrada **BSGatlas** como segunda especie curada (incluye DBTBS entre sus fuentes) |
| **BATTER-TPE** | **funciona en CPU**, con dos parches de compatibilidad | 100 % de sensibilidad sobre los curados de la región evaluada; no se versiona por falta de licencia |

---

## 1 · DBTBS

### Qué se probó y qué respondió

| petición | respuesta |
|---|---|
| `GET https://dbtbs.hgc.jp/` | `200`, **326 bytes** — solo un `<frame>`, sin datos |
| `GET https://dbtbs.hgc.jp/ver4/` | `200`, **291 bytes** — otro frame |
| `GET https://dbtbs.hgc.jp/download.html` | `200`, **2.583 bytes** — formulario |
| `GET https://dbtbs.hgc.jp/ver4/download.html` | `200`, 2.583 bytes (misma página) |
| `GET https://dbtbs.hgc.jp/robots.txt` | `200`, 34 bytes: `User-Agent: *` / `Disallow: /cgi-bin/` |
| `GET /sitemap.xml`, `/data/`, `/files/` | `404` (236 bytes) |
| `HEAD /cgi-bin/download.cgi` | `200`, `text/html` |

El formulario pide **Name\*, E-mail\*, Institution\*, Purpose\*** y su JavaScript
hace:

```js
frm.action = "/cgi-bin/download.cgi";
frm.target = "frame2";
frm.submit();
```

### Decisión

No se automatiza. Rellenar el formulario con datos inventados sería falsificar un
registro de uso, y `robots.txt` pide explícitamente no automatizar `/cgi-bin/`.
En lugar de forzarlo, se buscó una fuente abierta equivalente.

### Alternativa integrada: BSGatlas

**BSGatlas** (Geissler et al.) es la anotación integrada de *B. subtilis* 168 y
declara entre sus fuentes **DBTBS**, BsubCyc, Rfam, RefSeq y SubtiWiki. Se
descarga en GFF3 **sin registro**.

- Web: <https://rth.dk/resources/bsgatlas/download.php>
- DOI: <https://doi.org/10.5281/zenodo.4305869> (v1.1)
- Zenodo responde `403` ("unusual traffic") desde esta máquina; la web original
  (`rth.dk`) sí sirve los ficheros.

Ficheros descargados en `data/dbtbs/`:

| fichero | tamaño |
|---|---|
| `BSGatlas_v1.1.gff` | 7,1 MB |
| `genome.fna` | 4,3 MB |

### Lo que se integró

`scripts/validation/fetch_bsgatlas.py` — descargador + parser, idempotente,
con `--offline` para solo reparsar. Genera, en `data/dbtbs/`:

| salida | contenido |
|---|---|
| `gt_operons.tsv` | **2.266 operones**, 867 con ≥2 genes, 4.349 genes distintos; mismo esquema que el GT de RegulonDB |
| `gt_terminators.tsv` | **2.566 terminadores**, todos con energía libre |
| `gt_tss.tsv` | **3.390 TSS**, todos con factor sigma |

Dos detalles del GFF que costaron un intento cada uno:

1. Hay dos features `gene` por gen (dos fuentes) → se deduplica por coordenada.
2. El atributo `locus_tag` trae dos alias separados por coma
   (`locus_tag=BSU_00010,BSU00010`) → se toma el primero.

Además se dejó `scripts/validation/fetch_dbtbs.py`, que documenta el bloqueo y
convierte un fichero descargado a mano (`data/dbtbs/dbtbs_manual.tsv`) al mismo
esquema. Es el camino para el día que alguien complete el formulario.

---

## 2 · BATTER-TPE

### Disponibilidad

- Código: <https://github.com/lulab/BATTER> (32 MB, última actualización
  2026-02-16)
- Datos y pesos: <https://zenodo.org/records/18863501> (v3)
- Preprint: <https://doi.org/10.1101/2023.10.02.560326>

**Licencia: no declarada.** La API de GitHub devuelve `license: None` y el
repositorio no incluye `LICENSE`. Sin licencia los derechos quedan reservados:
por eso **no se versiona** y el script lo descarga a `data/terminators/`
(ignorado por git). Si el autor publica una licencia, hay que revisar esto.

### Viabilidad en esta máquina (sin GPU)

| requisito | estado |
|---|---|
| `torch` | ✅ 2.14.1+cpu |
| `transformers` | ✅ 4.46.3 |
| `pandas`, `numpy` | ✅ |
| `pyfaidx` | ⚠️ faltaba — instalado (0.9.2, dependencia pura Python) |
| GPU | no hace falta: el modelo son 32 MB |

### Dos incompatibilidades, dos parches

El modelo se entrenó con `transformers 4.18` y un torch antiguo:

1. **`load_state_dict` estricto** → `encountered unexpected key
   "encoder.embeddings.position_ids"`. Las versiones nuevas de transformers
   añaden ese buffer. Se resuelve con `strict=False`.
2. **`.view()` sobre tensores no contiguos** en la decodificación CRF →
   `view size is not compatible with input tensor's size and stride`. Se
   resuelve cambiando a `.reshape()` (6 apariciones).

Ambos parches los aplica `scripts/regulation/fetch_batter.py` sobre la
copia descargada, así que son reproducibles.

### Coste

~20 kb/min en un núcleo de CPU. Un genoma de 4,6 Mb tarda ~2,5 h; por eso la
evaluación se hizo sobre regiones, no sobre el genoma completo.

### Resultado

`scripts/regulation/compare_batter.py` compara contra
`data/regulondb/TerminatorSet.tsv`:

Región evaluada: primeros **300 kb** de K-12 (23 terminadores curados, 140
predicciones). Tolerancia ±25 pb.

| | curados detectados | predicciones |
|---|---|---|
| BATTER-TPE | **21/23 (91,3 %)** | 140 |
| horquillas propias, clase A | 6/23 (26 %) | 57 |
| horquillas propias, A + B | 19/23 (83 %) | 154 |
| **unión de ambos** | **23/23 (100 %)** | — |

Por tipo, BATTER acierta 20/22 rho-independientes y el único rho-dependiente
de la región. La primera región de 100 kb había dado 9/9 (100 %), consistente.

**Lo importante no es el 91 %, es la complementariedad**: hay 2 terminadores que
BATTER no ve y las horquillas sí, y 4 que las horquillas no ven y BATTER sí.
Juntos recuperan los 23. Es exactamente el argumento para tener dos detectores
independientes en vez de afinar uno solo.

El suelo de precisión de BATTER es 21,4 % (30/140). No es un mal número: con
RegulonDB como referencia incompleta, cualquier detector que encuentre
terminadores reales no curados baja artificialmente.

### Qué mide y qué no

- **Sensibilidad** es limpia: los curados son el denominador.
- **Precisión** es una **cota inferior**, porque RegulonDB no es exhaustivo:
  muchas predicciones "fuera" pueden ser terminadores reales no curados.

---

## 3 · Mapa de la integración

```
scripts/validation/fetch_bsgatlas.py    BSGatlas -> data/dbtbs/gt_*.tsv
scripts/validation/fetch_dbtbs.py       documenta el bloqueo + conversión manual
scripts/regulation/fetch_batter.py   descarga BATTER + parches (sin versionar)
scripts/regulation/compare_batter.py BATTER vs curados (sensibilidad/precisión)

data/dbtbs/          BSGatlas (GFF + GT extraído)      [ignorado por git]
data/terminators/    BATTER + predicciones             [ignorado por git]
```

Se añadieron `data/dbtbs/` y `data/terminators/` a `.gitignore` por tamaño y por
licencia.

---

## 4 · Pendiente y recomendación

1. **DBTBS directo**: solo si alguien completa el formulario una vez. El
   conversor ya está listo; el fichero va en `data/dbtbs/dbtbs_manual.tsv`.
2. **BSGatlas como validación**: es la pieza más valiosa que salió de aquí. Con
   2.266 operones y 867 multi-génicos se puede repetir el benchmark de operones
   fuera de *E. coli*, que es justo lo que faltaba para saber si el método
   generaliza o está sobreajustado a K-12.
3. **BATTER a escala**: si se quiere el genoma completo, hay que paralelizar por
   regiones (el modelo es independiente por ventana). 8 núcleos → ~20 min por
   genoma.
4. **Licencia de BATTER**: pedir al autor que declare una antes de integrarlo en
   la app o en cualquier distribución.

**Siguiente paso recomendado**: correr el benchmark de operones contra BSGatlas.
Es barato (los datos ya están), es la pregunta científica que importa
(¿generaliza fuera de K-12?) y no depende de nada pendiente.

## Actualización: DBTBS obtenido por las tablas públicas

El formulario de registro responde con una redirección a
`dbtbs_050922.xml.gz` (release 5, año 2005), pero esa URL **devuelve 404** en el
servidor de DBTBS. Comprobado con y sin `Referer`, y en varias rutas
(`/ver4/download/`, `/download/`, `/dbtbs/`): el archivo plano ya no existe.

Lo que sí funciona son las tablas HTML públicas, que `robots.txt` permite
(solo prohíbe `/cgi-bin/`):

| página | contenido | resultado |
|---|---|---|
| `promtable.html` | gen → operón, elementos reguladores | 2.030 genes, **1.093 operones** (410 multi-génicos) |
| `tfactable.html` | factores y sus sitios | no parseable con la heurística actual (0 filas) |

Mapeo a coordenadas contra el GFF de BSGatlas (misma cepa, `basu168`):
**1.852 de 2.030 genes (91%)**.

- Descargador: `scripts/validation/fetch_dbtbs_tables.py` (idempotente, 2 páginas).
- Salida: `data/dbtbs/gt_operons_dbtbs.tsv` (la tabla de factores no se pudo parsear; el fichero vacío se descartó).
- Créditos: entrada **DBTBS** en `webapp/credits.py` con la cita de Sierro et al.
  2008 y la nota de que no se redistribuye.
- Registro: se hizo con nombre, afiliación y propósito, como pide el sitio.

## ESMFold: servicio externo, no Modal

**Desde 2026-11-21, el plegamiento con ESMFold se hace vía neurosnap.ai**
(https://neurosnap.ai), no en Modal.

- La app `dna-esmfold` de Modal está **parada** y el script
  `scripts/structure/modal_esmfold.py` queda marcado como deprecado.
- Motivo: no consumir cuota de GPU (10 simultáneas por workspace) en una tarea
  que se resuelve con un servicio externo.
- El conector antiguo por API pública (`scripts/structure/26_esmfold_connect.py`,
  api.esmatlas.com) sigue como referencia; neurosnap es la vía vigente.
