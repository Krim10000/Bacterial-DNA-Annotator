# 04_network_trace.py — traza metabólica por genoma (v3.3)

## Qué hace

Cruza la red metabólica de un genoma (`data/metabolism/{genome}_network.json`,
generada por tu pipeline) contra TODA Rhea y clasifica cada metabolito del
grafo dirigido Sustrato→Producto:

- **📥 Fuentes**: consumidos sin productor conocido en Rhea → importables
  (nutrientes que entran por transporte/PTS, sustratos de daño ambiental...).
- **🧩 Brechas reales**: consumidos sin productor en el genoma pero con
  productor en Rhea → falta enzima. Cada brecha recibe una **tría** en 4
  cubetas accionables (v3.3, en consola y en el campo `triage` del reporte):

| Cubeta (`triage`) | Significado | Acción |
| --- | --- | --- |
| `closable_with_genome_enzyme` | la enzima YA está en la red; el hueco es de cableado/voto | revisar la dirección curada de esa reacción |
| `generic_class` | el metabolito es él mismo una clase genérica de Rhea ("an alcohol", ...) | artefacto probable; ignorar o instanciar |
| `generic_producer` | es instancia de una clase genérica (UTP de "a nucleoside triphosphate") o algún productor lo da vía una clase genérica | se cierra con inferencia is-a de ChEBI o `--resolve-generics` |
| `genuine_missing_enzyme` | sin productor específico en Rhea (prephenate, N-formyl-L-kynurenine) | candidata real de enzima faltante |

- **📤 Sinks** (producidos y nunca consumidos) y **⚠️ sinks con consumidores
  en Rhea** (posible brecha aguas abajo; no cuenta la reversa de la propia
  familia que lo produjo).

## Dirección (v3.1+)

Rhea es neutral en dirección (maestro UN + variantes LR/RL/BI). El cableado
usa el voto de `rhea2uniprot.tsv`: los ids LR/RL del genoma se cablean tal
cual; los UN/BI según el voto de las proteínas DEL GENOMA si la red contiene
ids UniProt (auto-detectados), o el voto global si no. Voto `{LR}` = sentido
natural, `{RL}` = invertido, vacío/mixto = reversible. Sin `rhea2uniprot.tsv`
todo se cuenta reversible y salen cientos de brechas/sinks falsos.

## Requisitos (data/metabolism/tsv/)

```
rhea.rdf.gz            # o su cache: rhea_participants.tsv + rhea_equations.tsv
chebiId_name.tsv
rhea2uniprot.tsv       # voto de dirección (¡imprescindible!)
```

```
curl -o data/metabolism/tsv/rhea.rdf.gz \
    https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz
curl -o data/metabolism/tsv/rhea2uniprot.tsv \
    https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv
```

## Uso

```
python scripts/metabolism/04_network_trace.py --genome GCF_000006985.1_ASM698v1
# modo batch (todas las redes del directorio; el RDF y el voto se cargan una vez):
python scripts/metabolism/04_network_trace.py
```

## Flags

| Flag | Efecto |
| --- | --- |
| `--genome ID` | una red; si se omite, batch sobre todos los `*_network.json` |
| `--top N` | ejemplos por categoría (10) |
| `--resolve-generics` | expandir clases genéricas de nucleótidos en el grafo Y el índice |
| `--all-reversible` | diagnóstico: todo bidireccional (debe dar 0 brechas) |
| `--no-reversible` | modo v2 (un sentido canónico): reproduce la patología vieja |
| `--keep-currency` | no filtrar ATP/H2O/NAD+/CoA... del grafo |
| `--force-parse` | re-parsear `rhea.rdf.gz` aunque exista cache |
| `--no-report` | no escribir el JSON de reporte |
| `--rhea-dir` / `--net-dir` | directorios alternativos |

## Reporte

`data/metabolism/{genome}_network_trace.json`: stats de la red, listas
completas de fuentes/brechas/sinks/brechas aguas abajo, y por brecha el
`triage`, `closable_in_genome`, `producers_total` y `producing_rhea`. En
consola, las familias productoras se muestran con su ecuación y una nota
cuando se cuenta el sentido inverso (familia reversible o voto RL).

La comparación entre genomas la hace `05_trace_compare.py`.
