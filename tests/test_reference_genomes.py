from pathlib import Path
import gzip
from webapp.reference_genomes import count_gff, load_manifest, compare_genome


def test_manifest():
    rows = load_manifest(Path('data'))
    assert len(rows) >= 30
    assert any(r['id'] == 'subtilis' for r in rows)


def test_count_gff(tmp_path):
    p = tmp_path / 'a.gff3'
    p.write_text('##gff-version 3\nc\ts\tCDS\t1\t9\t.\t+\t0\tID=a\nc\ts\ttRNA\t20\t30\t.\t-\t0\tID=b\n')
    assert count_gff(p) == {'CDS': 1, 'tRNA': 1}


def test_compare_available_subtilis():
    import pytest
    if not Path('results/subtilis/subtilis.summary.json').exists():
        pytest.skip('sin resultados de subtilis en esta instalación')
    c = compare_genome('subtilis', Path('data'), Path('results'))
    assert c['status'] == 'ok'
    assert {r['tipo'] for r in c['rows']} == {'CDS', 'tRNA', 'rRNA'}


def test_compare_planned_not_available(tmp_path):
    (tmp_path / 'reference_genomes.json').write_text('[{"id":"x","status":"planned"}]')
    c = compare_genome('x', tmp_path, tmp_path)
    assert c['status'] == 'not_available'
