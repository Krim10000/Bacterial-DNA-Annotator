from pathlib import Path
from webapp.ncrna_layer import _parse_tblout, available

def test_parse_tblout():
    p=Path('/tmp/ncrna_test.tbl')
    p.write_text("# header\ncontig1 - model - - - - 100 200 + no no 0.4 0.0 42.0 1e-10 ! desc\n")
    rows=_parse_tblout(p,'g')
    assert len(rows)==1
    assert rows[0]['start']==100 and rows[0]['end']==200
    assert rows[0]['family']=='model'

def test_available_without_panel(tmp_path):
    ok, why=available(tmp_path)
    assert not ok
    assert 'ncRNA' in why or 'cmsearch' in why
