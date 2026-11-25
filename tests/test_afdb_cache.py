"""Pruebas sintéticas de la caché AFDB/PAE.

No usan red: ``urllib.request.urlopen`` se reemplaza por un doble que sirve
respuestas en memoria y registra las URLs solicitadas.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import urllib.error
from pathlib import Path

from scripts.structure import fetch_afdb_cache as cache

PDB_TEXT = (
    "HEADER    TEST PROTEIN\n"
    "SEQRES   1 A    4  MET ALA LEU GLY\n"
    "ATOM      1  CA  MET A     1       0.000   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)
PDB_BYTES = PDB_TEXT.encode("ascii")
PAE_DATA = [{
    "predicted_aligned_error": [[0.0, 1.5], [1.5, 0.0]],
    "max_predicted_aligned_error": 1.5,
}]
PAE_BYTES = json.dumps(PAE_DATA).encode("utf-8")


def _meta(acc: str, sequence: str = "MALG") -> list[dict]:
    return [{
        "uniprotAccession": acc,
        "uniprotSequence": sequence,
        "pdbUrl": f"https://afdb.test/files/{acc}.pdb",
        "paeDocUrl": f"https://afdb.test/files/{acc}.json",
    }]


def _http_error(code: int, headers: dict | None = None):
    def _raise():
        raise urllib.error.HTTPError(
            "https://afdb.test/error", code, f"HTTP {code}",
            headers or {}, io.BytesIO(b""),
        )
    return _raise


class FakeResponse:
    def __init__(self, body: bytes):
        self.body = body

    def read(self) -> bytes:
        return self.body

    def __enter__(self) -> "FakeResponse":
        return self

    def __exit__(self, *exc_info) -> bool:
        return False


class FakeURL:
    """Doble de ``urlopen``: responde por URL y registra cada llamada."""

    def __init__(self, mapping: dict):
        self.mapping = {
            str(url): (list(value) if isinstance(value, list) else value)
            for url, value in mapping.items()
        }
        self.calls: list[str] = []

    def __call__(self, url, timeout=None):
        url = str(url)
        self.calls.append(url)
        if url not in self.mapping:
            raise AssertionError(f"URL inesperada: {url}")
        entry = self.mapping[url]
        if isinstance(entry, list):
            if not entry:
                raise AssertionError(f"Sin respuestas restantes para: {url}")
            entry = entry.pop(0)
        if callable(entry):
            return entry()
        if isinstance(entry, BaseException):
            raise entry
        if isinstance(entry, (bytes, bytearray)):
            return FakeResponse(bytes(entry))
        return entry


def _install(monkeypatch, mapping: dict) -> FakeURL:
    fake = FakeURL(mapping)
    monkeypatch.setattr(cache.urllib.request, "urlopen", fake)
    return fake


def _targets(tmp_path: Path) -> Path:
    path = tmp_path / "targets"
    path.mkdir()
    return path


def _write_genes(data_root: Path, accessions: list[str]) -> None:
    folder = data_root / "regulation" / "g1"
    folder.mkdir(parents=True)
    with (folder / "genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome", "gene_id", "acc"])
        for acc in accessions:
            writer.writerow(["g1", f"gene_{acc}", acc])


def _catalog_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Compresión y lectura transparente
# ---------------------------------------------------------------------------
def test_compression_roundtrip_for_pdb_and_pae(tmp_path):
    targets = _targets(tmp_path)

    pdb = cache.write_cached(targets / "A.pdb", PDB_BYTES, compress=True)
    pae = cache.write_cached(targets / "A.pae.json", PAE_BYTES, compress=True)

    assert pdb.name == "A.pdb.gz"
    assert pae.name == "A.pae.json.gz"
    assert pdb.read_bytes()[:2] == b"\x1f\x8b"
    assert gzip.decompress(pdb.read_bytes()) == PDB_BYTES
    assert cache.read_cached_bytes(pdb) == PDB_BYTES
    assert cache.read_cached_text(pdb) == PDB_TEXT
    assert json.loads(cache.read_cached_text(pae)) == PAE_DATA

    # Los archivos planos siguen funcionando igual que antes.
    plain = cache.write_cached(targets / "B.pdb", PDB_BYTES, compress=False)
    assert plain.name == "B.pdb"
    assert cache.read_cached_bytes(plain) == PDB_BYTES


def test_cached_path_finds_compressed_and_prefers_plain(tmp_path):
    targets = _targets(tmp_path)
    (targets / "A.pdb").write_bytes(PDB_BYTES)
    (targets / "A.pae.json.gz").write_bytes(gzip.compress(PAE_BYTES))
    (targets / "B.pdb.gz").write_bytes(gzip.compress(PDB_BYTES))

    assert cache.cached_path(targets, "A", "pdb") == targets / "A.pdb"
    assert cache.cached_path(targets, "A", "pae") == targets / "A.pae.json.gz"
    assert cache.cached_path(targets, "B", "pdb") == targets / "B.pdb.gz"
    assert cache.cached_path(targets, "Z", "pdb") is None
    assert cache.sequence_from_pdb(PDB_TEXT) == "MALG"
    assert cache.sequence_sha256("MALG") == hashlib.sha256(b"MALG").hexdigest()


# ---------------------------------------------------------------------------
# Catálogo: compatibilidad, deduplicación y campos nuevos
# ---------------------------------------------------------------------------
def test_catalog_legacy_roundtrip_keeps_new_fields(tmp_path):
    path = tmp_path / "catalog.tsv"
    path.write_text(
        "accession\tstatus\tpdb\tpae\tsha256\tpdb_bytes\tpae_bytes\t"
        "pdb_url\tpae_url\terror\tupdated\n"
        "P1\tok\t/x/P1.pdb\t\tdeadbeef\t10\t0\t\t\t\t2026-01-01T00:00:00\n",
        encoding="utf-8",
    )
    catalog = cache.load_catalog(path)
    row = catalog["P1"]
    assert row["compressed"] == "false"
    assert row["duplicate_of"] == ""
    assert row["sequence_sha256"] == ""
    assert row["source"] == ""

    row["duplicate_of"] = "P0"
    row["sequence_sha256"] = cache.sequence_sha256("MALG")
    row["source"] = cache.SOURCE_AFDB
    row["compressed"] = "true"
    cache.write_catalog(path, catalog)

    written = _catalog_rows(path)
    assert written[0]["accession"] == "P1"
    assert written[0]["duplicate_of"] == "P0"
    assert written[0]["source"] == "alphafold_db"
    assert written[0]["compressed"] == "true"
    assert set(cache.CAT_COLS).issubset(set(written[0]))


def test_accessions_deduplicates_by_accession(tmp_path):
    data_root = tmp_path / "data"
    folder = data_root / "regulation" / "g1"
    folder.mkdir(parents=True)
    with (folder / "genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["acc"])
        writer.writerow(["P2"])
        writer.writerow(["P1"])
        writer.writerow(["P2"])
        writer.writerow([""])
    (data_root / "regulation" / "g2").mkdir()
    (data_root / "regulation" / "g2" / "genes.tsv").write_text(
        "acc\nP1\n", encoding="utf-8")

    assert cache.accessions(data_root) == ["P1", "P2"]


def test_apply_dedup_by_sequence_hash_prefers_local_model():
    digest = cache.sequence_sha256("MALG")
    rows = {
        "A": {"accession": "A", "status": "error", "pdb": "",
              "sequence_sha256": digest},
        "B": {"accession": "B", "status": "ok", "pdb": "/x/B.pdb",
              "sequence_sha256": digest},
        "C": {"accession": "C", "status": "ok", "pdb": "/x/C.pdb",
              "sequence_sha256": cache.sequence_sha256("MALW")},
    }
    cache.apply_dedup(rows)
    assert rows["B"]["duplicate_of"] == ""
    assert rows["A"]["duplicate_of"] == "B"
    assert rows["C"]["duplicate_of"] == ""

    rows["A"]["sequence_sha256"] = cache.sequence_sha256("MALK")
    cache.apply_dedup(rows)
    assert rows["A"]["duplicate_of"] == ""
    assert rows["B"]["duplicate_of"] == ""


# ---------------------------------------------------------------------------
# Compresión en la descarga completa
# ---------------------------------------------------------------------------
def test_fetch_compressed_download_marks_catalog(tmp_path, monkeypatch):
    acc = "P12345"
    targets = _targets(tmp_path)
    meta = _meta(acc)
    fake = _install(monkeypatch, {
        cache.API.format(acc=acc): json.dumps(meta).encode(),
        meta[0]["pdbUrl"]: PDB_BYTES,
        meta[0]["paeDocUrl"]: PAE_BYTES,
    })

    row = cache.fetch(acc, targets, {}, timeout=1, retries=2, want_pae=True,
                      delay=0, compress=True, backoff_jitter=0)

    assert row["status"] == "ok"
    assert row["source"] == cache.SOURCE_AFDB
    assert row["compressed"] == "true"
    assert Path(row["pdb"]).name == f"{acc}.pdb.gz"
    assert Path(row["pae"]).name == f"{acc}.pae.json.gz"
    assert row["sha256"] == hashlib.sha256(PDB_BYTES).hexdigest()
    assert row["sequence_sha256"] == cache.sequence_sha256("MALG")
    assert cache.read_cached_bytes(row["pdb"]) == PDB_BYTES
    assert json.loads(cache.read_cached_text(row["pae"])) == PAE_DATA
    assert fake.calls == [
        cache.API.format(acc=acc), meta[0]["pdbUrl"], meta[0]["paeDocUrl"],
    ]


# ---------------------------------------------------------------------------
# Reanudación y compatibilidad con archivos existentes
# ---------------------------------------------------------------------------
def test_existing_uncompressed_is_not_rewritten_or_downloaded(
        tmp_path, monkeypatch):
    acc = "P00001"
    targets = _targets(tmp_path)
    existing = targets / f"{acc}.pdb"
    existing.write_bytes(PDB_BYTES)
    fake = _install(monkeypatch, {})  # cualquier llamada falla

    row = cache.fetch(acc, targets, {}, timeout=1, retries=2, delay=0,
                      compress=True, backoff_jitter=0)

    assert row["status"] == "already_present"
    assert row["compressed"] == "false"
    assert row["sha256"] == hashlib.sha256(PDB_BYTES).hexdigest()
    assert row["sequence_sha256"] == cache.sequence_sha256("MALG")
    assert existing.read_bytes() == PDB_BYTES
    assert not (targets / f"{acc}.pdb.gz").exists()
    assert fake.calls == []


def test_existing_compressed_is_reused_without_network(tmp_path, monkeypatch):
    acc = "P00003"
    targets = _targets(tmp_path)
    pdb = cache.write_cached(targets / f"{acc}.pdb", PDB_BYTES, compress=True)
    fake = _install(monkeypatch, {})

    row = cache.fetch(acc, targets, {}, timeout=1, retries=2, delay=0,
                      compress=True, backoff_jitter=0)

    assert row["status"] == "already_present"
    assert row["compressed"] == "true"
    assert Path(row["pdb"]) == pdb
    assert row["sequence_sha256"] == cache.sequence_sha256("MALG")
    assert fake.calls == []


def test_resume_only_fetches_missing_pae_with_known_url(tmp_path, monkeypatch):
    acc = "P00002"
    targets = _targets(tmp_path)
    pdb = targets / f"{acc}.pdb"
    pdb.write_bytes(PDB_BYTES)
    pae_url = f"https://afdb.test/files/{acc}.json"
    fake = _install(monkeypatch, {pae_url: PAE_BYTES})

    row = cache.fetch(
        acc, targets,
        {acc: {"accession": acc, "pdb": str(pdb), "pae_url": pae_url,
               "sha256": hashlib.sha256(PDB_BYTES).hexdigest()}},
        timeout=1, retries=2, want_pae=True, delay=0, backoff_jitter=0,
    )

    assert row["status"] == "ok"
    assert row["pae"].endswith(f"{acc}.pae.json")
    assert pdb.read_bytes() == PDB_BYTES
    assert fake.calls == [pae_url]


# ---------------------------------------------------------------------------
# 404 permanente
# ---------------------------------------------------------------------------
def test_404_is_permanent_and_not_retried(tmp_path, monkeypatch):
    acc = "P04040"
    targets = _targets(tmp_path)
    fake = _install(monkeypatch, {cache.API.format(acc=acc): _http_error(404)})

    row = cache.fetch(acc, targets, {}, timeout=1, retries=4, delay=0,
                      backoff_jitter=0)

    assert row["status"] == "not_found"
    assert "HTTP 404" in row["error"]
    assert row["error"].endswith("HTTP 404: HTTP 404")
    assert cache.is_permanent(row) is True
    assert fake.calls == [cache.API.format(acc=acc)]
    assert not (targets / f"{acc}.pdb").exists()

    # Un segundo intento con el catálogo ya marcado no vuelve a la API.
    row2 = cache.fetch(acc, targets,
                       {acc: row}, timeout=1, retries=4, delay=0)
    assert row2["status"] == "not_found"
    assert fake.calls == [cache.API.format(acc=acc)]


def test_404_on_pae_keeps_the_pdb_already_downloaded(tmp_path, monkeypatch):
    acc = "P04041"
    targets = _targets(tmp_path)
    meta = _meta(acc)
    _install(monkeypatch, {
        cache.API.format(acc=acc): json.dumps(meta).encode(),
        meta[0]["pdbUrl"]: PDB_BYTES,
        meta[0]["paeDocUrl"]: _http_error(404),
    })

    row = cache.fetch(acc, targets, {}, timeout=1, retries=2, want_pae=True,
                      delay=0, backoff_jitter=0)

    assert row["status"] == "not_found"
    assert row["pdb"].endswith(f"{acc}.pdb")
    assert Path(row["pdb"]).read_bytes() == PDB_BYTES
    assert row["pae"] == ""
    assert row["compressed"] == "false"


# ---------------------------------------------------------------------------
# 429: Retry-After, backoff y estado temporal
# ---------------------------------------------------------------------------
def test_429_honors_retry_after_and_then_succeeds(tmp_path, monkeypatch):
    acc = "P00429"
    targets = _targets(tmp_path)
    meta = _meta(acc)
    sleeps: list[float] = []
    monkeypatch.setattr(cache.time, "sleep", lambda seconds: sleeps.append(seconds))
    fake = _install(monkeypatch, {
        cache.API.format(acc=acc): [
            _http_error(429, {"Retry-After": "1"}),
            json.dumps(meta).encode(),
        ],
        meta[0]["pdbUrl"]: PDB_BYTES,
    })

    row = cache.fetch(acc, targets, {}, timeout=1, retries=3, delay=0,
                      backoff_base=2, backoff_jitter=0)

    assert row["status"] == "ok"
    assert fake.calls.count(cache.API.format(acc=acc)) == 2
    assert sleeps == [1.0]


def test_429_exhausted_is_rate_limited_not_permanent(tmp_path, monkeypatch):
    acc = "P00430"
    targets = _targets(tmp_path)
    sleeps: list[float] = []
    monkeypatch.setattr(cache.time, "sleep", lambda seconds: sleeps.append(seconds))
    fake = _install(monkeypatch, {
        cache.API.format(acc=acc): [
            _http_error(429), _http_error(429), _http_error(429),
        ],
    })

    row = cache.fetch(acc, targets, {}, timeout=1, retries=3, delay=0,
                      backoff_base=2, backoff_max=60, backoff_jitter=0)

    assert row["status"] == "rate_limited"
    assert "HTTP 429" in row["error"]
    assert cache.is_permanent(row) is False
    assert len(fake.calls) == 3
    assert sleeps == [2.0, 4.0]
    assert not (targets / f"{acc}.pdb").exists()


# ---------------------------------------------------------------------------
# main: límites, reanudación, 404 y dry-run
# ---------------------------------------------------------------------------
def test_main_limit_and_resume(tmp_path, monkeypatch):
    accs = ["P00001", "P00002", "P00003"]
    data_root = tmp_path / "data"
    _write_genes(data_root, accs)
    targets = tmp_path / "targets"
    catalog = tmp_path / "catalog.tsv"
    mapping = {}
    for acc in accs:
        meta = _meta(acc)
        mapping[cache.API.format(acc=acc)] = json.dumps(meta).encode()
        mapping[meta[0]["pdbUrl"]] = PDB_BYTES
    fake = _install(monkeypatch, mapping)
    common = ["--data-root", str(data_root), "--targets", str(targets),
              "--catalog", str(catalog), "--workers", "1", "--delay", "0",
              "--retries", "2", "--compress"]

    assert cache.main(common + ["--limit", "2"]) == 0
    rows = _catalog_rows(catalog)
    assert len(rows) == 2
    assert [r["status"] for r in rows] == ["ok", "ok"]
    assert all(r["compressed"] == "true" for r in rows)
    # La deduplicación por secuencia deja un canónico y un duplicado.
    assert sum(1 for r in rows if r["duplicate_of"]) == 1
    # Dos accesiones: API + PDB para cada una, sin repetir descargas.
    assert len(fake.calls) == 4
    first_calls = len(fake.calls)

    assert cache.main(common + ["--limit", "10"]) == 0
    rows = _catalog_rows(catalog)
    assert len(rows) == 3
    statuses = sorted(r["status"] for r in rows)
    assert statuses == ["already_present", "already_present", "ok"]
    # El tercero es el único pendiente; los dos previos no se re-descargan.
    assert len(fake.calls) == first_calls + 2
    assert (targets / "P00003.pdb.gz").exists()


def test_main_404_permanent_across_runs(tmp_path, monkeypatch):
    acc = "P04040"
    data_root = tmp_path / "data"
    _write_genes(data_root, [acc])
    targets = tmp_path / "targets"
    catalog = tmp_path / "catalog.tsv"
    fake = _install(monkeypatch, {cache.API.format(acc=acc): _http_error(404)})
    argv = ["--data-root", str(data_root), "--targets", str(targets),
            "--catalog", str(catalog), "--workers", "1", "--delay", "0",
            "--retries", "2", "--limit", "10"]

    assert cache.main(argv) == 0
    assert _catalog_rows(catalog)[0]["status"] == "not_found"
    assert len(fake.calls) == 1

    assert cache.main(argv) == 0
    assert len(fake.calls) == 1
    assert _catalog_rows(catalog)[0]["status"] == "not_found"


def test_main_dry_run_makes_no_network_calls_or_writes(tmp_path, monkeypatch):
    data_root = tmp_path / "data"
    _write_genes(data_root, ["P00001"])
    targets = tmp_path / "targets"
    catalog = tmp_path / "catalog.tsv"
    fake = _install(monkeypatch, {})

    assert cache.main(["--data-root", str(data_root), "--targets", str(targets),
                       "--catalog", str(catalog), "--limit", "10",
                       "--dry-run"]) == 0

    assert fake.calls == []
    assert not catalog.exists()
    assert not targets.exists()
