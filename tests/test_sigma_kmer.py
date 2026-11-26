import json
from pathlib import Path
from webapp.sigma_kmer import KmerModel, load_model, save_model


def test_kmer_model_roundtrip_and_score(tmp_path):
    model = KmerModel(
        genome="x", labels=("sigma70",), vocab={"aaaa": 0},
        scale=(1.0,), coef=((1.0,),), intercept=(0.0,), ngram=4,
    )
    p = model.probabilities("AAAA")
    assert 0 < p["sigma70"] < 1
    path = tmp_path / "m.json"
    save_model(model, path)
    again = load_model(path)
    assert again.probabilities("AAAA") == p


def test_kmer_model_rejects_bad_schema():
    import pytest
    with pytest.raises(ValueError):
        KmerModel.from_dict({"schema": "otro"})


def test_bundled_subtilis_model_is_usable():
    import pytest
    path = Path("webapp/sigma_models/subtilis_kmer.json")
    if not path.exists():
        pytest.skip("sin modelo incluido")
    model = load_model(path)
    probs = model.probabilities("TATAAT" + "A" * 30)
    assert set(probs) == set(model.labels)
    assert all(0.0 <= v <= 1.0 for v in probs.values())
