from webapp.sigma_inventory import inventory_from_rows, is_sigma_candidate, summary


def test_sigma_inventory_is_species_specific():
    rows = [
        {"gene_name": "rpoD", "product": "RNA polymerase sigma factor RpoD"},
        {"gene_name": "rpoN", "product": "RNA polymerase sigma factor RpoN"},
        {"gene_name": "sigB", "product": "RNA polymerase sigma factor SigB"},
        {"gene_name": "rseA", "product": "anti-sigma-E factor RseA"},
        {"gene_name": "rpoA", "product": "RNA polymerase subunit alpha"},
        {"gene_name": "x", "product": "sigma(54)-dependent regulator"},
    ]
    inv = inventory_from_rows("sp1", rows)
    assert {f.native_name for f in inv} == {"rpoD", "rpoN", "sigB"}
    assert summary(inv)["n_factors"] == 3
    assert summary(inv)["families"]["sigma54_like"] == 1


def test_does_not_force_names():
    assert is_sigma_candidate("sigZ", "")
    assert is_sigma_candidate("", "RNA polymerase sigma factor X")
    assert not is_sigma_candidate("rseA", "anti-sigma factor")
