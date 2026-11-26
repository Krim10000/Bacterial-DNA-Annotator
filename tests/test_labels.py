from webapp.labels import labels_for_text, is_regulator


def test_heat_and_flagellum_labels():
    assert "heat_shock" in labels_for_text("DnaK chaperone heat shock protein")
    assert "flagellum" in labels_for_text("flagellar hook protein FlgE")
    assert "nitrogen" in labels_for_text("nitrogen regulatory protein GlnA")


def test_sigma_regulator_is_not_its_own_target():
    assert is_regulator("rpoH", "RNA polymerase sigma factor RpoH")
    assert is_regulator("fliA", "RNA polymerase sigma factor FliA")
    assert not is_regulator("dnaK", "chaperone DnaK")
