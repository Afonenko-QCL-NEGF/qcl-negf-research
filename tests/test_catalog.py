"""Scientific input integrity and failure modes of catalog validation."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import shutil

from jsonschema import ValidationError
import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("research_validate", ROOT / "tools" / "validate.py")
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


@pytest.fixture
def catalog(tmp_path: Path) -> Path:
    shutil.copytree(ROOT / "config", tmp_path / "config")
    shutil.copytree(ROOT / "literature", tmp_path / "literature")
    return tmp_path


def test_catalog_is_closed_and_valid() -> None:
    counts = validator.validate(ROOT)
    assert counts == {"definitions": 16, "configurations": 30, "literature_files": 2}


def test_configurator_catalog_preserves_historical_candidate() -> None:
    # Independent literals: do not derive expected paths or science from the
    # configurator catalog or the Runner resolver being checked by consumers.
    catalog = json.loads((ROOT / "config/configurator-catalog.json").read_text())
    assert set(catalog) == {"schema", "templates"}
    assert catalog["schema"] == "qcl-negf-configurator-catalog-v1"
    assert len(catalog["templates"]) == 1
    template = catalog["templates"][0]
    assert set(template) == {"id", "name", "study", "model", "protocol", "description"}
    assert template["id"] == "equilibrium-0-70"
    assert template["study"] == "config/studies/equilibrium-0-70.yaml"
    for key, path in (
        ("model", "config/model/reference-2019-70k.yaml"),
        ("protocol", "config/policies/analytical-diagnostic.yaml"),
    ):
        assert set(template[key]) == {"name", "path"}
        assert template[key]["path"] == path
        assert isinstance(template[key]["name"], str) and template[key]["name"].strip()
        assert (ROOT / path).resolve().is_relative_to(ROOT.resolve())
        assert (ROOT / path).is_file()
    for key in ("name", "description"):
        assert isinstance(template[key], str) and template[key].strip()

    study = yaml.safe_load((ROOT / "config/studies/equilibrium-0-70.yaml").read_text())
    assert study["configuration"]["sources"] == [
        "../model/reference-2019-70k.yaml", "../policies/analytical-diagnostic.yaml"
    ]
    assert [
        (ROOT / template["study"]).parent.joinpath(source).resolve()
        for source in study["configuration"]["sources"]
    ] == [ROOT / template["model"]["path"], ROOT / template["protocol"]["path"]]
    assert study["axes"] == {
        "temperatures": {"values": [70], "unit": "K"},
        "voltages": {"values": [0], "unit": "mV"},
    }
    assert study["configuration"]["overrides"] == {
        "physical": {"voltage_per_period": "0 mV"},
        "solver": {"maximum_scba_iterations": 2048, "maximum_poisson_iterations": 128,
                   "self_energy_mixing": 0.2},
        "production": {"checkpoint_every_scba": 32},
    }
    assert study["variants"] == [{"id": "baseline", "method_id": "configured"}]
    assert study["policies"] == {
        "scba_to_poisson": "adaptive_working",
        "voltage": {"mode": "independent", "invalid_predecessor": "stop_branch"},
        "on_child_failure": "continue",
    }
    assert study["output"]["recovery"]["byte_budget"] == 137438953472
    model = yaml.safe_load((ROOT / "config/model/reference-2019-70k.yaml").read_text())
    policy = yaml.safe_load((ROOT / "config/policies/analytical-diagnostic.yaml").read_text())
    assert model["scattering"] == {
        "lo_phonon": True, "acoustic_phonon": True, "ionized_impurity": True,
        "interface_roughness": True, "alloy_disorder": False,
    }
    assert policy["scattering"] == {"alloy_disorder": True}
    assert "tolerances" not in policy["solver"]
    assert policy["solver"]["convergence"]["mode"] == "adaptive_working"


def test_unreferenced_fragment_is_an_error(catalog: Path) -> None:
    (catalog / "config" / "model" / "unused.yaml").write_text("run: {}\n")
    with pytest.raises(ValueError, match="unreachable"):
        validator.validate(catalog)


def test_variant_typo_cannot_bypass_schema(catalog: Path) -> None:
    path = catalog / "config/studies/grid-refinement-48.yaml"
    value = yaml.safe_load(path.read_text())
    value["variants"][0]["overrides"]["numerical"]["energy_node"] = 100
    path.write_text(yaml.safe_dump(value))
    with pytest.raises(ValidationError, match="energy_node"):
        validator.validate(catalog)


def test_ambiguous_yaml_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "ambiguous.yaml"
    path.write_text("physical: {}\nphysical: {}\n")
    with pytest.raises(ValueError, match="duplicate YAML key"):
        validator.read_mapping(path)


def test_inclusion_cycles_are_rejected(catalog: Path) -> None:
    path = catalog / "config/operators.yaml"
    value = yaml.safe_load(path.read_text())
    value["includes"].append("operators.yaml")
    path.write_text(yaml.safe_dump(value))
    with pytest.raises(ValueError, match="cyclic"):
        validator.validate(catalog)


def test_external_configuration_paths_are_rejected(catalog: Path) -> None:
    with pytest.raises(ValueError, match="escapes repository"):
        validator.resolve_member(catalog, catalog / "config", "../../outside.yaml")


def test_literature_mutation_requires_new_checksum(catalog: Path) -> None:
    path = catalog / "literature/reference-2019/negf_iv_digitized.csv"
    path.write_text(path.read_text().replace("3.432", "9.000"))
    with pytest.raises(ValueError, match="checksum mismatch"):
        validator.validate(catalog)


def test_unsupported_admission_is_rejected_by_schema(catalog: Path) -> None:
    path = catalog / "config/studies/reference-48.yaml"
    value = yaml.safe_load(path.read_text())
    value["admission"] = {"required_evidence": ["unscheduled-gate"]}
    path.write_text(yaml.safe_dump(value))
    with pytest.raises(ValidationError, match="admission"):
        validator.validate(catalog)


def test_reference_geometry_retains_literature_dimensions() -> None:
    model = validator.read_mapping(ROOT / "config/model/reference-2019-70k.yaml")
    with (ROOT / "literature/reference-2019/structure.csv").open() as stream:
        records = {row["quantity"]: row for row in csv.DictReader(stream)}
    quantities = ["barrier_1", "well_1", "barrier_2", "well_2", "doped_well_segment", "well_3_remainder"]
    for layer, quantity in zip(model["physical"]["layers"], quantities, strict=True):
        thickness, unit = layer["thickness"].split()
        assert unit == "nm"
        assert records[quantity]["unit"] == "angstrom"
        assert float(thickness) * 10 == pytest.approx(float(records[quantity]["value"]))
    assert model["physical"]["lattice_temperature"] == "70 K"
    assert float(records["published_operating_temperature"]["value"]) == 200
