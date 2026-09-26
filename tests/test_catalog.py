"""Scientific input integrity and failure modes of catalog validation."""

from __future__ import annotations

import csv
import importlib.util
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
    assert counts == {"definitions": 14, "configurations": 29, "literature_files": 2}


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
