#!/usr/bin/env python3
"""Validate the study graph, resolved configurations and literature provenance."""

from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from jsonschema import ValidationError
from qcl_negf_contracts import schema_validator
import yaml


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject ambiguous YAML mappings instead of accepting the last key."""


def _mapping(loader: UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False) -> dict:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise ValueError(f"duplicate YAML key: {key!r}")
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def read_mapping(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        data = yaml.load(stream, Loader=UniqueKeyLoader)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a mapping")
    return data


def merge(base: dict, patch: dict) -> dict:
    """Apply recursive mapping overrides; sequences are replaced as a whole."""
    result = deepcopy(base)
    for key, value in patch.items():
        result[key] = (
            merge(result[key], value)
            if isinstance(result.get(key), dict) and isinstance(value, dict)
            else deepcopy(value)
        )
    return result


def resolve_member(root: Path, parent: Path, relative: str) -> Path:
    path = (parent / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"configuration reference escapes repository: {relative}")
    if not path.is_file():
        raise ValueError(f"missing configuration source: {path}")
    return path


def validate(root: Path) -> dict[str, int]:
    root = root.resolve()
    definition_validator = schema_validator("scientific-definition.schema.json")
    run_validator = schema_validator("run.schema.json")
    definitions: dict[Path, dict] = {}
    ids: dict[str, Path] = {}
    reachable: set[Path] = set()
    active: set[Path] = set()
    configuration_count = 0

    def visit(path: Path) -> None:
        nonlocal configuration_count
        if path in active:
            raise ValueError(f"cyclic study inclusion: {path.relative_to(root)}")
        if path in definitions:
            return
        data = read_mapping(path)
        definition_validator.validate(data)
        if data["id"] in ids:
            raise ValueError(f"duplicate study id: {data['id']}")
        ids[data["id"]] = path
        reachable.add(path)
        active.add(path)
        if data["kind"] == "meta":
            for child in data["includes"]:
                if not isinstance(child, str):
                    raise ValueError("catalog includes must be explicit unmodified study paths")
                visit(resolve_member(root, path.parent, child))
        else:
            configuration = data.get("configuration", {})
            sources = configuration.get("sources", [])
            if not sources:
                raise ValueError(f"{path}: explicit configuration.sources are required")
            resolved: dict = {}
            for source in sources:
                source_path = resolve_member(root, path.parent, source)
                resolved = merge(resolved, read_mapping(source_path))
                reachable.add(source_path)
            resolved = merge(resolved, configuration.get("overrides", {}))
            for variant in data.get("variants", [{"id": "baseline"}]):
                run_validator.validate(merge(resolved, variant.get("overrides", {})))
                configuration_count += 1
        active.remove(path)
        definitions[path] = data

    roots = sorted((root / "config").glob("*.yaml"))
    if not roots:
        raise ValueError("no top-level scientific definitions")
    for path in roots:
        visit(path.resolve())
    orphaned = set((root / "config").rglob("*.yaml")) - reachable
    if orphaned:
        raise ValueError(f"unreachable configuration files: {sorted(str(p.relative_to(root)) for p in orphaned)}")

    data_files = 0
    for source in sorted((root / "literature").glob("*/source.yaml")):
        provenance = read_mapping(source)
        if not provenance.get("publication", {}).get("doi"):
            raise ValueError(f"{source}: missing publication DOI")
        listed = provenance.get("data_files", {})
        if not listed:
            raise ValueError(f"{source}: no literature data files")
        for name, record in listed.items():
            path = resolve_member(root, source.parent, name)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if record.get("sha256") != digest:
                raise ValueError(f"literature data checksum mismatch: {name}")
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
                raise ValueError(f"empty or malformed literature CSV: {name}")
            data_files += 1
        actual = {path.name for path in source.parent.glob("*.csv")}
        if actual != set(listed):
            raise ValueError(f"untracked literature CSV files: {sorted(actual ^ set(listed))}")
    if not data_files:
        raise ValueError("no literature data")
    return {"definitions": len(definitions), "configurations": configuration_count, "literature_files": data_files}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    arguments = parser.parse_args()
    try:
        summary = validate(arguments.root)
    except (OSError, ValueError, TypeError, KeyError, yaml.YAMLError, ValidationError) as error:
        parser.exit(1, f"Validation failed: {error}\n")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
