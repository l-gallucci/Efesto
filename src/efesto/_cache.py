# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Content-bound completion manifests for generated files."""

import hashlib
import json
from pathlib import Path
import tempfile


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cache_matches(manifest, inputs, outputs):
    try:
        record = json.loads(Path(manifest).read_text())
        return (isinstance(record, dict) and record.get("inputs") == inputs
                and record.get("outputs") == {
                    name: file_digest(path) for name, path in outputs.items()})
    except (OSError, ValueError, TypeError):
        return False


def write_manifest(manifest, inputs, outputs):
    """Publish the completion marker last; missing markers force recomputation."""
    manifest = Path(manifest)
    record = {"inputs": inputs,
              "outputs": {name: file_digest(path) for name, path in outputs.items()}}
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=manifest.parent,
                                         prefix=manifest.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(record, handle, sort_keys=True)
            handle.write("\n")
        temporary.replace(manifest)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
