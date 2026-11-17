# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Find bundled data in source checkouts and regular package installations."""
from pathlib import Path


def bundled_hmm_library():
    package = Path(__file__).resolve().parent
    for candidate in (package.parent.parent / "hmm_library", package / "data" / "hmm_library"):
        if (candidate / "hmm_registry.tsv").is_file():
            return candidate
    return None
