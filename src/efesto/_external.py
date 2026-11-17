# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Resolve Python tool scripts consistently with the environment check."""

from pathlib import Path
import shutil


def python_script_command(script):
    candidate = Path(script).expanduser()
    if not candidate.is_file():
        located = shutil.which(str(script))
        if located is None:
            raise FileNotFoundError(f"Python tool script not found: {script}")
        candidate = Path(located)
    # Preserve the interpreter selected by the active environment. Passing a
    # bare script name to Python would search the working directory, not PATH.
    return [shutil.which("python3") or "python3", str(candidate.resolve())]
