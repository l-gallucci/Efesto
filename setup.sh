#!/usr/bin/env bash
# Efesto setup — verifies the installation
# Not required if you used: conda env create -f environment.yml
# (the environment.yml already runs pip install -e . automatically)
set -e

echo "========================================"
echo "  Efesto  –  setup"
echo "========================================"

# ── Check conda environment is active ────────────────────────────────────────
if [[ -z "$CONDA_PREFIX" ]]; then
    echo "[ERROR] No conda environment is active."
    echo ""
    echo "  Please run:"
    echo "    conda env create -f environment.yml"
    echo "    conda activate efesto"
    exit 1
fi
echo "[INFO] Active environment : $CONDA_DEFAULT_ENV"

# ── Check required tools ──────────────────────────────────────────────────────
MISSING=()
for tool in python3 hmmsearch; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
# pyrodigal and pyrodigal-gv are Python packages — verify via import
python3 -c "import pyrodigal, pyrodigal_gv" 2>/dev/null || MISSING+=("pyrodigal/pyrodigal-gv")
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "[ERROR] Missing tools: ${MISSING[*]}"
    echo "        Run: conda env create -f environment.yml"
    exit 1
fi
echo "[INFO] All required tools found"

# ── Install package if not already installed ──────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v Efesto &>/dev/null; then
    echo "[INFO] Installing Efesto..."
    pip install -e "$REPO_DIR" --quiet
fi

# ── Verify ────────────────────────────────────────────────────────────────────
if command -v Efesto &>/dev/null; then
    echo "[INFO] Efesto installed successfully"
else
    echo "[ERROR] Installation failed"
    exit 1
fi

echo ""
echo "========================================"
echo "  Setup complete!"
echo ""
echo "  Usage:"
echo "    conda activate $CONDA_DEFAULT_ENV"
echo "    Efesto --help"
echo "========================================"
