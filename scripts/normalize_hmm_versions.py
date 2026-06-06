#!/usr/bin/env python3
"""
normalize_hmm_versions.py
=========================
Convert all HMM profiles in the MetalGenie-Evo library to the current
HMMER3/f format using hmmconvert.

Run once after adding new models, or whenever mixed-version warnings appear.

Usage:
    python scripts/normalize_hmm_versions.py [--hmm_dir hmm_library/] [--dry_run]
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

CURRENT_TAG = "HMMER3/f"


def hmm_version_tag(path):
    """Return the HMMER tag line from a profile file (first non-empty line)."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                return line
    return ""


def normalize_library(hmm_dir, dry_run=False):
    hmm_dir = Path(hmm_dir)
    all_hmms = sorted(hmm_dir.rglob("*.hmm"))

    if not all_hmms:
        print(f"[ERROR] No .hmm files found in {hmm_dir}")
        sys.exit(1)

    needs_convert = []
    already_current = 0
    for hmm in all_hmms:
        tag = hmm_version_tag(hmm)
        if tag.startswith(CURRENT_TAG):
            already_current += 1
        else:
            needs_convert.append((hmm, tag))

    print(f"[INFO] {len(all_hmms)} HMM files scanned")
    print(f"       {already_current} already HMMER3/f — no conversion needed")
    print(f"       {len(needs_convert)} need conversion")

    if not needs_convert:
        print("[INFO] Library already normalized.")
        return

    from collections import Counter
    tag_counts = Counter(tag for _, tag in needs_convert)
    for tag, n in sorted(tag_counts.items()):
        print(f"       {n:4d}  {tag}")

    if dry_run:
        print("[DRY RUN] No files modified.")
        return

    converted, failed = 0, []
    for hmm, old_tag in needs_convert:
        with tempfile.NamedTemporaryFile(suffix=".hmm", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        r = subprocess.run(
            ["hmmconvert", str(hmm)],
            capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            tmp_path.write_text(r.stdout)
            hmm.write_text(r.stdout)
            tmp_path.unlink(missing_ok=True)
            converted += 1
        else:
            tmp_path.unlink(missing_ok=True)
            failed.append((hmm, r.stderr.strip()))

    print(f"[INFO] Converted: {converted}  Failed: {len(failed)}")
    if failed:
        print("[WARN] Conversion failures:")
        for hmm, err in failed:
            print(f"       {hmm.name}: {err[:80]}")
    else:
        print("[INFO] All models now HMMER3/f format.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hmm_dir", default="hmm_library",
                   help="HMM library root directory (default: hmm_library/)")
    p.add_argument("--dry_run", action="store_true",
                   help="Report what would be converted without modifying files")
    args = p.parse_args()
    normalize_library(args.hmm_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
