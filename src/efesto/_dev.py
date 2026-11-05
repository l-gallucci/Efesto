"""Private benchmark adapters. Never imported by the normal command line."""

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from efesto.io import read_cutoffs
from efesto.categories import canonical_category


def add_dev_arguments(parser):
    """The caller registers these only with EFESTO_DEV_MODE=1."""
    parser.add_argument("--dev-fegenie-hmms", dest="_dev_fegenie_hmms",
                        type=Path, help=argparse.SUPPRESS)


def load_fegenie_catalog(directory):
    """Use original profiles and FeGenie's effective integer -T thresholds.

    FeGenie looks up the part of a filename preceding its FIRST dot, skips
    missing metadata, and formats hmmsearch -T with %d. Its two MtrB files
    differ only in DATE and share a tblout filename. Search that model once
    to avoid concurrent writers. Preserve both paths in the manifest.
    Efesto still performs its own downstream interpretation.
    """
    directory = Path(directory).expanduser().resolve()
    cutoff_file = directory / "HMM-bitcutoffs.txt"
    if not cutoff_file.is_file() or not (directory / "FeGenie-map.txt").is_file():
        raise ValueError("Expected original FeGenie hmms/iron directory with "
                         "HMM-bitcutoffs.txt and FeGenie-map.txt")
    metadata = read_cutoffs(cutoff_file)
    if not metadata:
        raise ValueError("FeGenie cutoff metadata is empty or invalid")
    categories, cutoffs, registry = defaultdict(list), {}, []
    selected, skipped, seen = [], [], {}
    files = sorted(f for d in directory.iterdir()
                   if d.is_dir() and not d.name.startswith((".", "_"))
                   for f in d.glob("*.hmm"))
    for profile in files:
        key = profile.name.split(".")[0]
        if key not in metadata:
            skipped.append(str(profile.relative_to(directory)))
            continue
        declared = metadata[key]
        if not math.isfinite(declared) or declared < 1:
            raise ValueError(f"Invalid positive FeGenie cutoff for {profile.name}")
        raw = profile.read_bytes()
        # DATE is non-scoring provenance, the only difference in upstream MtrB.
        scoring_signature = b"\n".join(
            line for line in raw.splitlines() if not line.startswith(b"DATE "))
        record = {"path": str(profile.relative_to(directory)),
                  "sha256": hashlib.sha256(raw).hexdigest(),
                  "source_category": profile.parent.name}
        if profile.stem in seen:
            previous, signature = seen[profile.stem]
            if scoring_signature != signature:
                raise ValueError(f"Different FeGenie models share stem {profile.stem}")
            previous["source_files"].append(record)
            continue
        category = canonical_category(profile.parent.name)
        effective = int(declared)  # Match upstream's -T %d, not a recalibration.
        nseq = next((line.split()[1] for line in raw.decode().splitlines()
                     if line.startswith("NSEQ ")), "")
        item = {"stem": profile.stem, "category": category,
                "declared_cutoff": declared, "effective_search_cutoff": effective,
                "source_files": [record]}
        seen[profile.stem] = item, scoring_signature
        selected.append(item)
        categories[category].append((profile.stem, profile))
        cutoffs[profile.stem] = effective
        registry.append({"stem": profile.stem, "source": "fegenie",
                         "status": "active", "category": category,
                         "cutoff": str(effective), "nseq": nseq})
    if not selected:
        raise ValueError("No FeGenie HMM files have usable cutoff metadata")
    manifest = {"mode": "developer_original_fegenie_library",
                "library": str(directory), "model_count": len(selected),
                "source_hmm_paths": len(files),
                "metadata_matched_paths": sum(len(x["source_files"]) for x in selected),
                "cutoff_file_sha256": hashlib.sha256(cutoff_file.read_bytes()).hexdigest(),
                "models": selected, "skipped_missing_metadata": skipped,
                "notes": ["No HMM profile files are modified",
                          "Integer -T reproduces FeGenie command formatting",
                          "Duplicate DATE-only variants searched once",
                          "Canonical category is first sorted source path",
                          "Efesto E-value, best-hit and context rules still apply",
                          "FeGenie built-in reference BLAST is not added to Efesto"]}
    return categories, cutoffs, registry, manifest


def write_manifest(path, manifest, args):
    manifest = dict(manifest, max_evalue=args.max_evalue,
                    zero_cutoff_min_bitscore=args.zero_cutoff_min_bitscore,
                    threads=args.threads, hmm_threads=args.hmm_threads,
                    all_results=args.all_results)
    Path(path).write_text(json.dumps(manifest, indent=2) + "\n")
