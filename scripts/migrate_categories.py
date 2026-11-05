#!/usr/bin/env python3
"""Rename library categories and apply documented curation; dry-run by default."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from efesto.categories import canonical_category


def plan_migration(library):
    library = Path(library).resolve()
    with (library / "hmm_registry.tsv").open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = list(reader.fieldnames)
        rows = list(reader)
    with (library / "category_curation.tsv").open(newline="") as fh:
        curations = {r["stem"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for field in ("category_original", "category_review_status", "category_note"):
        if field not in fields:
            fields.append(field)
    moves, changes, destinations = [], [], set()
    for row in rows:
        stem = row["stem"]
        old_cat, old_status = row["category"], row["status"]
        decision = curations.get(stem)
        new_cat = decision["category"] if decision else canonical_category(old_cat)
        if Path(new_cat).name != new_cat or new_cat in (".", ".."):
            raise ValueError(f"Category must be a single directory name: {new_cat}")
        old_rel = Path(row["hmm_file"])
        new_rel = Path(new_cat) / (stem + ".hmm")
        if old_rel.is_absolute() or ".." in old_rel.parts:
            raise ValueError(f"Non-local registry HMM path: {old_rel}")
        if new_rel in destinations:
            raise ValueError(f"Duplicate destination: {new_rel}")
        destinations.add(new_rel)
        source, target = library / old_rel, library / new_rel
        if not source.is_file():
            raise ValueError(f"Missing HMM: {source}")
        if target != source and target.exists():
            raise ValueError(f"Destination already exists: {target}")
        if target != source:
            moves.append((source, target, hashlib.sha256(source.read_bytes()).hexdigest()))
        row["category_original"] = row.get("category_original") or old_cat
        row["category"] = new_cat
        row["hmm_file"] = str(new_rel)
        if decision:
            row["category_review_status"] = decision["review_status"]
            row["category_note"] = decision["note"]
            if decision["status_override"]:
                row["status"] = decision["status_override"]
        else:
            row["category_review_status"] = row.get("category_review_status") or "not_individually_reviewed"
            row["category_note"] = row.get("category_note") or "Naming normalized; biological assignment not individually adjudicated."
        if old_cat != new_cat or old_status != row["status"]:
            changes.append({"stem": stem, "old_category": old_cat, "new_category": new_cat,
                            "old_status": old_status, "new_status": row["status"]})
    return library, fields, rows, moves, changes


def migrate(library, apply=False):
    library, fields, rows, moves, changes = plan_migration(library)
    summary = {"records": len(rows), "file_moves": len(moves),
               "runtime_enabled": sum(r["status"] in ("active", "experimental") for r in rows),
               "categories": sorted({r["category"] for r in rows}),
               "changes": changes, "applied": apply}
    if not apply:
        return summary
    # Validate the complete plan before moving anything; roll back on failure.
    done = []
    temp = library / "hmm_registry.tsv.category-migration.tmp"
    if temp.exists():
        raise ValueError(f"Migration temporary file already exists: {temp}")
    try:
        for source, target, digest in moves:
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
            done.append((source, target))
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError(f"HMM bytes changed during move: {target}")
        with temp.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        temp.replace(library / "hmm_registry.tsv")
    except Exception:
        for source, target in reversed(done):
            target.rename(source)
        if temp.exists():
            temp.unlink()
        raise
    for source, _ in done:
        if source.parent.exists() and not any(source.parent.iterdir()):
            source.parent.rmdir()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, default=Path(__file__).resolve().parents[1] / "hmm_library")
    parser.add_argument("--apply", action="store_true", help="Apply the validated migration plan")
    parser.add_argument("--manifest", type=Path, help="Write complete JSON migration plan/results")
    args = parser.parse_args()
    result = migrate(args.library, args.apply)
    if args.manifest:
        args.manifest.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{'Applied' if args.apply else 'Dry run'}: {result['records']} records, "
          f"{result['file_moves']} moves, {result['runtime_enabled']} enabled HMMs, "
          f"{len(result['categories'])} categories")


if __name__ == "__main__":
    main()
