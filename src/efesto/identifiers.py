# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Genome-scoped protein identities for combined tool inputs and exports."""

import csv
from itertools import chain
from collections import defaultdict
from urllib.parse import quote, unquote


def scoped_id(genome, orf):
    """Encode both parts reversibly, including literal colons and percent signs."""
    return quote(genome, safe="") + "::" + quote(orf, safe="")


def export_ids(rows):
    """Preserve existing labels if unique; otherwise scope every exported ID."""
    labels = {}
    collision = False
    for row in rows:
        key = (row["genome"], row["orf"])
        label = row.get("bakta_id") or row["orf"]
        if label in labels and labels[label] != key:
            collision = True
        labels[label] = key
    return {(r["genome"], r["orf"]): (
        scoped_id(r["genome"], r["orf"]) if collision
        else r.get("bakta_id") or r["orf"]) for r in rows}


class ProteinIdentifiers:
    def __init__(self, sequences):
        self.sequences = sequences

    def write_fasta(self, path, rows):
        seen = set()
        with open(path, "w") as handle:
            for row in rows:
                key = (row["genome"], row["orf"])
                seq = self.sequences.get(key[0], {}).get(key[1], "")
                if key in seen or not seq:
                    continue
                seen.add(key)
                handle.write(f">{scoped_id(*key)}\n{seq}\n")
        return len(seen)

    def group_hits(self, hits, required_rows, allow_raw=False):
        """Resolve query IDs; reject ambiguous legacy IDs rather than cross-assign."""
        required = {(r["genome"], r["orf"]) for r in required_rows}
        raw_needed = {orf for _, orf in required}
        owners = defaultdict(list)
        if allow_raw:
            # Check all input proteins, including ones without HMM hits.
            for genome, seqs in self.sequences.items():
                for orf in raw_needed.intersection(seqs):
                    owners[orf].append((genome, orf))
        grouped = defaultdict(dict)
        for query, hit in hits.items():
            matches = set()
            if "::" in query:
                parts = tuple(unquote(p) for p in query.split("::", 1))
                if scoped_id(*parts) == query and parts[1] in self.sequences.get(parts[0], {}):
                    matches.add(parts)
            if allow_raw:
                matches.update(owners.get(query, ()))
            if not matches.intersection(required):
                continue
            if len(matches) != 1:
                raise ValueError(
                    f"Ambiguous annotation query ID {query!r} occurs in multiple genomes. "
                    "Use genome-scoped query IDs from --export_flagged_faa; "
                    "see Efesto-protein-id-map.tsv.")
            genome, orf = next(iter(matches))
            grouped[genome][orf] = hit
        return dict(grouped)


def write_id_map(path, candidate_rows, final_rows):
    exported = export_ids(final_rows)
    seen = set()
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["query_id", "export_id", "genome", "orf", "bakta_id", "contig"])
        # Final rows first so mapped Bakta IDs take precedence over candidates.
        for row in chain(final_rows, candidate_rows):
            key = (row["genome"], row["orf"])
            if key in seen:
                continue
            seen.add(key)
            writer.writerow([scoped_id(*key), exported.get(key, ""), *key,
                             row.get("bakta_id", ""), row.get("contig", "")])
