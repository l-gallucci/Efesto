"""I/O helpers: FASTA, GFF, map/cutoff files, HMM catalog."""

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# ── HMM source provenance ─────────────────────────────────────────────────────

SOURCE_REFS = {
    "fegenie":  "doi:10.1038/s41396-019-0570-7",
    "tabuteau": "doi:10.1111/1462-2920.70218",
    "methmmdb": "doi:10.1101/2024.12.26.629440",
    "interpro": "doi:10.1093/nar/gkac993",
    "ncbifam":  "https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/",
    "curated":  "see registry reference column",
}

SOURCE_LABELS = {
    "fegenie":  "FeGenie (Garber et al. 2020, ISME J)",
    "tabuteau": "Tabuteau et al. 2025 (Environ Microbiol)",
    "methmmdb": "MetHMMDB (Kciuchcinski et al. 2025, bioRxiv)",
    "interpro": "InterPro (Paysan-Lafosse et al. 2023, NAR)",
    "ncbifam":  "NCBIfam (NCBI prokaryotic genome annotation)",
    "curated":  "manually curated models",
}

# ── Category catalog: user-facing token → internal HMM directory names ────────
# Tokens are used with --annotate. "all" / None means load everything.

_ANNOTATE_MAP = {
    # element-level selectors
    "Fe":        ["iron_acquisition-heme_oxygenase", "iron_acquisition-heme_transport",
                  "iron_acquisition-iron_transport", "iron_acquisition-siderophore_synthesis",
                  "iron_acquisition-siderophore_transport",
                  "iron_acquisition-siderophore_transport_potential",
                  "iron_gene_regulation", "iron_oxidation", "iron_reduction",
                  "iron_resistance", "iron_storage", "magnetosome_formation"],
    "Cu":        ["metal_resistance-copper"],
    "Zn":        ["metal_resistance-cobalt_zinc_cadmium"],
    "Mn":        ["metal_resistance-manganese"],
    "Ni":        ["metal_resistance-nickel"],
    "Co":        ["metal_resistance-cobalt_zinc_cadmium"],
    "Mo":        ["metal_resistance-molybdenum"],
    "As":        ["metal_resistance-arsenic"],
    "Hg":        ["metal_resistance-mercury"],
    "Cd":        ["metal_resistance-cobalt_zinc_cadmium"],
    "Cr":        ["metal_resistance-chromium"],
    "Ag":        ["metal_resistance-silver"],
    "Te":        ["metal_resistance-tellurite"],
    "Mg":        ["metal_resistance-magnesium"],
    "multimetal": ["metal_resistance-multimetal", "metal_resistance-non-specific"],
    # process-level selectors (Fe subcategories)
    "Fe-metabolism":  ["iron_acquisition-heme_oxygenase", "iron_acquisition-heme_transport",
                       "iron_acquisition-iron_transport", "iron_acquisition-siderophore_synthesis",
                       "iron_acquisition-siderophore_transport",
                       "iron_acquisition-siderophore_transport_potential",
                       "iron_oxidation", "iron_reduction", "iron_storage", "magnetosome_formation"],
    "Fe-acquisition": ["iron_acquisition-heme_oxygenase", "iron_acquisition-heme_transport",
                       "iron_acquisition-iron_transport", "iron_acquisition-siderophore_synthesis",
                       "iron_acquisition-siderophore_transport",
                       "iron_acquisition-siderophore_transport_potential"],
    "Fe-resistance":  ["iron_resistance"],
    "Fe-regulation":  ["iron_gene_regulation"],
    "Fe-reduction":   ["iron_reduction"],
    "Fe-oxidation":   ["iron_oxidation"],
    "Fe-storage":     ["iron_storage"],
    "Fe-magnetosome": ["magnetosome_formation"],
    "Cu-resistance":  ["metal_resistance-copper"],
    "Mn-resistance":  ["metal_resistance-manganese"],
    "Ni-resistance":  ["metal_resistance-nickel"],
    "Mo-resistance":  ["metal_resistance-molybdenum"],
    "As-resistance":  ["metal_resistance-arsenic"],
    "Hg-resistance":  ["metal_resistance-mercury"],
    "Co-Zn-Cd-resistance": ["metal_resistance-cobalt_zinc_cadmium"],
    "Cr-resistance":  ["metal_resistance-chromium"],
    "Ag-resistance":  ["metal_resistance-silver"],
    "Te-resistance":  ["metal_resistance-tellurite"],
    "Mg-resistance":  ["metal_resistance-magnesium"],
}

VALID_ANNOTATE_TOKENS = sorted(_ANNOTATE_MAP) + ["all"]


def read_fasta(path):
    seqs, header, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if header is not None:
        seqs[header] = "".join(parts)
    return seqs


def read_fasta_lengths(path):
    lengths, header, length = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = length
                header = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if header is not None:
        lengths[header] = length
    return lengths


def read_cutoffs(path):
    c = {}
    if not os.path.isfile(path):
        return c
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                try:
                    c[ls[0]] = float(ls[1])
                except ValueError:
                    pass
    return c


def read_map(path):
    m = {}
    if not os.path.isfile(path):
        return m
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                m[ls[0]] = ls[1]
    return m


def get_contig_lengths_from_gff(gff_path):
    lengths = defaultdict(int)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = re.search(r"##sequence-region\s+(\S+)\s+\d+\s+(\d+)", line)
                if m:
                    lengths[m.group(1)] = int(m.group(2))
                continue
            parts = line.rstrip().split("\t")
            if len(parts) >= 5 and parts[2] == "CDS":
                try:
                    end = int(parts[4])
                    if end > lengths[parts[0]]:
                        lengths[parts[0]] = end
                except ValueError:
                    pass
    return dict(lengths)


def build_contig_length_dict(faa_files, gff_dir=None, fna_dir=None, fna_ext="fna"):
    contig_lengths = {}
    for faa in faa_files:
        stem = faa.stem
        genome = faa.name
        found = None
        if fna_dir:
            for ext in [fna_ext, "fna", "fasta", "fa"]:
                c = Path(fna_dir) / f"{stem}.{ext}"
                if c.exists():
                    found = ("fna", str(c))
                    break
        if not found and gff_dir:
            for ext in [".gff", ".gff3", ".prodigal.gff"]:
                c = Path(gff_dir) / (stem + ext)
                if c.exists():
                    found = ("gff", str(c))
                    break
        if found:
            kind, path = found
            contig_lengths[genome] = (
                read_fasta_lengths(path) if kind == "fna"
                else get_contig_lengths_from_gff(path)
            )
    return contig_lengths


def load_prodigal_gff(gff_path):
    coords = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", parts[8])
            if m:
                coords[m.group(1).strip()] = {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
    return coords


def load_gff_dir(gff_dir, faa_files):
    gff_dir = Path(gff_dir)
    gc = {}
    for faa in faa_files:
        stem = faa.stem
        found = None
        for ext in (".gff", ".gff3", ".prodigal.gff"):
            c = gff_dir / (stem + ext)
            if c.exists():
                found = c
                break
        if found:
            gc[faa.name] = load_prodigal_gff(str(found))
        else:
            print(f"  [WARN] No GFF for {faa.name}, using index clustering",
                  file=sys.stderr)
    return gc


# ── Registry + catalog helpers ────────────────────────────────────────────────

def load_registry(hmm_dir):
    """Load hmm_registry.tsv; return list of row dicts or [] if absent."""
    path = Path(hmm_dir) / "hmm_registry.tsv"
    if not path.exists():
        return []
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def filter_categories(cat_hmms, annotate_tokens):
    """Return cat_hmms restricted to directories matching annotate_tokens.

    annotate_tokens: list of strings from --annotate, e.g. ['Fe', 'Cu-resistance'].
    Returns full cat_hmms unchanged when tokens are empty or contain 'all'.
    Unknown tokens are warned and skipped.
    """
    if not annotate_tokens or "all" in annotate_tokens:
        return cat_hmms
    keep = set()
    for tok in annotate_tokens:
        dirs = _ANNOTATE_MAP.get(tok)
        if dirs is None:
            print(f"[WARN] Unknown --annotate token '{tok}' — ignored. "
                  f"Valid tokens: {', '.join(VALID_ANNOTATE_TOKENS)}", file=sys.stderr)
            continue
        keep.update(dirs)
    filtered = {cat: hmms for cat, hmms in cat_hmms.items() if cat in keep}
    missing  = keep - set(cat_hmms)
    if missing:
        print(f"[WARN] --annotate requested categories not found in --hmm_dir: "
              f"{', '.join(sorted(missing))}", file=sys.stderr)
    return filtered


_NSEQ_WARN_THRESHOLD = 10


def build_nseq_map(registry):
    """Return {stem: int} for all registry rows that have a valid nseq value."""
    nseq = {}
    for r in registry:
        stem = r.get("stem", "")
        if not stem:
            continue
        try:
            nseq[stem] = int(float(r["nseq"]))
        except (KeyError, ValueError, TypeError):
            pass
    return nseq


def print_provenance(annotate_tokens, registry, cat_hmms):
    """Print HMM source/reference table for the active categories."""
    active_cats = set(cat_hmms)
    active_stems = {hf.stem for hmm_list in cat_hmms.values() for _, hf in hmm_list}
    active_rows = [r for r in registry if r.get("stem") in active_stems]
    if not active_rows:
        return

    source_counts = defaultdict(int)
    zero_cutoff   = defaultdict(int)
    for r in active_rows:
        src = r.get("source", "unknown")
        source_counts[src] += 1
        try:
            if float(r.get("cutoff", 1)) == 0.0:
                zero_cutoff[src] += 1
        except (ValueError, TypeError):
            pass

    cats_str = ", ".join(sorted(annotate_tokens)) if annotate_tokens else "all"
    print(f"[INFO] Annotating: {cats_str}")
    print(f"[INFO] HMM sources:")
    for src, count in sorted(source_counts.items(), key=lambda x: -x[1]):
        label = SOURCE_LABELS.get(src, src)
        ref   = SOURCE_REFS.get(src, "no reference")
        line  = f"       {src:<12}  {count:4d} models  {label}  ({ref})"
        zc = zero_cutoff.get(src, 0)
        if zc:
            line += f"  [!] {zc} with cutoff=0.0"
        print(line)

    no_ref = sum(1 for r in active_rows if not r.get("reference", "").strip())
    if no_ref:
        print(f"[WARN] {no_ref} active models lack a reference — treat hits cautiously",
              file=sys.stderr)

    low_by_src = defaultdict(list)
    for r in active_rows:
        try:
            if int(float(r.get("nseq", 100))) < _NSEQ_WARN_THRESHOLD:
                low_by_src[r.get("source", "unknown")].append(r["stem"])
        except (ValueError, TypeError):
            pass
    if low_by_src:
        total_low = sum(len(v) for v in low_by_src.values())
        src_summary = ", ".join(
            f"{src}={len(stems)}"
            for src, stems in sorted(low_by_src.items(), key=lambda x: -len(x[1]))
        )
        print(f"[WARN] {total_low} active models have nseq < {_NSEQ_WARN_THRESHOLD} "
              f"(limited training data — treat hits cautiously; "
              f"filter hmm_registry.tsv by nseq column for details): {src_summary}",
              file=sys.stderr)
