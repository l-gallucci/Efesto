"""I/O helpers: FASTA, GFF, map/cutoff files, HMM catalog."""

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

from efesto.categories import CANONICAL_CATEGORIES, canonical_category

# ── HMM source provenance ─────────────────────────────────────────────────────

SOURCE_REFS = {
    "fegenie":       "doi:10.3389/fmicb.2020.00037",
    "tabuteau":      "doi:10.1111/1462-2920.70218",
    "methmmdb":      "doi:10.1101/2024.12.26.629440",
    "interpro":      "doi:10.1093/nar/gkac993",
    "ncbifam":       "https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/",
    "mnoxgenetool":  "doi:10.1021/acs.est.5c01235",
    "curated":       "see registry reference column",
}

SOURCE_LABELS = {
    "fegenie":       "FeGenie (Garber et al. 2020, Front Microbiol)",
    "tabuteau":      "Tabuteau et al. 2025 (Environ Microbiol)",
    "methmmdb":      "MetHMMDB (Kciuchcinski et al. 2025, bioRxiv)",
    "interpro":      "InterPro (Paysan-Lafosse et al. 2023, NAR)",
    "ncbifam":       "NCBIfam (NCBI prokaryotic genome annotation)",
    "mnoxgenetool":  "MnOxGeneTool (Wang et al. 2025, Environ Sci Technol)",
    "curated":       "manually curated models",
}

# ── Category catalog: user-facing token → internal HMM directory names ────────
# Tokens are used with --annotate. "all" / None means load everything.

_ANNOTATE_MAP = {
    # element-level selectors
    "Fe":        ["Heme-Oxygenase", "Heme-Transport",
                  "IronAcquisition-Transport", "Siderophore-Synthesis",
                  "Siderophore-Transport",
                  "Siderophore-TransportPotential",
                  "Iron-Regulation", "Iron-Oxidation", "Iron-Reduction",
                  "Iron-Resistance", "Iron-Storage", "Magnetosome-Formation"],
    "Cu":        ["MetalResistance-Copper"],
    "Zn":        ["MetalResistance-CobaltZincCadmium"],
    "Mn":        ["Manganese-Acquisition", "Manganese-Oxidation", "Manganese-Regulation",
                  "MetalResistance-Manganese"],
    "Ni":        ["MetalResistance-Nickel"],
    "Co":        ["MetalResistance-CobaltZincCadmium"],
    "Mo":        ["MetalResistance-Molybdenum"],
    "As":        ["MetalResistance-Arsenic"],
    "Hg":        ["MetalResistance-Mercury"],
    "Cd":        ["MetalResistance-CobaltZincCadmium"],
    "Cr":        ["MetalResistance-Chromium"],
    "Ag":        ["MetalResistance-Silver"],
    "Te":        ["MetalResistance-Tellurite"],
    "Mg":        ["MetalResistance-Magnesium"],
    "multimetal": ["MetalResistance-MultiMetal", "MetalResistance-NonSpecific"],
    # process-level selectors (Fe subcategories)
    "Fe-metabolism":  ["Heme-Oxygenase", "Heme-Transport",
                       "IronAcquisition-Transport", "Siderophore-Synthesis",
                       "Siderophore-Transport",
                       "Siderophore-TransportPotential",
                       "Iron-Oxidation", "Iron-Reduction", "Iron-Storage", "Magnetosome-Formation"],
    "Fe-acquisition": ["Heme-Oxygenase", "Heme-Transport",
                       "IronAcquisition-Transport", "Siderophore-Synthesis",
                       "Siderophore-Transport",
                       "Siderophore-TransportPotential"],
    "Fe-resistance":  ["Iron-Resistance"],
    "Fe-regulation":  ["Iron-Regulation"],
    "Fe-reduction":   ["Iron-Reduction"],
    "Fe-oxidation":   ["Iron-Oxidation"],
    "Fe-storage":     ["Iron-Storage"],
    "Fe-magnetosome": ["Magnetosome-Formation"],
    "Cu-resistance":    ["MetalResistance-Copper"],
    "Mn-metabolism":    ["Manganese-Acquisition", "Manganese-Oxidation", "Manganese-Regulation",
                         "MetalResistance-Manganese"],
    "Mn-acquisition":   ["Manganese-Acquisition"],
    "Mn-oxidation":     ["Manganese-Oxidation"],
    "Mn-regulation":    ["Manganese-Regulation"],
    "Mn-resistance":    ["MetalResistance-Manganese"],
    "Ni-resistance":  ["MetalResistance-Nickel"],
    "Mo-resistance":  ["MetalResistance-Molybdenum"],
    "As-resistance":  ["MetalResistance-Arsenic"],
    "Hg-resistance":  ["MetalResistance-Mercury"],
    "Co-Zn-Cd-resistance": ["MetalResistance-CobaltZincCadmium"],
    "Cr-resistance":  ["MetalResistance-Chromium"],
    "Ag-resistance":  ["MetalResistance-Silver"],
    "Te-resistance":  ["MetalResistance-Tellurite"],
    "Mg-resistance":  ["MetalResistance-Magnesium"],
}

# Function-specific selectors distinguish acquisition/homeostasis from resistance.
_ANNOTATE_MAP.update({
    "Cu-acquisition": ["MetalAcquisition-Copper"],
    "Cu-homeostasis": ["MetalHomeostasis-Copper"],
    "Zn-acquisition": ["MetalAcquisition-Zinc", "MetalAcquisition-MultiMetal"],
    "Ni-acquisition": ["MetalAcquisition-Nickel", "MetalAcquisition-CobaltNickel"],
    "Co-acquisition": ["MetalAcquisition-CobaltNickel"],
    "Mo-acquisition": ["MetalAcquisition-Molybdenum"],
    "Mg-acquisition": ["MetalAcquisition-Magnesium"],
    "As-oxidation": ["Arsenic-Oxidation"],
    "Fe-stress": ["Iron-StressResponse"],
    "Fe-sulfur": ["IronSulfur-Assembly", "IronSulfur-Regulation"],
})
_ANNOTATE_MAP["Fe-acquisition"] += ["IronAcquisition-ComplexTransport", "MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Fe-metabolism"] += ["IronAcquisition-ComplexTransport", "IronSulfur-Assembly", "IronSulfur-Regulation", "Iron-StressResponse", "MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal", "MetalHomeostasis-IronZinc"]
_ANNOTATE_MAP["Fe"] += ["IronAcquisition-ComplexTransport", "IronSulfur-Assembly", "IronSulfur-Regulation", "Iron-StressResponse", "MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal", "MetalHomeostasis-IronZinc"]
_ANNOTATE_MAP["Cu"] += ["MetalAcquisition-Copper", "MetalHomeostasis-Copper"]
_ANNOTATE_MAP["Zn"] += ["MetalAcquisition-Zinc", "MetalAcquisition-MultiMetal", "MetalRegulation-Zinc", "MetalHomeostasis-IronZinc"]
_ANNOTATE_MAP["Ni"] += ["MetalAcquisition-Nickel", "MetalAcquisition-CobaltNickel", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Co"] += ["MetalAcquisition-CobaltNickel", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Mn"] += ["MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Mn-metabolism"] += ["MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Mn-acquisition"] += ["MetalAcquisition-IronManganese", "MetalAcquisition-MultiMetal"]
_ANNOTATE_MAP["Mo"] += ["MetalAcquisition-Molybdenum"]
_ANNOTATE_MAP["Mg"] += ["MetalAcquisition-Magnesium"]
_ANNOTATE_MAP["As"] += ["Arsenic-Oxidation"]
_ANNOTATE_MAP["multimetal"] += ["MetalAcquisition-MultiMetal", "MetalHomeostasis-MultiMetal", "MetalRegulation-MultiMetal", "Unclassified-MetalFunction"]
_ANNOTATE_MAP["Fe-regulation"] += ["IronSulfur-Regulation", "MetalRegulation-MultiMetal"]
for _element in ("Fe", "Mn", "Ni", "Co", "Zn"):
    _ANNOTATE_MAP[_element] += ["MetalRegulation-MultiMetal"]
_ANNOTATE_MAP.update({cat: [cat] for cat in CANONICAL_CATEGORIES})
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
            row["category"] = canonical_category(row.get("category", ""))
            rows.append(row)
    return rows


def filter_categories(cat_hmms, annotate_tokens):
    """Return cat_hmms restricted to directories matching annotate_tokens.

    annotate_tokens: list of strings from --annotate, e.g. ['Fe', 'Cu-resistance'].
    Returns full cat_hmms unchanged when tokens are empty or contain 'all'.
    Unknown tokens are warned and skipped.
    """
    normalized = defaultdict(list)
    for category, profiles in cat_hmms.items():
        normalized[canonical_category(category)].extend(profiles)
    cat_hmms = dict(normalized)
    if not annotate_tokens or "all" in annotate_tokens:
        return cat_hmms
    keep = set()
    for tok in annotate_tokens:
        dirs = _ANNOTATE_MAP.get(tok) or _ANNOTATE_MAP.get(canonical_category(tok))
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


def build_confirmation_map(registry):
    """
    Return {stem: {"gene_name": str, "kegg_ko": [str, ...]}} for every registry
    row with needs_confirmation == "yes".

    A stem with no kegg_ko value still gets an entry (empty ko list) — the
    eggNOG/Baktfold escalation can still run, it just won't have a KO-based
    confirm/contradict signal available and stays neutral for that part.
    kegg_ko may be a single ID or a comma-separated list (matches eggNOG-mapper's
    own KEGG_ko column convention of multiple IDs per hit).
    """
    flagged = {}
    for r in registry:
        stem = r.get("stem", "")
        if not stem or (r.get("needs_confirmation") or "").strip().lower() != "yes":
            continue
        ko_raw = (r.get("kegg_ko") or "").strip()
        kos = [k.strip() for k in ko_raw.split(",") if k.strip()]
        flagged[stem] = {"gene_name": r.get("gene_name", ""), "kegg_ko": kos}
    return flagged


def build_topology_expectation_map(registry):
    """
    Return {stem: {"expected_heme_c": int|None, "expected_topology": str|None}}
    for every registry row that sets expected_heme_c and/or expected_topology.

    Both fields are independently optional per row — a stem can set only one
    and stay neutral for the other, same no-data-no-penalty policy as
    build_confirmation_map. A row with neither field set is excluded entirely;
    topo_weight() never looks at it, so an unset stem is always neutral.
    """
    expectations = {}
    for r in registry:
        stem = r.get("stem", "")
        if not stem:
            continue
        heme_raw = (r.get("expected_heme_c") or "").strip()
        topo_raw = (r.get("expected_topology") or "").strip()
        if not heme_raw and not topo_raw:
            continue
        expected_heme_c = None
        if heme_raw:
            try:
                expected_heme_c = int(float(heme_raw))
            except ValueError:
                expected_heme_c = None
        expectations[stem] = {
            "expected_heme_c": expected_heme_c,
            "expected_topology": topo_raw or None,
        }
    return expectations


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
