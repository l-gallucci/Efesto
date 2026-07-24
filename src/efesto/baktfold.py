"""Baktfold tier-3 structural confirmation: Foldseek+ProstT5 running and comparison.

Used only for HMM hits on models flagged needs_confirmation in the registry
(narrow-margin models like sufA/sufD/sufE — see docs/hmm_library_curation.md),
same narrow scope as eggnog.py's tier-2. Two things kept intentionally small:

  - No separate ProstT5 Python step. Foldseek's `easy-search` accepts a plain
    amino-acid FASTA directly and translates it to the 3Di structural
    alphabet internally when given --prostt5-model (van Kempen et al. 2023,
    Nat Biotechnol; ProstT5 weights via `foldseek databases ProstT5`), so
    this module never imports torch/transformers itself.
  - The reference structure DB is a small curated set (e.g. MtrCAB, a real
    cryo-EM structure — Edwards et al. 2020), not the full AlphaFold DB/PDB.
    Users build and point at it themselves, same opt-in pattern as eggNOG's
    reference database.
"""

import subprocess
import sys
from pathlib import Path


def run_foldseek(flagged_faa_path, db_dir, out_dir, foldseek_path="foldseek",
                 prostt5_weights=None, gpu=False, threads=4):
    """
    Run `foldseek easy-search` on the (already small, pre-filtered)
    flagged-ORF FASTA against a small curated reference-structure DB.

    Never called automatically — only when the user passes
    --run_struct_confirmation explicitly, since Foldseek and the ProstT5
    weights are both opt-in installs (see docs/installation.md).

    Returns path to the produced .m8 results file, or None on failure
    (prints a [WARN], never raises — same convention as run_eggnog/run_uniop).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_path = out_dir / "efesto_struct.m8"
    tmp_dir = out_dir / "foldseek_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = [str(foldseek_path), "easy-search", str(flagged_faa_path), str(db_dir),
           str(results_path), str(tmp_dir), "--threads", str(threads)]
    if prostt5_weights:
        cmd += ["--prostt5-model", str(prostt5_weights)]
    if gpu:
        cmd += ["--gpu", "1"]

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except Exception as e:
        print(f"  [WARN] Could not run Foldseek: {e}", file=sys.stderr)
        return None

    if r.returncode != 0:
        print(f"  [WARN] Foldseek exited with code {r.returncode}:\n"
              f"{r.stderr[-2000:]}", file=sys.stderr)
        return None

    if not results_path.exists():
        print(f"  [WARN] Foldseek finished but no results file found at "
              f"{results_path}", file=sys.stderr)
        return None
    return str(results_path)


def _parse_foldseek_results(path):
    """
    Parse foldseek's default easy-search tabular output (BLAST-tab-like, 12
    columns, no header): query, target, fident, alnlen, mismatch, gapopen,
    qstart, qend, tstart, tend, evalue, bits.

    Keeps only the best (highest-bitscore) hit per query — a reference DB
    can contain multiple structures per gene family (e.g. several MtrCAB
    templates), and only the top hit is used for the confirm/contradict
    comparison.

    Returns {query_id: {"target": str, "bits": float, "evalue": float}}.
    Returns {} on any parse failure (never raises) — same policy as
    eggnog._parse_eggnog_annotations.
    """
    hits = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 12:
                    continue
                query, target = fields[0], fields[1]
                try:
                    evalue, bits = float(fields[10]), float(fields[11])
                except ValueError:
                    continue
                existing = hits.get(query)
                if existing is None or bits > existing["bits"]:
                    hits[query] = {"target": target, "bits": bits, "evalue": evalue}
    except Exception as e:
        print(f"  [WARN] Could not parse Foldseek results file '{path}': {e}",
              file=sys.stderr)
        return {}
    return hits


def confirm_struct_hit(expected_gene_name, foldseek_row):
    """
    Compare an HMM's guessed gene identity against Foldseek's independent
    structural hit for the same ORF.

    Args:
        expected_gene_name: gene_name from the registry row that produced the
            flagged HMM hit (e.g. "sufA")
        foldseek_row: dict from _parse_foldseek_results for this ORF's query
            ID, or None if Foldseek found no hit for this ORF at all

    Returns one of "confirmed", "contradicted", "neutral".

    Comparison is a loose substring match against the target structure's
    name — reference DB entries should be named so the expected gene name
    appears in the target ID (e.g. "MtrA_ref_1zzz"), the same naming
    convention used for curating HMM seeds. "neutral" (no evidence either
    way) when foldseek_row is None (no structural hit at all) or
    expected_gene_name is empty. Never penalizes absence of data — same
    policy as eggnog.confirm_hit / uniop_weight's fallback.
    """
    if foldseek_row is None or not expected_gene_name:
        return "neutral"
    target = (foldseek_row.get("target") or "").strip().lower()
    if not target:
        return "neutral"
    if expected_gene_name.strip().lower() in target:
        return "confirmed"
    return "contradicted"
