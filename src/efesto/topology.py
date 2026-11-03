"""Topology/motif consistency tier: Phobius parsing, running, and comparison logic.

Unlike the eggNOG (tier-2) and Baktfold (tier-3) escalation tiers — which are
scoped only to HMM hits on models flagged needs_confirmation in the registry —
this tier runs on every filtered HMM hit. It's cheap enough to afford that:
Phobius is a single lightweight batched call per genome (multi-FASTA input,
not one subprocess per ORF), and the heme-motif comparison reuses a count
(operon.count_heme) the pipeline already computes for every hit today.

Two independent signals, both optional per registry stem (see
io.build_topology_expectation_map):
  - expected_topology: transmembrane helix count / signal peptide prediction
    from Phobius, compared against a registry-declared expectation
    (e.g. MtrA/MtoA are lipoprotein-anchored decaheme cytochromes, expected
    to show a signal peptide and not a multi-pass membrane topology).
  - expected_heme_c: compared against the heme_motifs count operon.py
    already computes via count_heme() for every hit (e.g. MtrA/MtoA are
    decaheme — ~10 CXXCH-style motifs expected).

Same no-data-no-penalty policy as every other component: an unset
expectation, a missing Phobius binary, or an uninformative result all
resolve to neutral, never a penalty.
"""

import subprocess
import sys
from pathlib import Path


def run_phobius(faa_path, out_path, phobius_path="phobius.pl"):
    """
    Run Phobius (short-format output) on an entire filtered-hit FASTA for one
    genome, in a single batched call — Phobius accepts multi-FASTA input
    directly, so this stays one subprocess call per genome regardless of hit
    count.

    Returns path to the short-format output file, or None on failure (prints
    a [WARN], never raises — same convention as run_eggnog/run_uniop).
    """
    cmd = [str(phobius_path), "-short", str(faa_path)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except Exception as e:
        print(f"  [WARN] Could not run Phobius: {e}", file=sys.stderr)
        return None

    if r.returncode != 0:
        print(f"  [WARN] Phobius exited with code {r.returncode}:\n"
              f"{r.stderr[-2000:]}", file=sys.stderr)
        return None

    try:
        Path(out_path).write_text(r.stdout)
    except Exception as e:
        print(f"  [WARN] Could not write Phobius output to '{out_path}': {e}",
              file=sys.stderr)
        return None
    return str(out_path)


def _parse_phobius_output(path):
    """
    Parse Phobius' -short output format:

        SEQENCE_ID  TM  SP  PREDICTION

    TM = predicted transmembrane helix count, SP = "Y"/"0" signal-peptide
    flag. Whitespace-delimited (column widths vary with ID length).

    Returns {orf_id: {"signal_peptide": bool, "tm_helices": int}}.
    Returns {} on any parse failure (never raises) — same policy as
    eggnog._parse_eggnog_annotations.
    """
    hits = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.upper().startswith(("SEQENCE", "SEQUENCE")):
                    continue
                fields = line.split()
                if len(fields) < 3:
                    continue
                orf_id = fields[0]
                try:
                    tm = int(fields[1])
                except ValueError:
                    tm = 0
                sp = fields[2].strip().upper() == "Y"
                hits[orf_id] = {"signal_peptide": sp, "tm_helices": tm}
    except Exception as e:
        print(f"  [WARN] Could not parse Phobius output file '{path}': {e}",
              file=sys.stderr)
        return {}
    return hits


def topology_verdict(expected_topology, phobius_row):
    """
    Compare a registry's expected_topology for this stem against Phobius's
    prediction for this ORF.

    Conservative by design, mirrors eggnog.confirm_hit's neutral-biased
    policy: only contradicts on a clear architectural mismatch, stays
    neutral on borderline calls (a single predicted TM helix, or a
    lipoprotein signal Phobius' Sec-signal model may not catch cleanly).

    Returns one of "confirmed", "contradicted", "neutral".
    """
    if not expected_topology or phobius_row is None:
        return "neutral"
    tm = phobius_row.get("tm_helices", 0)
    sp = bool(phobius_row.get("signal_peptide"))

    if expected_topology == "signal_peptide":
        if sp:
            return "confirmed"
        if tm >= 2:
            return "contradicted"
        return "neutral"
    if expected_topology == "multi_pass_tm":
        if tm >= 2:
            return "confirmed"
        if tm == 0 and not sp:
            return "contradicted"
        return "neutral"
    if expected_topology == "cytoplasmic":
        if tm == 0 and not sp:
            return "confirmed"
        if sp or tm >= 2:
            return "contradicted"
        return "neutral"
    return "neutral"


def heme_verdict(expected_heme_c, observed_heme_c, tolerance=0.7):
    """
    Compare the observed heme-c motif count (operon.count_heme, already
    computed for every hit today) against the registry's expected_heme_c
    for this stem.

    Confirmed if observed >= tolerance * expected (allows minor regex
    undercount near sequence termini/overlapping matches). Contradicted only
    when the defining feature is essentially absent (observed == 0 while a
    substantial count, >=4, was expected) — a decaheme cytochrome hit with
    zero heme motifs is a strong signal; a partial undercount is not.
    Neutral otherwise, and always when expected_heme_c is unset.

    Returns one of "confirmed", "contradicted", "neutral".
    """
    if expected_heme_c is None or observed_heme_c is None:
        return "neutral"
    if observed_heme_c >= tolerance * expected_heme_c:
        return "confirmed"
    if observed_heme_c == 0 and expected_heme_c >= 4:
        return "contradicted"
    return "neutral"
