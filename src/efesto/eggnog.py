"""eggNOG-mapper tier-2 confirmation: parsing, running, and comparison logic.

Used only for HMM hits on models flagged needs_confirmation in the registry
(narrow-margin models like sufA/sufD/sufE — see docs/hmm_library_curation.md).
Two modes, both optional:
  - read-only: parse a .annotations file the user already generated elsewhere
  - internal:  shell out to emapper.py on a small flagged-ORF subset FASTA
    (never run automatically on a whole proteome — eggNOG's reference
    database is tens of GB, this is opt-in and scoped)
"""

import subprocess
import sys
from pathlib import Path


def _parse_eggnog_annotations(path):
    """
    Parse an eggNOG-mapper *.emapper.annotations TSV file.

    Format: several '##' comment/metadata lines, then a header line starting
    with '#query', then data rows, then a trailing '##' summary footer.
    Column of interest: query, Preferred_name, KEGG_ko, Description.

    Returns {query_id: {"preferred_name": str, "kegg_ko": [str, ...],
                         "description": str}}.
    Returns {} on any parse failure (never raises) — a missing/malformed
    annotations file should not crash the run, just leave tier 2 neutral.
    """
    hits = {}
    try:
        with open(path) as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").split("\t")
                    continue
                if header is None:
                    continue
                fields = line.split("\t")
                if len(fields) < len(header):
                    fields += [""] * (len(header) - len(fields))
                row = dict(zip(header, fields))
                query = row.get("query", "").strip()
                if not query:
                    continue
                ko_raw = row.get("KEGG_ko", "-").strip()
                kos = [k.strip().removeprefix("ko:") for k in ko_raw.split(",")
                       if k.strip() and k.strip() != "-"]
                # eggNOG-mapper uses "-" as its own placeholder for "no value"
                # in every column, not just KEGG_ko — normalize to "" so
                # downstream comparisons don't treat "-" as real data.
                pref = row.get("Preferred_name", "").strip()
                desc = row.get("Description", "").strip()
                hits[query] = {
                    "preferred_name": "" if pref == "-" else pref,
                    "kegg_ko": kos,
                    "description": "" if desc == "-" else desc,
                }
    except Exception as e:
        print(f"  [WARN] Could not parse eggNOG annotations file '{path}': {e}",
              file=sys.stderr)
        return {}
    return hits


def confirm_hit(expected_gene_name, expected_kos, eggnog_row):
    """
    Compare an HMM's guessed gene identity against eggNOG-mapper's independent
    annotation for the same ORF.

    Args:
        expected_gene_name: gene_name from the registry row that produced the
            flagged HMM hit (e.g. "sufA")
        expected_kos: kegg_ko list from build_confirmation_map for that stem
            (e.g. ["K05997"])
        eggnog_row: dict from _parse_eggnog_annotations for this ORF's query
            ID, or None if eggNOG has no entry for this ORF at all

    Returns one of "confirmed", "contradicted", "neutral".

    "neutral" (no evidence either way) when:
      - eggnog_row is None (no eggNOG hit for this ORF)
      - eggNOG assigned no KEGG_ko and no Preferred_name (uninformative hit)
      - expected_kos is empty (registry has no KO for this model — can't compare)
    Never penalizes absence of data — same policy as uniop_weight's fallback.
    """
    if eggnog_row is None:
        return "neutral"
    ko_hits = eggnog_row.get("kegg_ko") or []
    pref = (eggnog_row.get("preferred_name") or "").strip().lower()

    if not ko_hits and not pref:
        return "neutral"

    if expected_kos and ko_hits:
        if set(expected_kos) & set(ko_hits):
            return "confirmed"
        return "contradicted"

    # No KO on one side — fall back to a loose name comparison.
    if expected_gene_name and pref:
        if expected_gene_name.strip().lower() == pref:
            return "confirmed"
        return "contradicted"

    return "neutral"


def write_flagged_faa(flagged_orf_ids, faa_path, out_path):
    """
    Write a subset FASTA containing only the ORFs in flagged_orf_ids.

    Used both as input to run_eggnog()/Baktfold and as a standalone artifact
    (--export_flagged_faa) so a user can run their own structure/annotation
    tool of choice on exactly the same small candidate set efesto flagged.

    Returns the number of sequences written.
    """
    flagged_orf_ids = set(flagged_orf_ids)
    if not flagged_orf_ids:
        Path(out_path).write_text("")
        return 0
    written = 0
    write_this = False
    with open(faa_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                orf_id = line[1:].split()[0].strip()
                write_this = orf_id in flagged_orf_ids
                if write_this:
                    written += 1
            if write_this:
                fh_out.write(line)
    return written


def run_eggnog(flagged_faa_path, out_dir, eggnog_path="emapper.py",
               db_dir=None, cpu=4):
    """
    Optionally shell out to emapper.py on the (already small, pre-filtered)
    flagged-ORF FASTA. Never called automatically — only when the user passes
    --run_eggnog explicitly, since the eggNOG reference database is tens of
    GB and users manage that download themselves.

    Returns path to the produced .annotations file, or None on failure
    (prints a [WARN], never raises — same convention as run_uniop).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = "efesto_flagged"

    cmd = ["python3", str(eggnog_path),
           "-i", str(flagged_faa_path),
           "--itype", "proteins",
           "-o", prefix,
           "--output_dir", str(out_dir),
           "--cpu", str(cpu)]
    if db_dir:
        cmd += ["--data_dir", str(db_dir)]

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except Exception as e:
        print(f"  [WARN] Could not run eggNOG-mapper: {e}", file=sys.stderr)
        return None

    if r.returncode != 0:
        print(f"  [WARN] eggNOG-mapper exited with code {r.returncode}:\n"
              f"{r.stderr[-2000:]}", file=sys.stderr)
        return None

    annotations = out_dir / f"{prefix}.emapper.annotations"
    if not annotations.exists():
        print(f"  [WARN] eggNOG-mapper finished but no annotations file found "
              f"at {annotations}", file=sys.stderr)
        return None
    return str(annotations)
