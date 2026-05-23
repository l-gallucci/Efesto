"""Operon context filtering: FeGenie port, JSON rule engine, second-pass checks."""

import fnmatch
import json
import re
from pathlib import Path


# ── Default operon rules (used when no operon_rules.json is present) ──────────
_DEFAULT_OPERON_RULES = [
    {"name": "FLEET", "categories": ["iron_oxidation"],
     "genes": ["EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members"},
    {"name": "MAM", "categories": ["magnetosome_formation"],
     "genes": ["MamA", "MamB", "MamE", "MamK", "MamP", "MamM", "MamQ", "MamI", "MamL", "MamO"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members"},
    {"name": "FOXABC", "categories": ["iron_oxidation"],
     "genes": ["FoxA", "FoxB", "FoxC"],
     "rule": "require_n_of", "min_genes": 2, "on_fail": "passthrough_non_members"},
    {"name": "FOXEYZ", "categories": ["iron_oxidation"],
     "genes": ["FoxE", "FoxY", "FoxZ"],
     "rule": "require_anchor", "anchor": "FoxE", "on_fail": "passthrough_non_members"},
    {"name": "DFE1", "categories": ["iron_reduction", "probable_iron_reduction"],
     "genes": ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members"},
    {"name": "DFE2", "categories": ["iron_reduction", "probable_iron_reduction"],
     "genes": ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members"},
    {"name": "MtrMto",
     "categories": ["iron_oxidation", "iron_reduction",
                    "possible_iron_oxidation_and_possible_iron_reduction"],
     "genes": ["MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507", "MtoA", "CymA"],
     "rule": "mtr_disambiguation", "on_fail": "keep_all"},
    {"name": "SIDERO_TRANSPORT",
     "categories": ["iron_aquisition-siderophore_transport_potential",
                    "iron_aquisition-heme_transport",
                    "iron_aquisition-siderophore_transport"],
     "genes": [], "rule": "require_n_cat_or_lone_trusted", "min_genes": 2,
     "trusted_lone": ["FutA1-iron_ABC_transporter_iron-binding-rep",
                      "FutA2-iron_ABC_transporter_iron-binding-rep",
                      "FutC-iron_ABC_transporter_ATPase-rep",
                      "LbtU-LvtA-PiuA-PirA-RhtA", "LbtU-LbtB-legiobactin_receptor",
                      "LbtU_LbtB-legiobactin_receptor_2", "IroC-salmochelin_transport-rep"],
     "on_fail": "drop"},
    {"name": "SIDERO_SYNTH", "categories": ["iron_aquisition-siderophore_synthesis"],
     "genes": [], "rule": "require_n_cat", "min_genes": 3, "on_fail": "drop"},
    {"name": "IRON_TRANSPORT",
     "categories": ["iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"],
     "genes": [], "rule": "require_n_cat", "min_genes": 2, "on_fail": "drop"},
]

_REPORT_ALL_PATTERNS = ["metal_resistance-*", "iron_storage"]

# ── FeGenie gene sets ──────────────────────────────────────────────────────────
_FLEET_GENES  = {"EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"}
_MAM_GENES    = {"MamA", "MamB", "MamE", "MamK", "MamP", "MamM", "MamQ", "MamI", "MamL", "MamO"}
_FOXABC_GENES = {"FoxA", "FoxB", "FoxC"}
_FOXEYZ_GENES = {"FoxE", "FoxY", "FoxZ"}
_DFE1_GENES   = {"DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"}
_DFE2_GENES   = {"DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"}
_MTR_GENES    = {"MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507", "MtoA", "CymA"}
_TRUSTED_LONE = {"FutA1-iron_ABC_transporter_iron-binding-rep",
                 "FutA2-iron_ABC_transporter_iron-binding-rep",
                 "FutC-iron_ABC_transporter_ATPase-rep",
                 "LbtU-LvtA-PiuA-PirA-RhtA", "LbtU-LbtB-legiobactin_receptor",
                 "LbtU_LbtB-legiobactin_receptor_2", "IroC-salmochelin_transport-rep"}

_SIDERO_TRANS_CATS = {"iron_aquisition-siderophore_transport_potential",
                      "iron_aquisition-siderophore_transport",
                      "iron_aquisition-heme_transport"}
_SIDERO_SYNTH_CATS = {"iron_aquisition-siderophore_synthesis"}
_IRON_TRANS_CATS   = {"iron_aquisition-iron_transport", "iron_aquisition-heme_oxygenase"}
_IRON_ACQ_ALL      = _SIDERO_TRANS_CATS | _SIDERO_SYNTH_CATS | _IRON_TRANS_CATS

FE_REDOX = {"iron_reduction", "iron_oxidation"}


def load_operon_rules(hmm_dir):
    """Return (rules, report_all_pats, json_present).
    json_present=True → use JSON rule engine (model organisms).
    json_present=False → use FeGenie exact port (default).
    """
    p = Path(hmm_dir) / "operon_rules.json"
    if p.exists():
        with open(p) as fh:
            data = json.load(fh)
        rules = data.get("rules", [])
        pats  = data.get("report_all_categories", _REPORT_ALL_PATTERNS)
        print(f"[INFO] Loaded {len(rules)} operon rules from {p}")
        print(f"       Using JSON rule engine (model-organism mode)")
        return rules, pats, True
    print("[INFO] Using FeGenie exact operon logic (no operon_rules.json found)")
    return [], _REPORT_ALL_PATTERNS, False


def _cm(cat, pats):
    return any(fnmatch.fnmatch(cat, p) for p in pats)


def _mtr(rows):
    stems   = {r["hmm_stem"] for r in rows}
    updated = [dict(r) for r in rows]
    if "MtoA" in stems and "MtrB_TIGR03509" in stems and "MtrC_TIGR03507" not in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrB_TIGR03509", "MtoA", "CymA"}:
                r["cat"] = "iron_oxidation"
    elif "MtrA" in stems and "MtrB_TIGR03509" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "iron_reduction"
    elif "MtrC_TIGR03507" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "iron_reduction"
    return updated


def filter_cluster_fegenie(cluster_rows, report_all_pats, all_results=False,
                           catalog_mode=False):
    """
    Exact port of FeGenie's operon context filtering logic.

    catalog_mode=True: bypass co-occurrence count rules (FLEET ≥5, MAM ≥5,
    siderophore ≥2/3, iron_transport ≥2) but keep:
      - report_all bypass (metal_resistance-*, iron_storage)
      - Mtr/Mto disambiguation (category assignment, not count-based)
      - Cyc1/Cyc2 second-pass filters (handled separately in second_pass())
    This is appropriate for deduplicated gene catalogs where genomic context
    is unavailable by definition.
    """
    if all_results:
        return cluster_rows

    cats  = {r["cat"]      for r in cluster_rows}
    stems = {r["hmm_stem"] for r in cluster_rows}

    if all(_cm(c, report_all_pats) for c in cats):
        return cluster_rows

    if catalog_mode:
        if stems & _MTR_GENES:
            return _mtr(cluster_rows)
        return cluster_rows

    def n_unique_total():
        return len(stems)

    def n_unique_in_cats(cat_set):
        return len({r["hmm_stem"] for r in cluster_rows if r["cat"] in cat_set})

    def n_in_gene_set(gene_set):
        return len({r["hmm_stem"] for r in cluster_rows if r["hmm_stem"] in gene_set})

    if stems & _FLEET_GENES:
        n = n_in_gene_set(_FLEET_GENES)
        non_fleet = [r for r in cluster_rows if r["hmm_stem"] not in _FLEET_GENES]
        if n >= 5:
            return cluster_rows
        return non_fleet if non_fleet else []

    if stems & _MAM_GENES:
        n = n_in_gene_set(_MAM_GENES)
        non_mam = [r for r in cluster_rows if r["hmm_stem"] not in _MAM_GENES]
        if n >= 5:
            return cluster_rows
        return non_mam if non_mam else []

    if stems & _FOXABC_GENES:
        n = n_in_gene_set(_FOXABC_GENES)
        non_fox = [r for r in cluster_rows if r["hmm_stem"] not in _FOXABC_GENES]
        if n >= 2:
            return cluster_rows
        return non_fox if non_fox else []

    if stems & _FOXEYZ_GENES:
        non_fox = [r for r in cluster_rows if r["hmm_stem"] not in _FOXEYZ_GENES]
        if "FoxE" in stems:
            return cluster_rows
        return non_fox if non_fox else []

    if stems & _DFE1_GENES:
        n = n_in_gene_set(_DFE1_GENES)
        non_dfe = [r for r in cluster_rows if r["hmm_stem"] not in _DFE1_GENES]
        if n >= 3:
            return cluster_rows
        return non_dfe if non_dfe else []

    if stems & _DFE2_GENES:
        n = n_in_gene_set(_DFE2_GENES)
        non_dfe = [r for r in cluster_rows if r["hmm_stem"] not in _DFE2_GENES]
        if n >= 3:
            return cluster_rows
        return non_dfe if non_dfe else []

    # Cyc1 handled in second_pass, not here

    if "CymA" in stems:
        if stems & {"MtrA", "MtoA", "MtrB_TIGR03509", "MtrC_TIGR03507"}:
            return _mtr(cluster_rows)
        return [r for r in cluster_rows if r["hmm_stem"] != "CymA"] or []

    if stems & _MTR_GENES:
        return _mtr(cluster_rows)

    if cats & _IRON_ACQ_ALL:
        n_total = n_unique_total()
        kept = []
        skip_cluster = False

        for r in cluster_rows:
            cat = r["cat"]

            if cat in _SIDERO_TRANS_CATS:
                if n_total < 2:
                    skip_cluster = True
                    break
                if n_unique_in_cats(_SIDERO_TRANS_CATS) < 2:
                    pass
                else:
                    kept.append(r)

            elif cat in _SIDERO_SYNTH_CATS:
                if n_total < 3:
                    skip_cluster = True
                    break
                if n_unique_in_cats(_SIDERO_SYNTH_CATS) < 3:
                    pass
                else:
                    kept.append(r)

            elif cat in _IRON_TRANS_CATS:
                if n_total < 2:
                    skip_cluster = True
                    break
                if n_unique_in_cats(_IRON_TRANS_CATS) < 2:
                    pass
                else:
                    kept.append(r)

            else:
                kept.append(r)

        if skip_cluster:
            return []

        if not kept:
            if n_total > 1 or (stems & _TRUSTED_LONE):
                return cluster_rows
            return []

        return kept

    return cluster_rows


def filter_cluster_json(cluster_rows, operon_rules, report_all_pats,
                        all_results=False, contig_len=None, relaxed_threshold=10000):
    """
    JSON-rule engine (for operon_rules.json users — model organisms).
    Used when operon_rules.json is present in --hmm_dir.
    """
    if all_results:
        return cluster_rows
    cats = {r["cat"] for r in cluster_rows}
    if all(_cm(c, report_all_pats) for c in cats):
        return cluster_rows
    relaxed = (contig_len is not None and contig_len > 0
               and contig_len < relaxed_threshold)

    def _ugi(rows, gs):
        return len({r["hmm_stem"] for r in rows if r["hmm_stem"] in gs})

    def _ric(rows, cs):
        return [r for r in rows if r["cat"] in cs]

    def _apply(rows, rd):
        gs = set(rd.get("genes", []))
        cs = set(rd.get("categories", []))
        on_fail = rd.get("on_fail", "keep_all")
        rt = rd["rule"]
        sp = {r["hmm_stem"] for r in rows}
        cp = {r["cat"] for r in rows}
        if not ((gs and sp & gs) or (cs and cp & cs)):
            return rows
        raw_min = rd.get("min_genes", 1)
        min_n   = max(1, raw_min // 2) if relaxed else raw_min
        non_mbrs = [r for r in rows if r["hmm_stem"] not in gs]
        if rt == "require_n_of":
            if _ugi(rows, gs) >= min_n:
                return rows
            if on_fail == "passthrough_non_members" and non_mbrs:
                return non_mbrs
            return rows if on_fail == "keep_all" else []
        if rt == "require_anchor":
            anchor = rd.get("anchor", "")
            if anchor in {r["hmm_stem"] for r in rows if r["hmm_stem"] in gs}:
                return rows
            if on_fail == "passthrough_non_members" and non_mbrs:
                return non_mbrs
            return rows if on_fail == "keep_all" else []
        if rt == "require_n_cat":
            if len({r["hmm_stem"] for r in _ric(rows, cs)}) >= min_n:
                return rows
            return rows if on_fail == "keep_all" else []
        if rt == "require_n_cat_or_lone_trusted":
            tl    = set(rd.get("trusted_lone", []))
            cat_m = _ric(rows, cs)
            uq    = {r["hmm_stem"] for r in cat_m}
            if len(uq) > 1 or tl & sp or len(uq) >= min_n:
                return rows
            return rows if on_fail == "keep_all" else []
        if rt == "mtr_disambiguation":
            return _mtr(rows)
        return rows

    rows = cluster_rows
    for rd in operon_rules:
        rows = _apply(rows, rd)
        if not rows:
            return []
    return rows


def count_heme(seq):
    if not seq:
        return 0
    return sum(len(re.findall(p, seq)) for p in
               [r"C(..)CH", r"C(...)CH", r"C(....)CH",
                r"C(.{14})CH", r"C(.{15})CH"])


def second_pass(cluster_rows, g2c, seq_dict, all_results=False, catalog_mode=False):
    if all_results:
        return cluster_rows
    stems = [r["hmm_stem"] for r in cluster_rows]
    kept  = []
    for r in cluster_rows:
        stem = r["hmm_stem"]
        cat  = r["cat"]
        if stem == "Cyc1":
            if catalog_mode:
                continue
            if len({h for h in set(stems) if g2c.get(h, "") in FE_REDOX}) >= 2:
                kept.append(r)
            continue
        if re.match(r"Cyc2", stem):
            seq = seq_dict.get(r["genome"], {}).get(r["orf"], "")
            if len(seq) >= 365 and count_heme(seq) > 0:
                kept.append(r)
            continue
        if cat == "iron_gene_regulation":
            if catalog_mode or any("regulation" in g2c.get(h, "") for h in stems):
                kept.append(r)
            continue
        kept.append(r)
    return kept
