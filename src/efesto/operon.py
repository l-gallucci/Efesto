# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Operon context filtering: FeGenie-derived rules, JSON rule engine, second-pass checks."""

import fnmatch
import json
import re
from pathlib import Path

from efesto.categories import canonical_category


# ── Default operon rules (used when no operon_rules.json is present) ──────────
_DEFAULT_OPERON_RULES = [
    # max_bp_gap:     recommended max intergenic gap for this gene system (bp).
    #                 Used as metadata; the global --max_bp_gap is the detection window.
    {"name": "FLEET", "categories": ["Iron-Oxidation"],
     "genes": ["EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members",
     "max_bp_gap": 2000},
    {"name": "MAM", "categories": ["Magnetosome-Formation"],
     "genes": ["MamA", "MamB", "MamE", "MamK", "MamP", "MamM", "MamQ", "MamI", "MamL", "MamO"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members",
     "max_bp_gap": 3000},
    {"name": "FOXABC", "categories": ["Iron-Oxidation"],
     "genes": ["FoxA", "FoxB", "FoxC"],
     "rule": "require_n_of", "min_genes": 2, "on_fail": "passthrough_non_members",
     "max_bp_gap": 1000},
    {"name": "FOXEYZ", "categories": ["Iron-Oxidation"],
     "genes": ["FoxE", "FoxY", "FoxZ"],
     "rule": "require_anchor", "anchor": "FoxE", "on_fail": "passthrough_non_members",
     "max_bp_gap": 1000},
    {"name": "DFE1", "categories": ["Iron-Reduction", "Iron-ReductionPotential"],
     "genes": ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members",
     "max_bp_gap": 1000},
    {"name": "DFE2", "categories": ["Iron-Reduction", "Iron-ReductionPotential"],
     "genes": ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members",
     "max_bp_gap": 1000},
    {"name": "MtrMto",
     "categories": ["Iron-Oxidation", "Iron-Reduction",
                    "Iron-RedoxPotential"],
     # ImoA is excluded; native Fe(II) oxidation has not been demonstrated.
     "genes": ["MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507", "MtoA"],
     "rule": "mtr_disambiguation", "on_fail": "keep_all",
     "max_bp_gap": 2000},
    {"name": "SIDERO_TRANSPORT",
     "categories": ["Siderophore-TransportPotential",
                    "Heme-Transport",
                    "Siderophore-Transport", "IronAcquisition-ComplexTransport"],
     "genes": [], "rule": "require_n_cat_or_lone_trusted", "min_genes": 2,
     "trusted_lone": ["FutA1-iron_ABC_transporter_iron-binding-rep",
                      "FutA2-iron_ABC_transporter_iron-binding-rep",
                      "FutC-iron_ABC_transporter_ATPase-rep",
                      "LbtU-LvtA-PiuA-PirA-RhtA", "LbtU-LbtB-legiobactin_receptor",
                      "LbtU_LbtB-legiobactin_receptor_2", "IroC-salmochelin_transport-rep"],
     "on_fail": "drop", "max_bp_gap": 2000},
    {"name": "SIDERO_SYNTH", "categories": ["Siderophore-Synthesis"],
     "genes": [], "rule": "require_n_cat", "min_genes": 3, "on_fail": "drop",
     "max_bp_gap": 5000},
    {"name": "IRON_TRANSPORT",
     "categories": ["IronAcquisition-Transport", "Heme-Oxygenase"],
     "genes": [], "rule": "require_n_cat", "min_genes": 2, "on_fail": "drop",
     "max_bp_gap": 2000},
    {"name": "TAD_PILUS", "categories": ["Iron-Reduction"],
     "genes": ["PF04964-Flp_Fap_pilin", "cpaB_TIGR03177.1"],
     "rule": "require_n_of", "min_genes": 2, "on_fail": "passthrough_non_members",
     "max_bp_gap": 2000},
    {"name": "SUF_OPERON", "categories": ["IronSulfur-Assembly"],
     "genes": ["sufA_TIGR01997.1", "sufB_TIGR01980.1", "sufC_TIGR01978.1",
               "sufD_TIGR01981.1", "sufS_TIGR01979.1", "sufE_PF02657.22"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members",
     "max_bp_gap": 2000},
]

_DEFAULT_OPERON_RULES.extend(
[{'name': 'CNT_CORE_K23120',
  'categories': ['Metallophore-Synthesis'],
  'genes': ['K01780', 'K23120', 'K23122'],
  'rule': 'require_anchor',
  'anchor': 'K23120',
  'on_fail': 'passthrough_non_members',
  'max_bp_gap': 5000},
 {'name': 'CNT_CORE_K23122',
  'categories': ['Metallophore-Synthesis'],
  'genes': ['K01780', 'K23120', 'K23122'],
  'rule': 'require_anchor',
  'anchor': 'K23122',
  'on_fail': 'passthrough_non_members',
  'max_bp_gap': 5000}]
)

# Require all three specific families, including the short precursor/domain.
_DEFAULT_OPERON_RULES.extend([{'name': 'MBN_CORE_TIGR04071', 'categories': ['Metallophore-Synthesis'], 'genes': ['mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160'], 'rule': 'require_anchor', 'anchor': 'mbnA_TIGR04071', 'on_fail': 'passthrough_non_members', 'max_bp_gap': 5000}, {'name': 'MBN_CORE_TIGR04159', 'categories': ['Metallophore-Synthesis'], 'genes': ['mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160'], 'rule': 'require_anchor', 'anchor': 'mbnB_TIGR04159', 'on_fail': 'passthrough_non_members', 'max_bp_gap': 5000}, {'name': 'MBN_CORE_TIGR04160', 'categories': ['Metallophore-Synthesis'], 'genes': ['mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160'], 'rule': 'require_anchor', 'anchor': 'mbnC_TIGR04160', 'on_fail': 'passthrough_non_members', 'max_bp_gap': 5000}])

_REPORT_ALL_PATTERNS = ['MetalResistance-*', 'MetalAcquisition-*', 'MetalHomeostasis-*', 'MetalRegulation-*', 'Iron-Storage', 'Arsenic-Oxidation', 'Unclassified-MetalFunction']

# ── FeGenie gene sets ──────────────────────────────────────────────────────────
_FLEET_GENES  = {"EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"}
_MAM_GENES    = {"MamA", "MamB", "MamE", "MamK", "MamP", "MamM", "MamQ", "MamI", "MamL", "MamO"}
_FOXABC_GENES = {"FoxA", "FoxB", "FoxC"}
_FOXEYZ_GENES = {"FoxE", "FoxY", "FoxZ"}
_DFE1_GENES   = {"DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"}
_DFE2_GENES   = {"DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"}
_MTR_GENES    = {"MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507", "MtoA"}
_TAD_GENES    = {"PF04964-Flp_Fap_pilin", "cpaB_TIGR03177.1"}
_SUF_GENES    = {"sufA_TIGR01997.1", "sufB_TIGR01980.1", "sufC_TIGR01978.1",
                 "sufD_TIGR01981.1", "sufS_TIGR01979.1", "sufE_PF02657.22"}
_TRUSTED_LONE = {"FutA1-iron_ABC_transporter_iron-binding-rep",
                 "FutA2-iron_ABC_transporter_iron-binding-rep",
                 "FutC-iron_ABC_transporter_ATPase-rep",
                 "LbtU-LvtA-PiuA-PirA-RhtA", "LbtU-LbtB-legiobactin_receptor",
                 "LbtU_LbtB-legiobactin_receptor_2", "IroC-salmochelin_transport-rep"}

_SIDERO_TRANS_CATS = {"Siderophore-TransportPotential",
                      "Siderophore-Transport",
                      "IronAcquisition-ComplexTransport",
                      "Heme-Transport"}
_SIDERO_SYNTH_CATS = {"Siderophore-Synthesis"}
_CNT_SYNTH_GENES = {"K01780", "K23120", "K23122"}
_CNT_CORE_GENES = {"K23120", "K23122"}
_MBN_SYNTH_GENES = {'mbnA_TIGR04071', 'mbnB_TIGR04159', 'mbnC_TIGR04160'}
_IRON_TRANS_CATS   = {"IronAcquisition-Transport", "Heme-Oxygenase"}
_IRON_ACQ_ALL      = _SIDERO_TRANS_CATS | _SIDERO_SYNTH_CATS | _IRON_TRANS_CATS

FE_REDOX = {"Iron-Reduction", "Iron-Oxidation"}


def build_stem_gap_map(rules, gene_map):
    """
    Return {hmm_stem → max_bp_gap} from rules that define max_bp_gap for named genes.

    gene_map must be {hmm_stem → gene_name} (from read_map). Rules without a genes
    list (SIDERO_SYNTH, IRON_TRANSPORT, …) contribute their max_bp_gap as a fallback
    only for stems that belong to those categories — not implemented here since those
    rules have no named genes to anchor on. They continue using the global max_bp_gap.
    """
    name_to_gap = {}
    for rule in (rules if rules else _DEFAULT_OPERON_RULES):
        gap = rule.get("max_bp_gap")
        if gap:
            for g in rule.get("genes", []):
                name_to_gap[g] = gap
    result = {}
    for stem, name in (gene_map or {}).items():
        if name in name_to_gap:
            result[stem] = name_to_gap[name]
    return result


def load_operon_rules(hmm_dir):
    """Return (rules, report_all_pats, json_present).
    json_present=True → use JSON rule engine (model organisms).
    json_present=False → use the FeGenie-derived engine (default).
    """
    p = Path(hmm_dir) / "operon_rules.json"
    if p.exists():
        with open(p) as fh:
            data = json.load(fh)
        rules = [dict(rule, categories=[canonical_category(c) for c in rule.get("categories", [])])
                 for rule in data.get("rules", [])]
        pats  = [canonical_category(p) for p in data.get("report_all_categories", _REPORT_ALL_PATTERNS)]
        print(f"[INFO] Loaded {len(rules)} operon rules from {p}")
        print(f"       Using JSON rule engine (model-organism mode)")
        return rules, pats, True
    print("[INFO] Using FeGenie-derived operon rules (no operon_rules.json found)")
    return [], _REPORT_ALL_PATTERNS, False


def _cm(cat, pats):
    return any(fnmatch.fnmatch(canonical_category(cat), canonical_category(p)) for p in pats)


def _mtr(rows):
    stems   = {r["hmm_stem"] for r in rows}
    updated = [dict(r) for r in rows]
    if "MtoA" in stems and "MtrB_TIGR03509" in stems and "MtrC_TIGR03507" not in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrB_TIGR03509", "MtoA"}:
                r["cat"] = "Iron-Oxidation"
    elif "MtrA" in stems and "MtrB_TIGR03509" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "Iron-Reduction"
    elif "MtrC_TIGR03507" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "Iron-Reduction"
    return updated


def filter_cluster_fegenie(cluster_rows, report_all_pats, all_results=False,
                           catalog_mode=False):
    """
    Reimplementation of FeGenie's context rules with Efesto extensions.

    Based on Garber et al. (2020), doi:10.3389/fmicb.2020.00037.
    Efesto changes include catalog handling and additional gene systems;
    this is not a claim of exact behavior equivalence to every FeGenie mode.
    See Git history and CHANGELOG.md for modifications.

    catalog_mode=True: bypass co-occurrence count rules (FLEET ≥5, MAM ≥5,
    siderophore ≥2/3, iron_transport ≥2) but keep:
      - report_all bypass (MetalResistance-*, Iron-Storage)
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

    # Evaluate each system against its own members in this neighborhood.
    # A decision about one system must not short-circuit unrelated rules.
    rejected = set()
    for genes, minimum in ((_FLEET_GENES, 5), (_MAM_GENES, 5),
                           (_FOXABC_GENES, 2), (_DFE1_GENES, 3),
                           (_DFE2_GENES, 3), (_TAD_GENES, 2), (_SUF_GENES, 3)):
        if len(stems & genes) < minimum:
            rejected.update(stems & genes)
    for genes, anchors in ((_FOXEYZ_GENES, {"FoxE"}),
                           (_CNT_SYNTH_GENES, _CNT_CORE_GENES),
                           (_MBN_SYNTH_GENES, _MBN_SYNTH_GENES)):
        if not anchors <= stems:
            rejected.update(stems & genes)
    rows = [r for r in cluster_rows if r["hmm_stem"] not in rejected]
    # Count only distinct models in the relevant category group; another
    # system's members cannot satisfy a transport/synthesis requirement.
    for categories, minimum, trusted in ((_SIDERO_TRANS_CATS, 2, _TRUSTED_LONE),
                                          (_SIDERO_SYNTH_CATS, 3, set()),
                                          (_IRON_TRANS_CATS, 2, set())):
        members = {r["hmm_stem"] for r in rows
                   if canonical_category(r["cat"]) in categories}
        if len(members) < minimum and not members & trusted:
            rows = [r for r in rows if canonical_category(r["cat"]) not in categories]
    # Category disambiguation does not exempt other systems from filtering.
    return _mtr(rows) if any(r["hmm_stem"] in _MTR_GENES for r in rows) else rows


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
        cs = {canonical_category(c) for c in rd.get("categories", [])}
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
        if cat == "Iron-Regulation":
            if catalog_mode or any("regulation" in g2c.get(h, "").lower() for h in stems):
                kept.append(r)
            continue
        kept.append(r)
    return kept
