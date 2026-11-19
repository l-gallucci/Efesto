# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Operon context filtering: FeGenie-derived rules, JSON rule engine, second-pass checks."""

import fnmatch
import json
import re
from pathlib import Path

from efesto.categories import canonical_category, context_category, rule_category


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
    {"name": "SIDERO_TRANSPORT", "category_basis": "context",
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
    {"name": "SIDERO_SYNTH", "category_basis": "context", "categories": ["Siderophore-Synthesis"],
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


def _mtr(rows, decisions=None):
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
    _record_decision(decisions, [r for r in rows if r["hmm_stem"] in _MTR_GENES],
                     rule_id="MtrMto", rule_type="mtr_disambiguation",
                     status="evaluated", reason_code="category_disambiguation_only",
                     observed_models=sorted(stems & _MTR_GENES))
    return updated


def _record_decision(decisions, rows, *, stage="context", rule_id,
                     rule_type, status, reason_code, excluded=(), **observations):
    """Record an actual branch decision without affecting its result."""
    if decisions is not None and rows:
        decisions.append(dict(stage=stage, rule_id=rule_id, rule_type=rule_type,
                              status=status, reason_code=reason_code,
                              applies_to_orfs=[r["orf"] for r in rows],
                              excluded_orfs=list(excluded), **observations))


def filter_cluster_fegenie(cluster_rows, report_all_pats, all_results=False,
                           catalog_mode=False, decisions=None):
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
        _record_decision(decisions, cluster_rows, rule_id="ALL_RESULTS",
                         rule_type="bypass", status="bypassed", reason_code="all_results")
        return cluster_rows

    cats  = {r["cat"]      for r in cluster_rows}
    stems = {r["hmm_stem"] for r in cluster_rows}

    if all(_cm(c, report_all_pats) for c in cats):
        _record_decision(decisions, cluster_rows, rule_id="REPORT_ALL_CATEGORIES",
                         rule_type="bypass", status="bypassed", reason_code="report_all_categories")
        return cluster_rows

    if catalog_mode:
        _record_decision(decisions, cluster_rows, rule_id="CATALOG_MODE",
                         rule_type="bypass", status="bypassed", reason_code="catalog_context_bypass")
        if stems & _MTR_GENES:
            return _mtr(cluster_rows, decisions=decisions)
        return cluster_rows

    # Evaluate each system against its own members in this neighborhood.
    # A decision about one system must not short-circuit unrelated rules.
    rejected = set()
    for name, genes, minimum in (("FLEET", _FLEET_GENES, 5), ("MAM", _MAM_GENES, 5),
                                 ("FOXABC", _FOXABC_GENES, 2), ("DFE1", _DFE1_GENES, 3),
                                 ("DFE2", _DFE2_GENES, 3), ("TAD_PILUS", _TAD_GENES, 2),
                                 ("SUF_OPERON", _SUF_GENES, 3)):
        failed = len(stems & genes) < minimum
        if failed:
            rejected.update(stems & genes)
        _record_decision(decisions, [r for r in cluster_rows if r["hmm_stem"] in genes],
                         rule_id=name, rule_type="require_n_of",
                         status="failed" if failed else "passed",
                         reason_code="minimum_models_not_met" if failed else "minimum_models_met",
                         excluded=[r["orf"] for r in cluster_rows if failed and r["hmm_stem"] in genes],
                         observed_models=sorted(stems & genes), minimum_models=minimum)
    for name, genes, anchors in (("FOXEYZ", _FOXEYZ_GENES, {"FoxE"}),
                                 ("CNT_CORE", _CNT_SYNTH_GENES, _CNT_CORE_GENES),
                                 ("MBN_CORE", _MBN_SYNTH_GENES, _MBN_SYNTH_GENES)):
        failed = not anchors <= stems
        if failed:
            rejected.update(stems & genes)
        _record_decision(decisions, [r for r in cluster_rows if r["hmm_stem"] in genes],
                         rule_id=name, rule_type="require_all_anchors",
                         status="failed" if failed else "passed",
                         reason_code="required_anchors_not_met" if failed else "required_anchors_met",
                         excluded=[r["orf"] for r in cluster_rows if failed and r["hmm_stem"] in genes],
                         observed_models=sorted(stems & genes), required_anchors=sorted(anchors))
    rows = [r for r in cluster_rows if r["hmm_stem"] not in rejected]
    # Count only distinct models in the relevant category group; another
    # system's members cannot satisfy a transport/synthesis requirement.
    for name, categories, minimum, trusted in (("SIDERO_TRANSPORT", _SIDERO_TRANS_CATS, 2, _TRUSTED_LONE),
                                               ("SIDERO_SYNTH", _SIDERO_SYNTH_CATS, 3, set()),
                                               ("IRON_TRANSPORT", _IRON_TRANS_CATS, 2, set())):
        members = {r["hmm_stem"] for r in rows
                   if context_category(r) in categories}
        failed = len(members) < minimum and not members & trusted
        affected = [r for r in rows if context_category(r) in categories]
        _record_decision(decisions, affected, rule_id=name, rule_type="require_n_cat",
                         status="failed" if failed else "passed",
                         reason_code=("minimum_models_not_met" if failed else
                                      "trusted_lone_observed" if len(members) < minimum else "minimum_models_met"),
                         excluded=[r["orf"] for r in affected] if failed else [],
                         observed_models=sorted(members), minimum_models=minimum, category_basis="context",
                         trusted_models_observed=sorted(members & trusted))
        if failed:
            rows = [r for r in rows if context_category(r) not in categories]
    # Category disambiguation does not exempt other systems from filtering.
    return _mtr(rows, decisions=decisions) if any(r["hmm_stem"] in _MTR_GENES for r in rows) else rows


def filter_cluster_json(cluster_rows, operon_rules, report_all_pats,
                        all_results=False, contig_len=None, relaxed_threshold=10000,
                        decisions=None):
    """
    JSON-rule engine (for operon_rules.json users — model organisms).
    Used when operon_rules.json is present in --hmm_dir.
    """
    if all_results:
        _record_decision(decisions, cluster_rows, rule_id="ALL_RESULTS",
                         rule_type="bypass", status="bypassed", reason_code="all_results")
        return cluster_rows
    cats = {r["cat"] for r in cluster_rows}
    if all(_cm(c, report_all_pats) for c in cats):
        _record_decision(decisions, cluster_rows, rule_id="REPORT_ALL_CATEGORIES",
                         rule_type="bypass", status="bypassed", reason_code="report_all_categories")
        return cluster_rows
    relaxed = (contig_len is not None and contig_len > 0
               and contig_len < relaxed_threshold)

    def _ugi(rows, gs):
        return len({r["hmm_stem"] for r in rows if r["hmm_stem"] in gs})

    def _ric(rows, cs, rd):
        return [r for r in rows if rule_category(r, rd) in cs]

    def _apply(rows, rd):
        gs = set(rd.get("genes", []))
        cs = {canonical_category(c) for c in rd.get("categories", [])}
        on_fail = rd.get("on_fail", "keep_all")
        rt = rd["rule"]
        sp = {r["hmm_stem"] for r in rows}
        cp = {rule_category(r, rd) for r in rows}
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
            if len({r["hmm_stem"] for r in _ric(rows, cs, rd)}) >= min_n:
                return rows
            return rows if on_fail == "keep_all" else []
        if rt == "require_n_cat_or_lone_trusted":
            tl    = set(rd.get("trusted_lone", []))
            cat_m = _ric(rows, cs, rd)
            uq    = {r["hmm_stem"] for r in cat_m}
            if len(uq) > 1 or tl & sp or len(uq) >= min_n:
                return rows
            return rows if on_fail == "keep_all" else []
        if rt == "mtr_disambiguation":
            return _mtr(rows)
        return rows

    rows = cluster_rows
    for index, rd in enumerate(operon_rules):
        before = rows
        rows = _apply(rows, rd)
        gs = set(rd.get("genes", []))
        cs = {canonical_category(c) for c in rd.get("categories", [])}
        named = [r for r in before if r["hmm_stem"] in gs]
        categorized = [r for r in before if rule_category(r, rd) in cs]
        if named or categorized:
            kind = rd["rule"]
            observed = {r["hmm_stem"] for r in (categorized if kind.startswith("require_n_cat") else named)}
            minimum = max(1, rd.get("min_genes", 1) // 2) if relaxed else rd.get("min_genes", 1)
            trusted = set(rd.get("trusted_lone", [])) & {r["hmm_stem"] for r in before}
            if kind == "require_n_of" or kind == "require_n_cat":
                passed = len(observed) >= minimum
                reason = "minimum_models_met" if passed else "minimum_models_not_met"
            elif kind == "require_anchor":
                passed = rd.get("anchor", "") in {r["hmm_stem"] for r in named}
                reason = "required_anchor_met" if passed else "required_anchor_not_met"
            elif kind == "require_n_cat_or_lone_trusted":
                passed = len(observed) > 1 or bool(trusted) or len(observed) >= minimum
                reason = ("trusted_lone_observed" if trusted else "minimum_models_met" if len(observed) >= minimum
                          else "two_model_exception" if len(observed) > 1 else "minimum_models_not_met")
            else:
                passed = True
                reason = "category_disambiguation_only" if kind == "mtr_disambiguation" else "unrecognized_rule_passthrough"
            kept = {r["orf"] for r in rows}
            excluded = [r["orf"] for r in before if r["orf"] not in kept]
            members = (categorized or named) if kind.startswith("require_n_cat") else (named or categorized)
            affected_orfs = {r["orf"] for r in members} | set(excluded)
            affected = [r for r in before if r["orf"] in affected_orfs]
            _record_decision(decisions, affected, rule_id=rd.get("name", "rule_"+str(index)),
                             rule_type=kind, status="passed" if passed else "failed", reason_code=reason,
                             excluded=excluded, observed_models=sorted(observed),
                             minimum_models=minimum if kind.startswith("require_n") else None,
                             required_anchors=[rd.get("anchor", "")] if kind == "require_anchor" else [],
                             on_fail=rd.get("on_fail", "keep_all"), relaxed=relaxed,
                             category_basis=rd.get("category_basis", "annotation"))
        if not rows:
            return []
    return rows


def count_heme(seq):
    if not seq:
        return 0
    return sum(len(re.findall(p, seq)) for p in
               [r"C(..)CH", r"C(...)CH", r"C(....)CH",
                r"C(.{14})CH", r"C(.{15})CH"])


def second_pass(cluster_rows, g2c, seq_dict, all_results=False, catalog_mode=False,
                decisions=None):
    if all_results:
        _record_decision(decisions, cluster_rows, stage="sequence", rule_id="ALL_RESULTS",
                         rule_type="bypass", status="bypassed", reason_code="all_results")
        return cluster_rows
    stems = [r["hmm_stem"] for r in cluster_rows]
    kept  = []
    for r in cluster_rows:
        stem = r["hmm_stem"]
        cat  = r["cat"]
        if stem == "Cyc1":
            observed = {h for h in set(stems) if g2c.get(h, "") in FE_REDOX}
            passed = not catalog_mode and len(observed) >= 2
            _record_decision(decisions, [r], stage="sequence", rule_id="CYC1_REDOX_CONTEXT",
                             rule_type="require_n_redox_models", status="passed" if passed else "failed",
                             reason_code="cyc1_catalog_exclusion" if catalog_mode else
                                         "minimum_redox_models_met" if passed else "minimum_redox_models_not_met",
                             excluded=[] if passed else [r["orf"]], observed_models=sorted(observed), minimum_models=2)
            if catalog_mode:
                continue
            if len({h for h in set(stems) if g2c.get(h, "") in FE_REDOX}) >= 2:
                kept.append(r)
            continue
        if re.match(r"Cyc2", stem):
            seq = seq_dict.get(r["genome"], {}).get(r["orf"], "")
            passed = len(seq) >= 365 and count_heme(seq) > 0
            _record_decision(decisions, [r], stage="sequence", rule_id="CYC2_LENGTH_HEME",
                             rule_type="sequence_requirements", status="passed" if passed else "failed",
                             reason_code="cyc2_requirements_met" if passed else "cyc2_requirements_not_met",
                             excluded=[] if passed else [r["orf"]], protein_length=len(seq),
                             minimum_length=365, heme_motifs=count_heme(seq), minimum_heme_motifs=1)
            if len(seq) >= 365 and count_heme(seq) > 0:
                kept.append(r)
            continue
        if cat == "Iron-Regulation":
            passed = catalog_mode or any("regulation" in g2c.get(h, "").lower() for h in stems)
            _record_decision(decisions, [r], stage="sequence", rule_id="IRON_REGULATION_CONTEXT",
                             rule_type="category_requirement", status="passed" if passed else "failed",
                             reason_code="regulation_requirement_met" if passed else "regulation_requirement_not_met",
                             excluded=[] if passed else [r["orf"]])
            if catalog_mode or any("regulation" in g2c.get(h, "").lower() for h in stems):
                kept.append(r)
            continue
        kept.append(r)
    return kept
