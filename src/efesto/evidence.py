# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Observable hit and neighborhood evidence; no aggregate confidence score.

Context assessments describe retained hits, not a replay of the filtering
engine or a claim of pathway completeness. Model identifiers are used directly.
"""

from collections import Counter
import fnmatch
import hashlib
import json

from efesto.categories import canonical_category
from efesto.eggnog import confirm_hit
from efesto.topology import heme_verdict, topology_verdict

OUTPUT_SCHEMA = "separate-evidence-v1"


def rule_version(rules, patterns):
    payload = json.dumps({"rules": rules, "report_all_categories": patterns},
                         sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def coordinate_evidence(orf, coords, contig_lengths):
    """Distances are observations in bp, never a penalty for MAG fragmentation."""
    c = coords.get(orf)
    if not c:
        return {"coordinate_status": "unavailable", "distance_to_contig_edge_bp": ""}
    length = contig_lengths.get(c["contig"])
    return {"coordinate_status": "available",
            "distance_to_contig_edge_bp":
                max(0, min(c["start"] - 1, length - c["end"])) if length else ""}


def hit_evidence(row, *, coords, contig_lengths, confirmation_map, topology_map,
                 pending_keys, eggnog_hits, eggnog_status, foldseek_hits,
                 foldseek_status, topology_requested, phobius_hits,
                 phobius_status, bgc_regions, bgc_status, uniop_status,
                 operon_map):
    """Keep availability, comparison and raw observations in separate fields."""
    orf, stem = row["orf"], row["hmm_stem"]
    result = coordinate_evidence(orf, coords, contig_lengths)
    pending = (row["genome"], orf) in pending_keys
    info = confirmation_map.get(stem)
    egg = eggnog_hits.get(orf) if pending and info else None
    fold = foldseek_hits.get(orf) if pending and info else None
    result.update({
        "eggnog_status": eggnog_status if pending and info else "not_applicable",
        "eggnog_verdict": "not_evaluated",
        "eggnog_preferred_name": "", "eggnog_kos": "",
        "structure_status": foldseek_status if pending and info else "not_applicable",
        "structure_target": "", "structure_bitscore": "", "structure_evalue": "",
        "structure_label_match": "not_evaluated",
    })
    if egg:
        result.update(eggnog_preferred_name=egg.get("preferred_name", ""),
                      eggnog_kos=",".join(egg.get("kegg_ko") or []),
                      eggnog_verdict=confirm_hit(info["gene_name"], info["kegg_ko"], egg))
    elif pending and info and eggnog_status in {"provided", "completed"}:
        result["eggnog_verdict"] = "no_informative_result"
    if fold:
        target = fold.get("target", "")
        # A label comparison is not a structural/functional confirmation.
        name = (info.get("gene_name") or "").strip().lower()
        result.update(structure_target=target, structure_bitscore=fold.get("bits", ""),
                      structure_evalue=fold.get("evalue", ""),
                      structure_label_match=("matches_label" if name in target.lower()
                                             else "different_label") if name and target
                                            else "unresolved")
    elif pending and info and foldseek_status == "completed":
        result["structure_label_match"] = "no_informative_result"

    expectation = topology_map.get(stem, {})
    phobius = phobius_hits.get(orf)
    applicable = bool(expectation)
    result.update({
        "expected_heme_c": expectation.get("expected_heme_c", ""),
        "expected_topology": expectation.get("expected_topology", ""),
        "phobius_status": phobius_status if applicable else "not_applicable",
        "predicted_tm_helices": "", "predicted_signal_peptide": "",
        "topology_verdict": "not_evaluated" if expectation.get("expected_topology") else "not_applicable",
        "heme_verdict": "not_evaluated" if expectation.get("expected_heme_c") is not None else "not_applicable",
    })
    if topology_requested and applicable:
        if expectation.get("expected_heme_c") is not None:
            result["heme_verdict"] = heme_verdict(expectation["expected_heme_c"], row.get("heme_motifs"))
        if expectation.get("expected_topology"):
            result["topology_verdict"] = topology_verdict(expectation["expected_topology"], phobius)
            if not phobius:
                result["topology_verdict"] = "no_informative_result" if phobius_status == "completed" else "not_evaluated"
        if phobius:
            result["predicted_tm_helices"] = phobius["tm_helices"]
            result["predicted_signal_peptide"] = phobius["signal_peptide"]

    c = coords.get(orf)
    overlaps = [dict(contig=contig, start=start, end=end, product=product)
                for contig, start, end, product in (bgc_regions or [])
                if c and c["contig"] == contig and c["start"] <= end and c["end"] >= start]
    result.update(bgc_status=bgc_status, bgc_overlap="not_evaluated",
                  bgc_overlap_regions=json.dumps(overlaps, separators=(",", ":")),
                  uniop_status=uniop_status,
                  uniop_operon_id=operon_map.get(orf, ""),
                  uniop_context=row.get("uniop_context", "not_in_operon")
                    if uniop_status in {"completed", "cached"} else "not_evaluated")
    if bgc_status == "parsed":
        result["bgc_overlap"] = ("overlap" if overlaps else "no_overlap") if c else "coordinates_unavailable"
    return result


def context_evidence(rows, rules, patterns, *, mode, coords, relaxed=False,
                     json_mode=False, available_models=None, configuration_version=None):
    """Assess rule requirements on retained hits, without changing retention.

    Distinct model counts are not a completeness fraction. Bypass and missing
    genomic context are explicit. Version hashes identify the assessed rules.
    """
    if not rows:
        return []
    stems = {r["hmm_stem"] for r in rows}
    categories = {canonical_category(r["cat"]) for r in rows}
    copies = Counter(r["hmm_stem"] for r in rows)
    known = [coords[r["orf"]] for r in rows if r["orf"] in coords]
    coordinate_status = ("available" if len(known) == len(rows) else
                         "partial" if known else "unavailable")
    gaps = []
    for contig in sorted({c["contig"] for c in known}):
        ordered = sorted((c for c in known if c["contig"] == contig), key=lambda c: c["start"])
        gaps.extend(max(0, b["start"] - a["end"] - 1) for a, b in zip(ordered, ordered[1:]))
    base = dict(genome=rows[0]["genome"], cluster_id=rows[0]["cluster_id"],
                assessment_scope="retained_hits_not_filter_trace", mode=mode,
                rule_version=configuration_version or rule_version(rules, patterns),
                coordinate_status=coordinate_status,
                observed_models=sorted(stems), model_copy_counts=dict(sorted(copies.items())),
                inter_hit_gaps_bp=gaps)
    category_bypass = all(any(fnmatch.fnmatch(c, canonical_category(p)) for p in patterns)
                          for c in categories)
    records = []
    for rule in rules:
        named = set(rule.get("genes", []))
        cats = {canonical_category(c) for c in rule.get("categories", [])}
        members = stems & named if named else {r["hmm_stem"] for r in rows
                                               if canonical_category(r["cat"]) in cats}
        if not members:
            continue
        kind = rule.get("rule", "")
        minimum = rule.get("min_genes", 1) if kind in {"require_n_of", "require_n_cat", "require_n_cat_or_lone_trusted"} else ""
        if relaxed and json_mode and isinstance(minimum, int):
            minimum = max(1, minimum // 2)
        anchor = rule.get("anchor", "")
        not_observed = sorted(named - stems)
        not_searched = (sorted(named - available_models) if available_models is not None else [])
        bypass = mode == "all_results" or category_bypass
        if bypass:
            status = "bypassed"
        elif kind == "mtr_disambiguation":
            status = "category_disambiguation_only"
        elif mode == "catalog" and not json_mode:
            status = "bypassed"
        elif any(r.get("tier2_rescued") for r in rows):
            status = "not_evaluable_after_rescue"
        elif coordinate_status != "available" or mode == "catalog":
            status = "not_evaluable_without_coordinates"
        else:
            if kind in {"require_n_of", "require_n_cat"}:
                status = "satisfied_minimum" if len(members) >= minimum else "not_satisfied"
            elif kind == "require_anchor":
                status = "anchor_observed" if anchor in members else "not_satisfied"
            elif kind == "require_n_cat_or_lone_trusted":
                status = ("satisfied_minimum" if len(members) >= minimum else
                          "trusted_lone_observed" if stems & set(rule.get("trusted_lone", [])) else "not_satisfied")
            else:
                status = "unassessed_rule_type"
            if status == "not_satisfied" and not_searched:
                status = "not_evaluable_models_not_searched"
        records.append(dict(base, rule_id=rule.get("name", "unnamed"), rule_type=kind,
                            rule_status=status, rule_models=sorted(named),
                            observed_rule_models=sorted(members),
                            distinct_rule_models=len(members), minimum_models=minimum,
                            anchor=anchor, rule_models_not_observed=not_observed,
                            rule_models_not_searched=not_searched))
    if not records:
        records.append(dict(base, rule_id="", rule_type="", rule_status="not_applicable",
                            rule_models=[], observed_rule_models=[], distinct_rule_models=0,
                            minimum_models="", anchor="", rule_models_not_observed=[],
                            rule_models_not_searched=[]))
    return records
