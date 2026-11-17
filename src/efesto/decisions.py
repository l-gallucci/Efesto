# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Stream selected HMM candidates and actual filter decisions to a TSV audit."""

import csv
import json
import tempfile


DECISION_SCHEMA = "filter-decisions-v1"


class FilterDecisionLog:
    """Spool candidate records to disk; keep only one neighborhood in memory."""

    def __init__(self, directory=None):
        self._spool = tempfile.TemporaryFile(mode="w+t", encoding="utf-8", dir=directory)

    def close(self):
        self._spool.close()

    def record(self, rows, decisions, context_rows, sequence_rows, *, engine,
               version, coordinate_orfs=(), catalog_mode=False, known_contigs=(),
               min_contig_len=0):
        context = {r["orf"]: r for r in context_rows}
        sequence = {r["orf"] for r in sequence_rows}
        indexed_events = [(event, set(event["applies_to_orfs"]), set(event["excluded_orfs"]))
                          for event in decisions]
        for row in rows:
            orf = row["orf"]
            events = []
            for event, members, excluded in indexed_events:
                if orf not in members:
                    continue
                events.append({**{k: v for k, v in event.items()
                                  if k not in {"applies_to_orfs", "excluded_orfs"}},
                               "action": "excluded" if orf in excluded else "retained"})
            length_excluded = any(e["stage"] == "contig_length" and e["action"] == "excluded" for e in events)
            context_bypassed = any(e["stage"] == "context" and e["status"] == "bypassed" for e in events)
            sequence_bypassed = any(e["stage"] == "sequence" and e["status"] == "bypassed" for e in events)
            record = dict(schema_version=DECISION_SCHEMA, genome=row["genome"], contig=row["contig"],
                          original_orf=orf, hmm_stem=row["hmm_stem"],
                          category_before_filter=row["cat"],
                          category_after_context=context[orf]["cat"] if orf in context else "",
                          bitscore=row["bitscore"], evalue=row["evalue"], bitscore_cutoff=row["cutoff"],
                          search_bitscore_cutoff=row["search_bitscore_cutoff"],
                          threshold_type=row["threshold_type"], threshold_source=row["threshold_source"],
                          model_source=row["model_source"], raw_cluster_id=row["cluster_id"],
                          coordinate_status="available" if orf in coordinate_orfs else "unavailable",
                          clustering_basis="not_clustered" if length_excluded else "catalog_singleton" if catalog_mode
                                           else "coordinates" if orf in coordinate_orfs else "orf_index_or_coordinate_fallback",
                          contig_length_status="known" if row["contig"] in known_contigs else
                                               "lower_bound" if row["contig_len"] > 0 else "unavailable",
                          contig_len=row["contig_len"],
                          minimum_contig_length=min_contig_len,
                          context_engine=engine, rule_version=version,
                          contig_length_decision="excluded" if length_excluded else
                                                 "not_applied" if min_contig_len <= 0 else
                                                 "retained" if row["contig_len"] > 0 else "not_evaluable",
                          context_decision="not_evaluated" if length_excluded else "excluded" if orf not in context
                                           else "bypassed" if context_bypassed else "retained",
                          sequence_decision="not_evaluated" if orf not in context else "excluded" if orf not in sequence
                                            else "bypassed" if sequence_bypassed else "retained",
                          rule_decisions=events)
            self._spool.write(json.dumps(record, separators=(",", ":"))+"\n")

    def write(self, path, retained_rows, gene_map, bakta_map):
        retained = {(r["genome"], r["orf"]): r for r in retained_rows}
        fields = ["schema_version", "genome", "contig", "orf", "original_orf", "gene", "hmm_stem",
                  "category_before_filter", "category_after_context", "final_category", "bitscore", "evalue",
                  "bitscore_cutoff", "search_bitscore_cutoff", "threshold_type", "threshold_source", "model_source",
                  "raw_cluster_id", "final_cluster_id", "coordinate_status", "clustering_basis", "contig_length_status",
                  "contig_len", "minimum_contig_length", "context_engine", "rule_version", "contig_length_decision", "context_decision",
                  "sequence_decision", "final_decision", "retention_mode", "reason_codes", "rule_decisions"]
        self._spool.seek(0)
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for line in self._spool:
                record = json.loads(line)
                genome, orf = record["genome"], record["original_orf"]
                final = retained.get((genome, orf))
                rescue = final is not None and final.get("retention_mode") == "eggnog_rescue"
                record.update(orf=bakta_map.get(genome, {}).get(orf, orf),
                              gene=gene_map.get(record["hmm_stem"], record["hmm_stem"]),
                              final_category=final["cat"] if final else "",
                              final_cluster_id=final["cluster_id"] if final else "",
                              final_decision="rescued" if rescue else "retained" if final else "excluded",
                              retention_mode=final.get("retention_mode", "") if final else "not_retained")
                reasons = [e["rule_id"]+":"+e["reason_code"] for e in record["rule_decisions"] if e["action"] == "excluded"]
                if rescue:
                    reasons.append("EGGNOG_RESCUE:matching_annotation")
                record["reason_codes"] = json.dumps(reasons, separators=(",", ":"))
                record["rule_decisions"] = json.dumps(record["rule_decisions"], separators=(",", ":"))
                writer.writerow(record)
