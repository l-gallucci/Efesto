# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Output writers: CSV/TSV summary files and heatmaps."""

from efesto.identifiers import export_ids, scoped_id
from urllib.parse import quote

import csv
import json
import datetime
from collections import defaultdict
from pathlib import Path

from efesto.coverage import compute_coverage_per_million
from efesto.io import read_fasta

_NT_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def write_summary(path, rows):
    fields = ["category", "genome/assembly", "contig", "orf", "gene", "bitscore",
              "bitscore_cutoff", "threshold_type", "cluster_id", "heme_c_motifs",
              "protein_sequence"]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(fields)
        prev = None
        for row in rows:
            if prev is not None and row["cluster_id"] != prev:
                writer.writerow(["#"] * len(fields))
            writer.writerow([row["cat"], row["genome"], row["contig"],
                             row.get("bakta_id", row["orf"]), row["gene_name"],
                             f"{row['bitscore']:.1f}", row["cutoff"],
                             row.get("threshold_type", "generic"), row["cluster_id"],
                             row["heme_motifs"], row["sequence"]])
            prev = row["cluster_id"]


def write_gene_summary(path, rows):
    fields = ["process", "assembly", "contig", "orf", "gene", "bitscore", "cluster_id"]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(fields)
        prev = None
        for row in rows:
            if prev is not None and row["cluster_id"] != prev:
                writer.writerow(["#"] * len(fields))
            writer.writerow([row["cat"], row["genome"], row["contig"],
                             row.get("bakta_id", row["orf"]), row["gene_name"],
                             f"{row['bitscore']:.1f}", row["cluster_id"]])
            prev = row["cluster_id"]


HIT_EVIDENCE_FIELDS = [
    "threshold_source", "search_bitscore_cutoff", "model_source", "model_status", "retention_mode",
    "coordinate_status", "contig_length_status", "distance_to_contig_edge_bp", "expected_heme_c", "expected_topology",
    "eggnog_status", "eggnog_verdict", "eggnog_preferred_name", "eggnog_kos",
    "structure_status", "structure_target", "structure_bitscore", "structure_evalue",
    "structure_label_match", "phobius_status", "predicted_tm_helices",
    "predicted_signal_peptide", "topology_verdict", "heme_verdict",
    "bgc_status", "bgc_overlap", "bgc_overlap_regions", "uniop_status", "uniop_operon_id",
]


def write_long_format(path, rows):
    fields = ["category", "genome", "contig", "orf", "original_orf", "gene", "hmm_stem", "bitscore",
              "evalue", "bitscore_cutoff", "threshold_type", "cluster_id", "heme_c_motifs",
              "contig_len", "model_nseq", *HIT_EVIDENCE_FIELDS, "uniop_context"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            record = {field: row.get(field, "") for field in HIT_EVIDENCE_FIELDS}
            record["model_status"] = row.get("model_status") or "unknown"
            record.update(category=row["cat"], genome=row["genome"], contig=row["contig"],
                          orf=row.get("bakta_id", row["orf"]), original_orf=row["orf"], gene=row["gene_name"],
                          hmm_stem=row["hmm_stem"], bitscore=f"{row['bitscore']:.1f}",
                          evalue=row["evalue"], bitscore_cutoff=row["cutoff"],
                          threshold_type=row.get("threshold_type", "generic"),
                          cluster_id=row["cluster_id"], heme_c_motifs=row["heme_motifs"],
                          contig_len=row.get("contig_len", ""), model_nseq=row.get("model_nseq", ""),
                          uniop_context=row.get("uniop_context", "not_evaluated"))
            writer.writerow(record)


def write_context_evidence(path, records):
    fields = ["genome", "cluster_id", "assessment_scope", "mode", "rule_id", "rule_version",
              "rule_type", "rule_status", "category_basis", "coordinate_status", "observed_models",
              "model_copy_counts", "inter_hit_gaps_bp", "rule_models", "observed_rule_models",
              "distinct_rule_models", "minimum_models", "anchor", "rule_models_not_observed",
              "rule_models_not_searched"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow({key: json.dumps(value, separators=(",", ":"))
                             if isinstance(value, (list, dict)) else value
                             for key, value in record.items()})


def write_uniop_pair_evidence(path, rows, genome_pair_probs):
    """Available UniOP pair probabilities, including zero; no cluster summary."""
    identifiers = export_ids(rows)
    retained = {(row["genome"], row["orf"]): row for row in rows}
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome", "orf_a", "orf_b", "original_orf_a", "original_orf_b",
                         "cluster_a", "cluster_b", "operon_probability", "uniop_status"])
        for genome, pairs in sorted(genome_pair_probs.items()):
            for (a, b), probability in sorted(pairs.items()):
                ka, kb = (genome, a), (genome, b)
                if ka not in retained or kb not in retained:
                    continue
                writer.writerow([genome, identifiers[ka], identifiers[kb], a, b,
                                 retained[ka]["cluster_id"], retained[kb]["cluster_id"],
                                 probability, retained[ka].get("uniop_status", "")])


def write_heatmap(path, rows, all_genomes, norm_dict=None):
    all_cats = sorted({r["cat"] for r in rows})
    cm = defaultdict(lambda: defaultdict(set))
    for r in rows:
        cm[r["cat"]][r["genome"]].add(r["cluster_id"])
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["X", *all_genomes])
        for cat in all_cats:
            vals = []
            for g in all_genomes:
                raw = len(cm[cat].get(g, set()))
                vals.append(
                    f"{raw / norm_dict[g] * 1000:.4f}"
                    if norm_dict and norm_dict.get(g, 0) > 0
                    else str(raw)
                )
            writer.writerow([cat, *vals])


def write_coverage_heatmap(path, rows, all_genomes, genome_cov,
                           norm_coverage=False, contig_lengths=None):
    all_cats = sorted({r["cat"] for r in rows})
    normalised = ({genome: compute_coverage_per_million(genome_cov.get(genome, {}))
                   for genome in all_genomes} if norm_coverage else {})
    cm = defaultdict(lambda: defaultdict(float))
    for r in rows:
        info = genome_cov.get(r["genome"], {}).get(r["contig"], {})
        if not info:
            continue
        if norm_coverage:
            val = normalised.get(r["genome"], {}).get(r["contig"], 0.0)
        else:
            val = info.get("depth", 0.0)
        cm[r["cat"]][r["genome"]] += val
    label = "coverage_per_million_sum" if norm_coverage else "mean_depth_sum"
    with open(path, "w", newline="") as fh:
        fh.write(f"# coverage metric: {label}\n")
        writer = csv.writer(fh)
        writer.writerow(["X", *all_genomes])
        for cat in all_cats:
            writer.writerow([cat, *(f"{cm[cat].get(g, 0.0):.4f}" for g in all_genomes)])


def write_gff3(path, final_rows, genome_coords):
    """
    Write HMM hits as GFF3 features.

    Only rows with coordinate information (from genome_coords) are written.
    Rows without coordinates are silently skipped (no coordinate data without
    --gff_dir or --fna_dir).

    Attributes per feature:
        ID, gene, category, hmm_stem, cluster_id, threshold_type
    """
    identifiers = export_ids(final_rows)
    contig_owners = defaultdict(set)
    for genome, coords in genome_coords.items():
        for coord in coords.values():
            contig_owners[coord["contig"]].add(genome)
    scope_contigs = any(len(owners) > 1 for owners in contig_owners.values())
    n_written = 0
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"## Generated by Efesto  {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
        for r in final_rows:
            coords = genome_coords.get(r["genome"], {}).get(r["orf"])
            if coords is None:
                continue
            contig = coords["contig"]
            seqid = quote(scoped_id(r["genome"], contig) if scope_contigs else contig, safe="._:-")
            start  = coords["start"]
            end    = coords["end"]
            strand = coords.get("strand", ".")
            orf_id = identifiers[(r["genome"], r["orf"])]
            attrs  = (
                f"ID={quote(orf_id, safe='._:-')}"
                f";genome={quote(r['genome'], safe='')}"
                f";original_contig={quote(contig, safe='')}"
                f";original_orf={quote(r['orf'], safe='')}"
                f";bakta_id={quote(r.get('bakta_id') or r['orf'], safe='')}"
                f";gene={r['gene_name']}"
                f";category={r['cat']}"
                f";hmm_stem={r['hmm_stem']}"
                f";cluster_id={r['cluster_id']}"
                f";threshold_type={r.get('threshold_type', 'generic')}"
            )
            fh.write(
                f"{seqid}\tEfesto\tCDS\t{start}\t{end}"
                f"\t.\t{strand}\t.\t{attrs}\n"
            )
            n_written += 1
    return n_written


def write_summary_stats(path, final_rows, all_genomes, genome_coords=None,
                        runtime_s=None):
    """Write counts and evidence availability, without aggregate quality tiers."""
    clusters = {r["cluster_id"] for r in final_rows}
    cat_counts = defaultdict(int)
    genome_hits = defaultdict(int)
    threshold_counts = defaultdict(int)
    evidence_counts = defaultdict(int)
    for row in final_rows:
        cat_counts[row["cat"]] += 1
        genome_hits[row["genome"]] += 1
        threshold_counts[row.get("threshold_type", "generic")] += 1
        for field in ["eggnog_status", "structure_status", "phobius_status", "bgc_status", "uniop_status"]:
            evidence_counts[(field, row.get(field, "not_evaluated"))] += 1
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["section", "metric", "value", "notes"])
        def record(section, metric, value, note=""):
            writer.writerow([section, metric, value, note])
        record("RUN", "total_orf_hits", len(final_rows))
        record("RUN", "total_clusters", len(clusters))
        record("RUN", "genomes_with_hits", len(genome_hits))
        record("RUN", "genomes_zero_hits", len(all_genomes) - len(genome_hits))
        if runtime_s is not None:
            record("RUN", "runtime_seconds", f"{runtime_s:.1f}")
        for kind, count in sorted(threshold_counts.items()):
            record("THRESHOLD", kind, count, "ORF hits; threshold type is not validation status")
        for (field, status), count in sorted(evidence_counts.items()):
            record("EVIDENCE", field, count, status)
        for category, count in sorted(cat_counts.items()):
            record("CATEGORY", category, count, "ORF hits")
        for genome in sorted(all_genomes):
            record("GENOME", genome, genome_hits.get(genome, 0), "ORF hits")


def write_anvio_misc_data(path, final_rows, prodigal_to_bakta=None):
    """Raw numeric hit evidence; map exported labels to anvi'o integer gene IDs."""
    identifiers = export_ids(final_rows)
    fields = ["gene_callers_id", "bitscore", "evalue", "bitscore_cutoff",
              "search_bitscore_cutoff", "heme_c_motifs", "contig_len"]
    seen = set()
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(fields)
        for row in final_rows:
            key = (row["genome"], row["orf"])
            if key in seen:
                continue
            seen.add(key)
            writer.writerow([identifiers[key], row["bitscore"], row["evalue"], row["cutoff"],
                             row.get("search_bitscore_cutoff", ""), row["heme_motifs"],
                             row.get("contig_len", "")])


def write_anvio_functions(path, final_rows, prodigal_to_bakta=None):
    """
    Functions table for anvi-import-functions (tab-delimited):
      gene_callers_id    source    accession    function    e_value

    Exported labels follow the same collision-safe mapping as the evidence table.
    Map labels to the target database's integer IDs before anvi'o import.
    """
    identifiers = export_ids(final_rows)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_callers_id", "source", "accession", "function", "e_value"])
        seen = set()
        for r in final_rows:
            key = (r["genome"], r["orf"])
            if key in seen:
                continue
            seen.add(key)
            caller = identifiers[key]
            source = ("Efesto-generic-threshold"
                      if r.get("threshold_type") == "generic" else "Efesto")
            w.writerow([caller, source, r["hmm_stem"],
                        f"{r['gene_name']} [{r['cat']}]",
                        f"{r['evalue']:.2e}"])


def write_hit_faa(path, final_rows):
    """Write protein sequences for all HMM hits as FASTA with annotated headers.

    Header fields: gene, category, genome, contig, bitscore, cluster_id,
                   threshold_type
    """
    identifiers = export_ids(final_rows)
    n_written = 0
    with open(path, "w") as fh:
        for r in final_rows:
            seq = r.get("sequence", "")
            if not seq:
                continue
            orf_id = identifiers[(r["genome"], r["orf"])]
            fh.write(
                f">{orf_id} "
                f"gene={r['gene_name']} "
                f"category={r['cat']} "
                f"genome={r['genome']} "
                f"contig={r['contig']} "
                f"bitscore={r['bitscore']:.1f} "
                f"cluster_id={r['cluster_id']} "
                f"threshold_type={r.get('threshold_type', 'generic')}\n"
            )
            fh.write(seq + "\n")
            n_written += 1
    return n_written


def write_hit_fna(path, final_rows, genome_coords, fna_dir, fna_ext="fna"):
    """Write nucleotide sequences for all HMM hits as FASTA with annotated headers.

    Requires genome_coords (from GFF) and fna_dir (nucleotide assembly files).
    Hits without coordinate data or whose contig is absent from the FNA are skipped.
    Sequences on the minus strand are reverse-complemented.

    Header fields: gene, category, genome, contig, start, end, strand,
                   bitscore, cluster_id, threshold_type
    """
    identifiers = export_ids(final_rows)
    fna_dir = Path(fna_dir)
    by_genome = defaultdict(list)
    for r in final_rows:
        if r["genome"] in genome_coords:
            by_genome[r["genome"]].append(r)

    n_written = 0
    with open(path, "w") as fh:
        for genome_faa, rows in sorted(by_genome.items()):
            stem = Path(genome_faa).stem
            fna_path = None
            for ext in [fna_ext, "fna", "fasta", "fa"]:
                candidate = fna_dir / f"{stem}.{ext}"
                if candidate.exists():
                    fna_path = candidate
                    break
            if fna_path is None:
                continue
            contigs  = read_fasta(str(fna_path))
            coords_g = genome_coords[genome_faa]
            for r in rows:
                c = coords_g.get(r["orf"])
                if c is None:
                    continue
                contig_seq = contigs.get(c["contig"], "")
                if not contig_seq:
                    continue
                start  = c["start"] - 1   # GFF 1-based inclusive → 0-based
                end    = c["end"]
                subseq = contig_seq[start:end]
                if c.get("strand") == "-":
                    subseq = subseq.translate(_NT_COMP)[::-1]
                orf_id = identifiers[(r["genome"], r["orf"])]
                fh.write(
                    f">{orf_id} "
                    f"gene={r['gene_name']} "
                    f"category={r['cat']} "
                    f"genome={genome_faa} "
                    f"contig={c['contig']} "
                    f"start={c['start']} "
                    f"end={c['end']} "
                    f"strand={c.get('strand', '.')} "
                    f"bitscore={r['bitscore']:.1f} "
                    f"cluster_id={r['cluster_id']} "
                    f"threshold_type={r.get('threshold_type', 'generic')}\n"
                )
                fh.write(subseq + "\n")
                n_written += 1
    return n_written
