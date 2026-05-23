"""Output writers: CSV/TSV summary files and heatmaps."""

import csv
from collections import defaultdict

from metalgenie_evo.coverage import compute_tpm


def write_summary(path, rows):
    with open(path, "w") as fh:
        fh.write("category,genome/assembly,contig,orf,gene,bitscore,"
                 "bitscore_cutoff,cluster_id,heme_c_motifs,protein_sequence\n")
        prev = None
        for r in rows:
            if prev is not None and r["cluster_id"] != prev:
                fh.write("#,#,#,#,#,#,#,#,#\n")
            orf_id = r.get("bakta_id", r["orf"])
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{orf_id},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cutoff']},"
                     f"{r['cluster_id']},{r['heme_motifs']},{r['sequence']}\n")
            prev = r["cluster_id"]


def write_gene_summary(path, rows):
    with open(path, "w") as fh:
        fh.write("process,assembly,contig,orf,gene,bitscore,cluster_id\n")
        prev = None
        for r in rows:
            if prev is not None and r["cluster_id"] != prev:
                fh.write("#,#,#,#,#,#,#\n")
            orf_id = r.get("bakta_id", r["orf"])
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{orf_id},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cluster_id']}\n")
            prev = r["cluster_id"]


def write_long_format(path, rows):
    has_uniop = any("uniop_context" in r for r in rows)
    fields = ["category", "genome", "contig", "orf", "gene", "bitscore",
              "bitscore_cutoff", "confidence", "cluster_id", "heme_c_motifs",
              "contig_len"]
    if has_uniop:
        fields.append("uniop_context")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            row = {
                "category":       r["cat"],
                "genome":         r["genome"],
                "contig":         r["contig"],
                "orf":            r.get("bakta_id", r["orf"]),
                "gene":           r["gene_name"],
                "bitscore":       f"{r['bitscore']:.1f}",
                "bitscore_cutoff": r["cutoff"],
                "confidence":     r.get("confidence", "low_confidence"),
                "cluster_id":     r["cluster_id"],
                "heme_c_motifs":  r["heme_motifs"],
                "contig_len":     r.get("contig_len", ""),
            }
            if has_uniop:
                row["uniop_context"] = r.get("uniop_context", "not_in_operon")
            w.writerow(row)


def write_heatmap(path, rows, all_genomes, norm_dict=None):
    all_cats = sorted({r["cat"] for r in rows})
    cm = defaultdict(lambda: defaultdict(set))
    for r in rows:
        cm[r["cat"]][r["genome"]].add(r["cluster_id"])
    with open(path, "w") as fh:
        fh.write("X," + ",".join(all_genomes) + "\n")
        for cat in all_cats:
            vals = []
            for g in all_genomes:
                raw = len(cm[cat].get(g, set()))
                vals.append(
                    f"{raw / norm_dict[g] * 1000:.4f}"
                    if norm_dict and norm_dict.get(g, 0) > 0
                    else str(raw)
                )
            fh.write(cat + "," + ",".join(vals) + "\n")


def write_coverage_heatmap(path, rows, all_genomes, genome_cov,
                           norm_coverage=False, contig_lengths=None):
    all_cats = sorted({r["cat"] for r in rows})
    cm = defaultdict(lambda: defaultdict(float))
    for r in rows:
        info = genome_cov.get(r["genome"], {}).get(r["contig"], {})
        if not info:
            continue
        if norm_coverage and contig_lengths:
            tpm = compute_tpm(genome_cov.get(r["genome"], {}),
                              contig_lengths.get(r["genome"], {}))
            val = tpm.get(r["contig"], 0.0)
        else:
            val = info.get("depth", 0.0)
        cm[r["cat"]][r["genome"]] += val
    label = "TPM" if norm_coverage else "mean_depth_sum"
    with open(path, "w") as fh:
        fh.write(f"# coverage metric: {label}\n")
        fh.write("X," + ",".join(all_genomes) + "\n")
        for cat in all_cats:
            fh.write(cat + "," + ",".join(
                f"{cm[cat].get(g, 0.0):.4f}" for g in all_genomes) + "\n")


def write_anvio_functions(path, final_rows, prodigal_to_bakta=None):
    """
    Functions table for anvi-import-functions (tab-delimited):
      gene_callers_id    source    accession    function    e_value

    gene_callers_id is taken from r["bakta_id"] if set (populated by the
    main loop when --bakta_gff_dir is provided), otherwise falls back to
    r["orf"] (Prodigal name — requires manual mapping before Anvi'o import).
    """
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_callers_id", "source", "accession", "function", "e_value"])
        seen = set()
        for r in final_rows:
            orf = r["orf"]
            if orf in seen:
                continue
            seen.add(orf)
            caller = r.get("bakta_id") or orf
            w.writerow([caller, "MetalGenie-Evo", r["hmm_stem"],
                        f"{r['gene_name']} [{r['cat']}]",
                        f"{r['evalue']:.2e}"])
