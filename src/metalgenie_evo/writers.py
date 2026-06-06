"""Output writers: CSV/TSV summary files and heatmaps."""

import csv
import datetime
from collections import defaultdict
from pathlib import Path

from metalgenie_evo.coverage import compute_tpm
from metalgenie_evo.io import read_fasta

_NT_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _fmt_conf(r):
    c = r.get("cluster_confidence")
    return f"{c:.3f}" if c is not None else ""


def write_summary(path, rows):
    with open(path, "w") as fh:
        fh.write("category,genome/assembly,contig,orf,gene,bitscore,"
                 "bitscore_cutoff,cluster_id,cluster_confidence,"
                 "heme_c_motifs,protein_sequence\n")
        prev = None
        for r in rows:
            if prev is not None and r["cluster_id"] != prev:
                fh.write("#,#,#,#,#,#,#,#,#,#\n")
            orf_id = r.get("bakta_id", r["orf"])
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{orf_id},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cutoff']},"
                     f"{r['cluster_id']},{_fmt_conf(r)},"
                     f"{r['heme_motifs']},{r['sequence']}\n")
            prev = r["cluster_id"]


def write_gene_summary(path, rows):
    with open(path, "w") as fh:
        fh.write("process,assembly,contig,orf,gene,bitscore,cluster_id,"
                 "cluster_confidence\n")
        prev = None
        for r in rows:
            if prev is not None and r["cluster_id"] != prev:
                fh.write("#,#,#,#,#,#,#,#\n")
            orf_id = r.get("bakta_id", r["orf"])
            fh.write(f"{r['cat']},{r['genome']},{r['contig']},{orf_id},"
                     f"{r['gene_name']},{r['bitscore']:.1f},{r['cluster_id']},"
                     f"{_fmt_conf(r)}\n")
            prev = r["cluster_id"]


def write_long_format(path, rows):
    has_uniop = any("uniop_context" in r for r in rows)
    fields = ["category", "genome", "contig", "orf", "gene", "bitscore",
              "bitscore_cutoff", "confidence", "cluster_id", "heme_c_motifs",
              "contig_len", "cluster_confidence", "model_nseq"]
    if has_uniop:
        fields.append("uniop_context")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            row = {
                "category":          r["cat"],
                "genome":            r["genome"],
                "contig":            r["contig"],
                "orf":               r.get("bakta_id", r["orf"]),
                "gene":              r["gene_name"],
                "bitscore":          f"{r['bitscore']:.1f}",
                "bitscore_cutoff":   r["cutoff"],
                "confidence":        r.get("confidence", "low_confidence"),
                "cluster_id":        r["cluster_id"],
                "heme_c_motifs":     r["heme_motifs"],
                "contig_len":        r.get("contig_len", ""),
                "cluster_confidence": f"{r.get('cluster_confidence', ''):.3f}"
                                      if r.get("cluster_confidence") is not None
                                      else "",
                "model_nseq":        r.get("model_nseq", ""),
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


def write_gff3(path, final_rows, genome_coords):
    """
    Write HMM hits as GFF3 features.

    Only rows with coordinate information (from genome_coords) are written.
    Rows without coordinates are silently skipped (no coordinate data without
    --gff_dir or --fna_dir).

    Attributes per feature:
        ID, gene, category, hmm_stem, cluster_id, cluster_confidence, confidence
    """
    n_written = 0
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"## Generated by MetalGenie-Evo  {datetime.datetime.now().strftime('%Y-%m-%d')}\n")
        for r in final_rows:
            coords = genome_coords.get(r["genome"], {}).get(r["orf"])
            if coords is None:
                continue
            contig = coords["contig"]
            start  = coords["start"]
            end    = coords["end"]
            strand = coords.get("strand", ".")
            conf   = r.get("cluster_confidence")
            conf_s = f"{conf:.3f}" if conf is not None else "."
            orf_id = r.get("bakta_id", r["orf"])
            attrs  = (
                f"ID={orf_id}"
                f";gene={r['gene_name']}"
                f";category={r['cat']}"
                f";hmm_stem={r['hmm_stem']}"
                f";cluster_id={r['cluster_id']}"
                f";cluster_confidence={conf_s}"
                f";confidence={r.get('confidence', 'low_confidence')}"
            )
            fh.write(
                f"{contig}\tMetalGenie-Evo\tCDS\t{start}\t{end}"
                f"\t.\t{strand}\t.\t{attrs}\n"
            )
            n_written += 1
    return n_written


def write_summary_stats(path, final_rows, all_genomes, genome_coords=None,
                        runtime_s=None):
    """
    Write a human-readable summary statistics TSV.

    Sections:
        RUN        — overall counts and runtime
        CONFIDENCE — cluster confidence tier distribution
        CATEGORY   — hit counts per HMM category
        GENOME     — per-genome hit count and mean cluster_confidence
    """
    clusters  = defaultdict(list)
    for r in final_rows:
        clusters[r["cluster_id"]].append(r)

    confs = [r.get("cluster_confidence") for r in final_rows
             if r.get("cluster_confidence") is not None]
    by_cluster_conf = [grp[0].get("cluster_confidence")
                       for grp in clusters.values()
                       if grp[0].get("cluster_confidence") is not None]

    n_high   = sum(1 for c in by_cluster_conf if c >= 0.8)
    n_mid    = sum(1 for c in by_cluster_conf if 0.5 <= c < 0.8)
    n_low    = sum(1 for c in by_cluster_conf if c < 0.5)
    mean_all = sum(by_cluster_conf) / len(by_cluster_conf) if by_cluster_conf else 0.0

    cat_counts = defaultdict(int)
    for r in final_rows:
        cat_counts[r["cat"]] += 1

    genome_hits  = defaultdict(int)
    genome_confs = defaultdict(list)
    for r in final_rows:
        genome_hits[r["genome"]] += 1
        c = r.get("cluster_confidence")
        if c is not None:
            genome_confs[r["genome"]].append(c)

    n_genomes_hit  = len(genome_hits)
    n_genomes_zero = len(all_genomes) - n_genomes_hit

    with open(path, "w") as fh:
        def row(section, metric, value, note=""):
            fh.write(f"{section}\t{metric}\t{value}\t{note}\n")

        fh.write("# MetalGenie-Evo summary statistics\n")
        fh.write(f"# {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write("section\tmetric\tvalue\tnotes\n")

        row("RUN", "total_orf_hits",     len(final_rows))
        row("RUN", "total_clusters",     len(clusters))
        row("RUN", "genomes_with_hits",  n_genomes_hit)
        row("RUN", "genomes_zero_hits",  n_genomes_zero,
            "zero-hit genomes listed at bottom")
        if runtime_s is not None:
            row("RUN", "runtime_seconds", f"{runtime_s:.1f}")

        row("CONFIDENCE", "mean_cluster_confidence", f"{mean_all:.3f}",
            "per-cluster, not per-orf")
        row("CONFIDENCE", "high_confidence_clusters",   n_high,  "conf >= 0.8")
        row("CONFIDENCE", "medium_confidence_clusters", n_mid,   "0.5 <= conf < 0.8")
        row("CONFIDENCE", "low_confidence_clusters",    n_low,   "conf < 0.5")

        for cat, n in sorted(cat_counts.items()):
            row("CATEGORY", cat, n, "ORF hits")

        for g in sorted(all_genomes):
            n     = genome_hits.get(g, 0)
            confs_g = genome_confs.get(g, [])
            mean_g  = f"{sum(confs_g)/len(confs_g):.3f}" if confs_g else "NA"
            row("GENOME", g, n, f"mean_cluster_confidence={mean_g}")


def write_anvio_misc_data(path, final_rows, prodigal_to_bakta=None):
    """
    Per-gene numeric scores for anvi'o misc-data import.

    Format (tab-delimited):
        gene_callers_id  cluster_confidence  co_occ_score  hmm_weight
                         uniop_weight  bgc_boost

    gene_callers_id is the Bakta gene ID when --bakta_gff_dir was used, otherwise
    the Prodigal ORF name. For anvi'o import, integer gene_callers_ids are required:

        anvi-import-misc-data -c CONTIGS.db \\
            --target-data-table genes \\
            MetalGenie-Evo-anvio-gene-scores.tsv

    Prodigal users must first map ORF names to integer IDs:
        anvi-export-gene-calls -c CONTIGS.db -o gene_calls.tsv
    then join on the orf name in column 'contig'_'start'_'stop'_'direction'.
    Bakta users with external gene calls imported via anvi-import-external-gene-calls
    can use the Bakta locus tags directly.
    """
    fields = ["gene_callers_id", "cluster_confidence", "co_occ_score",
              "hmm_weight", "uniop_weight", "bgc_boost"]
    seen = set()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(fields)
        for r in final_rows:
            orf = r["orf"]
            if orf in seen:
                continue
            seen.add(orf)
            caller = r.get("bakta_id") or orf
            conf   = r.get("cluster_confidence")
            w.writerow([
                caller,
                f"{conf:.4f}"         if conf is not None else "",
                f"{r.get('co_occ_score', ''):.4f}"
                    if r.get("co_occ_score") is not None else "",
                f"{r.get('hmm_w', ''):.4f}"
                    if r.get("hmm_w") is not None else "",
                f"{r.get('uniop_w', ''):.4f}"
                    if r.get("uniop_w") is not None else "",
                f"{r.get('bgc_boost', 1.0):.4f}",
            ])


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
            source = ("MetalGenie-Evo-low_confidence"
                      if r.get("confidence") == "low_confidence"
                      else "MetalGenie-Evo")
            w.writerow([caller, source, r["hmm_stem"],
                        f"{r['gene_name']} [{r['cat']}]",
                        f"{r['evalue']:.2e}"])


def write_hit_faa(path, final_rows):
    """Write protein sequences for all HMM hits as FASTA with annotated headers.

    Header fields: gene, category, genome, contig, bitscore, cluster_id,
                   cluster_confidence, confidence
    """
    n_written = 0
    with open(path, "w") as fh:
        for r in final_rows:
            seq = r.get("sequence", "")
            if not seq:
                continue
            orf_id = r.get("bakta_id", r["orf"])
            conf   = r.get("cluster_confidence")
            conf_s = f"{conf:.3f}" if conf is not None else "NA"
            fh.write(
                f">{orf_id} "
                f"gene={r['gene_name']} "
                f"category={r['cat']} "
                f"genome={r['genome']} "
                f"contig={r['contig']} "
                f"bitscore={r['bitscore']:.1f} "
                f"cluster_id={r['cluster_id']} "
                f"cluster_confidence={conf_s} "
                f"confidence={r.get('confidence', 'low_confidence')}\n"
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
                   bitscore, cluster_id, cluster_confidence
    """
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
                orf_id = r.get("bakta_id", r["orf"])
                conf   = r.get("cluster_confidence")
                conf_s = f"{conf:.3f}" if conf is not None else "NA"
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
                    f"cluster_confidence={conf_s}\n"
                )
                fh.write(subseq + "\n")
                n_written += 1
    return n_written
