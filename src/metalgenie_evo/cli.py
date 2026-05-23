"""
MetalGenie-Evo  —  HMM-based annotation of iron cycling and metal resistance genes
===================================================================================
Built on FeGenie (Garber et al. 2020). See README for full documentation.

New in this version (metagenome support):
  7.  Integrated Prodigal  --fna_dir triggers internal ORF prediction
  8.  Contig length filter  --min_contig_len
  9.  Relaxed operon thresholds  --relaxed_operons
  10. TPM coverage normalisation  --norm_coverage
  11. Contig column in all outputs
  12. Long-format tidy output  MetalGenie-Evo-results-long.tsv
"""

import argparse
import datetime
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from metalgenie_evo.clustering import _orf_to_contig, build_clusters
from metalgenie_evo.coverage import build_contig_coverage, load_bams_tsv
from metalgenie_evo.gene_calling import run_prodigal
from metalgenie_evo.hmmer import collect_best_hits, run_all_hmmsearches
from metalgenie_evo.io import (build_contig_length_dict, load_gff_dir,
                                read_fasta, read_cutoffs, read_map)
from metalgenie_evo.operon import (count_heme, filter_cluster_fegenie,
                                    filter_cluster_json, load_operon_rules,
                                    second_pass)
from metalgenie_evo.uniop import (
    _uniop_context, build_prodigal_bakta_map, run_uniop, write_operon_structure)
from metalgenie_evo.writers import (
    write_anvio_functions, write_coverage_heatmap, write_gene_summary,
    write_heatmap, write_long_format, write_summary)


def main():
    p = argparse.ArgumentParser(
        prog="MetalGenie-Evo",
        description="HMM-based annotation of iron cycling and metal resistance genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ig = p.add_mutually_exclusive_group(required=True)
    ig.add_argument("--faa_dir", help="Directory of ORF .faa files (Prodigal output)")
    ig.add_argument("--fna_dir",
                    help="Directory of nucleotide assemblies — Prodigal run internally")
    p.add_argument("--faa_ext", default="faa")
    p.add_argument("--fna_ext", default="fna")
    p.add_argument("--meta", action="store_true",
                   help="Pyrodigal metagenomic mode (pre-trained models). "
                        "RECOMMENDED for MAGs and fragmented assemblies. "
                        "For complete/near-complete genomes, omit this flag. "
                        "Ref: Tang Li et al. 2022, Briefings in Bioinformatics.")
    p.add_argument("--gene_caller", default="pyrodigal-gv",
                   choices=["pyrodigal-gv", "pyrodigal"],
                   help="Gene prediction tool (default: pyrodigal-gv). "
                        "pyrodigal-gv: bacterial + viral models, best for mixed "
                        "or viral datasets. "
                        "pyrodigal: bacterial/archaeal models only, identical "
                        "results to Prodigal v2.6.3. "
                        "Both require the respective Python package to be installed.")
    p.add_argument("--gff_dir", help="Prodigal GFF files for bp clustering")
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--out", default="metalgenie_evo_out")
    p.add_argument("--min_contig_len", type=int, default=0,
                   help="Skip ORFs on contigs shorter than this (bp). 0=no filter")
    p.add_argument("--max_gap", type=int, default=5,
                   help="Max ORF-index gap (index mode)")
    p.add_argument("--max_bp_gap", type=int, default=5000,
                   help="Max bp gap (GFF mode)")
    p.add_argument("--strand_aware", action="store_true")
    p.add_argument("--relaxed_operons", action="store_true",
                   help="Halve operon min_genes for contigs < --relaxed_threshold")
    p.add_argument("--relaxed_threshold", type=int, default=10000,
                   help="Contig length (bp) threshold for relaxed operon rules")
    p.add_argument("--all_results", action="store_true",
                   help="Report all HMM hits; skip all operon-context filters")
    p.add_argument("--catalog_mode", action="store_true",
                   help="Gene catalog mode: bypass co-occurrence count rules "
                        "(FLEET ≥5, siderophore ≥2/3, iron_transport ≥2) but "
                        "keep bitscore cutoffs, Cyc2, and Mtr/Mto disambiguation. "
                        "Use when input is a deduplicated gene catalog where "
                        "genomic context is unavailable.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--hmm_threads", type=int, default=1)
    p.add_argument("--norm", action="store_true",
                   help="Normalise gene-count heatmap")
    p.add_argument("--bam", help="Single BAM file (requires samtools >=1.10)")
    p.add_argument("--bams", help="TSV: genome<TAB>bam_path")
    p.add_argument("--depth",
                   help="Pre-computed depth file (jgi/BBMap/samtools/plain)")
    p.add_argument("--depths", help="TSV: genome<TAB>depth_path")
    p.add_argument("--norm_coverage", action="store_true",
                   help="TPM-normalise coverage heatmap")
    p.add_argument("--keep_tblout", action="store_true")
    # ── Operon prediction (UniOP) ─────────────────────────────────────────────
    p.add_argument("--operon_prediction", action="store_true",
                   help="Run UniOP operon prediction and produce "
                        "OperonStructure.tsv. Requires --fna_dir "
                        "(or nucleotide files) and UniOP installed.")
    p.add_argument("--uniop_path", default="uniop",
                   help="Path to UniOP executable (default: 'uniop', assumed in PATH)")
    # ── Anvi'o output ─────────────────────────────────────────────────────────
    p.add_argument("--anvio", action="store_true",
                   help="Write MetalGenie-Evo-anvio-functions.tsv, compatible with "
                        "anvi-import-functions. gene_callers_id column contains ORF "
                        "names — map to integer IDs from your Anvi'o contigs database "
                        "before importing (see README).")
    p.add_argument("--bakta_gff_dir",
                   help="Directory of Bakta GFF3 files. When provided together with "
                        "--fna_dir, MetalGenie-Evo runs Prodigal internally for HMM "
                        "search and UniOP, then maps Prodigal ORF names back to Bakta "
                        "gene IDs via coordinate matching. The --anvio output will use "
                        "Bakta IDs directly, compatible with Anvi'o databases built "
                        "from Bakta external gene calls.")
    args = p.parse_args()

    out_dir    = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir = out_dir / "_tblout_cache"
    tblout_dir.mkdir(exist_ok=True)
    gff_dir_path = None

    if args.fna_dir:
        fna_dir   = Path(args.fna_dir)
        fna_files = sorted(fna_dir.glob(f"*.{args.fna_ext}"))
        if not fna_files:
            sys.exit(f"[ERROR] No .{args.fna_ext} in {fna_dir}")
        faa_dir, gff_dir_path = run_prodigal(
            fna_files, out_dir,
            meta_mode=args.meta,
            threads=args.threads,
            gene_caller=args.gene_caller)
        faa_ext = "faa"
    else:
        faa_dir = Path(args.faa_dir)
        faa_ext = args.faa_ext
        if args.gff_dir:
            gff_dir_path = Path(args.gff_dir)
        else:
            print("[WARN] --faa_dir used without --gff_dir: falling back to "
                  "ordinal-index clustering. If using Bakta FAA files, this "
                  "will produce incorrect clusters because Bakta locus tags "
                  "(e.g. CJMEHH_00001) don't encode contig boundaries. "
                  "Add --gff_dir pointing to your Bakta GFF3 files.",
                  file=sys.stderr)

    faa_files = sorted(faa_dir.glob(f"*.{faa_ext}"))
    if not faa_files:
        sys.exit(f"[ERROR] No .{faa_ext} in {faa_dir}")
    print(f"[INFO] {len(faa_files)} genome/bin FAA files")

    gene_map = read_map(str(Path(args.hmm_dir) / "MetalGenie-map.txt"))
    if not gene_map:
        gene_map = read_map(str(Path(args.hmm_dir) / "FeGenie-map.txt"))
    if not gene_map:
        print(f"[WARN] No MetalGenie-map.txt or FeGenie-map.txt found in {args.hmm_dir} — "
              f"gene names will show as raw HMM stems in all outputs.", file=sys.stderr)

    cutoffs  = read_cutoffs(str(Path(args.hmm_dir) / "HMM-bitcutoffs.txt"))
    cat_hmms = defaultdict(list)
    h2c      = {}
    for entry in sorted(Path(args.hmm_dir).iterdir()):
        if entry.is_dir() and not entry.name.startswith("."):
            for hf in sorted(entry.glob("*.hmm")):
                cat_hmms[entry.name].append((hf.stem, hf))
                h2c[hf.stem] = entry.name
    if not cat_hmms:
        sys.exit(f"[ERROR] No HMM dirs in {args.hmm_dir}")
    total_hmms = sum(len(v) for v in cat_hmms.values())
    print(f"[INFO] {total_hmms} HMMs across {len(cat_hmms)} categories")
    operon_rules, report_all_pats, json_mode = load_operon_rules(args.hmm_dir)

    contig_lengths = {}
    if args.min_contig_len > 0 or args.relaxed_operons or args.norm_coverage:
        print("[INFO] Loading contig lengths…")
        contig_lengths = build_contig_length_dict(
            faa_files,
            gff_dir=gff_dir_path,
            fna_dir=Path(args.fna_dir) if args.fna_dir else None,
            fna_ext=args.fna_ext if args.fna_dir else "fna")
        n_cl = sum(1 for f in faa_files if f.name in contig_lengths)
        print(f"       {n_cl}/{len(faa_files)} genomes have contig length data")

    genome_coords = {}
    if gff_dir_path:
        print(f"[INFO] Loading GFF from {gff_dir_path}…")
        genome_coords = load_gff_dir(gff_dir_path, faa_files)
        print(f"       {sum(1 for f in faa_files if f.name in genome_coords)}"
              f"/{len(faa_files)} with GFF")
    else:
        print("[INFO] No GFF: using ORF-index clustering (FeGenie-compatible)")

    print("[INFO] Loading protein sequences…")
    seq_dict = {faa.name: read_fasta(str(faa)) for faa in faa_files}
    print(f"[INFO] Launching hmmsearch ({args.threads} threads, "
          f"{args.hmm_threads} per job)…")
    run_all_hmmsearches(faa_files, cat_hmms, cutoffs, tblout_dir,
                        args.threads, args.hmm_threads)
    print("[INFO] Collecting best HMM hits…")
    best_hit = collect_best_hits(faa_files, cat_hmms, tblout_dir)

    print("[INFO] Clustering and filtering…")
    if args.min_contig_len > 0:
        print(f"       Skipping contigs < {args.min_contig_len} bp")
    if args.relaxed_operons:
        print(f"       Relaxed thresholds for contigs < {args.relaxed_threshold} bp")
    cluster_id = 0
    final_rows = []
    n_filt     = 0

    for faa in faa_files:
        genome    = faa.name
        orf_hits  = best_hit.get(genome, {})
        if not orf_hits:
            continue
        orf_coords = genome_coords.get(genome, {})
        clen       = contig_lengths.get(genome, {})
        if args.min_contig_len > 0 and clen:
            before   = len(orf_hits)
            orf_hits = {
                o: h for o, h in orf_hits.items()
                if clen.get(_orf_to_contig(o), args.min_contig_len + 1) >= args.min_contig_len
            }
            n_filt += before - len(orf_hits)
        raw_clusters = build_clusters(
            genome, orf_hits, orf_coords,
            args.max_gap, args.max_bp_gap, args.strand_aware)
        for orf_group in raw_clusters:
            cluster_rows = []
            for orf in orf_group:
                hit = orf_hits.get(orf)
                if hit is None:
                    continue
                contig = (orf_coords[orf]["contig"] if orf in orf_coords
                          else _orf_to_contig(orf))
                cluster_rows.append({
                    "cat":        hit["cat"],
                    "genome":     genome,
                    "contig":     contig,
                    "orf":        orf,
                    "hmm_stem":   hit["hmm_stem"],
                    "bitscore":   hit["bitscore"],
                    "cutoff":     cutoffs.get(hit["hmm_stem"], 0),
                    "evalue":     hit["evalue"],
                    "cluster_id": cluster_id,
                    "contig_len": clen.get(contig, 0),
                })
            if not cluster_rows:
                cluster_id += 1
                continue
            min_clen = (
                min(r["contig_len"] for r in cluster_rows if r["contig_len"] > 0)
                if any(r["contig_len"] > 0 for r in cluster_rows) else None)
            rel_thr = args.relaxed_threshold if args.relaxed_operons else 0
            if json_mode:
                filtered = filter_cluster_json(
                    cluster_rows, operon_rules, report_all_pats,
                    args.all_results, contig_len=min_clen,
                    relaxed_threshold=rel_thr)
            else:
                filtered = filter_cluster_fegenie(
                    cluster_rows, report_all_pats,
                    args.all_results, catalog_mode=args.catalog_mode)
            filtered = second_pass(filtered, h2c, seq_dict,
                                   args.all_results, catalog_mode=args.catalog_mode)
            for r in filtered:
                r["gene_name"]  = gene_map.get(r["hmm_stem"], r["hmm_stem"])
                r["sequence"]   = seq_dict.get(genome, {}).get(r["orf"], "")
                r["heme_motifs"] = count_heme(r["sequence"])
                final_rows.append(r)
            cluster_id += 1

    if n_filt:
        print(f"[INFO] {n_filt} ORFs removed (short contig)")
    final_rows.sort(key=lambda r: (r["cluster_id"], r["orf"]))
    norm_dict   = ({faa.name: len(seq_dict.get(faa.name, {})) for faa in faa_files}
                   if args.norm else None)
    all_genomes = sorted(f.name for f in faa_files)

    # ── Bakta ↔ Prodigal coordinate mapping ──────────────────────────────────
    # Built BEFORE writing outputs so all files use Bakta IDs when available.
    prodigal_to_bakta = {}
    if args.bakta_gff_dir:
        if gff_dir_path:
            print("[INFO] Building Prodigal↔Bakta coordinate map…")
            prodigal_faa_dir = gff_dir_path.parent / "faa"
            if not prodigal_faa_dir.exists():
                prodigal_faa_dir = out_dir / "_prodigal" / "faa"
            prodigal_to_bakta = build_prodigal_bakta_map(
                args.bakta_gff_dir, str(prodigal_faa_dir), faa_files)
            if not prodigal_to_bakta:
                print("[WARN] Prodigal↔Bakta mapping returned empty. "
                      "Check that --bakta_gff_dir contains .gff3 files whose "
                      "basenames match your genome FAA files.", file=sys.stderr)
        else:
            print("[WARN] --bakta_gff_dir requires --fna_dir so that Prodigal GFF "
                  "files are available for coordinate matching. "
                  "Bakta mapping skipped.", file=sys.stderr)

    if prodigal_to_bakta:
        n_mapped = 0
        for r in final_rows:
            b_map    = prodigal_to_bakta.get(r["genome"], {})
            bakta_id = b_map.get(r["orf"])
            if bakta_id:
                r["bakta_id"] = bakta_id
                n_mapped += 1
            else:
                r["bakta_id"] = r["orf"]
        n_total = len(final_rows)
        print(f"[INFO] Bakta IDs applied: {n_mapped}/{n_total} HMM hits mapped "
              f"({n_mapped / n_total * 100:.1f}%)")
        if n_mapped == 0:
            print(f"[WARN] No HMM hits could be mapped to Bakta IDs. "
                  f"Check that Prodigal GFF basenames match FAA basenames.",
                  file=sys.stderr)
            sample_orfs    = [r["orf"] for r in final_rows[:3]]
            sample_genome  = final_rows[0]["genome"] if final_rows else "?"
            sample_map_keys = list(prodigal_to_bakta.get(sample_genome, {}).keys())[:3]
            print(f"       Sample orf names  : {sample_orfs}", file=sys.stderr)
            print(f"       Sample map keys   : {sample_map_keys}", file=sys.stderr)
    else:
        for r in final_rows:
            r["bakta_id"] = r["orf"]
        print("[INFO] No Bakta mapping — orf column contains Prodigal ORF names")

    # summary and gene-summary don't carry uniop_context — write now
    for path, fn in [(out_dir / "MetalGenie-Evo-summary.csv", write_summary),
                     (out_dir / "MetalGenie-Evo-geneSummary-clusters.csv", write_gene_summary)]:
        print(f"[INFO] Writing {path.name}…")
        fn(str(path), final_rows)
    # long-format is written after UniOP so uniop_context column is populated

    print("[INFO] Writing MetalGenie-Evo-heatmap-data.csv…")
    write_heatmap(str(out_dir / "MetalGenie-Evo-heatmap-data.csv"),
                  final_rows, all_genomes, norm_dict)

    bam_map   = {}
    depth_map = {}
    if args.bams:
        bam_map = load_bams_tsv(args.bams)
    elif args.bam:
        for faa in faa_files:
            bam_map[faa.name] = args.bam
            bam_map[faa.stem] = args.bam
    if args.depths:
        depth_map = load_bams_tsv(args.depths)
    elif args.depth:
        for faa in faa_files:
            depth_map[faa.name] = args.depth
            depth_map[faa.stem] = args.depth
    if bam_map or depth_map:
        print("[INFO] Computing coverage…")
        gc = build_contig_coverage(faa_files, bam_map, depth_map, out_dir)
        if gc:
            print("[INFO] Writing MetalGenie-Evo-coverage-heatmap.csv…")
            write_coverage_heatmap(
                str(out_dir / "MetalGenie-Evo-coverage-heatmap.csv"),
                final_rows, all_genomes, gc,
                norm_coverage=args.norm_coverage,
                contig_lengths=contig_lengths if args.norm_coverage else None)
        else:
            print("[WARN] No coverage data loaded.", file=sys.stderr)

    if not args.keep_tblout:
        shutil.rmtree(tblout_dir, ignore_errors=True)
    else:
        print(f"[INFO] tblout cache at {tblout_dir}/")

    # ── UniOP operon prediction (optional) ────────────────────────────────────
    genome_operon_map = {}
    if args.operon_prediction:
        uniop_script = Path(args.uniop_path)
        if uniop_script.is_dir():
            candidate = uniop_script / "src" / "UniOP"
            if candidate.exists():
                print(f"[INFO] --uniop_path is a directory — auto-resolved to {candidate}")
                uniop_script = candidate
            else:
                print(f"[ERROR] --uniop_path '{args.uniop_path}' is a directory and "
                      f"UniOP script not found at {candidate}.\n"
                      f"        Pass the full path to the script, e.g.:\n"
                      f"        --uniop_path {args.uniop_path}/src/UniOP",
                      file=sys.stderr)
                uniop_script = None

        if uniop_script is not None:
            fna_dir_for_uniop = Path(args.fna_dir) if args.fna_dir else None
            print(f"[INFO] Running UniOP on {len(faa_files)} genomes…")
            print(f"       UniOP script: {uniop_script}")
            genome_operon_map = run_uniop(
                faa_files,
                fna_dir          = fna_dir_for_uniop,
                out_dir          = out_dir,
                uniop_path       = str(uniop_script),
                fna_ext          = args.fna_ext if args.fna_dir else "fna",
                prodigal_faa_dir = str(faa_dir) if args.fna_dir else None,
                bakta_gff_dir    = args.bakta_gff_dir,
            )
            if genome_operon_map:
                # Key by (genome, op_id) to prevent cross-genome operon ID collision.
                op_to_hits = defaultdict(list)
                for r in final_rows:
                    op_id = genome_operon_map.get(r["genome"], {}).get(r["orf"])
                    if op_id:
                        op_to_hits[(r["genome"], op_id)].append(r["orf"])
                for r in final_rows:
                    r["uniop_context"] = _uniop_context(
                        r["orf"], r["genome"], genome_operon_map, op_to_hits)

                op_path = out_dir / "MetalGenie-Evo-OperonStructure.tsv"
                print(f"[INFO] Writing {op_path.name}…")
                write_operon_structure(str(op_path), final_rows, genome_operon_map,
                                       prodigal_to_bakta=prodigal_to_bakta)
                ctx_counts = defaultdict(int)
                for r in final_rows:
                    ctx_counts[r.get("uniop_context", "not_in_operon")] += 1
                print(f"       in_operon_with_other_hits : "
                      f"{ctx_counts['in_operon_with_other_hits']}")
                print(f"       singleton_in_operon       : "
                      f"{ctx_counts['singleton_in_operon']}")
                print(f"       not_in_operon             : "
                      f"{ctx_counts['not_in_operon']}")
            else:
                print("[WARN] UniOP produced no predictions — "
                      "see MetalGenie-Evo-run.log for details.", file=sys.stderr)

    # Write long-format after UniOP so uniop_context column is included when available
    long_path = out_dir / "MetalGenie-Evo-results-long.tsv"
    print(f"[INFO] Writing {long_path.name}…")
    write_long_format(str(long_path), final_rows)

    # ── Anvi'o functions output (optional) ───────────────────────────────────
    if args.anvio:
        anvio_path = out_dir / "MetalGenie-Evo-anvio-functions.tsv"
        print(f"[INFO] Writing {anvio_path.name}…")
        write_anvio_functions(str(anvio_path), final_rows,
                              prodigal_to_bakta=prodigal_to_bakta)
        id_note = ("Bakta gene IDs" if prodigal_to_bakta
                   else "Prodigal ORF names — map to int IDs before import")
        print(f"       gene_callers_id: {id_note}")
        print(f"       Import with: anvi-import-functions -c CONTIGS.db "
              f"-i {anvio_path.name} -p MetalGenie-Evo")

    # ── Write run log ─────────────────────────────────────────────────────────
    log_path = out_dir / "MetalGenie-Evo-run.log"
    cc2 = defaultdict(int)
    for r in final_rows:
        cc2[r["cat"]] += 1
    with open(log_path, "w") as lf:
        lf.write("MetalGenie-Evo run log\n")
        lf.write(f"Date          : "
                 f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        lf.write(f"Command       : {' '.join(sys.argv)}\n")
        lf.write(f"Output dir    : {out_dir}\n")
        lf.write(f"Genomes input : {len(faa_files)}\n")
        lf.write(f"Genomes hit   : {len({r['genome'] for r in final_rows})}\n")
        lf.write(f"Total ORFs    : {len(final_rows)}\n")
        lf.write(
            f"Bakta mapping : "
            f"{'yes (' + str(sum(len(v) for v in prodigal_to_bakta.values())) + ' ORFs mapped)' if prodigal_to_bakta else 'no'}\n")
        lf.write(
            f"UniOP         : "
            f"{'yes (' + str(len(genome_operon_map)) + ' genomes)' if genome_operon_map else 'no/failed'}\n")
        lf.write("\nHits per category:\n")
        for cat, n in sorted(cc2.items()):
            lf.write(f"  {n:5d}  {cat}\n")
    print(f"[INFO] Run log → {log_path}")

    n_hit = len({r["genome"] for r in final_rows})
    cc    = defaultdict(int)
    for r in final_rows:
        cc[r["cat"]] += 1
    print(f"\n{'─' * 60}\n  MetalGenie-Evo  —  run complete")
    print(f"  {len(final_rows)} ORFs in {n_hit}/{len(faa_files)} genomes")
    print(f"\n  Hits per category:")
    for cat, n in sorted(cc.items()):
        print(f"    {n:5d}  {cat}")
    print(f"\n  Outputs  →  {out_dir}/\n{'─' * 60}")


if __name__ == "__main__":
    main()
