"""UniOP operon prediction: running, parsing, and output writing."""

import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def _parse_uniop_faa_index(faa_path):
    """
    Build dict: 1-based integer index → prodigal_orf_name
    from the FAA file produced by UniOP's internal Prodigal run.
    """
    idx_to_orf = {}
    idx = 0
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                idx += 1
                orf_id = line[1:].split()[0].strip()
                idx_to_orf[idx] = orf_id
    return idx_to_orf


def _parse_uniop_operon(operon_path, idx_to_orf):
    """
    Parse uniop.operon CSV:
      idx_genes,idx_op
      "[np.int64(1), np.int64(2), np.int64(3)]",0

    Returns dict: orf_name → operon_id (string)
    """
    orf_to_op = {}
    try:
        with open(operon_path) as fh:
            fh.readline()   # skip header
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                last_comma = line.rfind(",")
                if last_comma < 0:
                    continue
                idx_genes_str = line[:last_comma].strip().strip('"')
                idx_op_str    = line[last_comma + 1:].strip()
                # Extract gene indices from numpy array string.
                # e.g. "[np.int64(1), np.int64(2), np.int64(3)]" → [1, 2, 3]
                # Match \(digits\) to avoid capturing "64" from "np.int64".
                indices = [int(x) for x in re.findall(r'\((\d+)\)', idx_genes_str)]
                # Fallback: plain "[1, 2, 3]" format (no type prefix)
                if not indices:
                    indices = [int(x) for x in re.findall(r'\b(\d+)\b', idx_genes_str)]
                op_id = f"OP{int(idx_op_str):04d}"
                for idx in indices:
                    orf = idx_to_orf.get(idx)
                    if orf:
                        orf_to_op[orf] = op_id
    except Exception as e:
        print(f"  [WARN] Could not parse uniop.operon: {e}", file=sys.stderr)
    return orf_to_op


def _parse_uniop_pred(pred_path, idx_to_orf, threshold=0.5):
    """
    Parse uniop.pred (fallback when uniop.operon is empty):
      Gene A  Gene B  Prediction   (tab-separated, some rows have empty Prediction)

    Uses union-find to build connected components of pairs with prob > threshold.
    Returns dict: orf_name → operon_id
    """
    parent = {}

    def _find(x):
        root = x
        while parent.get(root, root) != root:
            root = parent.get(root, root)
        while x != root:
            parent[x], x = root, parent.get(x, x)
        return root

    def _union(a, b):
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[rb] = ra

    try:
        with open(pred_path) as fh:
            next(fh, None)  # skip header "Gene A  Gene B  Prediction"
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 3 or not parts[2].strip():
                    continue
                try:
                    ia   = int(parts[0].strip())
                    ib   = int(parts[1].strip())
                    prob = float(parts[2].strip())
                except ValueError:
                    continue
                if prob >= threshold:
                    oa = idx_to_orf.get(ia)
                    ob = idx_to_orf.get(ib)
                    if oa and ob:
                        parent.setdefault(oa, oa)
                        parent.setdefault(ob, ob)
                        _union(oa, ob)
    except Exception as e:
        print(f"  [WARN] Could not parse uniop.pred: {e}", file=sys.stderr)
        return {}

    comp_ids   = {}
    op_counter = 0
    orf_to_op  = {}
    for orf in list(parent.keys()):
        root = _find(orf)
        if root not in comp_ids:
            comp_ids[root] = f"OP{op_counter:04d}"
            op_counter += 1
        orf_to_op[orf] = comp_ids[root]
    return orf_to_op


def make_uniop_faa(bakta_faa_path, bakta_gff_path, out_faa_path):
    """
    Create a UniOP-compatible FAA from a Bakta FAA + Bakta GFF3.

    Bakta FAA headers:  >CJMEHH_00001 hypothetical protein
    UniOP needs:        >CJMEHH_00001 # start # end # strand # ID=...

    Reads coordinates from the Bakta GFF3 (locus_tag= attribute) and rewrites
    the FAA headers in Prodigal-compatible format. Sequences are unchanged.
    Returns True if successful, False if no coordinates could be matched.
    """
    coord_map = {}
    with open(bakta_gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue
            strand    = 1 if parts[6] == "+" else -1
            locus_tag = None
            for attr in parts[8].split(";"):
                attr = attr.strip()
                if attr.startswith("locus_tag="):
                    locus_tag = attr[10:].strip()
                    break
                if attr.startswith("ID=") and locus_tag is None:
                    locus_tag = attr[3:].strip()
            if locus_tag:
                coord_map[locus_tag] = (contig, start, end, strand)

    if not coord_map:
        return False

    matched = 0
    with open(bakta_faa_path) as fh_in, open(out_faa_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                current_id = line[1:].split()[0].strip()
                coords = coord_map.get(current_id)
                if coords:
                    contig, start, end, strand = coords
                    fh_out.write(
                        f">{current_id} # {start} # {end} # {strand} # "
                        f"ID={current_id};partial=00;start_type=ATG;"
                        f"rbs_motif=None;gc_cont=0.500\n"
                    )
                    matched += 1
                else:
                    fh_out.write(line)
            else:
                fh_out.write(line)

    return matched > 0


def run_uniop(faa_files, fna_dir, out_dir, uniop_path, fna_ext="fna",
              prodigal_faa_dir=None, bakta_gff_dir=None):
    """
    Run UniOP on each genome.

    Input mode priority:
      1. pyrodigal FAA (--fna_dir used): headers already Prodigal-compatible
         → pass directly to UniOP with -a
      2. Bakta FAA + GFF3 (--faa_dir + --bakta_gff_dir): reconstruct headers
         from GFF3 coordinates → write temp FAA → pass to UniOP with -a
      3. FNA fallback: UniOP runs Prodigal internally with -i

    Returns dict: genome_faa_name → {orf_id → operon_id}
    """
    uniop_dir = out_dir / "_uniop"
    uniop_dir.mkdir(exist_ok=True)
    genome_operon_map = {}

    for faa in faa_files:
        stem     = faa.stem
        work_dir = uniop_dir / stem
        work_dir.mkdir(exist_ok=True)

        operon_file = work_dir / "uniop.operon"
        pred_file   = work_dir / "uniop.pred"

        faa_for_uniop = None
        faa_for_index = None

        if prodigal_faa_dir:
            candidate = Path(prodigal_faa_dir) / f"{stem}.faa"
            if candidate.exists():
                faa_for_uniop = candidate
                faa_for_index = candidate

        if faa_for_uniop is None and bakta_gff_dir:
            bakta_gff = None
            for ext in (".gff3", ".gff"):
                c = Path(bakta_gff_dir) / (stem + ext)
                if c.exists():
                    bakta_gff = c
                    break
            if bakta_gff and faa.exists():
                temp_faa = work_dir / f"{stem}_uniop.faa"
                if not temp_faa.exists():
                    ok = make_uniop_faa(str(faa), str(bakta_gff), str(temp_faa))
                    if not ok:
                        print(f"  [WARN] UniOP: could not build temp FAA for {stem}",
                              file=sys.stderr)
                if temp_faa.exists():
                    faa_for_uniop = temp_faa
                    faa_for_index = temp_faa

        if operon_file.exists():
            print(f"  [INFO] UniOP cache hit for {stem}")
        else:
            if faa_for_uniop:
                cmd = ["python3", str(uniop_path),
                       "-a", str(faa_for_uniop.resolve()),
                       "-t", str(work_dir.resolve())]
            else:
                fna_path = None
                if fna_dir:
                    for ext in [fna_ext, "fna", "fasta", "fa"]:
                        candidate = Path(fna_dir) / f"{stem}.{ext}"
                        if candidate.exists():
                            fna_path = candidate
                            break
                if fna_path is None:
                    print(f"  [WARN] UniOP: no FAA or FNA found for {stem} — skipping.",
                          file=sys.stderr)
                    continue
                cmd = ["python3", str(uniop_path),
                       "-i", str(fna_path.resolve()),
                       "-t", str(work_dir.resolve())]

            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  [WARN] UniOP failed for {stem}:\n{r.stderr[:400]}",
                      file=sys.stderr)
                err_log = work_dir / "uniop_error.log"
                err_log.write_text(
                    f"Command: {' '.join(cmd)}\n\nSTDOUT:\n{r.stdout}\n\nSTDERR:\n{r.stderr}")
                print(f"         Full error → {err_log}", file=sys.stderr)
                continue

        if faa_for_index is None:
            faa_candidates = list(work_dir.glob("*.faa"))
            if not faa_candidates:
                print(f"  [WARN] UniOP: no FAA for index in {work_dir}",
                      file=sys.stderr)
                continue
            faa_for_index = faa_candidates[0]

        idx_to_orf = _parse_uniop_faa_index(str(faa_for_index))

        orf_to_op = {}
        if operon_file.exists():
            orf_to_op = _parse_uniop_operon(str(operon_file), idx_to_orf)
        if not orf_to_op and pred_file.exists():
            orf_to_op = _parse_uniop_pred(str(pred_file), idx_to_orf)

        n_operons = len(set(orf_to_op.values()))
        print(f"  [INFO] {stem}: {n_operons} operons, {len(orf_to_op)} genes assigned")

        genome_operon_map[faa.name] = orf_to_op

    return genome_operon_map


def parse_gff_coords(gff_path, source_hint=""):
    """
    Parse a GFF/GFF3 file and return:
      dict: contig → [(start, end, strand, gene_id)]

    Works with both Prodigal GFF (ID=contig_1_5) and Bakta GFF3
    (ID=AMXMAG_00053). source_hint is only used for logging.
    """
    coords = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue
            strand  = parts[6]
            gene_id = None
            for attr in parts[8].split(";"):
                attr = attr.strip()
                if attr.startswith("ID="):
                    gene_id = attr[3:].strip()
                    break
                if attr.startswith("locus_tag=") and gene_id is None:
                    gene_id = attr[10:].strip()
            if gene_id:
                coords[contig].append((start, end, strand, gene_id))
    return dict(coords)


def build_prodigal_bakta_map(bakta_gff_dir, prodigal_faa_dir, faa_files,
                              coord_tol=30):
    """
    Build mapping: genome → {prodigal_orf_id → bakta_locus_tag}

    Matches Prodigal ORFs to Bakta genes by (contig, start, end, strand).
    coord_tol=30bp covers alternative start codon choices between Prodigal
    and Pyrodigal (Bakta), where the end coordinate is usually identical
    but the start may differ by up to ~30bp.

    Matching priority:
      1. Exact (start, end, strand)
      2. Fuzzy: same strand, |start_diff| <= coord_tol, |end_diff| <= coord_tol
      3. End-anchor: same strand, same end, |start_diff| <= coord_tol
         (most common case: alternative start codon, same stop codon)
    """
    prodigal_to_bakta = {}
    stats = {"matched": 0, "unmatched": 0}

    for faa_file in faa_files:
        stem   = faa_file.stem
        genome = faa_file.name

        bakta_gff = None
        for ext in (".gff3", ".gff"):
            c = Path(bakta_gff_dir) / (stem + ext)
            if c.exists():
                bakta_gff = c
                break
        if bakta_gff is None:
            continue

        bakta_index = defaultdict(dict)
        with open(bakta_gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                contig = parts[0]
                try:
                    start = int(parts[3])
                    end   = int(parts[4])
                except ValueError:
                    continue
                strand    = parts[6]
                locus_tag = None
                id_val    = None
                for attr in parts[8].split(";"):
                    attr = attr.strip()
                    if attr.startswith("locus_tag="):
                        locus_tag = attr[10:].strip()
                        break
                    if attr.startswith("ID=") and id_val is None:
                        id_val = attr[3:].strip()
                if locus_tag is None:
                    locus_tag = id_val
                if locus_tag:
                    bakta_index[contig][(start, end, strand)] = locus_tag

        prodigal_faa = Path(prodigal_faa_dir) / f"{stem}.faa"
        if not prodigal_faa.exists():
            continue

        orf_map = {}
        with open(prodigal_faa) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                parts = line[1:].rstrip().split(" # ")
                if len(parts) < 4:
                    continue
                orf_id = parts[0].strip()
                try:
                    start      = int(parts[1].strip())
                    end        = int(parts[2].strip())
                    strand_int = int(parts[3].strip())
                except ValueError:
                    continue
                strand   = "+" if strand_int == 1 else "-"
                attr_str = parts[4].strip() if len(parts) > 4 else ""
                is_edge  = "start_type=Edge" in attr_str

                contig_parts = orf_id.rsplit("_", 1)
                contig = contig_parts[0] if len(contig_parts) == 2 else orf_id

                b_idx = bakta_index.get(contig, {})
                key   = (start, end, strand)
                if key in b_idx:
                    orf_map[orf_id] = b_idx[key]
                    stats["matched"] += 1
                    continue

                found = None
                for (bs, be, bst), bid in b_idx.items():
                    if (bst == strand and
                            abs(bs - start) <= coord_tol and
                            abs(be - end)   <= coord_tol):
                        found = bid
                        break

                if not found:
                    for (bs, be, bst), bid in b_idx.items():
                        if bst == strand and be == end and abs(bs - start) <= coord_tol:
                            found = bid
                            break

                if not found and is_edge:
                    for (bs, be, bst), bid in b_idx.items():
                        if bst == strand and be == end:
                            found = bid
                            break

                if found:
                    orf_map[orf_id] = found
                    stats["matched"] += 1
                else:
                    stats["unmatched"] += 1

        prodigal_to_bakta[genome] = orf_map

    total = stats["matched"] + stats["unmatched"]
    if total > 0:
        pct = stats["matched"] / total * 100
        print(f"[INFO] Prodigal↔Bakta mapping: "
              f"{stats['matched']}/{total} ORFs matched ({pct:.1f}%)")
        if stats["unmatched"] > 0:
            print(f"       {stats['unmatched']} Prodigal ORFs had no Bakta match "
                  f"(will keep Prodigal name)")
    return prodigal_to_bakta


def _uniop_context(orf, genome, genome_operon_map, op_to_orfs_with_hits):
    """
    Classify the UniOP context of an HMM hit:

    - 'in_operon_with_other_hits'  : in a UniOP operon AND at least one other
                                     ORF in that operon also has an HMM hit
    - 'singleton_in_operon'        : in a UniOP operon but no other HMM hits
                                     in the same operon
    - 'not_in_operon'              : not assigned to any UniOP operon

    Keys in op_to_orfs_with_hits must be (genome, op_id) tuples to avoid
    collision between identically-numbered operons from different genomes.
    """
    op_map = genome_operon_map.get(genome, {})
    op_id  = op_map.get(orf)
    if op_id is None:
        return "not_in_operon"
    other_hits = [o for o in op_to_orfs_with_hits.get((genome, op_id), []) if o != orf]
    return "in_operon_with_other_hits" if other_hits else "singleton_in_operon"


def write_operon_structure(path, final_rows, genome_operon_map,
                           prodigal_to_bakta=None):
    """
    OperonStructure.tsv — HMM hits linked to UniOP operon predictions.

    uniop_context values:
      in_operon_with_other_hits  — operon contains ≥1 other HMM-positive ORF
      singleton_in_operon        — operon contains no other HMM-positive ORFs
      not_in_operon              — ORF not assigned to any UniOP operon
    """
    # Key by (genome, op_id) to prevent collision across genomes.
    op_to_orfs_with_hits = defaultdict(list)
    for r in final_rows:
        op_map = genome_operon_map.get(r["genome"], {})
        op_id  = op_map.get(r["orf"])
        if op_id:
            op_to_orfs_with_hits[(r["genome"], op_id)].append(r["orf"])

    use_bakta = bool(prodigal_to_bakta)
    fields = ["operon_id", "genome", "contig", "orf", "gene", "category",
              "hmm_stem", "bitscore", "e_value",
              "uniop_context", "unioperon_members"]
    if use_bakta:
        fields.insert(fields.index("orf") + 1, "bakta_gene_id")

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in final_rows:
            genome = r["genome"]
            orf    = r["orf"]
            op_map = genome_operon_map.get(genome, {})
            op_id  = op_map.get(orf, "")
            members = [o for o in op_to_orfs_with_hits.get((genome, op_id), [])
                       if o != orf]
            ctx = _uniop_context(orf, genome, genome_operon_map,
                                 op_to_orfs_with_hits)
            row = {
                "operon_id":         op_id if op_id else "no_operon",
                "genome":            genome,
                "contig":            r["contig"],
                "orf":               orf,
                "gene":              r["gene_name"],
                "category":          r["cat"],
                "hmm_stem":          r["hmm_stem"],
                "bitscore":          f"{r['bitscore']:.1f}",
                "e_value":           f"{r['evalue']:.2e}",
                "uniop_context":     ctx,
                "unioperon_members": ",".join(members) if members else "",
            }
            if use_bakta:
                row["bakta_gene_id"] = prodigal_to_bakta.get(genome, {}).get(orf, orf)
            w.writerow(row)
