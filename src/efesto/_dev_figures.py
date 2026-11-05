"""Private export utility; deliberately separate from the public CLI."""
import argparse
import csv
import hashlib
import json
import math
import re
import textwrap
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from efesto.categories import canonical_category


def read_results(path):
    required = {"genome", "contig", "orf", "category", "gene", "cluster_id", "bitscore_cutoff"}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Missing result columns: " + ", ".join(sorted(missing)))
        seen, rows = {}, []
        for number, row in enumerate(reader, 2):
            if not row["orf"] or not row["genome"]:
                raise ValueError("Missing genome/ORF on row %d" % number)
            key = (row["genome"], row["contig"], row["orf"])
            if key in seen:
                if row != seen[key]:
                    raise ValueError("Conflicting assignments for %s; expected final best-hit output" % (key,))
                continue
            original = dict(row)
            try:
                cutoff = float(row["bitscore_cutoff"] or 0)
                if not math.isfinite(cutoff) or cutoff < 0:
                    raise ValueError("Invalid cutoff")
                row["_specific_cutoff"] = cutoff > 0
            except ValueError as exc:
                raise ValueError("Invalid cutoff on row %d" % number) from exc
            seen[key] = original
            row["category"] = canonical_category(row["category"])
            rows.append(row)
    return rows


def read_gff(path):
    coords = defaultdict(set)
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = dict(item.split("=", 1) for item in parts[8].split(";") if "=" in item)
            identity = unquote(attrs.get("ID", ""))
            if not identity:
                continue
            start, end = int(parts[3]), int(parts[4])
            if start < 1 or end < start:
                raise ValueError("Invalid CDS coordinates in " + str(path))
            value = (start, end, parts[6])
            coords[(parts[0], identity)].add(value)
            # Prodigal names FASTA records contig_index and uses seqnum_index
            # in GFF ID. Restrict this alias to its explicit source/ID schema.
            if parts[1].lower().startswith(("prodigal", "pyrodigal")) and re.fullmatch(r"\d+_\d+", identity):
                coords[(parts[0], parts[0] + "_" + identity.rsplit("_", 1)[1])].add(value)
    return coords


def attach_coordinates(rows, results_dir, gff_dir=None):
    used, warnings, cache = [], [], {}
    global_gff = Path(results_dir) / "Efesto-hits.gff3"
    global_coords = read_gff(global_gff) if global_gff.exists() else {}
    if global_gff.exists():
        used.append(global_gff)
    roots = [Path(gff_dir)] if gff_dir else [Path(results_dir) / "_prodigal" / "gff"]
    for row in rows:
        genome = row["genome"]
        if genome not in cache:
            cache[genome] = {}
            for root in roots:
                for suffix in (".gff", ".gff3", ".prodigal.gff"):
                    path = root / (Path(genome).stem + suffix)
                    if path.is_file():
                        cache[genome] = read_gff(path)
                        used.append(path)
                        break
                if cache[genome]:
                    break
        key = (row["contig"], row["orf"])
        values = cache[genome].get(key) or global_coords.get(key, set())
        if len(values) == 1:
            row["_coordinates"] = next(iter(values))
        elif len(values) > 1:
            warnings.append("Ambiguous GFF coordinates: %s/%s" % (genome, row["orf"]))
    return sorted(set(used)), warnings


def aggregate(rows):
    counts = Counter((r["genome"], r["category"]) for r in rows)
    by_category = Counter(r["category"] for r in rows)
    by_genome = Counter(r["genome"] for r in rows)
    cutoff = Counter("Model-specific cutoff" if r["_specific_cutoff"] else "Generic cutoff" for r in rows)
    return counts, by_category, by_genome, cutoff


def select_clusters(rows, limit):
    groups = defaultdict(list)
    for row in rows:
        if row["cluster_id"] not in ("", "."):
            groups[(row["genome"], row["contig"], row["cluster_id"])].append(row)
    candidates = [(key, group) for key, group in groups.items() if len(group) >= 2]
    marker = re.compile(r"\b(?:MtrA|MtoA|CbcL|DA_402|Cyc2|ImoA|CymA)\b", re.I)
    def rank(item):
        key, group = item
        has_coords = all("_coordinates" in r for r in group)
        redox = any(marker.search(r["gene"]) for r in group)
        return (-int(has_coords), -int(redox), -min(len(group), 14), key)
    return sorted(candidates, key=rank)[:limit], len(candidates)


def _label(value):
    return Path(value).stem.replace("_", " ")


def export_figures(results_dir, out_dir, *, gff_dir=None, genome=None, max_clusters=6,
                   top_categories=24, top_genomes=30, dpi=300, title="Efesto"):
    # Import plotting libraries only in this private command, never on a normal run.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrow, Patch
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError as exc:
        raise RuntimeError("Optional plotting dependency missing. Install matplotlib in this environment.") from exc
    source = Path(results_dir) / "Efesto-results-long.tsv"
    rows = read_results(source)
    if genome:
        rows = [r for r in rows if r["genome"] == genome or Path(r["genome"]).stem == genome]
        if not rows:
            raise ValueError("No result rows match genome " + genome)
    out = Path(out_dir)
    if out.resolve() == Path(results_dir).resolve():
        raise ValueError("Choose a separate figure output directory")
    out.mkdir(parents=True, exist_ok=True)
    gffs, warnings = attach_coordinates(rows, results_dir, gff_dir)
    counts, categories, genomes, cutoff = aggregate(rows)
    heatmap_source = Path(results_dir) / "Efesto-heatmap-data.csv"
    if heatmap_source.is_file() and not genome:
        with heatmap_source.open(newline="", encoding="utf-8-sig") as handle:
            header = next(csv.reader(handle), [])
        for sample in header[1:]:
            if sample:
                genomes.setdefault(sample, 0)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "axes.titlesize": 20,
        "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white"})
    navy, blue, pale, gray = "#103E62", "#2486B5", "#D9EEF6", "#667580"
    palette = sum((list(plt.get_cmap(name).colors) for name in ("tab20", "tab20b", "tab20c")), [])
    category_colors = {c: palette[i % len(palette)] for i, c in enumerate(sorted(categories))}
    exports, selected_rows = [], []
    def save(fig, name):
        for extension in ("png", "svg", "pdf"):
            path = out / (name + "." + extension)
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            exports.append(path.name)
        plt.close(fig)
    def footer(fig, text="Predicted genomic potential; annotations do not demonstrate biological activity."):
        fig.text(.02, .015, text, fontsize=9, color=gray)
    top_c = sorted(categories, key=lambda c: (-categories[c], c))[:top_categories]
    top_g = sorted(genomes, key=lambda g: (-genomes[g], g))[:top_genomes]

    fig, ax = plt.subplots(figsize=(13.5, max(5.5, .33 * len(top_c) + 1.8)))
    if rows:
        ordered = list(reversed(top_c))
        ax.barh(range(len(ordered)), [categories[c] for c in ordered], color=blue)
        ax.set_yticks(range(len(ordered)), [c.replace("-", " · ", 1) for c in ordered])
        for y, c in enumerate(ordered):
            ax.text(categories[c], y, "  %d" % categories[c], va="center", fontsize=10)
        ax.margins(x=.13)
    else:
        ax.text(.5, .5, "No annotations in this result set", transform=ax.transAxes, ha="center")
    ax.set_xlabel("Annotated ORFs (raw counts)")
    ax.set_title(title + " | Functional annotations", loc="left", color=navy, pad=28)
    ax.text(0, 1.025, "%d ORFs · %d sample%s · showing %d of %d categories" %
            (len(rows), len(genomes), "" if len(genomes)==1 else "s", len(top_c), len(categories)), transform=ax.transAxes, color=gray, fontsize=11)
    ax.grid(axis="x", alpha=.15);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.39, bottom=.13, top=.85, right=.93)
    footer(fig);save(fig, "01_functional_overview")

    fig, ax = plt.subplots(figsize=(max(9.5, .40 * len(top_g) + 6), max(5.5, .30 * len(top_c) + 2)))
    if top_c and top_g:
        matrix = [[math.log2(1 + counts[(g,c)]) for g in top_g] for c in top_c]
        cmap = LinearSegmentedColormap.from_list("efesto_blue", ["#FFFFFF", pale, blue, navy])
        plot=ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0)
        ax.set_xticks(range(len(top_g)), [_label(g) for g in top_g], rotation=45, ha="right", fontsize=10)
        ax.set_yticks(range(len(top_c)), top_c, fontsize=10)
        fig.colorbar(plot, ax=ax, fraction=.025, pad=.025, label="log₂(1 + annotated ORF count)")
        if len(top_c) * len(top_g) <= 180:
            maximum=max(max(r) for r in matrix) or 1
            for i,c in enumerate(top_c):
                for j,g in enumerate(top_g):
                    n=counts[(g,c)]
                    if n:ax.text(j,i,str(n),ha="center",va="center",fontsize=9,color="white" if matrix[i][j]>maximum*.6 else navy)
    else:
        ax.text(.5,.5,"No annotations to display",transform=ax.transAxes,ha="center")
    ax.set_title(title + " | Category profiles",loc="left",color=navy,pad=30)
    ax.text(0,1.025,"Top %d of %d samples · %d of %d categories · raw counts, no abundance normalization" %
            (len(top_g),len(genomes),len(top_c),len(categories)),transform=ax.transAxes,fontsize=10,color=gray)
    fig.subplots_adjust(left=.30,bottom=.28,top=.86,right=.91)
    footer(fig);save(fig,"02_category_heatmap")

    fig,ax=plt.subplots(figsize=(10,5))
    labels=["Model-specific cutoff","Generic cutoff"]
    ax.barh(labels,[cutoff[x] for x in labels],color=[blue,"#E3AA50"])
    for i,label in enumerate(labels):ax.text(cutoff[label],i,"  %d" % cutoff[label],va="center")
    ax.set_xlabel("Annotated ORFs");ax.margins(x=.15)
    ax.set_title(title+" | Threshold provenance",loc="left",color=navy,pad=22)
    ax.grid(axis="x",alpha=.15);ax.set_axisbelow(True)
    fig.subplots_adjust(left=.25,bottom=.22,top=.78,right=.93)
    footer(fig,"Cutoff provenance is not a measurement of classifier accuracy or a validation probability.")
    save(fig,"03_threshold_provenance")

    clusters,total_clusters=select_clusters(rows,max_clusters)
    cluster_manifest=[]
    for number,(key,group) in enumerate(clusters,1):
        genome_name,contig,cluster=key
        with_coords=all("_coordinates" in row for row in group)
        ordered=sorted(group,key=lambda r:r["_coordinates"][0]) if with_coords else list(group)
        selected=ordered[:14]
        fig,ax=plt.subplots(figsize=(15.5,5.5))
        if with_coords:
            origin=selected[0]["_coordinates"][0];finish=max(r["_coordinates"][1] for r in selected);span=max(finish-origin+1,1)
            ax.axhline(0,color="#CBD5DC",lw=1,zorder=0)
            for i,row in enumerate(selected):
                start,end,strand=row["_coordinates"];length=end-start+1
                direction=-1 if strand=="-" else 1
                x=start-origin if direction==1 else end-origin+1
                ax.add_patch(FancyArrow(x,0,direction*length,0,width=.24,head_width=.43,
                    head_length=min(length*.4,span*.025),length_includes_head=True,
                    facecolor=category_colors[row["category"]],edgecolor="white",linewidth=.6))
                if strand not in ("+","-"):
                    # Unknown-strand features must not be assigned an invented direction.
                    ax.patches[-1].remove()
                    ax.plot([start-origin,end-origin],[0,0],lw=12,color=category_colors[row["category"]],solid_capstyle="butt")
                label="\n".join(textwrap.wrap(row["gene"],22)[:3])
                ax.text((start+end)/2-origin,.45+(i%3)*.48,label,ha="center",va="bottom",fontsize=9)
            ax.set_xlim(-span*.04,span*1.04);ax.set_xlabel("Position relative to first displayed CDS (bp)")
            note="Genomic coordinates · only annotated genes shown"
        else:
            for i,row in enumerate(selected):
                ax.scatter(i,0,s=400,marker="s",color=category_colors[row["category"]])
                ax.text(i,.45+(i%3)*.48,"\n".join(textwrap.wrap(row["gene"],20)[:3]),ha="center",va="bottom",fontsize=9)
            ax.set_xlim(-.7,len(selected)-.3);ax.set_xticks([])
            note="Cluster membership only · coordinates unavailable; spacing and table order are schematic"
        ax.set_ylim(-1.3,2.7);ax.set_yticks([])
        for side in ("top","right","left"):ax.spines[side].set_visible(False)
        if not with_coords:ax.spines["bottom"].set_visible(False)
        ax.set_title("%s | Cluster %s" % (_label(genome_name),cluster),loc="left",color=navy,pad=24)
        ax.text(0,1.01,contig+" · "+note,transform=ax.transAxes,fontsize=10,color=gray)
        cats=sorted({r["category"] for r in selected})
        ax.legend(handles=[Patch(facecolor=category_colors[c],label=c) for c in cats],
            loc="lower left",frameon=False,ncol=min(3,len(cats)),fontsize=9)
        fig.subplots_adjust(left=.04,right=.98,bottom=.18,top=.80)
        footer(fig,"%d of %d annotated members shown. Neighbourhood does not establish a functional operon." % (len(selected),len(group)))
        save(fig,"cluster_%02d" % number)
        cluster_manifest.append({'figure':"cluster_%02d" % number,'genome':genome_name,'contig':contig,
            'cluster_id':cluster,'members_total':len(group),'members_displayed':len(selected),'genomic_coordinates':with_coords})
        for row in selected:
            start,end,strand=row.get("_coordinates",("","",""))
            selected_rows.append({'figure':"cluster_%02d" % number,'genome':genome_name,'contig':contig,
                'orf':row['orf'],'gene':row['gene'],'category':row['category'],'start':start,'end':end,'strand':strand})
    if total_clusters == 0:warnings.append("No multi-gene clusters found; no cluster example figures generated.")
    with (out/"figure_category_counts.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(["genome","category","annotated_orfs"])
        w.writerows((g,c,n) for (g,c),n in sorted(counts.items()))
    with (out/"figure_cluster_members.tsv").open("w",newline="") as f:
        fields=["figure","genome","contig","orf","gene","category","start","end","strand"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(selected_rows)
    manifest={'source':str(source.resolve()),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'gff_sources':[{'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in gffs],
        'sample_inventory_source':{'path':str(heatmap_source.resolve()),'sha256':hashlib.sha256(heatmap_source.read_bytes()).hexdigest()} if heatmap_source.is_file() and not genome else None,
        'orf_count':len(rows),'total_samples':len(genomes),'samples_with_hits':sum(n>0 for n in genomes.values()),'category_count':len(categories),
        'selected_genomes':top_g,'selected_categories':top_c,'clusters_available':total_clusters,
        'cluster_selection':cluster_manifest,'warnings':sorted(set(warnings)),
        'category_handling':'Spelling aliases only; no biological reassignment of historical results.',
        'interpretation':'Raw annotation counts; genomic potential, not expression, activity, abundance, or validation accuracy.',
        'outputs':exports,'dpi':dpi,'formats':['png','svg','pdf']}
    (out/"figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest


def main(argv=None):
    parser=argparse.ArgumentParser(description="Private figure export from an existing Efesto run")
    parser.add_argument("--results-dir",required=True,type=Path)
    parser.add_argument("--out-dir",required=True,type=Path)
    parser.add_argument("--gff-dir",type=Path)
    parser.add_argument("--genome")
    parser.add_argument("--title",default="Efesto")
    parser.add_argument("--max-clusters",type=int,default=6)
    parser.add_argument("--top-categories",type=int,default=24)
    parser.add_argument("--top-genomes",type=int,default=30)
    parser.add_argument("--dpi",type=int,default=300)
    args=parser.parse_args(argv)
    if min(args.top_categories,args.top_genomes,args.dpi)<=0 or args.max_clusters<0:
        parser.error("Counts and DPI must be positive; max-clusters may be zero")
    try:
        result=export_figures(args.results_dir,args.out_dir,gff_dir=args.gff_dir,genome=args.genome,
            max_clusters=args.max_clusters,top_categories=args.top_categories,
            top_genomes=args.top_genomes,dpi=args.dpi,title=args.title)
    except (ValueError,RuntimeError,OSError) as exc:
        parser.exit(2,"Figure export failed: %s\n" % exc)
    print("Exported %d image files to %s" % (len(result['outputs']),args.out_dir))
    for warning in result['warnings']:print("[WARN] "+warning)


if __name__ == "__main__":
    main()
