# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Gene prediction via pyrodigal / pyrodigal-gv."""

# pyrodigal-gv (Camargo et al. 2023) extends pyrodigal with viral/giant-virus
# metagenomic models. For bacterial sequences in meta mode it selects the
# bacterial model → results identical to pyrodigal. For viral sequences it
# uses virus-specific models. pyrodigal produces results identical to Prodigal
# v2.6.3 (verified by Julian Hahnfeld). Prodigal binary is no longer needed.
#
# Mode recommendation (Tang Li et al. 2022, Briefings Bioinformatics):
#   meta   → MAGs, metagenomes, gene catalogs (fragmented assemblies)
#   single → complete/near-complete isolated genomes

import sys
import tempfile
from importlib.metadata import version, PackageNotFoundError
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from efesto._cache import cache_matches, file_digest, write_manifest


def _pyrodigal_job(args_tuple):
    """Picklable gene prediction job for ProcessPoolExecutor."""
    fna_path, faa_out, gff_out, meta_mode, gene_caller = args_tuple
    stem = Path(fna_path).stem

    faa_out, gff_out = Path(faa_out), Path(gff_out)
    manifest = faa_out.with_suffix(".cache.json")
    try:
        if gene_caller == "pyrodigal-gv":
            import pyrodigal_gv as library
            finder_class = library.ViralGeneFinder
        else:
            import pyrodigal as library
            finder_class = library.GeneFinder
        try:
            caller_version = version(gene_caller)
        except PackageNotFoundError:
            caller_version = getattr(library, "__version__", "unavailable")
        inputs = {"schema": 1, "fna_sha256": file_digest(fna_path),
                  "meta": bool(meta_mode), "caller": gene_caller,
                  "caller_version": caller_version}
        outputs = {"faa": faa_out, "gff": gff_out}
        if cache_matches(manifest, inputs, outputs):
            return stem, True, ""
        manifest.unlink(missing_ok=True)
    except (ImportError, OSError) as error:
        return stem, False, f"Could not prepare {gene_caller}: {error}"

    sequences = {}
    header, seq_parts = None, []
    try:
        with open(fna_path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header:
                        sequences[header] = "".join(seq_parts)
                    header = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
        if header:
            sequences[header] = "".join(seq_parts)
    except Exception as e:
        return stem, False, f"Could not read FASTA: {e}"

    if not sequences:
        return stem, False, "No sequences found"

    try:
        finder = finder_class(meta=True) if meta_mode else finder_class()
        if not meta_mode:
            finder.train(*sequences.values())
    except Exception as error:
        return stem, False, f"Could not initialize/train {gene_caller}: {error}"

    try:
        with tempfile.TemporaryDirectory(dir=faa_out.parent, prefix=stem + ".") as staging:
            temp_faa, temp_gff = Path(staging) / "genes.faa", Path(staging) / "genes.gff"
            with temp_faa.open("w") as faa_fh, temp_gff.open("w") as gff_fh:
                gff_fh.write("##gff-version 3\n")
                gene_idx = 0
                for seq_id, seq in sequences.items():
                    for gene in finder.find_genes(seq):
                        gene_idx  += 1
                        orf_id     = f"{seq_id}_{gene_idx}"
                        start      = gene.begin
                        end        = gene.end
                        strand     = "+" if gene.strand == 1 else "-"
                        start_type = gene.start_type or "Edge"
                        p_start    = start_type == "Edge"
                        p_end      = getattr(gene, "partial_end", False)
                        partial_flag = f"{'1' if p_start else '0'}{'1' if p_end else '0'}"
                        gc = sum(c in "GCgc" for c in seq[start - 1:end]) / max(end - start + 1, 1)
                        aa_seq = gene.translate()

                        faa_fh.write(
                            f">{orf_id} # {start} # {end} # {gene.strand} # "
                            f"ID={gene_idx};partial={partial_flag};"
                            f"start_type={start_type};"
                            f"rbs_motif={gene.rbs_motif or 'None'};"
                            f"gc_cont={gc:.3f}\n"
                        )
                        for i in range(0, len(aa_seq), 60):
                            faa_fh.write(aa_seq[i:i + 60] + "\n")

                        gff_fh.write(
                            f"{seq_id}\t{gene_caller}\tCDS\t{start}\t{end}\t"
                            f".\t{strand}\t0\t"
                            f"ID={orf_id};partial={partial_flag};"
                            f"start_type={start_type}\n"
                        )
            if file_digest(fna_path) != inputs["fna_sha256"]:
                raise RuntimeError("Assembly changed during gene prediction; rerun with stable inputs")
            temp_faa.replace(faa_out)
            temp_gff.replace(gff_out)
            write_manifest(manifest, inputs, outputs)
        return stem, True, ""
    except Exception as e:
        return stem, False, str(e)


def run_prodigal(fna_files, out_dir, meta_mode=False, threads=1,
                 gene_caller="pyrodigal-gv"):
    prodigal_dir = out_dir / "_prodigal"
    faa_dir = prodigal_dir / "faa"
    gff_dir = prodigal_dir / "gff"
    faa_dir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(exist_ok=True)

    if len({Path(fna).stem for fna in fna_files}) != len(fna_files):
        raise RuntimeError("Assembly filenames must have unique stems")
    failures = []
    jobs = [
        (str(fna), str(faa_dir / f"{fna.stem}.faa"),
         str(gff_dir / f"{fna.stem}.gff"), meta_mode, gene_caller)
        for fna in fna_files
    ]

    print(f"[INFO] Running {gene_caller} on {len(fna_files)} assemblies "
          f"({'meta' if meta_mode else 'single'} mode)…")

    n_workers = min(threads, len(fna_files))
    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_pyrodigal_job, j): j for j in jobs}
            done = errors = 0
            for fut in as_completed(futures):
                done += 1
                try:
                    stem, ok, err = fut.result()
                except Exception as error:
                    stem, ok, err = Path(futures[fut][0]).stem, False, str(error)
                if not ok:
                    errors += 1
                    failures.append(f"{stem}: {err}")
                    print(f"\n  [WARN] Gene prediction failed {stem}: {err[:80]}",
                          file=sys.stderr)
                sys.stdout.write(f"\r  {done}/{len(fna_files)}  ({errors} errors)  ")
                sys.stdout.flush()
        print()
    else:
        for i, job in enumerate(jobs, 1):
            stem, ok, err = _pyrodigal_job(job)
            sys.stdout.write(f"\r  {i}/{len(fna_files)}  ")
            sys.stdout.flush()
            if not ok:
                failures.append(f"{stem}: {err}")
                print(f"\n  [WARN] Gene prediction failed {stem}: {err[:80]}",
                      file=sys.stderr)
        print()

    if failures:
        raise RuntimeError("Gene prediction failed; annotation stopped: " + "; ".join(failures))
    return faa_dir, gff_dir
