"""I/O helpers: FASTA, GFF, map/cutoff files."""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    seqs, header, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if header is not None:
        seqs[header] = "".join(parts)
    return seqs


def read_fasta_lengths(path):
    lengths, header, length = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = length
                header = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if header is not None:
        lengths[header] = length
    return lengths


def read_cutoffs(path):
    c = {}
    if not os.path.isfile(path):
        return c
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                try:
                    c[ls[0]] = float(ls[1])
                except ValueError:
                    pass
    return c


def read_map(path):
    m = {}
    if not os.path.isfile(path):
        return m
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                m[ls[0]] = ls[1]
    return m


def get_contig_lengths_from_gff(gff_path):
    lengths = defaultdict(int)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = re.search(r"##sequence-region\s+(\S+)\s+\d+\s+(\d+)", line)
                if m:
                    lengths[m.group(1)] = int(m.group(2))
                continue
            parts = line.rstrip().split("\t")
            if len(parts) >= 5 and parts[2] == "CDS":
                try:
                    end = int(parts[4])
                    if end > lengths[parts[0]]:
                        lengths[parts[0]] = end
                except ValueError:
                    pass
    return dict(lengths)


def build_contig_length_dict(faa_files, gff_dir=None, fna_dir=None, fna_ext="fna"):
    contig_lengths = {}
    for faa in faa_files:
        stem = faa.stem
        genome = faa.name
        found = None
        if fna_dir:
            for ext in [fna_ext, "fna", "fasta", "fa"]:
                c = Path(fna_dir) / f"{stem}.{ext}"
                if c.exists():
                    found = ("fna", str(c))
                    break
        if not found and gff_dir:
            for ext in [".gff", ".gff3", ".prodigal.gff"]:
                c = Path(gff_dir) / (stem + ext)
                if c.exists():
                    found = ("gff", str(c))
                    break
        if found:
            kind, path = found
            contig_lengths[genome] = (
                read_fasta_lengths(path) if kind == "fna"
                else get_contig_lengths_from_gff(path)
            )
    return contig_lengths


def load_prodigal_gff(gff_path):
    coords = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", parts[8])
            if m:
                coords[m.group(1).strip()] = {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
    return coords


def load_gff_dir(gff_dir, faa_files):
    gff_dir = Path(gff_dir)
    gc = {}
    for faa in faa_files:
        stem = faa.stem
        found = None
        for ext in (".gff", ".gff3", ".prodigal.gff"):
            c = gff_dir / (stem + ext)
            if c.exists():
                found = c
                break
        if found:
            gc[faa.name] = load_prodigal_gff(str(found))
        else:
            print(f"  [WARN] No GFF for {faa.name}, using index clustering",
                  file=sys.stderr)
    return gc
