# Workflows and examples

## Genomes and MAGs

Provide one nucleotide FASTA per genome or bin. Extensions must match the input;
the default is `.fna`, so specify `--fna_ext fasta` for `.fasta` files.

```bash
efesto --fna_dir genomes/ --fna_ext fasta --meta \
  --out efesto_out --threads 8 --hmm_threads 1 --strand_aware
```

Use `--meta` for fragmented assemblies and MAGs. For complete genomes, the default
single-genome mode may be appropriate. `--gene_caller pyrodigal` selects the
alternative bacterial/archaeal caller; the default is `pyrodigal-gv`.

Gene-calling outputs are stored under `_prodigal/` for compatibility; the directory
name does not mean that Efesto executed the original Prodigal program.

## Existing proteins with coordinates

```bash
efesto --faa_dir proteins/ --gff_dir gff/ --out efesto_out --threads 4
```

Match file stems, contig IDs and protein/CDS IDs. For example, `sample.faa` and
`sample.gff` should describe the same gene calls. A GFF from a different caller
may have incompatible identifiers or coordinates.

Without GFF coordinates, Efesto falls back to suffix/index-based clustering.
That fallback assumes Prodigal-style identifiers encoding contig and gene order;
it is inappropriate for arbitrary Bakta locus tags or shuffled catalogs.

## Gene catalogs

```bash
efesto --faa_dir catalog/ --catalog_mode --out catalog_out --threads 4
```

Each ORF becomes a singleton, because catalog numbering does not establish genomic
adjacency. Co-occurrence filters are bypassed, while score gates and applicable
per-protein checks remain. Interpret these as individual-family annotations.

## Restrict the annotation scope

```bash
efesto --faa_dir proteins/ --gff_dir gff/ --annotate Fe Cu Zn --out metal_out
```

Other selectors include `Ln`, `W`, `metallophores` and canonical category names.
`Ln` selects LanM/LanD handling profiles. `W` includes Tup and shared Mo/W
families. See [library scope](hmm_library.md).

## Reuse Bakta labels and external annotations

```bash
efesto --fna_dir genomes/ --fna_ext fasta --meta \
  --bakta_gff_dir bakta_gff/ \
  --bgc_dir antismash_gff/ \
  --eggnog_annotations annotations.emapper.annotations \
  --out efesto_out --threads 8 --hmm_threads 1 \
  --max_bp_gap 5000 --strand_aware --keep_tblout
```

Efesto still calls proteins for its search, then maps compatible calls to Bakta
labels by coordinates. Matching allows alternative starts; unmatched ORFs retain
their original identifiers. This does not import all Bakta annotations or replace
the HMM results with Bakta calls. See [anvi'o identifiers](anvio_integration.md).

## Reproducible comparisons

Use a **new output directory** after changing the input, library, thresholds or
software version. Existing cached searches or gene calls can otherwise be reused.
Record the command, source revision, library and tool versions together with
`Efesto-run.log`. `--keep_tblout` retains raw HMM search tables for inspection.
