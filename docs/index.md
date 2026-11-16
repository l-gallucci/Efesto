# Efesto

![Efesto logo](img/logo.png){ .project-logo }

**Profile HMM annotation of metal cycling, acquisition, homeostasis and resistance genes.**

Efesto searches proteins from genomes, MAGs and gene catalogs, then combines
family matches with genomic context and optional annotation evidence. It extends
[FeGenie](https://github.com/Arkadiy-Garber/FeGenie) with an expanded library,
coordinate-based clustering and integrations with UniOP, antiSMASH and anvi'o.

## Start here

1. [Install Efesto](installation.md).
2. [Choose a workflow](workflows.md) for assemblies, existing proteins or a gene catalog.
3. Read the [pipeline overview](pipeline.md) and [output reference](output_formats.md).
4. Check [interpretation and limitations](interpretation.md) before making biological claims.

## HMM library

{{catalog_summary}}

Registered `active` and `experimental` profiles are enabled in ordinary runs;
deprecated profiles are skipped. This runtime status is distinct from independent
biological validation. See the [current catalog](hmm_library.md).

## Quick start

```bash
efesto --fna_dir genomes/ --fna_ext fasta --meta --out efesto_out --threads 4
```

For proteins that already have gene coordinates:

```bash
efesto --faa_dir proteins/ --gff_dir gff/ --out efesto_out --threads 4
```

The default nucleotide gene caller is Pyrodigal-GV. The HMM library is located
relative to the installed source checkout; use `--hmm_dir` for a custom library.
See [installation](installation.md) and [workflows](workflows.md).

## Interpretation

A family hit supports genetic potential. It does not demonstrate expression,
metal specificity for every homolog, a complete pathway or a resistance phenotype.
Cluster confidence is a heuristic evidence score, not a probability of correctness.

[References and citations](references.md) ·
[Source code](https://github.com/l-gallucci/Efesto) ·
[MIT license](https://github.com/l-gallucci/Efesto/blob/main/LICENSE)

Logo: [Elena Panariello](https://www.instagram.com/elena_panariello/).
