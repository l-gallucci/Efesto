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

The default nucleotide gene caller is Pyrodigal-GV. Release wheels include the
HMM library; editable installs use the source checkout. Use `--hmm_dir` for a
custom library.
See [installation](installation.md) and [workflows](workflows.md).

## Read the results

| Question | Start with |
|---|---|
| Which protein family matched? | `Efesto-results-long.tsv`: HMM, bitscore, E-value and threshold type |
| What was observed nearby? | `Efesto-context-evidence.tsv`: distinct models, copies, distances and rule observations |
| Was a rejected candidate recovered? | `Efesto-confirmation-candidates.tsv`: eggNOG comparisons and retention outcome |
| What did UniOP predict? | Operon membership and `Efesto-UniOP-pair-evidence.tsv` |
| How many retained clusters or hits are present? | Cluster-count heatmap or per-ORF results, respectively |

The [pipeline](pipeline.md) explains which checks decide retention. The
[evidence guide](evidence.md) shows how to read individual observations, with
MtrA/LanM and SUF examples. Tool results and missing data remain separate.

A family hit supports genetic potential. It does not demonstrate expression,
metal specificity for every homolog, a complete pathway or a resistance phenotype.
See [interpretation](interpretation.md) for the biological limits.

[References and citations](references.md) ·
[Source code](https://github.com/l-gallucci/Efesto) ·
[AGPLv3 code license](https://github.com/l-gallucci/Efesto/blob/main/LICENSE) ·
[Third-party terms](third-party-notices.md)

Logo: [Elena Panariello](https://www.instagram.com/elena_panariello/).
