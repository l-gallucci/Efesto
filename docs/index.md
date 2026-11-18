# Efesto

![Efesto logo](img/logo.png){ .project-logo }

**Profile HMM annotation of metal cycling, acquisition, homeostasis and resistance genes.**

Efesto integrates profile HMMs from multiple databases and published collections
to annotate metal-related genes in genomes, metagenomes and gene catalogs. It
combines protein-family searches with genomic-context filtering and reports
optional annotation, operon, BGC, structure and topology evidence separately.

## Start here

1. [Install Efesto](installation.md).
2. [Choose a workflow](workflows.md) for assemblies, existing proteins or a gene catalog.
3. Read the [pipeline overview](pipeline.md) and [output reference](output_formats.md).
4. Check [interpretation and limitations](interpretation.md) before making biological claims.

## HMM library

Sources include **FeGenie, MetHMMDB, MnOxGeneTool, Tabuteau et al., Pfam and NCBI
protein family models**, together with profiles built for Efesto. Coverage spans
iron, manganese, copper, zinc, nickel, cobalt, tungsten, lanthanides and other
metal/metalloid functions. The
[catalog and source records](hmm_library.md#sources) distinguish import channels
from original providers and retain model-specific provenance and cutoffs.
FeGenie also provides the basis of the default context engine, with
Efesto-specific changes; see [pipeline](pipeline.md) and [references](references.md).

{{catalog_summary}}

Registered `active` and `experimental` profiles are enabled in ordinary runs;
deprecated profiles are skipped. See the [current catalog](hmm_library.md) for
model-specific scope and limitations.

## Scope and model changes compared with FeGenie

The [source and workflow comparison](comparison.md) describes model changes relative to FeGenie,
contributions from other sources, Efesto-built profiles, categories, thresholds
and performance comparisons.

A matched-library benchmark (**3–4 November 2026, before v1.0.0**) measured
**38.0 s / 562.8 MiB** sampled peak RSS for Efesto and **182.7 s / 886.4 MiB**
for FeGenie using the same proteins and 191 original FeGenie HMMs. The paired
HMMER tables contained identical hit rows. Each program was measured once,
including FeGenie's built-in BLAST searches and waits. See [benchmark methods
and results](comparison.md#performance-comparisons).

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
| Why was a candidate retained, excluded or rescued? | `Efesto-filter-decisions.tsv`: actual filter actions and final outcome |
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
