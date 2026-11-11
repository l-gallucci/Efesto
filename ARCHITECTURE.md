# Architecture

Efesto is a Python command-line pipeline. Core modules use standard-library
parsers and data structures; external tools and optional Python packages are
required only for their corresponding execution paths.

## Execution stages

| Stage | Main modules | Responsibility |
|---|---|---|
| CLI and orchestration | `cli.py` | Inputs, options, stage ordering, reporting |
| Gene prediction | `gene_calling.py` | Pyrodigal/Pyrodigal-GV, Prodigal-compatible FAA/GFF |
| Library/input loading | `io.py`, `categories.py` | Registry, category aliases/selectors, FASTA/GFF and metadata |
| Profile search | `hmmer.py` | Bounded HMMER job queue, search tables, best accepted hit per ORF |
| Neighborhoods | `clustering.py` | Coordinate/strand grouping, ordinal fallback, catalog singletons |
| Context filtering | `operon.py` | Default and JSON rules, Mtr/Mto context and applicable checks |
| Evidence weighting | `scoring.py` | Heuristic confidence factors |
| Optional context | `uniop.py`, `bgc.py` | Operon probabilities and qualifying BGC overlap |
| Optional confirmation | `eggnog.py`, `baktfold.py`, `topology.py` | Applicable annotation, structure and topology evidence |
| Coverage | `coverage.py` | Coverage inputs and optional normalization |
| Export | `writers.py` | Summary, tidy TSV, GFF3, sequence and anvi'o exports |

## Library contract

The HMM filename stem is the external lookup key for registry, cutoff and gene
maps. The internal HMM name is separate. Deprecated/inactive registry entries
are excluded before search; active/experimental entries remain eligible.
One accepted HMM assignment is retained per ORF by highest raw bitscore.

Custom JSON context rules are loaded from the selected library directory.
The repository-root example does not automatically override the default engine.
Canonical category aliases preserve spelling compatibility without certifying
biological equivalence.

## Concurrency and caching

HMMER runs as external processes, coordinated with a bounded Python thread pool.
At most two searches overlap by default, subject to the supplied CPU budget.
Gene prediction uses a separate process-pool allocation. Sequence loading is
ordered to avoid overlap with the HMMER memory peak where possible.

Raw searches and gene calls can be reused from existing output directories.
The cache is not a complete content-hash provenance system: changing inputs or
parameters requires a fresh output directory for a reproducible comparison.

## Confidence and evidence

The composite score combines HMM threshold status, co-occurrence, available UniOP
pairs, qualifying BGC overlap and applicable optional confirmation. Missing optional
evidence is neutral. The score is capped at 1.2 and is not a correctness probability.
There is no global confidence≥0.5 retention gate in ordinary runs.

See [the user pipeline](docs/pipeline.md) and
[confidence scoring](docs/cluster_confidence_scoring.md).
