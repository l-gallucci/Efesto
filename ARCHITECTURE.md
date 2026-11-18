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
| Evidence reporting | `evidence.py` | Per-hit observations, availability and retained-neighborhood rule assessments |
| Filter decisions | `decisions.py`, `operon.py` | Actual decisions for selected candidates, disk spooling and final rescue status |
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

HMMER, gene-call and UniOP caches verify recorded input/configuration and completed
output hashes before reuse. Use a fresh output directory for comparisons so stale
final reports cannot be mistaken for new output. Record the source revision and
external-tool versions alongside the run log.

## Separate evidence

`evidence.py` records search gates, available annotation/structure/topology
observations and retained-neighborhood rule assessments without a combined score.
UniOP probabilities remain pair-level observations; qualifying BGC regions are
reported per overlapping ORF. Missing data and bypassed rules are explicit.
Rule satisfaction does not imply completeness or physiological activity.

`decisions.py` spools selected HMM candidates to disk with decisions emitted by
the executing context and second-pass filters. The final TSV preserves excluded
candidates and their initial rule failures even after eggNOG rescue. It does not
alter retained annotations or duplicate the full candidate dataset in RAM.

See [the user pipeline](docs/pipeline.md) and [separate evidence](docs/evidence.md).
