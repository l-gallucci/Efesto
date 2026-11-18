# Performance and memory

Efesto limits concurrent HMM searches and uses a bounded job queue to reduce
peak memory use. This page explains the default scheduling behavior and how
thread settings affect resource use.

For measured wall time and memory, see the [benchmarks](comparison.md#performance-comparisons):
matched-library and full-library comparisons with FeGenie, and a comparison of
Efesto search-scheduling configurations.

## HMM search concurrency

Efesto runs at most **two HMM searches simultaneously by default**. The
CPU budget also applies: `--threads 1` or a larger `--hmm_threads` value can
reduce concurrency further.

For example:

```bash
efesto --fna_dir genomes --out efesto_out --threads 4 --hmm_threads 1
```

The limit reduces overlap between HMMER working matrices. Individual
profile/target combinations can use hundreds of MiB per search. Fewer
simultaneous searches can increase wall time. This is a concurrency limit,
**not a hard RAM ceiling**: individual searches may still need substantial
memory, and higher numbers of HMMER workers per search can raise it.

External HMMER processes use a bounded job queue. Protein sequences load after
HMM searches finish, reducing overlap between sequence storage and search memory.
Concurrency settings affect scheduling rather than annotation thresholds or rules.

Filter-decision reporting spools selected candidates to a temporary file in the
output directory and writes the final TSV after retention and rescue decisions.
It avoids keeping a second full candidate table in RAM, at the cost of temporary
disk space and additional I/O.

The HMM search limit does not change Pyrodigal/Pyrodigal-GV worker allocation.
`--threads` continues to control that stage. Lower it as well if gene prediction
is the memory bottleneck on a larger dataset.

Use a fresh output directory when comparing resources: existing HMM cache
files can otherwise bypass searches. Check gene-level output equivalence as
well as memory and runtime.
