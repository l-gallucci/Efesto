# Memory use and HMM search concurrency

Efesto runs at most **two HMM searches simultaneously by default**. The existing
CPU budget still applies: `--threads 1` or a larger `--hmm_threads` value can
reduce concurrency further. No new command-line option is required.

A typical command remains:

```bash
efesto --fna_dir genomes --out efesto_out --threads 4 --hmm_threads 1
```

The limit reduces overlap between HMMER working matrices. Difficult
profile/target combinations can use hundreds of MiB per search, so reducing
Python interpreter overhead alone does not solve peak memory use. Fewer
simultaneous searches can increase wall time. This is a concurrency limit,
**not a hard RAM ceiling**: individual searches may still need substantial
memory, and higher numbers of HMMER workers per search can raise it.

Python threads now coordinate external HMMER processes, with at most two
submitted jobs per active search worker. Protein sequences load after HMM
searches finish. These changes avoid extra coordinator interpreters, bound
queue size and avoid overlapping the sequence dictionary with the HMMER peak.
HMMER arguments, score/E-value cutoffs and annotation rules remain unchanged.

The HMM search limit does not change Pyrodigal/Pyrodigal-GV worker allocation.
`--threads` continues to control that stage. Lower it as well if gene prediction
is the memory bottleneck on a larger dataset.

The runtime uses the standard-library CSV/TSV readers and dictionaries, not
pandas. Replacing these with Polars would add a dependency without addressing
the measured HMMER working-matrix peak. Reconsider columnar processing if
future profiling identifies large result-table operations as the bottleneck.

Use a fresh output directory when comparing resources: existing HMM cache
files can otherwise bypass searches. Check gene-level output equivalence as
well as memory and runtime.
