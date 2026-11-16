# Changelog

## 1.0.0

First tagged-release preparation. Earlier development checkouts used internal
1.0.x version strings; no earlier release tags were published.

- Integrate the curated metal-function HMM catalog, canonical categories,
  profile provenance and model-specific interpretation guidance.
- Support genome and metagenome inputs, gene-neighbourhood filtering and
  optional UniOP, eggNOG, topology, structure and BGC evidence.
- Distinguish repeated protein IDs across genomes in external-tool queries and
  combined exports; provide an identifier crosswalk.
- Resolve Bakta ORFs to GFF contigs before short-contig filtering.
- Stop on HMMER failures; reuse only complete, input-matched search caches.
- Parse samtools depth position tables correctly and stream per-base coverage.
- Read UniOP pair probabilities alongside operon membership, including zero.
- Write all confidence-bearing tables after final evidence calculation.
- Replace double length-normalisation with coverage-per-million scaling;
  coverage normalisation does not alter annotation confidence.
- Include the model library and provenance metadata in wheel/source packages.
- Add software citation metadata and separate third-party source conditions.
