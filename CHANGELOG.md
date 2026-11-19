# Changelog

## Unreleased

- Export registry model status per hit and selected candidate, including experimental
  profiles and an explicit unknown status for libraries without that metadata.
- Preserve provenance, exclusions, confirmation expectations and source notes when
  reconstructing libraries from existing models or reviewed reports.

- Separate siderophore synthesis, export and utilization labels from the stable
  neighborhood groups used for candidate retention; include the roles in Fe selectors.
- Distinguish nutritional ferric-chelate reduction from respiratory iron reduction.
- Report the shared FhuF/IucA/IucC domain PF06276 as function unresolved; mark
  PF06276 and the unresolved PHYLLO assignment experimental, retaining both in searches.
- Clarify IroB/K23725 functional overlap, substrate-ambiguous transport families,
  vibriobactin VibH, VabS and shared NRPS accessories; preserve HMMs and cutoffs.
- Add a complete siderophore/metallophore scope table and reference scope regressions.

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
- Stop on HMMER and gene-prediction failures; reuse only complete, input-matched caches.
- Bind gene-call caches to assembly/caller/mode and restrict collection to current inputs.
- Bind UniOP caches to input content, script/interpreter and the exact saved FAA index;
  reconstruct Bakta headers only when all proteins have matching coordinates.
- Evaluate all default context rules in mixed neighborhoods without early-return bypasses;
  apply category minima independently and honor declared trusted-lone exceptions.
- Export every selected HMM candidate with its actual contig/context/second-pass
  decisions and final retention/rescue status in `Efesto-filter-decisions.tsv`;
  spool records to disk without changing thresholds or retained annotations.
- Exclude local HMM staging experiments from source archives as well as wheels.
- Preserve registry cutoff precision in the runtime table and reject inconsistent packages.
- Resolve UniOP/eggNOG script names through PATH and quote count/coverage heatmap CSVs.
- Parse samtools depth position tables correctly and stream per-base coverage.
- Read UniOP pair probabilities alongside operon membership, including zero.
- Replace the composite confidence score and weights with separate hit, context,
  candidate and UniOP-pair evidence; distinguish missing checks and bypassed rules.
- Replace double length-normalisation with coverage-per-million scaling;
  coverage normalisation does not alter match or context evidence.
- Include the model library and provenance metadata in wheel/source packages.
- Add software citation metadata and separate third-party source conditions.
- License all Efesto-owned code, standalone tooling, tests and documentation
  under AGPLv3, and record source-specific HMM terms and unresolved permissions.
- Reconstruct all eight CbcL seed records and alignment from public sequences;
  reproduce the existing model exactly without changing its matrix or cutoff.
- Clarify that the optional R script uses separately installed pheatmap;
  the package implementation is not redistributed by Efesto.
- Include upstream license texts, per-model source notices and the curated
  seed alignments in the appropriate release archives.
