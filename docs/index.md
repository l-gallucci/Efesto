# Efesto Documentation

Efesto is a profile HMM annotation pipeline for iron-cycling and metal resistance
genes in genomes and metagenomes. It extends [FeGenie](https://github.com/Arkadiy-Garber/FeGenie)
with coordinate-aware operon clustering, a cluster confidence scoring system, an expanded
curated HMM library, and integration with external tools (UniOP, antiSMASH, Anvi'o).

---

## Contents

### Getting started
- [[Installation|installation]]
- [Usage — quick start](https://github.com/l-gallucci/Efesto#usage)

### HMM library
- [[HMM library overview|hmm_library]] — categories, sources, model counts, registry format
- [[Iron-sulfur assembly models|iron_sulfur_assembly]] — SUF / ISC biology and HMM details
- [[Library curation log|hmm_library_curation]] — deduplication decisions, deprecation rationale
- [[Biological rationale for new models|hmm_expansion_biological_rationale]] — all new categories

### Pipeline logic
- [[Operon rules|operon_rules]] — rule engine, JSON schema, TonB/ExbBD energizer-guard
- [[Cluster confidence scoring|cluster_confidence_scoring]] — formula and all four components
- [[AntiSMASH BGC integration|antismash_bgc]] — siderophore BGC boost

### Outputs and integration
- [[Output formats|output_formats]] — all files, all columns
- [[Anvi'o integration|anvio_integration]] — functions import, gene-scores import
- [[Visualisation|visualisation]] — R heatmap script

### Development
- [[Extending the HMM library|extending_the_library]] — adding models, curation workflow
- [[Publication plan|publication_plan]] — manuscript roadmap

---

## Quick orientation

```
Input: FAA/FNA + GFF  →  hmmsearch  →  genomic clustering  →  operon rules
  →  confidence scoring  →  GFF3 + TSV + summary stats + optional Anvi'o TSVs
```

Key design decisions:
- **One hit per ORF** (highest bitscore wins across all HMMs).
- **Operon clustering** uses bp coordinates (GFF mode) or Prodigal ordinal index (fallback).
- **Confidence score** is multiplicative: `hmm_weight × co_occ_score × uniop_pair_score × bgc_boost`.
- **TonB-ExbBD** are energy-transducer proteins, not substrate-specific; a cluster of only
  Ton-motor proteins is dropped by the `SIDERO_TRANSPORT` rule (see [[Operon rules|operon_rules]]).
- **Zero-cutoff models** (MetHMMDB + 14 FeGenie siderophore models) use `--zero_cutoff_min_bitscore`
  as a bitscore floor instead of a calibrated threshold.

---

## Citing Efesto

If you use Efesto, cite all sources that apply to your run. See [README — Citations](https://github.com/l-gallucci/Efesto#citations).
