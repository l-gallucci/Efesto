# HMM Library Curation — Decision Log

This document records the rationale for every model excluded from the active
MetalGenie-Evo HMM library. All decisions are reflected in `hmm_library/hmm_registry.tsv`
via the `status` column.

---

## Layer A — Minimum-Coverage Filter

**Status tag:** `deprecated_nseq_insufficient`  
**Criterion:** Models built from fewer sequences than the source-specific threshold
(fegenie ≤ 5 seqs, methmmdb ≤ 10 seqs, others ≤ 3 seqs).

**Rationale:** HMM profiles built from very few sequences do not generalise reliably
beyond the training sequences. They produce high false-positive rates when applied
to diverse metagenomic datasets because the profile has not seen enough sequence
variation to define the family boundaries accurately.

**Count:** 60 models deprecated.

---

## Layer A — Name-Based Deduplication

**Status tags:** `deprecated_exact_duplicate`, `deprecated_dedup_lower_coverage`  
**Criterion:** Models sharing an identical stem name within a source, where one model
has lower sequence coverage of the alignment or fewer training sequences.

**Rationale:** Several sources (primarily methmmdb) contain organism-specific variants
of the same gene (e.g. `Cu_efflux_copA_1` through `copA_12`). Where the underlying
protein family is the same, retaining multiple narrow models inflates hit counts and
makes output interpretation harder. The representative with the highest nseq and
alignment coverage was retained; organism-specific duplicates were deprecated.

**Note on FeGenie same-name groups:** FeGenie contains groups of models that share
a gene name but have distinct calibrated bitscore cutoffs (e.g. iron reduction models).
These were inspected manually and retained as intentional subfamilies — FeGenie used
separate models to distinguish functional variants within a gene family.

**Count:** 88 models deprecated (1 exact duplicate + 87 lower-coverage duplicates).

---

## Layer B — Sequence-Level Deduplication

**Tool:** MMseqs2 `easy-cluster`, applied to hmmemit consensus sequences.  
**Parameters:** 70% sequence identity, 80% bidirectional coverage, sensitivity 7.5.  
**Script:** `scripts/layer_b_dedup_all.py`

### B1 — Same-category redundant pairs

**Status tag:** `deprecated_dedup_sequence_cluster`

Models within the same functional category whose consensus sequences cluster at ≥ 70%
identity with ≥ 80% bidirectional coverage. The representative was selected by:
1. Calibrated cutoff (> 0) preferred over uncalibrated (= 0)
2. Higher `nseq` (more training sequences → better-defined family boundaries)
3. Source tier: interpro > fegenie/tabuteau > methmmdb

| Deprecated model | Category | Kept representative | Reason |
|---|---|---|---|
| `IsdG-heme-oxygenase_fam2-rep` | iron_aquisition-heme_oxygenase | `IsdG-heme-oxygenase_fam1-rep` | fam1 has higher nseq (15 vs 8); near-identical cutoffs (102 vs 101.8) |
| `Sid_MbtC_Mycobactin_biosynthesis_…` | iron_aquisition-siderophore_synthesis | `K04790` (tabuteau) | tabuteau model calibrated (626.5), nseq=128 vs fegenie nseq=25, cutoff=0 |
| `Sid_AsbD_Petrobactin_biosynthesis_…` | iron_aquisition-siderophore_synthesis | `K24111` (tabuteau) | tabuteau model calibrated (97.4), nseq=39 vs fegenie nseq=25 |
| `FpuC-FhuC-FpuD-YusV-CbrD-FepC-PvuE-FecE-HatC-ATPase` | iron_aquisition-siderophore_transport | `FpuCD-YusV-CbrD-FepC-FepE-PvuE-ABC_transporters_ATPase` | near-identical ABC ATPase models (nseq 51 vs 50); retained the one with slightly higher nseq |
| `FecE-PvuE-FhuC-CbrD-FepC-YusV-FpuD-FpuC-PiuA-PirA-FepA-RhtA` | iron_aquisition-siderophore_transport | `PvuE-FecB-vibrioferrin_transport` | kept higher-cutoff model (350.0 vs 200.1) |
| `CdCoZn_efflux_czcA_5` | metal_resistance-multimetal | `CdCoZn_efflux_czcA_1` | both methmmdb czcA variants; czcA_1 has more than double the nseq (50 vs 23) |
| `Ag_efflux_silA_1` | metal_resistance-silver | `CuAg_cusA_1` (metal_resistance-multimetal) | CusA and SilA are paralogous RND efflux pumps (Cu/Ag vs Ag); the broader CuAg model with higher nseq (75 vs 25) covers both families |
| `Ag_lipoprotein_silC_1` | metal_resistance-silver | `CuAg_cusC_1` (metal_resistance-multimetal) | CusC and SilC are paralogous OMF components of the same RND system; CuAg model retained (nseq 101 vs 24) |

### B2 — Category-mismatch deprecations

**Status tag:** `deprecated_category_mismatch`

Models that cluster with well-characterised iron **uptake** genes but are catalogued
under `iron_resistance` in their source database. Iron homeostasis genes are sometimes
annotated as "resistance" in broad databases because excess or misregulated iron can
be toxic, but from a functional genomics perspective they are uptake/transport genes.
Retaining them in the resistance category would inflate false-positive resistance calls.

| Deprecated model | Source category | Clusters with | Evidence |
|---|---|---|---|
| `Fe_binding_yfeA_1` (methmmdb) | iron_resistance | FeGenie `Iron_uptake_YfeA_…` [iron_aquisition-iron_transport] | YfeA is the periplasmic substrate-binding protein of the Yfe chelated-iron ABC transporter (Yersinia pestis). It is an iron acquisition gene, not a resistance determinant. |
| `Fe_transport_1` (methmmdb) | iron_resistance | FeGenie `Iron_uptake_YfeB_…` [iron_aquisition-iron_transport] | YfeB is the membrane permease subunit of the same Yfe system. Same rationale as above. |

---

## Models kept despite cross-category clustering

The Layer B clustering (70% id / 80% cov) also identified pairs spanning unrelated
functional categories (e.g. siderophore synthesis vs siderophore transport, or iron
regulation vs transport). These were **not deprecated** for the following reasons:

1. **Low-nseq hmmemit artifact.** Several FeGenie models were built from ≤ 25
   sequences. The `hmmemit -c` consensus for such models is degenerate — it collapses
   to the most frequent residue at each alignment column and does not represent a real
   protein sequence. Consensus sequences from different structural families can
   spuriously cluster at 70% identity when both are degenerate.

2. **Distinct structural families.** Pairs flagged include NRPS adenylation domains
   (siderophore synthesis) clustering with TonB-dependent outer-membrane receptors
   (siderophore transport) — these are unrelated folds. Deprecating one would silently
   remove a whole protein class.

3. **FutA1 / FutA2 (iron_aquisition-iron_transport).** These cyanobacterial
   iron-binding ABC transporter substrate-binding proteins are paralogous (both
   annotated as FutA), share high sequence identity, and cluster at 70%. However
   FutA1 and FutA2 have been shown to have distinct expression patterns and iron
   affinities in *Synechocystis* sp. PCC 6803. Both models are retained to allow
   separate annotation.

Pairs retained:

| Pair | Reason kept |
|---|---|
| `K11604` [synthesis] + FeGenie `yfeA` [transport] | Different folds; yfeA consensus degenerate (nseq=25) |
| `K02363` [transport] + `DhbE-PchD-…` [synthesis] | TonB receptor vs NRPS; DhbE nseq=9 (degenerate) |
| `K04780` [transport] + `DhbF-PvdDIJL-VabF` [synthesis] | Transport permease vs NRPS; DhbF nseq=10 |
| `K23181` [synthesis] + `FecB-PvuB-vibrioferrin` [transport] | Synthesis domain vs periplasmic binding protein |
| `put_FhuA_V3_corr` [transport] + `PRK10044.1` [synthesis] | Outer membrane beta-barrel vs synthesis enzyme |
| `Sid_FpvI_regulator` [regulation] + `Sid_FpvH_iron_release` [transport] | ECF sigma factor vs periplasmic iron-release protein |
| `FutA1` + `FutA2` [both iron_transport] | Functionally distinct cyanobacterial paralogs |

---

## Summary of all deprecations

| Status | Count | Reason |
|---|---|---|
| `deprecated_nseq_insufficient` | 60 | Too few training sequences |
| `deprecated_dedup_lower_coverage` | 87 | Name-based duplicate with lower alignment coverage |
| `deprecated_exact_duplicate` | 1 | Byte-identical to another active model |
| `deprecated_dedup_sequence_cluster` | 12 | Layer B: redundant at 70% id / 80% cov |
| `deprecated_category_mismatch` | 2 | Iron uptake genes misclassified as resistance |
| **Total deprecated** | **162** | |
| **Active** | **456** | |
| **Library total** | **618** | |
