# HMM Library Curation — Decision Log

This document records the rationale for every model excluded from the active
Efesto HMM library. All decisions are reflected in `hmm_library/hmm_registry.tsv`
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
| `IsdG-heme-oxygenase_fam2-rep` | iron_acquisition-heme_oxygenase | `IsdG-heme-oxygenase_fam1-rep` | fam1 has higher nseq (15 vs 8); near-identical cutoffs (102 vs 101.8) |
| `Sid_MbtC_Mycobactin_biosynthesis_…` | iron_acquisition-siderophore_synthesis | `K04790` (tabuteau) | tabuteau model calibrated (626.5), nseq=128 vs fegenie nseq=25, cutoff=0 |
| `Sid_AsbD_Petrobactin_biosynthesis_…` | iron_acquisition-siderophore_synthesis | `K24111` (tabuteau) | tabuteau model calibrated (97.4), nseq=39 vs fegenie nseq=25 |
| `FpuC-FhuC-FpuD-YusV-CbrD-FepC-PvuE-FecE-HatC-ATPase` | iron_acquisition-siderophore_transport | `FpuCD-YusV-CbrD-FepC-FepE-PvuE-ABC_transporters_ATPase` | near-identical ABC ATPase models (nseq 51 vs 50); retained the one with slightly higher nseq |
| `FecE-PvuE-FhuC-CbrD-FepC-YusV-FpuD-FpuC-PiuA-PirA-FepA-RhtA` | iron_acquisition-siderophore_transport | `PvuE-FecB-vibrioferrin_transport` | kept higher-cutoff model (350.0 vs 200.1) |
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
| `Fe_binding_yfeA_1` (methmmdb) | iron_resistance | FeGenie `Iron_uptake_YfeA_…` [iron_acquisition-iron_transport] | YfeA is the periplasmic substrate-binding protein of the Yfe chelated-iron ABC transporter (Yersinia pestis). It is an iron acquisition gene, not a resistance determinant. |
| `Fe_transport_1` (methmmdb) | iron_resistance | FeGenie `Iron_uptake_YfeB_…` [iron_acquisition-iron_transport] | YfeB is the membrane permease subunit of the same Yfe system. Same rationale as above. |
| `Cation_symport_actP_1` | metal_resistance-non-specific | — | ActP is a cation/acetate symporter involved in acetate uptake, not metal homeostasis. No homology to any metal resistance or transport function. |
| `Multidrug_efflux_yfmO_2` | metal_resistance-non-specific | — | YfmO is a subunit of the MdtABC multidrug efflux complex. Primary function is antibiotic resistance; no evidence for metal-specific efflux. |
| `Multidrug_resistance_mdtA_1` | metal_resistance-non-specific | — | MdtA: membrane fusion protein of MdtABC multidrug efflux. Not metal-specific. |
| `Multidrug_resistance_mdtB_1` | metal_resistance-non-specific | — | MdtB: RND permease subunit of MdtABC. Not metal-specific. |
| `Multidrug_resistance_mdtC_1` | metal_resistance-non-specific | — | MdtC: RND permease subunit of MdtABC. Not metal-specific. |

### Category reclassifications (2026-05-23)

Models whose source category was incorrect but which are retained as active under the correct category.

| Model | Old category | New category | Reason |
|---|---|---|---|
| `Fe_periplasmic_1`, `Fe_permease_efeU_1`, `Fe_transport_2/3`, `Fe_transport_fbpB/C_1`, `Fe_transport_fecD/E_1` | `iron_resistance` | `iron_acquisition-iron_transport` | These are iron uptake transporters (EfeU, FbpBC, FecDE, YfeCD); erroneously placed in resistance by FeGenie. [doi:10.3389/fcimb.2013.00090](https://doi.org/10.3389/fcimb.2013.00090) |
| `Fe_ferripyoverdine_receptor_1`, `Fe_pyochelin_1` | `iron_resistance` | `iron_acquisition-siderophore_transport` | TonB-dependent outer membrane receptors for siderophore-iron complexes; acquisition, not resistance. |
| `Transferrin_TbpB_binding_protein_Haemophilus_influenzae_P44971` | `iron_storage` | `iron_acquisition-iron_transport` | TbpB captures host transferrin for iron extraction — acquisition from host, not storage. |
| `FMN_reductase_arsH_2`, `Quinone_reductase_arsH_1` | `metal_resistance-non-specific` | `metal_resistance-arsenic` | ArsH is the auxiliary NADPH-oxidoreductase of the *ars* operon; function is arsenic detoxification. |
| `Metal_binding_zinT_1` | `metal_resistance-non-specific` | `metal_resistance-cobalt_zinc_cadmium` | ZinT is a periplasmic zinc-binding chaperone that assists ZnuABC zinc uptake. |
| `Cation_efflux_fieF_1` | `metal_resistance-non-specific` | `metal_resistance-cobalt_zinc_cadmium` | FieF (YiiP) is a CDF-family Fe²⁺/Zn²⁺ efflux transporter; more specific to Zn/Fe homeostasis than "non-specific". |
| `flavodoxin_long`, `flavodoxin_short` | `iron_gene_regulation` | `iron_stress` | Flavodoxins are electron carriers expressed under iron limitation as ferredoxin substitutes — biomarkers of iron stress, not transcriptional regulators. [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0) |
| `Zur` | `iron_gene_regulation` | `metal_resistance-cobalt_zinc_cadmium` | Zur is a zinc-sensing Fur-family repressor (Zn²⁺, not Fe²⁺). Placed with ZnuABC and other zinc homeostasis models. [doi:10.3389/fcimb.2013.00059](https://doi.org/10.3389/fcimb.2013.00059) |

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

3. **FutA1 / FutA2 (iron_acquisition-iron_transport).** These cyanobacterial
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

## MtrA / MtoA Disambiguation — Calibrated (2026-05-23)

**Status:** Curated seeds deployed; TC/GA/NC calibrated against 3794-sequence
universe. GA cutoffs deployed. Inherent Shewanella MtrD twilight zone documented.

### Problem (FeGenie-derived models, nseq=9/10)

FeGenie MtrA and MtoA HMMs cross-hit each other at ~2× the calibrated cutoff:

```
MtrA HMM (FeGenie, nseq=9)  vs MtoA consensus: 265.8 bits  (cutoff=140) — 1.9× cutoff
MtoA HMM (FeGenie, nseq=10) vs MtrA consensus: 289.8 bits  (cutoff=140) — 2.1× cutoff
```

NCBIfam/TIGRFAM survey: no model discriminates the two subfamilies.
TIGR03508 ("DmsE family decaheme c-type cytochrome") covers BOTH families
across 2531 RefSeq hits — too broad. No dedicated MtrA or MtoA model exists.

### Resolution — HMM rebuild from curated seeds

Script: `scripts/build_mtr_mto_subfamily_hmms.py`
Seed file: `data/seeds/mtr_mto_seeds.tsv`

**Seeds used (final set — 15 MtrA, 6 MtoA):**

| Accession | Subfamily | Organism | Len | Role in seed set |
|---|---|---|---|---|
| WP_011071900.1 | MtrA | *Shewanella oneidensis* MR-1 | 333 | Canonical MtrA |
| AAN54830.1     | MtrA | *Shewanella oneidensis* MR-1 (SO_1777) | 333 | Original locus |
| CAQ9398997.1   | MtrA | *Shewanella xiamenensis* (MAG) | 333 | Shewanella diversity |
| WP_482097394.1 | MtrA | *Shewanella* sp. | 333 | Shewanella diversity |
| WP_432404861.1 | MtrA | *Shewanella baltica* | 333 | Shewanella diversity |
| QYX66249.1     | MtrA | *Shewanella putrefaciens* | 333 | Shewanella diversity |
| WP_011623263.1 | MtrA | *Shewanella* sp. (unclassified) | 329 | Shewanella diversity |
| BGZ64319.1     | MtrA | *Shewanella algae* | 329 | Shewanella diversity |
| BGZ78722.1     | MtrA | *Shewanella marina* | 325 | Shewanella diversity |
| XAG79646.1     | MtrA | *bacterium* 19NY03SH02 (Shewanellales) | 329 | Unclassified Shewanellales |
| WP_233132600.1 | MtrA | *Paraferrimonas haliotis* | ~330 | Non-Shewanella iron reducer |
| WP_095506167.1 | MtrA | *Paraferrimonas sedimenticola* | ~330 | Non-Shewanella iron reducer |
| WP_095505338.1 | MtrA | *Paraferrimonas sedimenticola* (paralog) | ~330 | Non-Shewanella iron reducer |
| WP_222466948.1 | MtrA | *Ferrimonas balearica* | ~330 | Non-Shewanella iron reducer |
| WP_013344866.1 | MtrA | *Ferrimonas balearica* DSM 9799 | ~330 | Non-Shewanella iron reducer |
| ADE12722.1     | MtoA | *Sideroxydans lithotrophicus* ES-1 (Slit_2497) | 355 | Canonical MtoA |
| WP_013030620.1 | MtoA | *Sideroxydans lithotrophicus* | 355 | Sideroxydans diversity |
| WP_283743965.1 | MtoA | *Sideroxydans* sp. CL21 | 359 | Sideroxydans diversity |
| WP_283743089.1 | MtoA | *Sideroxydans* sp. CL21 | 354 | Sideroxydans diversity |
| ADL56010.1     | MtoA | *Gallionella capsiferriformans* ES-2 | 343 | Genus diversity |
| WP_013293942.1 | MtoA | *Gallionella capsiferriformans* | 343 | Genus diversity |

**HMM build statistics:**

| Model | NSEQ | EFFN | LENG |
|---|---|---|---|
| MtrA | 15 | ~0.43 | 333 |
| MtoA | 6  | 0.630 | 344 |

### Full calibration against 3794-sequence universe

The calibration universe = all NCBI proteins titled "DmsE family decaheme c-type
cytochrome" (NCBI query: `decaheme[Title] AND "DmsE"[Title]`, 2026-05-23).
This represents the TIGR03508 superfamily — the broadest boundary within which
MtrA and MtoA must be distinguished.

**Key calibration finding: Shewanella MtrD is the critical NC**

NCBI does not distinguish MtrA from MtrD (or DmsE) in protein titles — all are
called "DmsE family decaheme c-type cytochrome". The true discriminating boundary
was found by scoring confirmed MtrD sequences (DMSO/alternate-respiration periplasmic
decaheme paralog in *Shewanella*; locus SO_1780 family) directly:

| Protein | Organism | Score vs MtrA HMM |
|---|---|---|
| MtrD — GIU44267.1 | *Shewanella* sp. | 551.7 |
| MtrD — GIU17583.1 | *Shewanella* sp. | 549.4 |
| MtrD — GLD76960.1 | *Shewanella* sp. | 547.2 |
| MtrD — GIU07609.1 | *Shewanella* sp. | 546.9 |
| MtrD — AAN54835.2 | *S. oneidensis* MR-1 | 539.9 |
| MtrD — GCF87789.1 | *Shewanella* sp. | 529.9 |

These are the highest-scoring confirmed NON-MtrA sequences in the Shewanellales.

**Score landscape for MtrA HMM across the full universe:**

```
MtoA (Gallionellaceae iron oxidizers)      :  223 – 262 bits  [correctly excluded]
MtrD (Shewanella DMSO/alternate respiration):  530 – 552 bits  ← NC
────────── twilight zone: 552 – 603 bits ─────────────────────── GA = 580
Ferrimonas / Paraferrimonas (iron reducers) :  603 – 637 bits  ← TC
Shewanella MtrA                            :  603 – 706 bits  ← TC (same range)
Unclassified Shewanellales (MAG)           :  690 – 724 bits  
```

**Score landscape for MtoA HMM:**

```
Shewanella MtrA / MtrD (iron reducers)     :  255 – 303 bits  [correctly excluded]
betaproteobacteria DmsE-family (unknown)   :  460 – 505 bits  ← NC (uncertain)
────────── twilight zone: 505 – 536 bits ─────────────────────── GA = 520
Gallionella capsiferriformans ES-2         :  612 bits         ← TC
Sideroxydans sp. CL21 (most distant)       :  536 bits         ← TC (lowest)
Sideroxydans lithotrophicus ES-1           :  646 bits         
```

### TC / GA / NC values

**MtrA HMM:**

| Value | Bits | Basis |
|---|---|---|
| NC (Noise Cutoff) | 552 | Max MtrD score (*Shewanella* SO_1780 family) |
| **GA (Gathering, deployed)** | **580** | Midpoint NC–TC; 28 bits above NC |
| TC (Trusted Cutoff) | 603 | Min confirmed MtrA (*Ferrimonas balearica*) |
| Twilight zone | 552–603 | Gene-context required for sequences in this range |

**MtoA HMM:**

| Value | Bits | Basis |
|---|---|---|
| NC (Noise Cutoff) | 505 | Max score for betaproteobacteria DmsE-family (*Rhodoferax*, *Polynucleobacter*) |
| **GA (Gathering, deployed)** | **520** | Conservative; 15 bits above NC |
| TC (Trusted Cutoff) | 536 | Min confirmed MtoA (*Sideroxydans* sp. CL21) |
| Twilight zone | 505–536 | Narrow; gene-context required |

### Metagenome annotation guidance

The following table describes expected annotation behavior for real metagenomes:

| Organism type | MtrA score | MtoA score | Expected output | Notes |
|---|---|---|---|---|
| *Shewanella* MtrA | 600–706 | 255–292 | iron_reduction | Correctly assigned |
| *Shewanella* MtrD (DMSO) | 530–552 | ~260 | **not annotated** | Below GA (580); gene context confirms |
| *Ferrimonas* / *Paraferrimonas* | 603–637 | 263–273 | iron_reduction | Correctly assigned |
| Unclassified Shewanellales MAG | 690–724 | ~292 | iron_reduction | Correctly assigned |
| *Geobacter* (uses OmcB/Z, not MtrA) | <50 | <50 | **not annotated** | Different protein family entirely |
| *Gallionella* / *Sideroxydans* MtoA | 223–262 | 536–646 | iron_oxidation | Correctly assigned |
| betaproteobacteria DmsE-family | ~200–400 | 460–505 | **not annotated** | Below both GAs |
| Sequences scoring 552–603 (MtrA zone) | — | — | **low_confidence** | Twilight zone; flag for review |
| Sequences scoring 505–536 (MtoA zone) | — | — | **low_confidence** | Twilight zone; narrow margin |

**MtrA/MtoA cannot be assigned to Geobacter**, which uses outer-membrane
multiheme cytochromes OmcB, OmcS, OmcZ (different protein family, captured
by separate HMMs `OmcS`, `OmcZ` already in the library).

### Known remaining limitations

1. **Shewanella MtrD twilight zone (552–603 bits)** is inherent to the protein
   family — MtrA and MtrD share the decaheme c-type fold. Gene context (presence
   of MtrB + MtrC in the same cluster) is the definitive discriminator for sequences
   in this score range. Efesto's UniOP/operon module provides this context.

2. **MtoA NC margin is narrow (30.8 bits)**: The betaproteobacteria DmsE-family
   sequences at 460–505 bits (NC = 505.5, *Limnohabitans sp.*) are **confirmed
   non-iron-oxidizers** (validated 2026-05-23): *Limnohabitans* = freshwater
   heterotroph; *Rhodoferax* = iron **reducer** (scores 486–497, below GA);
   *Polynucleobacter* = obligate aerobe, tiny genome; *Sulfurimicrobium* = sulfur
   oxidizer; *Methylococcales/Methylophilaceae* = methanotrophs. No action needed.
   NC = 505 confirmed. GA = 520 is appropriately placed (14.5 bits above NC,
   16.3 bits below TC).

3. **Geobacter and Desulfuromonas** use completely different MHC proteins for EET
   (OmcB, OmcS, OmcZ, DmsE-family with different topology). These are already
   in the library as separate HMMs and do not interact with MtrA/MtoA annotation.

4. **MtoA database sparsity**: Only 6 unique MtoA sequences in all of NCBI (2026-05-23).
   The model will improve significantly as Gallionellaceae MAGs accumulate.

### Calibration procedure used

```bash
# Universe: all "DmsE family decaheme c-type cytochrome" sequences (3794 total)
python scripts/calibrate_mtr_mto_cutoffs.py
# Output: hmm_library/_calibration/mtr_mto/calibration_report.tsv
# Output: hmm_library/_calibration/mtr_mto/calibration_summary.txt

# NC validation (MtrD scoring):
# Searched MtrD sequences manually from NCBI:
# '"MtrD" AND Shewanella[Organism] AND decaheme[All Fields]'
# → GIU44267.1, GLD76960.1, AAN54835.2, GCF87789.1 etc.
# → scored 530–552 bits → sets NC at 552
```

---

## New categories added (2026-05-23)

| Category | Models | Rationale | Key references |
|---|---|---|---|
| `iron_stress` | `flavodoxin_long` (TIGR01752.1), `flavodoxin_short` (TIGR01753.1) | Iron-starvation response biomarkers. Flavodoxins substitute for ferredoxin when iron is limiting. Detecting them in metagenomes indicates Fe-stress capacity. NOT transcriptional regulators; separated from `iron_gene_regulation`. | [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0) |
| `iron_sulfur_assembly` | `sufC` (TIGR01978.1), `sufB` (TIGR01980.1), `sufS` (TIGR01979.1), `IscS` (TIGR02006.1) | Fe-S cluster biosynthesis machinery. SUF system (oxidative-stress / iron-limitation pathway); ISC system (constitutive). Allows distinguishing organisms relying on stress-tolerant SUF vs housekeeping ISC assembly. | SUF: [doi:10.1074/jbc.M308004200](https://doi.org/10.1074/jbc.M308004200); IscS: [doi:10.1074/jbc.M401261200](https://doi.org/10.1074/jbc.M401261200) |

---

## Exceptions to the Layer A minimum-coverage rule (2026-07-16)

**Status:** Kept `active` despite violating the ≤3-seq ("others" tier)
deprecation threshold documented above.

| Model | nseq | Category | Reason kept active |
|---|---|---|---|
| `mofA_MnOxGeneTool` | 3 | `manganese_oxidation` | At the threshold, not below it. MofA is a genuinely rare, narrowly-distributed Mn(II)-oxidizing multicopper oxidase clade — MnOxGeneTool's own reference database only contains 3 representative sequences for this clade. No broader alternative source exists. |
| `mopA_A_MnOxGeneTool` | 1 | `manganese_oxidation` | Single-sequence model (Alphaproteobacteria MopA clade). Same reason: MnOxGeneTool's curated reference set has only 1 sequence for this specific clade split. |
| `mopA_E_MnOxGeneTool` | 1 | `manganese_oxidation` | Single-sequence model (Epsilonproteobacteria MopA clade). Same reason. |

**Rationale for the exception:** unlike the FeGenie/methmmdb low-nseq models
deprecated above (which had *more diverse* alternative models available or were
simply under-sampled duplicates), these three MopA/MofA clade splits reflect
genuine sparsity in the manganese-oxidation literature and reference database —
MnOxGeneTool (Wang et al. 2025, [doi:10.1021/acs.est.5c01235](https://doi.org/10.1021/acs.est.5c01235))
deliberately keeps clade-specific subfamily models (matching the precedent of
FeGenie's own same-name subfamily groups, see Layer A note above) rather than
collapsing them into one over-broad MopA model. There is currently no
alternative, better-sampled source for these three specific clades.

**Known limitation:** single/near-single-sequence profiles generalise poorly
to divergent members of the same clade not represented in the training set —
expect these three models to under-call (false negatives), not over-call, on
novel lineages. Revisit if MnOxGeneTool publishes an updated reference set
with more representatives for these clades.

---

## SUF operon deployment — `sufA`/`sufD`/`sufE` flagged for confirmation (2026-07-19)

**Status:** Deployed active, but marked `needs_confirmation` in the registry
(new column) — a distinct, weaker trust tier from ordinary `active` models.

`sufA` (TIGR01997.1), `sufD` (TIGR01981.1), and `sufE` (Pfam PF02657.22) were
added to complete the `sufABCDSE` operon alongside the pre-existing `sufB`/
`sufC`/`sufS`/`IscS`. Unlike those four, all three have a near-zero GA–NC gap
(2.4, 7, and ~0.2 bits respectively) — confirmed empirically, not just from the
HMM's own embedded statistics: tested `sufE` against the real *E. coli* K-12
proteome and it cross-hits `csdE` (a related but functionally distinct,
non-SUF-operon paralog — per TIGRFAM's own annotation of `csdE`, "not found
next to other such genes as are its paralogs from the SUF... systems") at
138.2 bits, essentially the same range as the true `sufE` hit (154.5 bits).

**Why deployed anyway instead of left out:** `sufA`/`sufD` are true `sufB`/
`sufC`/`sufS`-scale accessory proteins (deployed to complete the annotatable
operon) and the weak standalone signal is compensated by the new
`SUF_OPERON` operon co-occurrence rule (`src/efesto/operon.py`,
`operon_rules.json`, docs in `operon_rules.md`) — a lone `needs_confirmation`
hit unsupported by ≥ 3 total SUF genes (in practice, the calibrated anchors)
is dropped from the reported cluster.

**Why no custom HMM rebuild was attempted:** investigated first. Confirmed
seed scarcity is a genuine data-availability limit, not a fixable curation
mistake — Swiss-Prot review coverage for `sufA`/`sufD` is essentially just
*E. coli* + one uncertain *B. subtilis* entry ("Uncharacterized protein
SufA"), with no real cross-taxon diversity to build a tighter model from.
`sufE` has more reviewed sequences but they're almost all Enterobacterales,
which would only shift the same problem to a different sequence space.
Synteny-bootstrapping candidate seeds (finding neighbors of confirmed
`sufB`/`sufC`/`sufS` and assuming they're `sufD`/`sufA`/`sufE`) was
considered and explicitly rejected — deliberate policy decision — because
operon gene order isn't universally conserved across distant taxa and using
genomic position to build training data for a rule that itself checks
genomic position risks circularity.

**Escalation pipeline (2026-07-20) — tier 2 implemented, tier 3 pending:** a
general low-confidence-hit escalation, not SUF-specific — any hit landing on
a `needs_confirmation`-flagged model (currently just `sufA`/`sufD`/`sufE`,
but the mechanism is general-purpose) that gets dropped by tier-1 operon
filtering (`SUF_OPERON`, etc.) is collected into a small pending list and
checked in a second pass, scoped only to that flagged subset — never the
whole proteome:

- **Tier 2 (eggNOG-mapper) — implemented.** `src/efesto/eggnog.py` parses a
  `*.emapper.annotations` file (`--eggnog_annotations`, read-only — Efesto
  never modifies it) or optionally runs `emapper.py` itself on just the
  flagged subset (`--run_eggnog`, requires `--eggnog_db_dir`; never
  auto-invoked — the reference database is tens of GB, users manage that
  download themselves). Comparison is by KEGG Ortholog (KO) number where
  both sides have one, falling back to a loose `Preferred_name` match
  otherwise. A confirmed hit is re-injected into the final output as a
  rescued single-gene cluster; a contradicted or uninformative one stays
  dropped. New `annot_weight` factor in `scoring.py`
  (`cluster_confidence = ... × annot_weight × struct_weight`), same
  neutral-on-no-data policy as `uniop_weight`. `--export_flagged_faa` writes
  the same flagged-subset FASTA independent of eggNOG, for running any other
  tool of choice on exactly that candidate set. Verified end-to-end: a
  synthetic isolated `sufA` hit (no operon neighbors) is correctly dropped
  with no eggNOG data, correctly rescued when a matching eggNOG confirmation
  is supplied, and correctly stays dropped when eggNOG contradicts it.
- **Tier 3 (Baktfold, ProstT5→Foldseek structural search) — designed, not
  yet implemented.** `struct_weight` exists as a placeholder parameter in
  `cluster_confidence` (always 1.0) for forward compatibility.

KO numbers verified via the KEGG REST API for all six SUF genes: `sufA`
K05997, `sufB` K09014, `sufC` K09013, `sufD` K09015, `sufS` K11717, `sufE`
K02426 — added as a new `kegg_ko` registry column rather than a separate
mapping file, to avoid fragmenting model metadata across files.

---

## Summary of all deprecations

| Status | Count | Reason |
|---|---|---|
| `deprecated_nseq_insufficient` | 60 | Too few training sequences |
| `deprecated_dedup_lower_coverage` | 87 | Name-based duplicate with lower alignment coverage |
| `deprecated_dedup_sequence_cluster` | 12 | Layer B: redundant at 70% id / 80% cov |
| `deprecated_category_mismatch` | 7 | Acquisition genes in resistance; non-metal genes (ActP, MdtABC, YfmO) |
| **Total deprecated** | **166** | |
| `experimental` | 10 | Lower-trust tier, not yet fully validated |
| **Active** | **478** | |
| **Library total** | **654** | |
