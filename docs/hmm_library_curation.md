# HMM Library Curation — Decision Log

This document records the rationale for every model excluded from the active
Efesto HMM library. The current registry is authoritative; historical narrative below is not a complete validation of every retained assignment.

## Category cleanup — 5 November 2026

All category directories and runtime references now use the canonical vocabulary in `category_aliases.tsv`. Explicit per-model decisions and evidence are in `category_curation.tsv`; current counts and review limitations are described in [HMM library](hmm_library.md).

The previous audit identified 47 Tabuteau import/biosynthesis conflicts, seven live InterPro conflicts, ten live MetHMMDB iron-acquisition conflicts and the FeGenie TbpB storage assignment. Their categories have been corrected; ComGD and three glyceraldehyde dehydrogenase models are now deprecated out of scope. Additional source-supported uptake, homeostasis and stress-response assignments have been separated from resistance/iron uptake. Unresolved mechanisms remain explicitly unclassified.

The Fe uptake models were imported through **MetHMMDB**, not FeGenie. Earlier text saying these moves had already been applied was inconsistent with the actual library; the migration now synchronizes directories and registry. TbpB's storage category was inherited from FeGenie. Historical deprecation thresholds and retained-model counts below describe earlier curation steps, not a fresh measurement of generalization or classifier accuracy. Low NSEQ alone does not establish an empirical false-positive rate.

Current status: 652 registered, 473 active, 12 experimental, 167 deprecated. The three CutABC profiles were physically removed after the sequence/literature audit; the mislabeled duplicate FpvI was quarantined. HMM bytes, cutoffs and existing experimental statuses were preserved. `source_supported` review labels are source evidence, not independent validation.


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
| `IsdG-heme-oxygenase_fam2-rep` | Heme-Oxygenase | `IsdG-heme-oxygenase_fam1-rep` | fam1 has higher nseq (15 vs 8); near-identical cutoffs (102 vs 101.8) |
| `Sid_MbtC_Mycobactin_biosynthesis_…` | Siderophore-Synthesis | `K04790` (tabuteau) | tabuteau model calibrated (626.5), nseq=128 vs fegenie nseq=25, cutoff=0 |
| `Sid_AsbD_Petrobactin_biosynthesis_…` | Siderophore-Synthesis | `K24111` (tabuteau) | tabuteau model calibrated (97.4), nseq=39 vs fegenie nseq=25 |
| `FpuC-FhuC-FpuD-YusV-CbrD-FepC-PvuE-FecE-HatC-ATPase` | Siderophore-Transport | `FpuCD-YusV-CbrD-FepC-FepE-PvuE-ABC_transporters_ATPase` | near-identical ABC ATPase models (nseq 51 vs 50); retained the one with slightly higher nseq |
| `FecE-PvuE-FhuC-CbrD-FepC-YusV-FpuD-FpuC-PiuA-PirA-FepA-RhtA` | Siderophore-Transport | `PvuE-FecB-vibrioferrin_transport` | kept higher-cutoff model (350.0 vs 200.1) |
| `CdCoZn_efflux_czcA_5` | MetalResistance-MultiMetal | `CdCoZn_efflux_czcA_1` | both methmmdb czcA variants; czcA_1 has more than double the nseq (50 vs 23) |
| `Ag_efflux_silA_1` | MetalResistance-Silver | `CuAg_cusA_1` (MetalResistance-MultiMetal) | CusA and SilA are paralogous RND efflux pumps (Cu/Ag vs Ag); the broader CuAg model with higher nseq (75 vs 25) covers both families |
| `Ag_lipoprotein_silC_1` | MetalResistance-Silver | `CuAg_cusC_1` (MetalResistance-MultiMetal) | CusC and SilC are paralogous OMF components of the same RND system; CuAg model retained (nseq 101 vs 24) |

### B2 — Category-mismatch deprecations

**Status tag:** `deprecated_category_mismatch`

Models that cluster with well-characterised iron **uptake** genes but are catalogued
under `Iron-Resistance` in their source database. Iron homeostasis genes are sometimes
annotated as "resistance" in broad databases because excess or misregulated iron can
be toxic, but from a functional genomics perspective they are uptake/transport genes.
Retaining them in the resistance category would inflate false-positive resistance calls.

| Deprecated model | Source category | Clusters with | Evidence |
|---|---|---|---|
| `Fe_binding_yfeA_1` (methmmdb) | Iron-Resistance | FeGenie `Iron_uptake_YfeA_…` [IronAcquisition-Transport] | YfeA is the periplasmic substrate-binding protein of the Yfe chelated-iron ABC transporter (Yersinia pestis). It is an iron acquisition gene, not a resistance determinant. |
| `Fe_transport_1` (methmmdb) | Iron-Resistance | FeGenie `Iron_uptake_YfeB_…` [IronAcquisition-Transport] | YfeB is the membrane permease subunit of the same Yfe system. Same rationale as above. |
| `Cation_symport_actP_1` | MetalResistance-NonSpecific | — | ActP is a cation/acetate symporter involved in acetate uptake, not metal homeostasis. No homology to any metal resistance or transport function. |
| `Multidrug_efflux_yfmO_2` | MetalResistance-NonSpecific | — | YfmO is a subunit of the MdtABC multidrug efflux complex. Primary function is antibiotic resistance; no evidence for metal-specific efflux. |
| `Multidrug_resistance_mdtA_1` | MetalResistance-NonSpecific | — | MdtA: membrane fusion protein of MdtABC multidrug efflux. Not metal-specific. |
| `Multidrug_resistance_mdtB_1` | MetalResistance-NonSpecific | — | MdtB: RND permease subunit of MdtABC. Not metal-specific. |
| `Multidrug_resistance_mdtC_1` | MetalResistance-NonSpecific | — | MdtC: RND permease subunit of MdtABC. Not metal-specific. |

### Category reclassifications (2026-05-23)

Models whose source category was incorrect but which are retained as active under the correct category.

| Model | Old category | New category | Reason |
|---|---|---|---|
| `Fe_periplasmic_1`, `Fe_permease_efeU_1`, `Fe_transport_2/3`, `Fe_transport_fbpB/C_1`, `Fe_transport_fecD/E_1` | `Iron-Resistance` | `IronAcquisition-Transport` | These are iron uptake transporters (EfeU, FbpBC, FecDE, YfeCD); previously assigned to resistance during the MetHMMDB import. [doi:10.3389/fcimb.2013.00090](https://doi.org/10.3389/fcimb.2013.00090) |
| `Fe_ferripyoverdine_receptor_1`, `Fe_pyochelin_1` | `Iron-Resistance` | `Siderophore-Transport` | TonB-dependent outer membrane receptors for siderophore-iron complexes; acquisition, not resistance. |
| `Transferrin_TbpB_binding_protein_Haemophilus_influenzae_P44971` | `Iron-Storage` | `IronAcquisition-Transport` | TbpB captures host transferrin for iron extraction — acquisition from host, not storage. |
| `FMN_reductase_arsH_2`, `Quinone_reductase_arsH_1` | `MetalResistance-NonSpecific` | `MetalResistance-Arsenic` | ArsH is the auxiliary NADPH-oxidoreductase of the *ars* operon; function is arsenic detoxification. |
| `Metal_binding_zinT_1` | `MetalResistance-NonSpecific` | `MetalResistance-CobaltZincCadmium` | ZinT is a periplasmic zinc-binding chaperone that assists ZnuABC zinc uptake. |
| `Cation_efflux_fieF_1` | `MetalResistance-NonSpecific` | `MetalResistance-CobaltZincCadmium` | FieF (YiiP) is a CDF-family Fe²⁺/Zn²⁺ efflux transporter; more specific to Zn/Fe homeostasis than "non-specific". |
| `flavodoxin_long`, `flavodoxin_short` | `Iron-Regulation` | `Iron-StressResponse` | Flavodoxins are electron carriers expressed under iron limitation as ferredoxin substitutes — biomarkers of iron stress, not transcriptional regulators. [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0) |
| `Zur` | `Iron-Regulation` | `MetalResistance-CobaltZincCadmium` | Zur is a zinc-sensing Fur-family repressor (Zn²⁺, not Fe²⁺). Placed with ZnuABC and other zinc homeostasis models. [doi:10.3389/fcimb.2013.00059](https://doi.org/10.3389/fcimb.2013.00059) |

---

## Reassessment of apparent cross-category redundancy (2026-11-05)

The former explanation of these pairs as unrelated folds with “degenerate consensuses” is withdrawn. Several category labels in that table were wrong: K02363/DhbE are related ligases, K04780/DhbF are NRPS profiles, K11604/YfeA are related binding proteins, K23181/FecB are uptake proteins, and the two FhuA profiles are receptors. Their related consensus sequences do not establish interchangeable biological scope.

A fresh audit of all 655 pre-removal profiles, including experimental models, identified eight enabled pairs at >=70% consensus identity and >=80% coverage of both sequences. The other seven pairs remain enabled pending independent evaluation; none is removed solely for a consensus match. On ten local proteomes (37,633 proteins), several pairs retain additional hits under their own thresholds. Broad ATPase-family co-hits are a specificity/threshold question, not proof of removable duplicates.

The exception is `Sid_FpvI_regulator_PA2387_Paeruginosa_PAO1_180620`: its emissions, transitions, STATS and checksum are identical to the FpvH profile. Both score 425.0 bits against official PA2404/NP_251094.1; neither detects the true FpvI PA2387/NP_251077.1 above cutoff 50. The duplicate exists in the pinned upstream FeGenie snapshot. FpvH is retained; the mislabeled FpvI copy is archived with `deprecated_misassigned_duplicate`. No replacement for the true FpvI is claimed. [FpvI primary study](https://doi.org/10.1128/jb.185.4.1261-1265.2003).

FutA1/FutA2 are biological paralogs, but neither high consensus similarity nor one shared hit establishes that their models discriminate them. Independent cyanobacterial references are required. NSEQ and positive cutoff values alone are not a ranking of model quality.

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
| AAN54830.1 | MtrA | *Shewanella oneidensis* MR-1 (SO_1777) | 333 | Identical-sequence duplicate; not additional diversity |
| CAQ9398997.1 | MtrA | *Shewanella xiamenensis* (MAG) | 333 | Shewanella diversity |
| WP_482097394.1 | MtrA | *Shewanella* sp. | 333 | Shewanella diversity |
| WP_432404861.1 | MtrA | *Shewanella baltica* | 333 | Shewanella diversity |
| QYX66249.1 | MtrA | *Shewanella putrefaciens* | 333 | Shewanella diversity |
| WP_011623263.1 | MtrA | *Shewanella* sp. (unclassified) | 329 | Shewanella diversity |
| BGZ64319.1 | MtrA | *Shewanella algae* | 329 | Shewanella diversity |
| BGZ78722.1 | MtrA | *Shewanella marina* | 325 | Shewanella diversity |
| XAG79646.1 | MtrA | *bacterium* 19NY03SH02 (Shewanellales) | 329 | DmsE-family homolog; specific MtrA function inferred |
| WP_233132600.1 | MtrA | *Paraferrimonas haliotis* | 312 | DmsE-family homolog; specific MtrA function inferred |
| WP_095506167.1 | MtrA | *Paraferrimonas sedimenticola* | 318 | DmsE-family homolog; specific MtrA function inferred |
| WP_095505338.1 | MtrA | *Paraferrimonas sedimenticola* (paralog) | 315 | DmsE-family homolog; specific MtrA function inferred |
| WP_222466948.1 | MtrA | *Ferrimonas balearica* | 316 | DmsE-family homolog; specific MtrA function inferred |
| WP_013344866.1 | MtrA | *Ferrimonas balearica* DSM 9799 | 316 | DmsE-family homolog; specific MtrA function inferred |
| ADE12722.1 | MtoA | *Sideroxydans lithotrophicus* ES-1 (Slit_2497) | 355 | Canonical MtoA |
| WP_013030620.1 | MtoA | *Sideroxydans lithotrophicus* | 355 | Identical-sequence duplicate; not additional diversity |
| WP_283743965.1 | MtoA | *Sideroxydans* sp. CL21 | 359 | Sideroxydans diversity |
| WP_283743089.1 | MtoA | *Sideroxydans* sp. CL21 | 354 | Sideroxydans diversity |
| ADL56010.1 | MtoA | *Gallionella capsiferriformans* ES-2 | 343 | Genus diversity |
| WP_013293942.1 | MtoA | *Gallionella capsiferriformans* | 343 | Identical-sequence duplicate; not additional diversity |

**HMM build statistics:**

| Model | NSEQ | EFFN | LENG |
|---|---|---|---|
| MtrA | 15 | 0.487061 | 326 |
| MtoA | 6 | 0.629883 | 344 |

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
| *Shewanella* MtrA | 600–706 | 255–292 | Iron-Reduction | Correctly assigned |
| *Shewanella* MtrD (DMSO) | 530–552 | ~260 | **not annotated** | Below GA (580); gene context confirms |
| *Ferrimonas* / *Paraferrimonas* | 603–637 | 263–273 | Iron-Reduction | Correctly assigned |
| Unclassified Shewanellales MAG | 690–724 | ~292 | Iron-Reduction | Correctly assigned |
| *Geobacter* (uses OmcB/Z, not MtrA) | <50 | <50 | **not annotated** | Different protein family entirely |
| *Gallionella* / *Sideroxydans* MtoA | 223–262 | 536–646 | Iron-Oxidation | Correctly assigned |
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

## ImoA / CymA correction — MtrMto oxidation-side gene fix (2026-07-24)

### Problem

The `MtrMto` rule's oxidation-direction gene list was `MtoA, MtrB_TIGR03509, CymA`. This
"CymA" traces back to Liu et al. 2012 (*Front. Microbiol.*, doi:10.3389/fmicb.2012.00037),
who originally called the third gene in the *Sideroxydans lithotrophicus* ES-1 `mtoAB`
locus "CymA-like" based on loose sequence similarity to *Shewanella*'s CymA — both are
members of the broad NapC/NirT cytochrome c superfamily (Pfam PF03264), which also
includes NirT, NapC, and TorC across unrelated respiratory pathways (nitrate/nitrite/TMAO
respiration) in completely different organisms.

A follow-up study (Jain et al. 2022; [doi:10.1128/mbio.02150-22](https://doi.org/10.1128/mbio.02150-22)) renamed Slit_2495 as ImoA and supported its distinction from CymA using phylogeny and sequence comparison. Its direct functional result is **quinol oxidation** and complementation in a Shewanella CymA background. This does not demonstrate direct Fe(II) oxidation by ImoA or establish a single physiological electron-flow direction for the Mto pathway. The Iron-Oxidation label denotes pathway association.

So the registered "CymA" model was never validated against a real CymA ortholog for the
oxidation side of this rule — it was misnamed from the start. Real CymA (Myers & Myers
1997, *J. Bacteriol.*, doi:10.1128/jb.179.4.1143-1152.1997) is itself a real, separate,
well-characterized gene: a tetraheme cytochrome c required for *Shewanella*'s reduction of
Fe(III), nitrate, and fumarate (TnphoA knockout mutant `CMTn-1` lost all three activities;
complementation with the cloned gene restored Fe(III) reduction to 47.1 µmol/day vs.
wild-type's 46.1, from 0.890 in the mutant — Table 2 of that paper). It's mechanistically
upstream of MtrA (feeds electrons from the menaquinone pool toward MtrA), not itself part
of the periplasmic MtrA-MtrB-MtrC trio, and not the correct gene for the `MtrMto` rule's
oxidation side.

The existing registry `CymA` model itself (`hmm_library/Iron-Reduction/CymA.hmm`) was
FeGenie-derived with only 4 training sequences and had been deprecated as
`deprecated_nseq_insufficient` — so in practice the "CymA" entry in the `MtrMto` gene list
had *no active model backing it at all*, meaning the oxidation side of this rule could
never reach 3/3 completeness regardless of what was misnamed.

### Resolution — two separate models, correctly named and calibrated

**ImoA** (`hmm_library/Iron-Oxidation/ImoA.hmm`) — new model for the real oxidation-side
gene. Seed data availability mirrors the MtrA/MtoA situation: the documented build uses two distinct sequences; this is not an exhaustive inventory of ImoA-containing genomes in NCBI.

| Accession | Organism | Len | Role |
|---|---|---|---|
| ADE12720.1 | *Sideroxydans lithotrophicus* ES-1 | 198 | Canonical (Slit_2495, the paper's own subject) |
| VVC82587.1 | *Sideroxydans* sp. CL21 | 198 | Genus diversity |

**CymA** (`hmm_library/Iron-Reduction/CymA.hmm`) — rebuilt from real, verified CymA
sequences (previous FeGenie-derived model deprecated for insufficient seed count).

| Accession | Organism | Len | Role |
|---|---|---|---|
| Q8E8S0.1 | *Shewanella oneidensis* MR-1 | 187 | Canonical — directly cross-referenced to GenBank U75974, the exact accession from Myers & Myers 1997 |
| CAQ9471801.1 | *Shewanella xiamenensis* (MAG) | 187 | Species diversity |
| CAQ8020556.1 | *Shewanella* sp. (MAG) | 187 | Species diversity |
| BGZ77281.1 | *Shewanella marina* | 188 | Species diversity |
| BGZ62728.1 | *Shewanella algae* | 187 | Species diversity |

Both built with the same `mafft --auto` → `hmmbuild` pipeline used for the MtrA/MtoA
rebuild. Checked against each other and against a limited negative-reference set (real
NirT, NapC, TorC — P24038, P0ABL5, P33226) to verify genuine separation, not just
self-consistency:

```
ImoA HMM   vs self (n=2):            473.5, 475.1 bits
           vs best confounder (NirT): 300.4 bits   → ~173-bit gap
           vs CymA seeds:            94.8–106.8 bits (clearly non-members)

CymA HMM   vs self (n=5):            393.6–431.0 bits
           vs best confounder (TorC): 119.3 bits   → ~274-bit gap
           vs ImoA seeds:            106.2–110.2 bits (clearly non-members)
```

These checks show separation on the documented training sequences and selected confounders. They are not independent cross-validation or proof of specificity across the entire superfamily.

### TC values

Following the same "TC = minimum confirmed true positive" convention as MtrA/MtoA:

| Model | TC (bits) | Rationale |
|---|---|---|
| **ImoA** | 473 | Minimum documented ImoA seed score (CL21 homolog, 473.5; inferred member) |
| **CymA** | 393 | Minimum documented CymA seed score (*S. marina*, 393.6; inferred member) |

### Known limitations

Both models are registered `status=experimental` (nseq < 10, same threshold/warning
convention as every other sparse model in this registry) — `n=2` for ImoA and `n=5` for
CymA are genuinely small seed sets, a real data-availability limit like MtoA's own n=6, not
a curation shortcut. The score gaps above apply to these selected references; neither
model has been validated against as broad a "universe" scan as the 3794-sequence one used
for MtrA/MtoA — treat hits with the same caution the registry's `experimental` status flag
already signals throughout the pipeline.

---

## Topology pilot extended to all experimental models (2026-07-24)

`expected_topology` (Component 7, `docs/cluster_confidence_scoring.md`) was originally
piloted on MtrA/MtoA only. Extended to the remaining 10 `status=experimental` models —
every one independently verified against real literature/structural data before being
registered, not assumed from gene-family name alone (a naive "P-type ATPase = multi-pass
TM" blanket rule would have been wrong for at least two of these ten: `Cu_resistance_K_1`
and the RND-system components, none of which are simple polytopic membrane proteins).
None of the ten are cytochromes, so `expected_heme_c` stays unset for all — only
`expected_topology` applies.

| Stem | Real identity | `expected_topology` | Evidence |
|---|---|---|---|
| `Cr_transport_2` | ChrA (CHR superfamily chromate efflux) | `multi_pass_tm` | 13-TMS topology in *P. aeruginosa* ChrA (Díaz-Pérez et al. 2007, *FEBS J.*, doi:10.1111/j.1742-4658.2007.06141.x); consistent across TCDB family 2.A.51 and multiple independent papers. |
| `Zn_transport_troB_1` | TroB, ABC transporter ATPase subunit | `cytoplasmic` | Cluster-9 ABC transporter architecture: TroA (periplasmic solute-binding) / TroB (ATPase) / TroC+TroD (membrane permease) — TroB is the cytoplasmic nucleotide-binding domain, not membrane-integral (Hazlett et al., *JBC*; Agrobacterium TroCBA work, *AEM* doi:10.1128/aem.00299-16). |
| `Zn_transport_troC_1` | TroC, ABC transporter membrane permease | `multi_pass_tm` | Same system as above — TroC(+TroD) form the "heterodimeric cytoplasmic membrane permease" component. |
| `Cu_ATPase_1` | Generic P1B-type Cu-exporting ATPase | `multi_pass_tm` | P1B-ATPases have 8 TM helices (2 N-terminal MA/MB + 6 core M1-M6) — Structural basis of ion uptake in Cu-ATPases, *Nat. Commun.* doi:10.1038/s41467-022-32751-w. |
| `Cu_ATPase_5` | HMA5 | `multi_pass_tm` | Same P1B-ATPase structural family. |
| `Cu_ATPase_pacS_2` | PacS | `multi_pass_tm` | Same P1B-ATPase structural family. |
| `Cu_resistance_K_1` | CopK | `signal_peptide` | Solved structure (PDB 2K0Q, 2LEL) — small periplasmic all-β protein, secreted, zero TM helices (Bianco et al., *J. Mol. Biol.* doi via ScienceDirect S0022283608005718). **Not** membrane-integral — would have been wrong as `multi_pass_tm`. |
| `CuNi_ATPase_1` | CtpD-like Co/Ni-exporting P1B ATPase | `multi_pass_tm` | Same P1B-ATPase family, 6-8 TM helices (*M. tuberculosis* CtpD, UniProt P9WPT3, as a directly comparable named example). |
| `NiCoCd_resistance_nccB_1` | NccB, RND membrane-fusion protein (CzcB/AcrA family) | `signal_peptide` | **Lower confidence.** MFPs are anchored via N-terminal lipidation (palmitoylated Cys), not a classic Sec-cleaved signal sequence or a true multi-pass TM bundle (MexA/AcrA structural work, *Nat. Commun.* doi:10.1038/ncomms10731) — mapped to the closest of our three categories, but whether Phobius reliably flags a lipoprotein anchor as signal-peptide-positive in practice is genuinely uncertain; the conservative `topology_verdict` design means a misclassification here degrades to neutral, not a wrong penalty. |
| `NiCo_resistance_cnrC_3` | CnrC, RND outer-membrane factor (CzcC/TolC family) | `signal_peptide` | OMF proteins are β-barrels, Sec-secreted with a cleaved signal peptide for OM insertion (Kittleson et al./RND efflux architecture review). Phobius models α-helical TM topology only, not β-barrels, so this presents as signal-peptide-positive with ~0 detected TM helices — genuinely `signal_peptide`, not `multi_pass_tm`, despite being a membrane protein. |

All 10 verified against real *Cupriavidus metallidurans* CH34 genome data (the classic
multi-metal-resistance reference strain, in `benchmark/genomes/`): confirmed real HMM hits
for chromate transport, TroB, NccB, PacS, and Cu_ATPase_1 fire correctly with the new
`expected_topology` registered and flowing through to `topo_weight` (neutral without
Phobius installed, as expected — the registrations are wired and will activate the moment
Phobius runs).

**CymA and ImoA** (the two cytochromes among the 12 experimental models) were registered
separately with `expected_heme_c=4` — both confirmed tetraheme via
`operon.count_heme()` run on the real curated sequences (Q8E8S0.1 for CymA per the 1997
paper's own count; ADE12720.1 for ImoA, computed directly with the pipeline's own
function rather than counted by hand). `expected_topology` is deliberately left **unset**
for both: Myers & Myers 1997 describes the whole NapC/NirT family's shared fold as a
single *uncleaved* N-terminal hydrophobic membrane anchor with the hydrophilic
heme-binding majority facing the periplasm — an architecture that doesn't cleanly fit any
of the three `expected_topology` categories (`signal_peptide` implies a cleaved Sec
signal; `multi_pass_tm` implies ≥2 helices; `cytoplasmic` implies neither). Forcing a label
here would be less accurate than leaving it unset — `heme_verdict` alone already provides
real signal for these two (as demonstrated in the earlier real-genome validation: both
fired `topo_weight=1.2` from heme count alone, without Phobius).

---

## New categories added (2026-05-23)

| Category | Models | Rationale | Key references |
|---|---|---|---|
| `Iron-StressResponse` | `flavodoxin_long` (TIGR01752.1), `flavodoxin_short` (TIGR01753.1) | Iron-starvation response biomarkers. Flavodoxins substitute for ferredoxin when iron is limiting. Detecting them in metagenomes indicates Fe-stress capacity. NOT transcriptional regulators; separated from `Iron-Regulation`. | [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0) |
| `IronSulfur-Assembly` | `sufC` (TIGR01978.1), `sufB` (TIGR01980.1), `sufS` (TIGR01979.1), `IscS` (TIGR02006.1) | Fe-S cluster biosynthesis machinery. SUF system (oxidative-stress / iron-limitation pathway); ISC system (constitutive). Allows distinguishing organisms relying on stress-tolerant SUF vs housekeeping ISC assembly. | SUF: [doi:10.1074/jbc.M308004200](https://doi.org/10.1074/jbc.M308004200); IscS: [doi:10.1074/jbc.M401261200](https://doi.org/10.1074/jbc.M401261200) |

---

## Exceptions to the Layer A minimum-coverage rule (2026-07-16)

**Status:** Kept `active` despite violating the ≤3-seq ("others" tier)
deprecation threshold documented above.

| Model | nseq | Category | Reason kept active |
|---|---|---|---|
| `mofA_MnOxGeneTool` | 3 | `Manganese-Oxidation` | At the threshold, not below it. MofA is a genuinely rare, narrowly-distributed Mn(II)-oxidizing multicopper oxidase clade — MnOxGeneTool's own reference database only contains 3 representative sequences for this clade. No broader alternative source exists. |
| `mopA_A_MnOxGeneTool` | 1 | `Manganese-Oxidation` | Single-sequence model (Alphaproteobacteria MopA clade). Same reason: MnOxGeneTool's curated reference set has only 1 sequence for this specific clade split. |
| `mopA_E_MnOxGeneTool` | 1 | `Manganese-Oxidation` | Single-sequence model (Epsilonproteobacteria MopA clade). Same reason. |

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

**Escalation pipeline (2026-07-20, tier 3 added 2026-07-24) — all three tiers
implemented:** a
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
- **Tier 3 (Baktfold, Foldseek+ProstT5 structural search) — implemented.**
  `src/efesto/baktfold.py` runs `foldseek easy-search` on the same flagged
  subset (`--run_struct_confirmation`, requires `--foldseek_db_dir` and
  Foldseek installed; never auto-invoked, same opt-in policy as tier 2).
  Foldseek translates the amino-acid sequences to the 3Di structural
  alphabet internally via ProstT5 (`--prostt5-model`) — no separate
  structure-prediction step or torch/transformers dependency. Comparison is
  a loose substring match against the reference structure's name (same
  naming convention as HMM seeds, e.g. `MtrA_ref_1zzz`). New `struct_weight`
  factor in `scoring.py`, same weakest-link/neutral-on-no-data policy as
  `annot_weight`. The reference structure DB is intentionally small and
  curated — not the full AlphaFold DB/PDB, same narrow-scope reasoning as
  the synteny-bootstrapping rejection above. Recommended first pilot target:
  **MtrCAB**, a real experimental structure (crystal structure, Edwards
  et al. 2020, *Cell* 181(3):665-673, doi:10.1016/j.cell.2020.03.032) —
  a stronger starting point than
  sufA/D/E, whose structural data is sparser. Populating that DB is separate
  curatorial follow-up work, not part of the `struct_weight` code path.
  See `docs/cluster_confidence_scoring.md`, Component 6.

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

## Sequence and primary-literature audit (2026-11-05)

Inputs and a per-sequence evidence manifest are now preserved under `data/seeds/curated/`. Five current models (MtrA, MtoA, CymA, ImoA, DA_402) were rebuilt with identical scoring matrices and HMM statistical parameters. This confirms build provenance, not independent predictive accuracy. MtrA's actual LENG is 326; its 15 input records contain 14 unique sequences. MtoA's six records contain four unique sequences. The original record multiplicities were preserved; deduplication would require rebuilding and reassessing thresholds.

Six MtrA seeds are officially described only as DmsE-family proteins. Their specific MtrA orthology remains to be tested with genomic context and phylogeny. A source organism's iron-reducing phenotype is not per-protein evidence. CymA supports several respiratory acceptors. ImoA has demonstrated quinol oxidation activity; the current category records pathway association, not direct Fe(II)-oxidase activity.

The CbcL anchor is confirmed: Q74GH2 is identical to GSU0274/AAR33608.1, explicitly identified by [Antunes et al. 2022](https://doi.org/10.3389/fmicb.2022.898015). It scores 1331.3 bits against CbcL (cutoff 200). Seven of the eight original seed accessions and the original alignment remain unavailable, so the full CbcL build is not yet reproducible.

The 65 earlier category conflicts were reviewed against primary literature, with qualifications retained in `category_curation.tsv`. NrpY/NrpG have locus-level biosynthetic support rather than individual enzyme assays; the original Fit study supports iron acquisition but does not establish coprogen as the substrate; YfeC/D belong to an Fe/Mn uptake system. Broad KO families must not be interpreted as ligand-specific experimental identifications.

NF041018/19/20 are NCBIfam CutA/B/C glyceraldehyde dehydrogenase profiles imported via the `interpro` channel. Swiss-Prot cross-references and real reference-sequence scores (1295.7/436.4/313.1 bits, respective runtime cutoffs 790/350/290) confirm their identity. Their exclusion from iron transport remains appropriate. [Kardinahl et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10095793/).

## CutABC physical removal (2026-11-05)

NF041018, NF041019 and NF041020 were removed from the shipped library, registry and map/cutoff files after reference-sequence and primary-literature verification. These are carbon-metabolism glyceraldehyde dehydrogenase subunits. Metal cofactors provide no basis to retain them as metal transport/cycling markers. Their curation records remain as removal tombstones so rebuilding from upstream or a reviewed import report does not restore them. This removal does not remove an enabled metal-cycling model: all three were already excluded from runtime.
