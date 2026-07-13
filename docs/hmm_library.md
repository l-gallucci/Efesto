# HMM Library

Efesto ships with a curated library of profile HMMs in `hmm_library/`.
All models are tracked in a versioned registry (`hmm_library/hmm_registry.tsv`)
with provenance, bitscore cutoffs, and curation status.

---

## Sources

| Source | Models (active) | Description |
|--------|----------------|-------------|
| FeGenie | 196 | Original iron-cycling profiles from Garber et al. 2020 |
| Tabuteau et al. 2025 | 130 | Iron acquisition profiles from KOfam, FeGenie, and NCBI NF* |
| MetHMMDB | 115 | Metal mobility resistance gene HMMs (broad metal scope) |
| NCBIfam / TIGRFAM | 16 | Curated iron oxidation, regulation, sulfur assembly, molybdenum resistance, and Tad-pilus (CpaB) models |
| Pfam | 1 | Flp/Fap pilin (Tad pilus), sourced directly from Pfam-A |
| Curated (custom-built) | 2 (+2 pre-existing: MtrA, MtoA) | CbcL and DA_402 (MISO extracellular MHC) — built from literature-confirmed seed sequences, background-tested against the rest of the library |
| **Total active** | **470** | After deduplication (167 deprecated; 637 total) |

---

## Categories

### Iron cycling

| Category | Description | Notable models |
|----------|-------------|----------------|
| `iron_oxidation` | Enzymatic oxidation of Fe²⁺ → Fe³⁺ | MtoA, MtrA, FoxABC, Cyc2, rusticyanin, cytochrome579 |
| `iron_reduction` | Enzymatic reduction of Fe³⁺ → Fe²⁺ | MtrA, MtrC, OmcS, OmcZ, DFE operons, CbcL, DA_402, Flp/CpaB (Tad pilus) |
| `probable_iron_reduction` | Fe reduction genes with lower specificity | CymA, omcE |
| `possible_iron_oxidation_and_possible_iron_reduction` | Dual-assignment before Mtr/Mto disambiguation | See MtrMto operon rule |
| `iron_storage` | Cellular iron storage proteins | Ferritin, bacterioferritin, Dps |
| `iron_gene_regulation` | Fe-responsive transcriptional regulators | Fur, DtxR, IscR, NsrR, PerR, SoxR |
| `iron_stress` | Iron-starvation biomarkers | Flavodoxin long, flavodoxin short |
| `iron_sulfur_assembly` | Fe-S cluster biosynthesis machinery | SufB, SufC, SufS, IscS |
| `magnetosome_formation` | Magnetosome island genes | MamA/B/E/K/P/M/Q/I/L/O |

### Iron acquisition (siderophore-centred)

| Category | Description |
|----------|-------------|
| `iron_acquisition-siderophore_synthesis` | NRPS and related enzymes; NRP/PK siderophore biosynthesis |
| `iron_acquisition-siderophore_transport` | Outer membrane TonB-dependent receptors for Fe-siderophore complexes |
| `iron_acquisition-siderophore_transport_potential` | TonB-ExbBD motor proteins; ABC transporter permease/ATPase subunits |
| `iron_acquisition-heme_oxygenase` | Enzymes releasing iron from heme |
| `iron_acquisition-heme_transport` | Outer membrane heme receptors, ABC transporters |
| `iron_acquisition-iron_transport` | ABC iron transporters (FbpABC, FeoABC, EfeU, FutABC, etc.) |

### Metal resistance (MetHMMDB)

| Category | Description |
|----------|-------------|
| `metal_resistance-arsenic` | Arsenate reductase, arsenite efflux (ArsABC, ArsH) |
| `metal_resistance-chromium` | Chromate efflux (ChrA) |
| `metal_resistance-cobalt_zinc_cadmium` | CzcABC RND efflux, ZnuABC uptake, ZinT, FieF, Zur |
| `metal_resistance-copper` | CopABCD, CusSR, CopY |
| `metal_resistance-mercury` | MerABCDEPRT |
| `metal_resistance-molybdenum` | ModABC ABC transporter, ModC ATPase (TIGRFAM) |
| `metal_resistance-multimetal` | CuAg/CdCoZn RND systems, broad-spectrum |
| `metal_resistance-nickel` | NikABCDE, NicT, CznABC |
| `metal_resistance-silver` | CuAg_CusA/C (overlapping CopA) |
| `metal_resistance-non-specific` | Broad-spectrum and non-metal genes (after curation) |
| `metal_resistance-tellurium` | TerABCDEZ, TehAB |

---

## Registry format (`hmm_library/hmm_registry.tsv`)

Tab-separated, one row per model. Columns:

| Column | Description |
|--------|-------------|
| `stem` | HMM file stem; must exactly match filename without `.hmm` |
| `name` | Human-readable gene name (shown in outputs) |
| `accession` | Pfam/TIGRFAM accession or empty |
| `category` | Functional category (directory name) |
| `hmm_file` | Relative path from `hmm_library/` |
| `nseq` | Training sequences used to build the HMM |
| `cutoff` | Calibrated bitscore cutoff (0 = zero-cutoff fallback) |
| `date_added` | ISO date |
| `status` | `active` or `deprecated_*` |
| `reference` | DOI or citation |

> **Registry rule:** `stem` must exactly match the HMM file's `NAME` field and the
> filename without `.hmm`. Mismatches silently break cutoff lookup and gene-name mapping.
> Run `python scripts/curate_hmm_library.py --verify hmm_library/` to detect drift.

---

## Bitscore thresholds

`hmm_library/HMM-bitcutoffs.txt` — two-column TSV: `stem<TAB>bitscore`.

| Models | Threshold type |
|--------|----------------|
| 337 calibrated models | Per-HMM gathering cutoff (GA) from TIGRFAM/Pfam or manual calibration |
| 129 zero-cutoff models | No calibrated threshold; use `--zero_cutoff_min_bitscore` (default 30.0) |

Zero-cutoff models: MetHMMDB (115) + 14 FeGenie siderophore models. Raise
`--zero_cutoff_min_bitscore` to reduce false positives from these models.

---

## Gene name mapping

`hmm_library/FeGenie-map.txt` — two-column TSV: `stem<TAB>gene_name`.

Gene names from this file appear in the `gene` column of all output files.
When a model has no entry, the raw stem is used as the gene name.

---

## Deduplication layers

See [[Library curation log|hmm_library_curation]] for full decision log.

| Layer | Method | Deprecated |
|-------|--------|-----------|
| A — Minimum nseq | Source-specific thresholds (fegenie ≤ 5, methmmdb ≤ 10) | 60 |
| A — Name dedup | Same stem, keep higher nseq/coverage | 88 |
| B — Sequence dedup | MMseqs2 70% id / 80% cov on hmmemit consensus | 12 |
| Category mismatch | Manual review of cross-category hits | 7 |

---

## Startup provenance output

At startup, Efesto prints a provenance table listing active sources,
model counts, and any models with `nseq < 10` (flagged as `low_training_data`):

```
Source          Models  Calibrated  Zero-cutoff  nseq<10
FeGenie            196         182           14        0
Tabuteau           130          98           32        0
NCBIfam             15          15            0        0
MetHMMDB           115           0          115       10
────────────────────────────────────────────────────────
Total              456         295          161       10

WARNING: 10 models have nseq < 10 (limited training data):
  [methmmdb] CdCoZn_efflux_czcD_1 (nseq=7), ...
```

Models with `nseq < 10` are included in runs but appear as `low_confidence` in
output. Identify them via the `model_nseq` column in `results-long.tsv`.

---

## Normalising legacy HMM formats

All models in the shipped library are in HMMER3/f format. If you add models built
with HMMER < 3.1, run once with `--normalize_hmms` to convert them in-place.
