# MetalGenie-Evo — Session Checkpoint
**Date:** 2026-05-24  
**Target:** mSystems publication  
**Test suite:** 79 tests, all passing  
**Version:** pyproject.toml v1.0.6

---

## Project overview

HMM-based annotation of iron cycling and metal resistance genes in genomes and metagenomes. Built on FeGenie (Garber et al. 2020, ISME J). Adds metagenome support (Pyrodigal ORF calling), cluster confidence scoring, operon prediction (UniOP), antiSMASH BGC integration, and rich output formats.

Entry point: `MetalGenie-Evo` CLI → `src/metalgenie_evo/cli.py:main()`

---

## Module map

| File | Purpose |
|---|---|
| `cli.py` (624 lines) | Main pipeline: arg parsing, FAA loading, hmmsearch, clustering, scoring, output |
| `hmmer.py` | hmmsearch execution (parallel), tblout parsing, HMM library normalization |
| `io.py` | FASTA/GFF readers, registry loader, cutoff/map readers, provenance printing, nseq map |
| `operon.py` | Operon rule engine (FeGenie port + JSON), `build_stem_gap_map`, `build_canonical_size_map`, `second_pass` |
| `clustering.py` | `cluster_by_coordinates` (bp-gap, per-stem gap override), `cluster_by_index` fallback |
| `scoring.py` | `distance_decay`, `co_occurrence_score`, `uniop_pair_score`, `hmm_weight`, `cluster_confidence` |
| `uniop.py` | UniOP binary wrapper, FAA index parsing, pair probability extraction |
| `bgc.py` | antiSMASH GFF3 parser, siderophore region detection, BGC boost |
| `writers.py` | All output writers: CSV/TSV/GFF3/heatmap/summary-stats/anvi'o |
| `coverage.py` | BAM/depth coverage loading, TPM normalization |
| `gene_calling.py` | Pyrodigal + Pyrodigal-GV ORF prediction |

---

## Cluster confidence framework (scoring.py)

Formula: `min(1.2, hmm_w × co_occ_w × uniop_w × bgc_boost)`

**Components:**

| Component | Formula | Notes |
|---|---|---|
| `hmm_weight` | mean bitscore confidence per cluster | calibrated=1.0, low_confidence=0.5 |
| `co_occurrence_score` | completeness × distance_factor × edge_factor | completeness = n_obs / canonical_size (capped 1.0) |
| `uniop_pair_score` | min pairwise UniOP probability | 1.0 if no UniOP data |
| `bgc_boost` | 1.2 if overlaps siderophore antiSMASH region | 1.0 otherwise |

Distance decay: `exp(-ln(2) × gap_bp / 500)` (half-life = 500 bp, Zhang et al. 2006)  
Edge penalty: 0.7× if any ORF within 3000 bp of contig edge

---

## HMM library (hmm_library/)

**Registry:** `hmm_registry.tsv` — 634 total rows, 466 active  
**Columns:** stem, name, acc, category, gene_name, source, hmm_file, nseq, cutoff, added_date, status, reference, validated_in

**Source breakdown (active):**
- fegenie: 196 (doi:10.1038/s41396-019-0570-7)
- tabuteau: 130 (doi:10.1111/1462-2920.70218)
- methmmdb: 115 (doi:10.1101/2024.12.26.629440)
- ncbifam: 15
- interpro: 8
- curated: 2

**Cutoff status:**
- 337 models have calibrated bitscores in `HMM-bitcutoffs.txt`
- 129 models have cutoff=0: methmmdb (115) + fegenie (14)
  - Current behavior: fallback to `--zero_cutoff_min_bitscore` (default 30.0) + E-value < 0.1

**Zero-cutoff fegenie models (14):** HemO/HmoB/HupZ heme oxygenases + 11 siderophore synthesis/transport models. FeGenie shipped 40.0 as a blanket cutoff — apply this for publication.

**nseq < 10 warning:** fires at startup grouped by source (153 total: methmmdb=84, fegenie=45, tabuteau=9, ncbifam=8, interpro=6, curated=1)

**Files:**
- `HMM-bitcutoffs.txt`: `{file_stem}\t{bitscore}` — 15 ncbifam entries added (2026-05-24)
- `FeGenie-map.txt`: `{file_stem}\t{gene_name}` — 15 ncbifam entries added (2026-05-24)
- `hmm_registry.tsv` stem column: fixed for 15 ncbifam/curated models (file stem now matches HMM filename stem)

**Important:** `stem` in registry must equal `Path(hmm_file).stem`. This is what `cat_hmms` uses (via `hf.stem` filesystem scan). If they diverge, cutoffs and gene_names silently fail.

---

## Key operon rules (operon_rules.json + _DEFAULT_OPERON_RULES)

| Rule | canonical_size | max_bp_gap | Key genes |
|---|---|---|---|
| FLEET (MtrABC+Cyc) | 8 | 2000 bp | MtrA/B/C, OmcS/Z, CymA, Cyc1/2 |
| MAM (magnetosome) | 10 | 3000 bp | MamABEKMOPQ |
| FOXABC (iron oxidation) | 3 | 1000 bp | FoxA/B/C |
| FOXEYZ | 3 | 1000 bp | FoxE/Y/Z |
| DFE1/DFE2 | 4–5 | 1000 bp | DFE genes |
| MtrMto | 3 | 2000 bp | MtrA, MtoA, MtoB |
| SIDERO_TRANSPORT | — | 2000 bp | — |
| SIDERO_SYNTH | — | 5000 bp | — |
| IRON_TRANSPORT | — | 2000 bp | — |

Per-category clustering: `build_stem_gap_map()` → `effective_gap = min(gap_a, gap_b)` for each adjacent ORF pair.

---

## Output files

| File | Format | Notes |
|---|---|---|
| `MetalGenie-Evo-results.csv` | CSV | Main output, `#` separator between clusters |
| `MetalGenie-Evo-results-long.tsv` | TSV | Tidy format, includes `cluster_confidence`, `model_nseq`, `uniop_context` |
| `MetalGenie-Evo-gene-summary.csv` | CSV | Per-gene, no sequences |
| `MetalGenie-Evo-heatmap.csv` | CSV | Category × genome presence matrix |
| `MetalGenie-Evo-results.gff3` | GFF3 | Written when genome_coords available |
| `MetalGenie-Evo-summary-stats.tsv` | TSV | RUN/CONFIDENCE/CATEGORY/GENOME sections |
| `MetalGenie-Evo-anvio-functions.tsv` | TSV | anvi-import-functions compatible |
| `MetalGenie-Evo-anvio-gene-scores.tsv` | TSV | anvi-import-misc-data compatible (--anvio) |

---

## Pending tasks (priority order)

### Before submission (required)

1. **14 fegenie zero-cutoff models → 40.0** (30 min)  
   Models: HemO/HmoB/HupZ heme oxygenases + 11 Sid_ siderophore models  
   Action: add entries to `HMM-bitcutoffs.txt` with bitscore=40.0, update registry cutoff column

2. **Swiss-Prot calibration for all fegenie models** (1–2 days)  
   Run `hmmsearch` against Swiss-Prot, find bitscore gap between annotated iron-function hits and non-hits. Required for defensible per-model cutoffs in mSystems methods section.

3. **Reference genome benchmark** (publication-critical, ~1 week)  
   Test MetalGenie-Evo against genomes with known iron cycling / metal resistance phenotypes. Report precision/recall per category. Compare against FeGenie on same dataset.

4. **MetHMMDB developer contact** (see questions list below)

5. **Fill `validated_in` column** in registry after benchmark runs

6. **MtrA/MtoA false positive table** from `_calibration/calibration_report.tsv`

### Before submission (recommended)

7. **Expand nseq** for rusticyanin (4), cytochrome579 (4), MtoA (6), Fox genes (5–6)  
   Options: add sequences from RefSeq/UniProt for these well-characterized genes. Would change HMM file and require re-running calibration.

8. **Environmental case study** on SRA metagenomes  
   Demonstrates metagenome applicability for mSystems audience.

9. **docs/cluster_confidence_scoring.md** — already written, verify max_bp_gap values against literature (flagged as estimates, not extracted from papers)

### Nice-to-have

10. Metal resistance HMM calibration (copper/arsenic/nickel/etc. — all methmmdb, cutoff=0)  
    Needs phenotype-linked reference genome set. Contact MetHMMDB developers first.

11. Use `--cut_ga` flag for ncbifam/interpro HMMs instead of `-T` (domain- vs sequence-level difference is minor; not worth complexity)

---

## MetHMMDB developer contact — questions

1. Do per-model recommended bitscore cutoffs exist, or is E-value < 0.1 their design intent?
2. Are low-nseq models (84 of 115 have nseq < 10) considered reliable or experimental?
3. Do they have a reference genome benchmark with known metal resistance phenotypes (shareable)?
4. Which model families have known false positive issues (especially transporters vs. resistance genes)?
5. Is journal publication of the bioRxiv preprint (Dec 2024) imminent? (methods citation)
6. Courtesy notice: their library is integrated in MetalGenie-Evo, targeting mSystems.

---

## Registry facts (critical for future work)

- `stem` column = HMM file stem (filename without `.hmm`) — must match or cutoffs/gene_names silently fail
- `hmm_file` format: `category/stem.hmm` (relative to hmm_library root, no leading `hmm_library/`)
- `status != "active"` → model skipped at runtime
- `cutoff = 0` → uses `fallback_bitscore` (CLI arg `--zero_cutoff_min_bitscore`, default 30.0)
- `gene_name` column → display name in outputs (via FeGenie-map.txt lookup by file stem)
- 58 HMM files have embedded GA/TC/NC — registry cutoffs already match GA values, no action needed

---

## Code constraints

- NEVER run git commit/push (user handles all git operations)
- Never hardcode paths in scripts
- `scripts/` = developer/curation tools only
- No comments unless WHY is non-obvious
- Good code formatting always

---

## Source references

| Source | DOI/URL | Label |
|---|---|---|
| fegenie | doi:10.1038/s41396-019-0570-7 | Garber et al. 2020, ISME J |
| tabuteau | doi:10.1111/1462-2920.70218 | Tabuteau et al. 2025, Environ Microbiol |
| methmmdb | doi:10.1101/2024.12.26.629440 | Kciuchcinski et al. 2025, bioRxiv |
| ncbifam | NCBI prokaryotic annotation | NCBIfam |
| interpro | doi:10.1093/nar/gkac993 | Paysan-Lafosse et al. 2023, NAR |
