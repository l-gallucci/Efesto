# Efesto — Session Checkpoint
**Date:** 2026-06-10  
**Target:** mSystems publication  
**Test suite:** 91 tests, all passing  
**Version:** pyproject.toml v1.0.6  
**Active branch:** manganese-metabolism

---

## Project overview

HMM-based annotation of iron cycling and metal resistance genes in genomes and metagenomes. Built on FeGenie (Garber et al. 2020, ISME J). Adds metagenome support (Pyrodigal ORF calling), cluster confidence scoring, operon prediction (UniOP), antiSMASH BGC integration, and rich output formats.

Entry point: `Efesto` CLI → `src/efesto/cli.py:main()`

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

**Registry:** `hmm_registry.tsv` — 649 total rows  
**Columns:** stem, name, acc, category, gene_name, source, hmm_file, nseq, cutoff, added_date, status, reference, validated_in

**Source breakdown (active):**
- fegenie: 196 (doi:10.1038/s41396-019-0570-7)
- tabuteau: 130 (doi:10.1111/1462-2920.70218)
- methmmdb: 115 (doi:10.1101/2024.12.26.629440)
- mnoxgenetool: 12 (doi:10.1021/acs.est.5c01235)  ← NEW
- ncbifam: 18  ← updated (+3 Mn models)
- interpro: 8
- curated: 2

**Cutoff status:**
- All ncbifam and mnoxgenetool models have calibrated bitscores in `HMM-bitcutoffs.txt`
- 129 models have cutoff=0: methmmdb (115) + fegenie (14)
  - Current behavior: fallback to `--zero_cutoff_min_bitscore` (default 30.0) + E-value < 0.1

**Zero-cutoff fegenie models (14):** HemO/HmoB/HupZ heme oxygenases + 11 siderophore synthesis/transport models. FeGenie shipped 40.0 as a blanket cutoff — apply this for publication.

**nseq < 10 warning:** fires at startup grouped by source. Includes mofA (3), mopA_A (1), mopA_E (1), mopA_R (15) in mnoxgenetool — expected, treat hits cautiously.

**Files:**
- `HMM-bitcutoffs.txt`: `{file_stem}\t{bitscore}` — 12 mnoxgenetool + 3 ncbifam entries added 2026-06-10
- `FeGenie-map.txt`: `{file_stem}\t{gene_name}` — 12 mnoxgenetool + 3 ncbifam entries added 2026-06-10
- `hmm_registry.tsv` stem column: must equal `Path(hmm_file).stem` — verified for all new rows

**Important:** `stem` in registry must equal `Path(hmm_file).stem`. This is what `cat_hmms` uses (via `hf.stem` filesystem scan). If they diverge, cutoffs and gene_names silently fail.

---

## Manganese pathway — current state (branch: manganese-metabolism)

### Categories and HMM counts

| Category | HMMs | Status |
|---|---|---|
| `manganese_acquisition` | 2 | Complete (MntH + MntA SBP) |
| `mn_gene_regulation` | 1 | Complete (MntR) |
| `metal_resistance-manganese` | 2 | Complete (MntP efflux + MntATPase) |
| `manganese_oxidation` | 12 | Complete (MnOxGeneTool models) |

**`manganese_acquisition/`:**
- `mntH_NF001923.hmm` — NRAMP Mn²⁺ importer, ncbifam, GA=430.0, nseq=284
- `Mn_binding_mntA_1.hmm` — MntA SBP (ABC importer periplasmic component), methmmdb, cutoff=0, nseq=25

**`mn_gene_regulation/`:**
- `mntR_NF003025.hmm` — DtxR-family Mn²⁺ sensor/repressor, ncbifam, GA=138.0, nseq=52

**`metal_resistance-manganese/`:**
- `mntP_NF008546.hmm` — CDF-family Mn²⁺ efflux pump, ncbifam, GA=214.0, nseq=44
- `Mn_transport_ATPase_1.hmm` — Mn-exporting P-type ATPase, methmmdb, cutoff=0, nseq=69

**`manganese_oxidation/`** — all from MnOxGeneTool (Wang et al. 2025, doi:10.1021/acs.est.5c01235):

| Stem | Gene | nseq | Cutoff | Notes |
|------|------|------|--------|-------|
| `mnxG_B_MnOxGeneTool` | mnxG | 195 | 668.1 | Bacillus-clade MCO |
| `mnxG_P_MnOxGeneTool` | mnxG | 1071 | 994.4 | Proteobacteria-clade MCO |
| `mcoA_MnOxGeneTool` | mcoA | 1336 | 587.8 | P. putida MCO |
| `mofA_MnOxGeneTool` | mofA | 3 | 906.9 | Leptothrix MCO ⚠️ nseq=3 |
| `moxA_MnOxGeneTool` | moxA | 6716 | 384.0 | Pedomicrobium peroxidase |
| `mopA_P_MnOxGeneTool` | mopA | 514 | 1058.5 | Pseudomonadota heme peroxidase |
| `mopA_R_MnOxGeneTool` | mopA | 15 | 858.3 | Roseobacter heme peroxidase |
| `mopA_A_MnOxGeneTool` | mopA | 1 | 1568.0 | Alphaproteo ⚠️ nseq=1 |
| `mopA_E_MnOxGeneTool` | mopA | 1 | 1181.1 | Epsilonproteo ⚠️ nseq=1 |
| `cotA_MnOxGeneTool` | cotA | 730 | 461.5 | Bacillus laccase |
| `boxA_MnOxGeneTool` | boxA | 247 | 392.9 | Arthrobacter-type |
| `mokA_MnOxGeneTool` | mokA | 192 | 390.1 | Novel Mn oxidase |

**Excluded from manganese_oxidation intentionally:**
- `cueO` — already in `metal_resistance-copper`
- `copA` — already in `metal_resistance-copper`
- `katG` — universal stress gene, too nonspecific for Mn oxidation annotation

**NOT implemented (future):**
- `manganese_reduction` — uses same MtrCAB machinery as iron reduction; no separate category needed
- MnxE/MnxF accessory protein HMMs — unique to Mnx complex (no homologs outside Mn oxidation), good for operon disambiguation; not in any public database, need curated seeds from Butterfield et al. 2013/2015 sequence data
- Multi-label annotation — would require changes to `clustering.py` and `writers.py`

### io.py _ANNOTATE_MAP Mn selectors

```python
"Mn":           ["manganese_acquisition", "manganese_oxidation", "mn_gene_regulation",
                  "metal_resistance-manganese"],
"Mn-metabolism":["manganese_acquisition", "manganese_oxidation", "mn_gene_regulation",
                  "metal_resistance-manganese"],
"Mn-acquisition":["manganese_acquisition"],
"Mn-oxidation":  ["manganese_oxidation"],
"Mn-regulation": ["mn_gene_regulation"],
"Mn-resistance": ["metal_resistance-manganese"],
```

### Overlap decisions (MCO fold problem)

MCO genes (mnxG, mcoA, cotA, moxA) share the MCO fold with copper MCOs (CueO, CopA). Disambiguation strategy:
- **Option A (implemented):** MnOxGeneTool uses phylogenetically curated subfamily HMMs with validated bitsccore thresholds — `mnxG_B` and `mnxG_P` are distinct from CueO at the sequence level with high cutoffs
- **Option B (future):** Operon-context rules: MnxG always co-occurs with MnxE/MnxF in the mnxDEFG operon; need to add a MNX_COMPLEX operon rule and build MnxE/MnxF HMMs from seeds

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
| `Efesto-results.csv` | CSV | Main output, `#` separator between clusters |
| `Efesto-results-long.tsv` | TSV | Tidy format, includes `cluster_confidence`, `model_nseq`, `uniop_context` |
| `Efesto-gene-summary.csv` | CSV | Per-gene, no sequences |
| `Efesto-heatmap.csv` | CSV | Category × genome presence matrix |
| `Efesto-results.gff3` | GFF3 | Written when genome_coords available |
| `Efesto-summary-stats.tsv` | TSV | RUN/CONFIDENCE/CATEGORY/GENOME sections |
| `Efesto-anvio-functions.tsv` | TSV | anvi-import-functions compatible |
| `Efesto-anvio-gene-scores.tsv` | TSV | anvi-import-misc-data compatible (--anvio) |

---

## Pending tasks (priority order)

### Manganese pathway — next steps

1. **Layer B dedup — MntH redundancy** (30 min)  
   `mntH_NF001923` (ncbifam, 284 seqs) and `Cation_transport_mntH_2` (methmmdb) both model the same NRAMP protein family.  
   Action: run `scripts/layer_b_dedup_all.py`, confirm they hit same sequences, then deprecate the methmmdb model in registry (`status=deprecated`).

2. **MnxE/MnxF HMMs for operon disambiguation** (Option B, ~1–2 hours)  
   MnxE/MnxF are unique to the Mnx Mn-oxidation complex (no homologs elsewhere). Literature sources: Butterfield et al. 2013 (PNAS), Butterfield et al. 2015 (BBA). Sequences in Bacillus sp. PL-12 and SG-1.  
   Action: fetch MnxE/MnxF protein sequences from NCBI, build HMMs with `hmmbuild`, add to `manganese_oxidation/`, add operon rule `MNX_COMPLEX` to `operon_rules.json`.

3. **SitABCD — Fe/Mn dual transporter decision** (15 min discussion)  
   SitABCD is an ABC transporter used for both Fe²⁺ and Mn²⁺ in pathogens (Salmonella, Streptococcus). Pending decision: add to `manganese_acquisition` only, `iron_acquisition-iron_transport` only, or both (multi-label, not yet implemented).  
   Action: decide, then look for NCBIfam equivalog.

4. **Mur — Fur-family Mn regulator** (deferred, no NCBIfam equivalog)  
   Mur is found in Rhizobiaceae; no validated NCBIfam equivalog exists. TIGR02783 was confirmed wrong (plasmid conjugation gene). Would require curated HMM from seeds.

### Before submission (required)

5. **14 fegenie zero-cutoff models → 40.0** (30 min)  
   Models: HemO/HmoB/HupZ heme oxygenases + 11 Sid_ siderophore models  
   Action: add entries to `HMM-bitcutoffs.txt` with bitscore=40.0, update registry cutoff column

6. **Swiss-Prot calibration for all fegenie models** (1–2 days)  
   Run `hmmsearch` against Swiss-Prot, find bitscore gap between annotated iron-function hits and non-hits. Required for defensible per-model cutoffs in mSystems methods section.

7. **Reference genome benchmark** (publication-critical, ~1 week)  
   Test Efesto against genomes with known iron cycling / metal resistance phenotypes. Report precision/recall per category. Compare against FeGenie on same dataset.

8. **MetHMMDB developer contact** (see questions list below)

9. **Fill `validated_in` column** in registry after benchmark runs

10. **MtrA/MtoA false positive table** from `_calibration/calibration_report.tsv`

### Before submission (recommended)

11. **Expand nseq** for rusticyanin (4), cytochrome579 (4), MtoA (6), Fox genes (5–6)  
    Options: add sequences from RefSeq/UniProt for these well-characterized genes.

12. **Environmental case study** on SRA metagenomes  
    Demonstrates metagenome applicability for mSystems audience.

13. **docs/cluster_confidence_scoring.md** — verify max_bp_gap values against literature

14. **docs/hmm_expansion_biological_rationale.md** — add manganese section documenting rationale for all 4 Mn categories and model choices

### Nice-to-have

15. Metal resistance HMM calibration (copper/arsenic/nickel/etc. — all methmmdb, cutoff=0)

16. Multi-label annotation implementation (clustering.py + writers.py) — needed for SitABCD (Fe+Mn), MtrCAB (Fe+Mn reduction)

17. MnxE/MnxF HMMs → operon rule for Mn oxidation complex disambiguation from Cu MCOs

---

## MetHMMDB developer contact — questions

1. Do per-model recommended bitscore cutoffs exist, or is E-value < 0.1 their design intent?
2. Are low-nseq models (84 of 115 have nseq < 10) considered reliable or experimental?
3. Do they have a reference genome benchmark with known metal resistance phenotypes (shareable)?
4. Which model families have known false positive issues (especially transporters vs. resistance genes)?
5. Is journal publication of the bioRxiv preprint (Dec 2024) imminent? (methods citation)
6. Courtesy notice: their library is integrated in Efesto, targeting mSystems.

---

## Registry facts (critical for future work)

- `stem` column = HMM file stem (filename without `.hmm`) — must match or cutoffs/gene_names silently fail
- `hmm_file` format: `category/stem.hmm` (relative to hmm_library root, no leading `hmm_library/`)
- `status != "active"` → model skipped at runtime
- `cutoff = 0` → uses `fallback_bitscore` (CLI arg `--zero_cutoff_min_bitscore`, default 30.0)
- `gene_name` column → display name in outputs (via FeGenie-map.txt lookup by file stem)
- 58 HMM files have embedded GA/TC/NC — registry cutoffs already match GA values, no action needed

---

## Test runner

```bash
conda run -n metalgenie-evo python3 -m pytest tests/ -x -q --import-mode=importlib
```

Conda env: `metalgenie-evo` (NOT `efesto`)

---

## Code constraints

- NEVER run git commit/push (user handles all git operations)
- Never hardcode paths in scripts
- `scripts/` = developer/curation tools only
- No comments unless WHY is non-obvious
- Good code formatting always

---

## Source references

| Source key | DOI/URL | Label |
|---|---|---|
| fegenie | doi:10.1038/s41396-019-0570-7 | Garber et al. 2020, ISME J |
| tabuteau | doi:10.1111/1462-2920.70218 | Tabuteau et al. 2025, Environ Microbiol |
| methmmdb | doi:10.1101/2024.12.26.629440 | Kciuchcinski et al. 2025, bioRxiv |
| ncbifam | NCBI prokaryotic annotation | NCBIfam |
| interpro | doi:10.1093/nar/gkac993 | Paysan-Lafosse et al. 2023, NAR |
| mnoxgenetool | doi:10.1021/acs.est.5c01235 | Wang et al. 2025, Environ Sci Technol |
