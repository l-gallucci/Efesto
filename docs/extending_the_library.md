# Extending the HMM Library

---

## Adding individual HMMs manually

1. Place `.hmm` files in `hmm_library/<category>/`
2. Add `<stem><TAB><bitscore>` to `hmm_library/HMM-bitcutoffs.txt`
   (use `0` for zero-cutoff fallback)
3. Add `<stem><TAB><gene_name>` to `hmm_library/FeGenie-map.txt`
4. Add a row to `hmm_library/hmm_registry.tsv`
   - `stem` must exactly match the HMM filename without `.hmm`
   - `status` = `active`
5. Optionally add or update a rule in `operon_rules.json`
6. Run `python scripts/curate_hmm_library.py --verify hmm_library/` to check consistency

> **Registry rule:** `stem` must exactly match the HMM `NAME` field AND the filename
> without `.hmm`. Mismatches silently break cutoff and gene-name lookup.

---

## Normalising HMM format

Models built with HMMER < 3.1 must be converted to HMMER3/f format:

```bash
efesto --normalize_hmms --faa_dir orfs/ --hmm_dir hmm_library/ --out results/
```

Or manually:

```bash
python scripts/normalize_hmm_versions.py hmm_library/
```

Safe to run repeatedly (already-current files are skipped).

---

## Full library rebuild from sources

```bash
python scripts/curate_hmm_library.py \
    --fegenie_dir   /path/to/FeGenie/hmms/iron/ \
    --flat_dir      /path/to/new_iron_hmms  iron_acquisition  tabuteau \
    --methmmdb_json /path/to/MetHMMDB/metadata.json \
    --methmmdb_dir  /path/to/MetHMMDB/ \
    --out_dir       hmm_library/ \
    --log           curation_report.tsv

# Verify after rebuild
python scripts/curate_hmm_library.py --verify hmm_library/
```

---

## Selective annotation with `--annotate`

`--annotate` restricts the run to a subset of functional categories without
modifying the library itself:

```bash
# Annotate only iron metabolism
efesto --annotate Fe-metabolism --faa_dir orfs/ --hmm_dir hmm_library/ --out results/

# Multiple tokens
efesto --annotate Fe-metabolism Cu Zn --faa_dir orfs/ --hmm_dir hmm_library/ --out results/

# All categories (default)
efesto --annotate all --faa_dir orfs/ --hmm_dir hmm_library/ --out results/
```

Accepted tokens: element-level (`Fe`, `Cu`, `Zn`, `Mn`, `Ni`, `Co`, `Mo`, `As`, `Hg`,
`Cd`, `Cr`, `Ag`, `Te`, `Mg`, `multimetal`) or process-level
(`Fe-metabolism`, `Fe-resistance`, `Fe-acquisition`), or `all`.

---

## Layer B deduplication (cross-source)

To run sequence-level deduplication after adding new models:

```bash
python scripts/layer_b_dedup_all.py \
    --hmm_dir hmm_library/ \
    --registry hmm_library/hmm_registry.tsv \
    --out hmm_library/_dedup_work/
```

This runs MMseqs2 `easy-cluster` at 70% identity / 80% coverage on `hmmemit`
consensus sequences across all active models. Review the cluster TSV before
deprecating — low-nseq consensus sequences are degenerate and can cluster
spuriously (see [[Library curation log|hmm_library_curation]] for details).

---

## Building subfamily HMMs (MtrA/MtoA example)

When an existing broad model cross-hits multiple subfamilies (as FeGenie MtrA/MtoA
did), build subfamily HMMs from curated seeds:

```bash
# Build HMMs from seed alignment file
python scripts/build_mtr_mto_subfamily_hmms.py \
    --seeds data/seeds/mtr_mto_seeds.tsv \
    --out   hmm_library/

# Calibrate against a representative universe
python scripts/calibrate_mtr_mto_cutoffs.py \
    --hmm_dir hmm_library/ \
    --universe hmm_library/_calibration/mtr_mto/calibration_universe.faa \
    --out      hmm_library/_calibration/mtr_mto/
```

Calibration reports are stored in `hmm_library/_calibration/<name>/`:
- `calibration_report.tsv` — per-sequence scores
- `calibration_summary.txt` — TC/GA/NC derivation

See [[Library curation log — MtrA/MtoA|hmm_library_curation]]
for the full calibration rationale.

---

## Adding a new category

1. Create `hmm_library/<new_category>/` and place HMMs there
2. Add models to registry, bitcutoffs, and gene map (as above)
3. Decide whether the category should bypass rules:
   - Single-gene hits informative (e.g. `iron_stress`, `iron_sulfur_assembly`) →
     add to `report_all_categories` in `operon_rules.json`
   - Co-occurrence required (e.g. new transport system) → add a new rule object
4. If co-occurrence is required, define `canonical_size` for confidence scoring
5. Run `pytest tests/` — add a test if the new category has custom logic
