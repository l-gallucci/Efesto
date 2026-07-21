# Anvi'o Integration

Efesto produces two Anvi'o-compatible output files when run with `--anvio`.

---

## Output files

### `Efesto-anvio-functions.tsv`

Import with `anvi-import-functions`. Provides gene functional annotations
(gene name + functional category) that appear in the Anvi'o interactive interface
under the Efesto annotation source.

**Columns:** `gene_callers_id`, `source`, `accession`, `function`, `e_value`

```bash
anvi-import-functions \
    -c CONTIGS.db \
    -i results/Efesto-anvio-functions.tsv \
    -p Efesto
```

### `Efesto-anvio-gene-scores.tsv`

Import with `anvi-import-misc-data` (target: genes). Provides numeric per-gene
scores for visualisation in the Anvi'o interactive interface.

**Columns:** `gene_callers_id`, `cluster_confidence`, `co_occ_score`,
`hmm_weight`, `uniop_weight`, `bgc_boost`

```bash
anvi-import-misc-data \
    -c CONTIGS.db \
    --target-data-table genes \
    results/Efesto-anvio-gene-scores.tsv
```

---

## Without Bakta mapping (simple workflow)

```bash
# 1. Run Efesto with FAA input (Prodigal ORFs)
efesto \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --anvio

# 2. Import into Anvi'o (gene_callers_id = Prodigal ORF integer IDs)
anvi-import-functions \
    -c CONTIGS.db \
    -i results/Efesto-anvio-functions.tsv \
    -p Efesto
```

---

## With Bakta mapping (`--bakta_gff_dir`) — recommended

When using `--fna_dir` + `--bakta_gff_dir`, Efesto runs Prodigal internally,
then matches each Prodigal ORF to its Bakta gene ID via coordinate overlap (±3 bp
tolerance). The `gene_callers_id` column in both TSVs will contain Bakta IDs
(e.g. `AMXMAG_00053`), mapping directly to Anvi'o integer IDs when the contigs
database was built from Bakta external gene calls.

### Step-by-step

```bash
# Step 1 — Build Anvi'o contigs database from Bakta external gene calls
anvi-script-process-genbank \
    --genbank-file  genome.gbff \
    --output-db     CONTIGS.db \
    --gene-calls    bakta_gene_calls.tsv \
    --annotation    bakta_annotation.tsv

anvi-gen-contigs-database -f genome.fna -o CONTIGS.db \
    --external-gene-calls bakta_gene_calls.tsv \
    --annotation          bakta_annotation.tsv

# Step 2 — Run Efesto with Bakta GFF3 for ID mapping
efesto \
    --fna_dir        assemblies/ \
    --bakta_gff_dir  bakta_output/ \
    --hmm_dir        hmm_library/ \
    --out            results/ \
    --anvio

# Step 3 — Import both files
anvi-import-functions \
    -c CONTIGS.db \
    -i results/Efesto-anvio-functions.tsv \
    -p Efesto

anvi-import-misc-data \
    -c CONTIGS.db \
    --target-data-table genes \
    results/Efesto-anvio-gene-scores.tsv
```

---

## What you get in Anvi'o

After import, in the interactive interface:

- **Functions layer:** gene annotations visible in the gene-detail panel; search
  by gene name or category (`anvi-search-functions`)
- **Gene scores layer:** `cluster_confidence`, `hmm_weight`, `uniop_weight`,
  `co_occ_score`, `bgc_boost` as numeric tracks per gene. Colour genes by
  confidence to visually identify high-quality iron cycling loci on contigs.

---

## Architecture note

Anvi'o's `anvi-import-misc-data` with `--target-data-table genes` requires a
TSV where the first column is `gene_callers_id` (integer) and all other columns
are numeric. Efesto's gene-scores TSV is formatted to satisfy this
requirement. See the [Anvi'o misc-data documentation](https://anvio.org/help/main/programs/anvi-import-misc-data/)
for details on the expected format.

The `confidence` string column (`calibrated` / `low_confidence`) is intentionally
omitted from the gene-scores TSV (non-numeric); it appears in
`anvio-functions.tsv` as part of the function annotation instead.
