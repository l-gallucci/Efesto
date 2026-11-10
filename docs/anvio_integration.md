# Anvi'o integration

`--anvio` exports functional annotations and numeric gene scores:

```bash
efesto --faa_dir proteins/ --gff_dir gff/ --out efesto_out --anvio
```

| File | Columns |
|---|---|
| `Efesto-anvio-functions.tsv` | `gene_callers_id`, `source`, `accession`, `function`, `e_value` |
| `Efesto-anvio-gene-scores.tsv` | `gene_callers_id`, `cluster_confidence`, `co_occ_score`, `hmm_weight`, `uniop_weight`, `bgc_boost` |

## Match gene identifiers first

The column name `gene_callers_id` does **not** mean that Efesto has converted
protein labels to the integer IDs in an anvi'o contigs database. It normally
contains the original ORF label, or a mapped Bakta locus tag when available.

Export the gene calls from the actual target database and build a mapping using
matching contigs, coordinates and strand:

```bash
anvi-export-gene-calls -c CONTIGS.db -o gene_calls.tsv
```

Rewrite the first column of both Efesto exports with the corresponding integer
IDs, checking for unmatched and duplicate assignments. Bakta label mapping alone
does not guarantee the required anvi'o integer IDs. If every label is already the
correct integer ID for that database, a further conversion is unnecessary.

## Import mapped files

```bash
anvi-import-functions -c CONTIGS.db -i mapped-functions.tsv
anvi-import-misc-data mapped-gene-scores.tsv \
  -c CONTIGS.db --target-data-table genes
```

Efesto can use `--fna_dir` with `--bakta_gff_dir` to recover Bakta labels by
coordinate matching. It still searches its internally called proteins and can
retain unmapped original IDs. Inspect the mapping summary before import.

Keep anvi'o and its databases in an appropriate separate environment. Consult
the official [functions import](https://anvio.org/help/main/programs/anvi-import-functions/)
and [misc-data import](https://anvio.org/help/main/programs/anvi-import-misc-data/)
references for your installed version.
