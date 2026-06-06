# AntiSMASH BGC Integration

Efesto can integrate antiSMASH biosynthetic gene cluster (BGC) predictions
to boost the confidence score of clusters that co-localise with siderophore or
metallophore biosynthetic regions.

---

## Usage

```bash
Efesto \
    --faa_dir  orfs/ \
    --gff_dir  gff/ \
    --hmm_dir  hmm_library/ \
    --out      results/ \
    --bgc_dir  antismash_output/
```

`--bgc_dir` accepts a directory of antiSMASH output. Efesto searches for
GFF3/GFF files matching each genome's stem (e.g. `genome.gff3`, `genome/genome.gff3`).

---

## What triggers the boost

Only siderophore-type regions trigger the 1.2× `bgc_boost`:

- GFF3 features with `product=siderophore` or `type=siderophore`
- GFF3 features with `product=metallophore` or `type=metallophore`

Matching is case-insensitive. Generic NRPS/PKS clusters without explicit siderophore
annotation do **not** trigger the boost.

For a cluster to receive the boost, at least one ORF in the cluster must overlap a
siderophore BGC region on the same contig.

---

## Effect on confidence score

```
cluster_confidence = min(1.2,  hmm_weight  ×  co_occ_score  ×  uniop_pair_score  ×  bgc_boost)
```

`bgc_boost` is 1.2× when overlap detected, 1.0 otherwise. This can push a
near-perfect cluster (score ~1.0) to 1.2 — the only mechanism that allows
`cluster_confidence > 1.0`.

The boost reflects that co-location with an antiSMASH-annotated siderophore BGC
is strong independent genomic evidence for iron acquisition function, beyond what
HMM hits and operon structure alone can confirm.

---

## antiSMASH GFF3 file discovery

Efesto resolves GFF files from `--bgc_dir` using the genome stem:

1. `<bgc_dir>/<stem>.gff3`
2. `<bgc_dir>/<stem>.gff`
3. `<bgc_dir>/<stem>/<stem>.gff3`
4. `<bgc_dir>/<stem>/<stem>.gff`

Where `stem` is the genome filename without extension (e.g. `MAG_001.faa` → stem `MAG_001`).

If no matching file is found for a genome, that genome's clusters receive `bgc_boost = 1.0`
(no boost, no error).

---

## Interpreting results

`bgc_boost` appears as a column in:
- `Efesto-results-long.tsv`
- `Efesto-anvio-gene-scores.tsv`

A boost of 1.2 in the gene-scores TSV identifies specific ORFs overlapping a
siderophore BGC. Use this to cross-validate HMM-based siderophore synthesis calls
with antiSMASH NRPS/PKS predictions.

---

## Limitations

- antiSMASH must be run separately before Efesto; GFF3 output is expected.
- Boost is applied at the cluster level — if any ORF overlaps the BGC, all ORFs
  in that cluster receive the boost.
- The boost only applies to siderophore-annotated BGCs. Iron reduction, oxidation,
  and resistance clusters are unaffected by `--bgc_dir`.
- Coordinate overlap is evaluated on the same contig only.
