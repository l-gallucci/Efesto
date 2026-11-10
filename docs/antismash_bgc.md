# antiSMASH integration

Supply previously generated antiSMASH GFF3 annotations:

```bash
efesto --faa_dir proteins/ --gff_dir gff/ \
  --bgc_dir antismash_gff/ --out efesto_out
```

Efesto matches each genome stem to `<stem>.gff3` or `<stem>.gff`, either directly
in the directory or below `<stem>/`. Contig names and coordinates must refer to
the same assembly used by Efesto.

Only GFF features of type `region` or `cluster` are considered. The `product`
attribute, or `type` when no product is available, must contain a case-insensitive
`siderophore` or `metallophore` keyword. Generic NRPS/PKS annotations do not
qualify by themselves.

An overlapping retained ORF gives the cluster `bgc_boost = 1.2`; otherwise the
factor is 1.0. This is one component of the full [confidence formula](cluster_confidence_scoring.md),
which is capped at 1.2. Other optional components can also boost the score.

BGC overlap adds related annotation context. It does not independently establish
product chemistry, expression, Fe specificity or metal uptake. Efesto does not
run antiSMASH itself or convert a generic NRPS region into a validated metallophore.
