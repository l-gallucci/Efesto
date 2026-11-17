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

Matching region coordinates and product labels are reported for each overlapping
retained ORF in the [evidence table](evidence.md). A neighboring gene outside
the region does not inherit that overlap. Missing GFF input, unavailable
coordinates and no qualifying overlap are reported separately. Overlap does
not change which HMM hits are retained.

BGC overlap adds related annotation context. It does not independently establish
product chemistry, expression, Fe specificity or metal uptake. Efesto does not
run antiSMASH itself or convert a generic NRPS region into a validated metallophore.
