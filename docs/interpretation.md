# Interpreting results

## What a hit supports

An accepted match identifies a protein-family relationship under the configured
search gates. Functional specificity depends on the profile and the evidence
behind its members. A source-supported label is not an independent accuracy test.

- Acquisition, resistance, homeostasis and regulation describe different roles.
- Generic transporter components do not determine substrate specificity.
- A shared domain can occur in paralogs with different genomic roles.
- DNA presence does not demonstrate expression, stress induction or activity.
- A missing hit can reflect incomplete assemblies, uncalled small genes or limited
  model coverage; it does not establish biological absence.

## Confidence and completeness

`cluster_confidence` combines selected evidence using heuristic weights. It is
not a probability or a directly comparable measure of biological activity across
gene systems. Optional evidence can be missing and neutral. A high score is not
a guarantee that every component is present or correctly assigned.

The normal output does not globally discard scores below 0.5. Use the
[scoring reference](cluster_confidence_scoring.md) to decide and document your own
analysis thresholds.

## Specific limitations

| Case | Interpretation |
|---|---|
| Fe-S assembly profiles | Cofactor-biogenesis potential, not proof of iron respiration or stress |
| SufE-related domain | May match CsdE; inspect orthology and genomic context |
| CymA | Multifunctional respiratory component; not an Fe-specific terminal reductase |
| ImoA | Deprecated in this library; native Fe(II) oxidation not established by the cited evidence |
| TupA NF041772 | Known characterized references fall below GA; negative detection is not absence of W uptake |
| LanM / LanD | Lanthanide binding and handling, not a complete uptake system |
| Cnt families | Opine-type metallophore potential; not a unique product or metal assignment |
| MbnA/B/C | Predicted biosynthetic potential; the short precursor may be missed by gene calling |
| MntX, Lut/LanA, Vup candidates | Development results, not enabled runtime classifiers |

## Comparisons and publication

Compare the same proteins, library versions and filtering conditions. Distinguish
raw family assignment from the effects of genomic-context rules. Changes in
category naming or retention can change summaries without changing HMM scores.

A metagenomic case study demonstrates applicability; it is not automatically a
labelled accuracy benchmark. The [validation page](validation.md) separates
software checks, targeted reference panels and missing independent validation.
