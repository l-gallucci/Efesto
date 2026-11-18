# Interpreting results

## What a hit supports

An accepted match identifies a protein-family relationship under the configured
search gates. Functional specificity depends on the profile and the evidence
behind its members. Consult the model's scope and functional references.

- Acquisition, resistance, homeostasis and regulation describe different roles.
- Generic transporter components do not determine substrate specificity.
- A shared domain can occur in paralogs with different genomic roles.
- DNA presence does not demonstrate expression, stress induction or activity.
- A missing hit can reflect incomplete assemblies, uncalled small genes or limited
  model coverage; it does not establish biological absence.

## Evidence and completeness

Match statistics, threshold type, neighborhood rules and optional checks are
reported separately. A model-specific cutoff controls search acceptance;
functional interpretation also depends on the profile and genomic context.
There is no combined score or global confidence threshold.

A neighborhood rule can require a minimum number of distinct model families;
passing it does not show that a complete or active pathway exists. Repeated
copies and genes from other systems do not become missing functional components.
See [separate evidence](evidence.md) for the availability, scope and meaning of
each observation.

## Practical reading order

1. Check the HMM identifier, model source, bitscore, E-value and effective threshold.
2. Inspect `retention_mode` and any applicable context-rule observations. Use
   `Efesto-filter-decisions.tsv` to trace why a selected HMM candidate was retained,
   excluded or rescued; see the [decision guide](evidence.md#why-was-a-candidate-excluded).
3. Read optional-tool availability before interpreting a comparison or missing value.
4. Keep discordant observations visible; another result does not cancel them out.
5. Distinguish per-ORF counts from the cluster-count matrix before comparing samples.

For a minimum SUF rule, `satisfied_minimum` means that the configured minimum of
distinct model families was observed. It does not mean every SUF component exists.
For MtrA plus LanM hits, LanM does not supply a missing Mtr component. See
[worked examples](evidence.md#worked-examples).

## Specific limitations

| Case | Interpretation |
|---|---|
| Fe-S assembly profiles | Cofactor-biogenesis potential, not proof of iron respiration or stress |
| SufE-related domain | May match CsdE; inspect orthology and genomic context |
| CbcL | Inner-membrane electron transfer associated with low-potential Fe(III) respiration in the reference organism; not a terminal reductase |
| MtrA / MtoA | Related multiheme families; use partner proteins and context to interpret pathway role and direction |
| CymA | Multifunctional respiratory component; not an Fe-specific terminal reductase |
| DA_402 | Candidate electron-transfer family; functional scope remains provisional |
| ImoA | Deprecated in this library; native Fe(II) oxidation not established by the cited evidence |
| TupA NF041772 | Known characterized references fall below GA; negative detection is not absence of W uptake |
| LanM / LanD | Lanthanide binding and handling, not a complete uptake system |
| Cnt families | Opine-type metallophore potential; not a unique product or metal assignment |
| MbnA/B/C | Predicted biosynthetic potential; the short precursor may be missed by gene calling |

## Comparisons and publication

Compare the same proteins, library versions and filtering conditions. Distinguish
raw family assignment from the effects of genomic-context rules. Changes in
category naming or retention can change summaries without changing HMM scores.

Use independently labelled reference data to estimate annotation accuracy.
Metagenomic examples and runtime comparisons illustrate use and resource
requirements; they do not provide sensitivity or specificity estimates.
