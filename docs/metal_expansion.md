# Lanthanides, tungsten and metallophores

## Imported profiles

The eleven profiles below are unchanged official imports. They were not trained
from new Efesto seed collections. Full-sequence GA values are inherited, not
recalibrated from the small Efesto reference panels.

{{expansion_table}}

Source URLs, versions, hashes and reference expectations are recorded under
`data/validation/metal_expansion/` and `data/validation/methanobactin/` in the
[source repository](https://github.com/l-gallucci/Efesto/tree/main/data/validation).

## Lanthanides and tungsten

LanM and LanD describe binding/handling. LanD interaction with LanM is supported;
a proposed physiological sorting/export mechanism is not established by a gene hit.
No generic EF-hand, TonB or ABC match is substituted for Ln-specific uptake.

Wtp components are assigned to shared Mo/W acquisition. Tup profiles retain their
W-associated source scope, but TupA NF041772 has known coverage gaps: characterized
D. alaskensis Q316W1 and P. acidaminophilum Q93KD6 fall below GA400. The model is
experimental and its negative result must not be read as absence of W acquisition.
The cutoff was not lowered to recover these few references.

## Opine-type metallophores

CntK/L/M are `Metallophore-Synthesis`, not Fe-specific siderophore synthesis.
Normal genomic-context filtering requires CntL+CntM. CntK is optional across
products, so its absence does not by itself exclude the core pathway. Individual
hits do not resolve staphylopine, pseudopaline or another opine product.

## Methanobactin

Normal genomic-context output requires the three specific MbnA/B/C families in
the same detected cluster. Generic MNIO/DUF692 membership alone is insufficient.
Catalog/all-results modes retain individual-family evidence.

The canonical MbnA precursor is about 30 residues and may escape gene calling.
Thus a real MbnBC locus with an uncalled precursor can be absent from ordinary
context-filtered output. Efesto does not yet rescue small ORFs. TIGR04071 also
matches MbnA-like domains fused to some sigma factors: `mbnA-domain` is not an
automatic claim of a standalone precursor. Co-occurrence remains predicted
biosynthetic potential, not demonstrated product, expression or copper uptake.

MbnT is not assigned from a generic TonB profile.
