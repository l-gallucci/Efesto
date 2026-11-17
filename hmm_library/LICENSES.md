# HMM license scope

Imported models retain the terms of their original providers. Read
`model_licenses.tsv` together with the model registry and
[third-party notices](../THIRD_PARTY_NOTICES.md). Import channel, original source,
license notice and uncertainty are recorded separately; an empty or unresolved
license field does not grant permission.

## Efesto rebuilt profiles

Copyright (c) 2026 Efesto contributors.

Efesto-owned copyrightable contributions to the following six HMM profiles
are licensed under **Creative Commons Attribution-ShareAlike 4.0 International**:

- `Iron-Reduction/MtrA.hmm`
- `Iron-Oxidation/MtoA.hmm`
- `Iron-Reduction/CymA.hmm`
- `Unclassified-MetalFunction/ImoA.hmm`
- `Iron-Reduction/DA_402.hmm`
- `Iron-Reduction/CbcL.hmm`

License and disclaimer: https://creativecommons.org/licenses/by-sa/4.0/legalcode.
Attribute Efesto contributors and retain the original sequence-source and
literature references. Indicate modifications and preserve the applicable
ShareAlike conditions. This grant does not replace third-party rights in the
training inputs or extend to other imported HMMs.

These models were rebuilt from the documented inputs; their build/provenance
records are in `data/seeds/curated/` in the source repository. The stored NCBI
and UniProt input records retain their original conditions. CbcL's eight public
input sequences are now documented and cross-verified against NCBI; its freshly
reconstructed alignment reproduces the deployed model exactly, apart from its
build date. The historical input files were not recovered; this is a documented
reconstruction, not a claim of independent functional validation.

## Changes to imported material

Category/filename organization, selected cutoffs and interpretation metadata
have been curated by Efesto. Seventeen Pfam profiles imported through FeGenie
were converted from HMMER3/b to HMMER3/f. Imported HMM matrices are not assigned
a new Efesto license by this reorganization or conversion. Per-model source and
modification notes are in `model_licenses.tsv`; source notices remain applicable.
