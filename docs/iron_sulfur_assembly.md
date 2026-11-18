# Fe-S assembly

Fe-S biogenesis profiles describe machinery that builds or transfers iron-sulfur
cofactors used in electron transfer, metabolism and regulation. They do not detect
all Fe-S-containing proteins or specifically identify iron-respiring organisms.

## Shipped assembly profiles

{{fes_table}}

These are imported NCBIfam/TIGRFAM and Pfam profiles. The cutoffs in the table
are taken directly from the **full-sequence GA (gathering threshold)** in each
source HMM; they were not estimated de novo by Efesto. Where an HMM provides
separate sequence and domain GA values, Efesto uses the sequence value. NSEQ
is the number of training records reported by the source model.

Source versions and the corresponding runtime thresholds are preserved in
the registry and cutoff table. See [cutoff provenance](hmm_library_curation.md#cutoff-provenance).

## Interpretation

- SUF can be stress-responsive or constitutive/sole machinery depending on the
  organism. DNA presence does not establish current oxidative stress or Fe limitation.
- IscS is a sulfur-mobilizing family; one IscS hit is not an inventory of a complete ISC system.
- The SufE-related profile can recognize CsdE. A shared domain does not by itself
  assign the gene to a SUF locus.
- Related SufA/IscA/ErpA and SufB/SufD families warrant orthology and context checks.

The SUF contextual rule requires multiple distinct components in a detected cluster;
this supports a local assembly system without demonstrating physiological activity.
Some profiles carry `needs_confirmation` metadata for optional annotation/structure
checks. Their availability and comparison results are reported separately; a
missing result does not confirm the assignment. Passing the minimum context
requirement does not establish a complete SUF variant.

See [curation](hmm_library_curation.md) and [separate evidence](evidence.md).
