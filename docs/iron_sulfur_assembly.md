# Fe-S assembly

Fe-S biogenesis profiles describe machinery that builds or transfers iron-sulfur
cofactors used in electron transfer, metabolism and regulation. They do not detect
all Fe-S-containing proteins or specifically identify iron-respiring organisms.

## Shipped assembly profiles

{{fes_table}}

These are official NCBIfam/Pfam imports selected for Efesto, not newly trained
Efesto models. Exact source versions and cutoffs are retained in the registry.

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
checks. Missing confirmation is neutral, not evidence of correctness.

See [curation](hmm_library_curation.md) and [confidence scoring](cluster_confidence_scoring.md).
