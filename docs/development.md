# Research and development

!!! warning "Development status"
    The candidates below are not enabled runtime classifiers. A file or a research
    result in `data/validation/` does not add a model to the ordinary Efesto run.

Status reflects the source revision used to build this site. Current user-facing
capabilities are documented in [workflows](workflows.md) and the [HMM catalog](hmm_library.md).

| Target | Completed work | Remaining requirement |
|---|---|---|
| Marine MntX | Three-anchor candidate, grouped nested pilot, official-profile comparison, exploratory phylogeny and fragment diagnostics | Broader functional labels, independent external positives and a demonstrated advantage over the source profile |
| LutH / LutH-like | AM1, PA1 and OB3b reference mapping; exploratory TonB phylogeny | More independent groups and related functional controls; AM1/PA1 are nearly identical |
| LanA | Characterized 5GB1C reference and published-query correspondence | Additional independent characterized groups and calibrated discrimination |
| LutAEF | Verified AM1 reference locus | Additional loci, component-specific confounders and threshold calibration |
| VupABC | Historical-to-current accession mapping and functional component reconstruction | Multiple independent loci and separation from Tup, or a justified shared W/V family scope |
| Methylolanthanin / other metallophores | Literature and locus-level candidate assessment | Product-specific context rules and validation before adding metal-specific labels |

## What has been implemented

The development workflow separates curated evidence, similarity groups, calibration
and evaluation. It rejects duplicate/group leakage, rebuilds alignments and HMMs
from training data and preserves outputs for review. The current MntX candidate
has three seeds; it has not replaced an official runtime profile.

## What has not been completed

No proGenomes3 expansion or EPA-ng placement pipeline has been executed. Existing
exploratory trees and broad UniProt retrievals do not constitute that workflow.
No independent external positive test is available for the MntX pilot.

The older curated set also has remaining provenance/specificity limits, including
incomplete original CbcL inputs. See [curation](hmm_library_curation.md).

Read the [methodology](hmm_methodology.md) and [validation results](validation.md).
