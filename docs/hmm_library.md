# HMM library

Efesto combines profiles from FeGenie, MetHMMDB, MnOxGeneTool, Tabuteau et al.,
Pfam and NCBI protein family models, alongside profiles built for Efesto.
Shared categories and runtime tables allow these sources to be searched in one
workflow. Per-model records preserve their origins, references and search cutoffs;
a shared category does not make the profiles functionally interchangeable.
See [source citations](references.md) and [curation and provenance](hmm_library_curation.md).

{{catalog_summary}}

Counts below are generated from `hmm_registry.tsv` during the website build.
They describe the source revision used to build this site, which may differ from
an older installed checkout.

## Runtime status

| Registry status | Ordinary run |
|---|---|
| `active` | Enabled unless removed by category selection |
| `experimental` | Enabled, with model-specific interpretation limitations |
| `deprecated*` / `inactive*` | Excluded before search |

## Categories

Names distinguish acquisition, homeostasis, regulation, resistance and other
processes. Shared-substrate and unresolved categories retain that uncertainty.
A category name does not establish a measured phenotype.

{{category_table}}

## Sources

{{source_table}}

Source labels record the import channel as well as the original family database.
Check each registry row and its reference; for example, an InterPro import may
contain a Pfam or NCBIfam model. Functional scope and references are recorded
for each profile.

The [per-model license/source table](../hmm_library/model_licenses.tsv) records
the original provider, identified terms and unresolved conditions for every
registered profile. See [HMM license scope](../hmm_library/LICENSES.md) and
[third-party notices](third-party-notices.md); runtime status does not
change the redistribution terms of bundled files.

## Selectors

Examples: `--annotate Fe Cu Zn`, `--annotate Ln W`, or an exact canonical category
such as `Metallophore-Synthesis`. `all` is the default. See the complete valid
list in the [command reference](command-reference.md).

- `Ln` selects LanM/LanD homeostasis profiles.
- `W` includes Tup and shared Mo/W families; `Mo` also includes shared Mo/W.
- `metallophores` includes metallophore and siderophore synthesis categories.
- Resistance-specific selectors preserve their resistance meaning.

See [siderophores and metallophores](siderophores.md) for the complete profile
table, shared-family limits and the distinction between synthesis, uptake, export
and utilization.

## Registry and search thresholds

The registry records filename stem, model name/accession, category, source,
status, cutoff, NSEQ, references and interpretation notes. Filename stems key the
runtime maps; the HMM `NAME` can legitimately differ from that stem.
`category_review_status` describes the supporting source or an unresolved
category assignment.

The runtime reads full-sequence cutoffs from `HMM-bitcutoffs.txt`. A missing or
nonpositive cutoff uses the generic bitscore floor, default 30; an E-value gate
also applies. Keep registry, cutoff table, gene map and files synchronized.
NSEQ is a count of training records, not independent experimental positives.

See [curation](hmm_library_curation.md), [functional scope](hmm_expansion_biological_rationale.md)
and [custom libraries](extending_the_library.md).

## Cutoff provenance

Direct imports of official Pfam and NCBIfam/TIGRFAM profiles, including enabled
member-database profiles imported through InterPro, use the full-sequence GA
value from the source HMM. These thresholds are inherited from the provider,
not estimated de novo in Efesto. Profiles imported through another collection
retain that collection’s operating cutoff, which can differ from the HMM GA;
for example, FeGenie supplies 25.2 bits for the bundled PF06276.7 profile.

Other published collections can supply cutoffs in external metadata.
Efesto-built profiles use Efesto-selected operating thresholds.
Profiles without a positive model-specific cutoff use the generic bitscore
floor. See [curation and construction methods](hmm_library_curation.md#cutoff-provenance)
for the distinction between these threshold sources.

## Cutoff consistency

The maintained `hmm_registry.tsv` records the authoritative curated model cutoffs.
`HMM-bitcutoffs.txt` is the runtime table and must contain the same numeric values,
without rounding or duplicate entries. Both files preserve the full numeric precision of model cutoffs. The
cutoff source and model-specific interpretation notes are recorded in the registry.
