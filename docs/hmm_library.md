# HMM library

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
contain a Pfam or NCBIfam model. Source presence is not a per-profile validation.

## Selectors

Examples: `--annotate Fe Cu Zn`, `--annotate Ln W`, or an exact canonical category
such as `Metallophore-Synthesis`. `all` is the default. See the complete valid
list in the [command reference](command-reference.md).

- `Ln` selects LanM/LanD homeostasis profiles.
- `W` includes Tup and shared Mo/W families; `Mo` also includes shared Mo/W.
- `metallophores` includes metallophore and siderophore synthesis categories.
- Resistance-specific selectors preserve their resistance meaning.

## Registry and search thresholds

The registry records filename stem, model name/accession, category, source,
status, cutoff, NSEQ, references and review metadata. Filename stems key the
runtime maps; the HMM `NAME` can legitimately differ from that stem.
`category_review_status` records source support or unresolved review, not a
classifier accuracy estimate.

The runtime reads full-sequence cutoffs from `HMM-bitcutoffs.txt`. A missing or
nonpositive cutoff uses the generic bitscore floor, default 30; an E-value gate
also applies. Keep registry, cutoff table, gene map and files synchronized.
NSEQ is a count of training records, not independent experimental positives.

See [curation](hmm_library_curation.md), [functional scope](hmm_expansion_biological_rationale.md)
and [custom libraries](extending_the_library.md).
