# Pipeline

```text
FNA ── Pyrodigal/Pyrodigal-GV ── FAA + coordinates
                                │
FAA + optional matching GFF ─────┘
       ↓
enabled profile HMMs → score and E-value gates → best hit per ORF
       ↓
genomic clustering → context rules → applicable confirmation checks
       ↓
confidence annotation → summaries, TSV, matched sequences and GFF3
```

## Search and assignment

Profiles whose registry status starts with `deprecated` or `inactive` are excluded.
Active and experimental profiles are searched unless category selection removes
them. A positive model-specific full-sequence cutoff is used where available;
otherwise the generic bitscore floor defaults to 30. The E-value limit defaults
to `1e-5` and applies to all accepted hits.

One best HMM hit is retained per ORF using the highest raw bitscore across accepted
matches. This rule can conceal competing-family matches; inspect retained tblout
files when investigating paralog assignments.

## Clusters and context

With coordinates, hits are grouped on a contig by base-pair separation and,
optionally, strand. Named systems may use tighter rule-specific gaps than the
global `--max_bp_gap`. Without coordinates, the fallback uses Prodigal-style
ordinal identifiers. A gene catalog instead uses singleton clusters.

The default context engine derives from FeGenie, with Efesto-specific additions.
A library-local `operon_rules.json` selects the JSON engine. Context can change
retention and redox-category assignment; it does not demonstrate physiological
activity. See [operon rules](operon_rules.md).

## Additional evidence

UniOP, antiSMASH, eggNOG-mapper, Foldseek and Phobius are optional. Their
applicability differs: annotation and structural confirmation target models
flagged `needs_confirmation`; topology checks target models with explicit
expectations. These checks are not applied uniformly to every protein.
See [confidence scoring](cluster_confidence_scoring.md).

## Resource use

At most two HMM searches overlap by default. `--threads` and `--hmm_threads`
can reduce that concurrency; they do not define a hard RAM limit. Gene prediction
has its own worker allocation. See [performance](performance.md).
