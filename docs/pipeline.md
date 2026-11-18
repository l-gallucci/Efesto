# How Efesto processes a dataset

Efesto first decides which protein-family matches to retain, then reports the
available evidence for those matches. HMMER bitscores remain part of the search;
there is no combined confidence score or final confidence cutoff.

```text
Assemblies (FNA) → Pyrodigal/Pyrodigal-GV → proteins + coordinates
Existing proteins (FAA) + optional matching GFF ────────────┘
       ↓
Enabled HMM profiles → bitscore and E-value gates → best accepted hit per ORF
       ↓
Genomic grouping → context/per-protein filters
       ↓
Optional eggNOG rescue of flagged, first-pass rejected candidates
       ↓
Retained annotations + separate context, tool and UniOP-pair evidence
       ↓
Tables, cluster-count/coverage matrices, matched sequences and GFF3
```

## 1. Prepare proteins and coordinates

With `--fna_dir`, Efesto predicts genes using Pyrodigal-GV by default; Pyrodigal
can be selected explicitly. With `--faa_dir`, it searches the supplied proteins.
Gene-prediction caches are checked against assembly content, caller/version and
meta/single mode. Only complete, unchanged FAA/GFF outputs are reused; a prediction
failure stops annotation. Only assemblies selected for the current run are searched.
Matching GFF coordinates allow grouping by contig, genomic separation and,
optionally, strand. See [input workflows](workflows.md).

Without coordinates, the fallback assumes Prodigal-style ordinal identifiers.
A gene catalog instead uses singleton clusters; catalog numbering does not
establish genomic proximity.

## 2. Search the enabled profiles

Active and experimental profiles are enabled unless category selection removes
them. Profiles whose registry status starts with `deprecated` or `inactive` are
excluded. Experimental profiles have model-specific interpretation notes in the
[catalog and scope guide](hmm_library_curation.md).

Each search uses a positive model-specific full-sequence cutoff where available.
Otherwise the generic bitscore floor defaults to **30**. The E-value limit defaults
to **1e-5** and applies to every accepted hit. `threshold_type` records which
cutoff path was used. The maintained registry records the authoritative model
cutoffs; the bundled runtime table preserves those values without rounding.

Efesto retains one accepted HMM per ORF using the highest raw bitscore across
accepted matches. This can conceal competing-family matches: inspect tblout
files retained with `--keep_tblout` when investigating paralogs. Bitscores across
models are not a common probability scale.

## 3. Group and filter family matches

With coordinates, neighboring hits are grouped by base-pair separation and,
optionally, strand. Named systems can use tighter rule-specific gaps than the
global `--max_bp_gap`.

The default context engine derives from FeGenie with Efesto additions. A
library-local `operon_rules.json` selects the configurable engine. These rules
can remove matches or alter Mtr/Mto category assignments. In the default engine,
each system is evaluated independently: one supported system does not exempt a
neighboring incomplete system, and unrelated models do not supply category minima.
The declared trusted-lone exceptions remain available for siderophore transport.
Existing per-protein
checks, including Cyc2 length/motif requirements, still apply.

The [context reference](operon_rules.md) describes modes and requirements.
A minimum count refers to distinct models and does not establish a complete or
active biological system.

## 4. Optionally recover flagged candidates

Only hits on `needs_confirmation` models that failed the first-pass filters enter
the eggNOG rescue set. Informative matching name/KO evidence can recover a
candidate as a singleton. Candidates without matching informative evidence remain
excluded. Supplying a whole-proteome annotation file does not check every retained
hit automatically.

`Efesto-confirmation-candidates.tsv` reports both recovered and unretained
candidates. It is not an additional annotation count table.

## 5. Report optional evidence separately

| Analysis | Added observations | Changes retention? |
|---|---|---|
| UniOP | Operon membership and available pair probabilities | No |
| antiSMASH inputs | Qualifying siderophore/metallophore region overlaps | No |
| Foldseek | Candidate reference target, match statistics and label comparison | No |
| Phobius and optional motif comparison | Predictions and consistency with model expectations | No |

Tool applicability and availability are explicit. A missing result is not support
or a demonstrated contradiction. Optional Phobius/motif reporting is separate
from sequence checks already used in the filtering stage.

See [evidence and decisions](evidence.md) for worked MtrA/LanM, SUF and singleton
examples, and [output files](output_formats.md) for the exact tables.
`Efesto-filter-decisions.tsv` records the actual filter actions on selected HMM
candidates, including exclusions and the final outcome after any eggNOG rescue.
It starts after HMM acceptance and best-hit selection, so subthreshold matches
and competing models for the same ORF are not included.
`Efesto-context-evidence.tsv` describes retained neighborhoods and their rule
observations; it is not the execution trace.

## Resource use

At most two HMM searches overlap by default. `--threads` and `--hmm_threads` can
reduce concurrency; they do not define a hard RAM limit. Gene prediction has its
own worker allocation. See [performance](performance.md).
