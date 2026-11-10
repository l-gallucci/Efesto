# HMM development methodology

The proposed hybrid protocol combines reference-anchored sequence/phylogeny
curation, as used in the supplied [Mankowski preprint](https://doi.org/10.64898/2026.09.02.748891),
with nested calibration and independent testing inspired by
[HMS-S-S](https://doi.org/10.1111/1755-0998.13642). It is an Efesto adaptation,
not a claim to reproduce every step of either study.

## Evidence and retrieval

Collect experimentally linked anchors and related sequences from family searches,
retaining accession/version, sequence hash, source and inclusion rationale.
Distinguish experimental anchors, inferred homologs, related controls and unknowns.
A Pfam/KO family match or gene-name similarity is not a functional label by itself.

The current initial collection uses UniProt and official NCBI/Pfam records.
Mankowski's later proGenomes3 expansion and EPA-ng confirmation have not been
performed here. Newly retrieved homologs remain candidates until curatorial
evidence supports their assignment.

## Group before splitting

Collapse exact aliases and keep near homologs in the same partition. The current
pilot uses 80% identity and 80% bilateral coverage with conservative genus grouping
for background controls. Separate analyses recompute links at 70% and 90% identity.
Clustering is approximate and does not establish functional equivalence.

Freeze an external test before model selection when enough independent examples
exist. Do not manufacture external positives by relabeling HMM hits.

## Build and calibrate

MAFFT `--auto` aligns training positives and HMMER 3.4 `hmmbuild` produces the
profiles. Alignment and build are repeated in each fold. The MntX pilot uses
three outer and two inner folds because only three independent positive groups
are available; that is an exploratory stability exercise.

Cutoffs maximize MCC on inner out-of-fold predictions, using full-sequence scores
and an E-value gate. Ties favor fewer false positives, then the lower cutoff.
Outer folds evaluate the resulting assignments. The final candidate is trained
on development data only; its cutoff is not chosen from an external test.
HMS-S-S originally optimizes F1; adopting MCC here is an Efesto choice.

Do not confuse these quantities:

- Full-sequence and domain thresholds are different.
- Official GA values are retained for unchanged source imports.
- A development cutoff is not an independently validated production threshold.
- NSEQ counts records, not independent experimental examples.

## Phylogeny and evaluation

Current exploratory trees use MAFFT, a declared gap-column filter and IQ-TREE
LG+G4 with 1,000 SH-aLRT/UFBoot replicates. Branch support measures a split;
it does not prove substrate or transport direction for every member.

Assess related-family confusion, model competition, genomic-context effects and
fragmentation separately. Report TP/FP/FN/TN, precision, recall, F1, MCC and
uncertainty, including failed folds and evidence gaps. Promote a candidate only
with a justified scope and sufficient validation. See [current results](validation.md).
