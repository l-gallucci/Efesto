# Siderophore model scope examples

The four public reference proteins illustrate a known domain-specificity limit.
With the bundled PF06276.7 profile and the unchanged runtime cutoff (25.2 bits),
AcsD, IucA and FhuF all exceed the threshold. PF06276 therefore cannot distinguish
NIS biosynthesis from ferric-siderophore reduction by itself. The YqjH profile
recognizes the reference YqjH in this panel.

`reference_sources.json` records sequence sources, retrieval date, checksums and
functional papers. `scores.tsv` records full-sequence HMMER 3.4 results on the
four-protein FASTA. `PF06276.json` is the provider's entry description, which
explicitly includes both FhuF and IucA/IucC proteins. `profile_integrity.json`
records the pre-reclassification profile hashes, cutoffs, status and categories.
The scope tests preserve these search parameters; they do not estimate precision,
recall, coverage of sequence diversity or a false-positive rate.

The IroB/KO K23725 pair has the same functional target but different HMMs and
training sizes. Both are retained; they are not independent evidence of two
biological pathway steps. Distinct-model neighborhood minima are contextual
heuristics and do not establish a complete biosynthetic pathway.
