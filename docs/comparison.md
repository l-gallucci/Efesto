# Efesto and FeGenie: scope, model changes and comparisons

Efesto integrates multiple profile-HMM collections to annotate metal cycling,
acquisition, homeostasis and resistance. FeGenie contributes iron profiles and
the basis of the default context engine. MetHMMDB, MnOxGeneTool, Tabuteau et al.,
public protein-family databases and Efesto-built profiles contribute other
families and functional layers.

This page recaps how Efesto differs from the pinned FeGenie workflow: shared
profiles retained unchanged, rebuilt redox models, additional source collections,
annotation changes and completed local comparisons. It describes the distributed
library ahead of a dedicated methods publication. Cite the
software version using [CITATION.cff](../CITATION.cff), together with the
[original model and tool sources](references.md) applicable to the analysis.

## Current library at a glance

{{catalog_summary}}

Runtime status describes which profiles are searched; it does not certify
independent biological validation. Source labels below distinguish import
channels, which may differ from the original database behind a model.

## Functional scope beyond iron

The combined library includes manganese, copper, zinc, nickel, cobalt, cadmium,
chromium, mercury, silver, tungsten, lanthanide and metal/metalloid-associated
families. Acquisition, homeostasis, regulation, redox and resistance remain
separate functional categories. The following are examples of enabled profiles,
not counts of experimentally established functions in a sample.

| Scope | Example profile families | Source contribution |
|---|---|---|
| Iron redox, acquisition, heme, storage and magnetosomes | Cyc2, Mtr/Mto, heme-handling and siderophore families | FeGenie, Tabuteau et al., public family databases and Efesto-built profiles |
| Manganese oxidation and acquisition | BoxA, MnxG, MopA; MntA and MntH | MnOxGeneTool, MetHMMDB and NCBI protein family models |
| Copper homeostasis and resistance | CopZ, copper-transporting ATPases, CueO | MetHMMDB |
| Other metal and metalloid functions | ZitB, cadmium-associated ATPases, mercuric reductase, ArsB/Acr3; Ni/Cr/Ag/Te-associated families | MetHMMDB, with source-specific category review |
| Lanthanide binding and handling | LanM, LanD | NCBI protein family models and Pfam |
| Tungsten and shared Mo/W acquisition | TupA/B/C, WtpA/B/C | NCBI protein family models, including TIGRFAM-derived profiles |
| Metallophore biosynthesis | CntK/L/M and MbnA/B/C | KOfam through Tabuteau et al.; NCBI protein family models |
| Fe–S cofactor biogenesis | SufA/B/C/D/E/S and IscS | NCBI protein family models and Pfam |

**Additional profiles are not all newly trained models.** The eleven LanM/LanD,
Wtp/Tup and MbnA/B/C profiles are unchanged official imports with inherited GA
cutoffs. Current counts, exclusions and experimental status are documented in the
[catalog](hmm_library.md),
[curation records](hmm_library_curation.md) and
[lanthanide, tungsten and metallophore guide](metal_expansion.md).

A family match supports genetic potential. It does not by itself establish metal
specificity, pathway completeness, activity or a resistance phenotype. In
particular, Fe–S assembly is a cofactor-biogenesis annotation, LanM/LanD do not
identify a complete uptake pathway, and TupA has a documented coverage gap.

## Model changes relative to FeGenie

Efesto retains shared source profiles, rebuilds selected redox families and brings
in other collections. Changes also include selecting among overlapping entries,
reviewing functional assignments and standardizing runtime metadata. These are different operations: importing a model
or renaming its category does not retrain its HMM.

### Inclusions and exclusions by source

Counts and exclusion reasons below come from the current registry, as in the
[catalog](hmm_library.md#sources). Enabled includes active and experimental
profiles. These counts concern the subset imported into Efesto, not every model
available in the upstream collections.

{{source_decision_table}}

Excluded registry entries remain bundled for provenance and are not searched.
The three glyceraldehyde-dehydrogenase imports were removed entirely and are
absent from these counts. Low seed count and redundancy labels record curation
choices, not measured error rates. In particular, MamG remains excluded by a
historical seed-count rule despite its relevant magnetosome scope; BimA also lies
outside metal scope although its recorded exclusion status is low seed count.

These are not counts of models searched by the original FeGenie workflow.
Compare model stems/accessions and runtime status rather than treating a
difference in totals as newly trained models.

### Source contributions and modifications

| Source label | Contribution and changes in Efesto |
|---|---|
| `fegenie` | Iron-related source profiles were selected and reviewed for sparse training sets, redundancy and misassignment. The archived audit found unchanged scoring parameters for all 211 imported profiles after format conversion; current files match those audited copies. The mislabeled FpvI/FpvH duplicate is excluded. Rebuilt MtrA/MtoA/CymA belong to the separate `curated` group. |
| `methmmdb` | Adds metal acquisition, homeostasis and resistance families beyond iron. Overlapping variants, low-seed entries and category mismatches were excluded; reviewed uptake and homeostasis proteins were reassigned from resistance where justified. Enabled profiles lack positive model-specific cutoffs and use the generic bitscore/E-value gates. |
| `tabuteau` | Adds acquisition and siderophore/metallophore biosynthetic profiles, including KOfam-derived models. Source cutoffs are retained; selected sparse or redundant entries are excluded. CntK/L/M are classified as metallophore synthesis. |
| `mnoxgenetool` | Adds 12 Mn-oxidation profile imports with published external cutoffs. Sparse source subclades include documented exceptions to seed-count exclusions. These were imported, not retrained by Efesto. |
| `ncbifam` | Adds official prokaryotic family profiles, including Fe–S machinery, regulation, stress, lanthanide handling, Mo/W uptake and methanobactin-related families. Source models and thresholds are retained, with model-specific cautions such as TupA coverage. |
| `pfam` | Adds direct domain/family imports, including SufE, LanD and a pilin-related profile. Context and interpretation guidance distinguish broad family detection from metal-specific function. |
| `interpro` | Records an import channel for member-database profiles. Functional assignments were reviewed: the MHB domain belongs to heme handling, and MamH/MamJ to magnetosome formation. Three glyceraldehyde-dehydrogenase profiles were removed from the metal inventory; unrelated competence/virulence entries remain excluded. |
| `curated` | Six redox profiles were built or rebuilt from documented sequence inputs. Five remain enabled; ImoA is excluded. See their individual status and rationale below. |

Exclusion reasons and source support are recorded per model. Low seed count and
sequence similarity are curation heuristics, not measured error rates. The
[source table](../hmm_library/model_licenses.tsv),
[category decisions](../hmm_library/category_curation.tsv) and
[curation guide](hmm_library_curation.md) retain the detailed provenance.

### An updated imported model is not a rebuild

The original FeGenie `PF16525-MHB` file contains **Pfam PF16525.5**. Efesto uses
**PF16525.11**, imported through InterPro, under `Heme-Transport`. It has the
newer source model's parameters and GA of 27.4 bits, versus 27.3 in the old file.
This is a source-version replacement and annotation review, not an Efesto-trained
HMM or evidence of improved accuracy by itself.

### Profiles built or rebuilt for Efesto

MtrA, MtoA and CymA replace the corresponding original FeGenie HMMs in Efesto's
ordinary library. CbcL and DA_402 add profile-based targets; ImoA was built but
subsequently excluded. FeGenie also detects Geobacter/Thermincola reference
proteins through built-in BLAST searches, so these profile additions do not mean
that FeGenie lacked all detection for those organisms or electron-transfer systems.

In the pinned FeGenie metadata, MtrA and MtoA each used a 140-bit full-sequence
cutoff, and CymA used 237 bits. The rebuilt Efesto profiles use the registry
cutoffs below. Their HMM parameters also changed, so the different score scales
must not be interpreted as a direct probability or accuracy gain. The matched
runtime experiment later on this page used the original profiles and thresholds,
not these rebuilt ones.

The following table is generated from registry rows labelled `curated`. NSEQ
counts training records, including any duplicate or closely related sequences;
it is not a count of independent functional experiments. Cutoffs are
full-sequence bitscores.

{{curated_table}}

| Profile | What changed and how to interpret it |
|---|---|
| MtrA | Rebuilt to separate the MtrA-side family from MtoA-related matches, with a manually selected cutoff. Some seeds have broad DmsE-family annotations and still require orthology/context interpretation. |
| MtoA | Rebuilt with a distinct seed set and cutoff for the MtoA-side family. Six training records represent four unique sequences; electron-flow direction is not established by a homolog hit alone. |
| CymA | Rebuilt separately from ImoA. Retained as experimental because its characterized Fe(III)-respiration role coexists with other respiratory functions; it is not an Fe-specific terminal reductase. |
| ImoA | A separate profile was built, then excluded from ordinary runs because the cited experiments establish quinol oxidation/complementation, not direct native Fe(II) oxidation. It is retained for provenance. |
| CbcL | Adds a restricted inner-membrane electron-transfer family linked to low-potential Fe(III) respiration in the reference organism. Eight public inputs and a reconstructed alignment reproduce the existing HMM exactly, excluding the build date; this does not independently calibrate the classifier. |
| DA_402 | Adds a candidate electron-transfer family using corrected sequence provenance. Four same-genome paralogs and one external homolog form the five-record training set; broader functional generalization remains limited. |

Recovered inputs, alignments and reproduction records are distributed under
[`data/seeds/curated/`](https://github.com/l-gallucci/Efesto/tree/main/data/seeds/curated).
Seed recovery and build reproducibility document how a model was generated;
they do not substitute for held-out positive and near-paralog negative testing.

### Categories and thresholds

- Category names were standardized across sources, for example
  `iron_acquisition-iron_transport` → `IronAcquisition-Transport`.
  [Aliases](../hmm_library/category_aliases.tsv) preserve historical spellings.
- Functional reassignment is recorded separately from spelling changes.
  CntK/L/M describe `Metallophore-Synthesis`; the glyceraldehyde-dehydrogenase
  imports were removed rather than given another metal-function label.
- Official imported cutoffs and the locally selected redox cutoffs retain their
  distinct provenance. The registry and runtime cutoff table now preserve the
  same numeric precision.
- A missing or nonpositive model-specific cutoff uses the ordinary generic floor
  of 30 bits; the default E-value gate is 1e-5. These generic settings are not
  individually calibrated GA thresholds. The archived same-HMM control below
  deliberately used different settings to match the pinned FeGenie workflow.

## Annotation and reporting changes

| Area | Current behavior |
|---|---|
| Genomic context | Coordinate-based grouping, optional strand awareness and system-specific rules; distinct systems in a mixed neighborhood are evaluated independently. |
| Retention | HMM gates, applicable sequence/context filters and restricted eggNOG rescue determine which candidates are retained. |
| Evidence | HMM statistics, context, UniOP pair probabilities, BGC overlap and optional structure/topology observations remain separate; the composite confidence score has been removed. |
| Traceability | `Efesto-filter-decisions.tsv` records actual decisions and final retention/rescue outcomes for selected HMM candidates, including exclusions. |
| Resource use | At most two HMM searches overlap by default; a bounded queue and disk-spooled decision audit reduce avoidable memory use without setting a hard RAM ceiling. |
| Input handling | Genome-scoped identifiers, Bakta-aware contig resolution, checked external-tool failures and input-bound caches protect multi-genome runs. |

See [pipeline](pipeline.md), [outputs](output_formats.md),
[evidence examples](evidence.md) and [performance](performance.md) for practical use.

## What the completed comparisons establish

The local experiments below answer different questions. The **matched-library
control** compares the invoked workflows with the same proteins and original
FeGenie HMMs. The **full-library comparison** measures workflows with different
functional scope. Neither is an independent biological accuracy benchmark.

These are archived measurements from **3–4 November 2026**, before subsequent
library and release corrections. They have not been rerun on the current release.
A [downloadable summary](data/comparison_20261104.json) records exact measurements,
input checksums, software provenance and source-snapshot hashes. Full raw outputs
and sampling traces are retained locally and are not included in this documentation.

## Same proteins and HMMs

The control used the same **37,633 Prodigal-called proteins from 10 reference
genomes**, the same HMMER 3.4 executable and **191 distinct searched HMMs** from
unmodified FeGenie commit `8efc74cabfd09680050e35b4deafc5db33e381fb`.
Gene prediction was performed once and excluded from these timings.

| Program | Runs | Wall time (s) | Sampled process-tree peak RSS (MiB) | ORFs annotated through HMMs |
|---|---:|---:|---:|---:|
| Efesto | 1 | 38.0 | 562.8 | 1,304 |
| FeGenie | 1 | 182.7 | 886.4 | 1,304 |

For this paired pilot, the wall-time ratio was **4.82**, and Efesto's sampled
peak RSS was **36.5% lower**. These are observations from one run per program,
not replicated estimates or a general speed/memory guarantee.

### Agreement of the search results

All **1,910 paired HMMER tables** (191 profiles × 10 genomes) had identical
non-comment hit rows after sorting: **3,235 raw hit rows**. Both workflows
annotated the same **1,304 ORFs through HMMs**. Raw-search agreement does not mean
that all downstream labels or annotation methods were identical:

- Seven ORFs had different category labels: six involved the MtrB
  oxidation/reduction category alias, and one involved an equal-bitscore
  PchH/FpvD assignment.
- FeGenie's built-in reference BLAST searches added **12 ORFs**. These are
  additional BLAST results, not HMM matches missing from the paired search tables.

### What was controlled

- Both programs received the same FAA files; Efesto also received matching
  Prodigal GFF coordinates. Both used all-results mode to bypass context filters;
  their other sequence checks, best-hit selection and output logic remained.
- The control reproduced upstream FeGenie's integer-formatted search thresholds.
  Efesto used an E-value limit of `<=0.1`, versus upstream `<0.1`; this boundary
  difference did not change the observed results. These are control settings,
  not Efesto's ordinary default thresholds.
- Both had a nominal four-CPU budget. Efesto ran at most two concurrent searches,
  with HMMER `--cpu 1`; FeGenie ran searches serially with `--cpu 4`.
  This controls the nominal budget, not the scheduling strategy or actual CPU use.
- The original library contained 216 category paths. The upstream metadata lookup
  enabled 192 paths representing 191 distinct searched names. A duplicate MtrB
  profile occurred in two categories; 24 paths lacked cutoff metadata. Thus the
  control concerns the 191 searched models, not every file in the source library.

The measured wall time covers each complete invoked protein-input workflow.
It includes **FeGenie's built-in BLAST searches and deliberate waits** in the
upstream code. Optional external evidence and R plotting were disabled.
The ratio therefore describes these workflows; it does not isolate HMMER speed.

## Full libraries: different scope

An earlier protein-input experiment used the same proteins but each program's
then-enabled library, before Efesto's memory optimization. The three repetitions
per program gave the following medians:

| Program | HMMs searched per genome | Runs | Median wall time (s) | Median sampled peak RSS (MiB) |
|---|---:|---:|---:|---:|
| Efesto | 490 | 3 | 73.3 | 1707.9 |
| FeGenie | 191 | 3 | 187.9 | 884.8 |

These libraries were unequal: Efesto searched 490 profiles, versus 191 in the
pinned FeGenie workflow. Efesto also received GFF coordinates; FeGenie's protein
workflow lacked a contig-name map and used singleton neighborhoods. Consequently,
annotation-count differences combine library scope and context behavior. They
cannot be interpreted as sensitivity or specificity differences.

The historical 490-profile count is not the current library inventory. Use the
[registry-backed catalog](hmm_library.md) for the installed source revision's scope.

## Memory optimization within Efesto

A separate before/after pilot used the same full-library snapshot, the same
FAA/GFF inputs and unchanged command arguments. It tested the change to at most
two overlapping HMM searches, bounded job submission and later sequence loading.

| Efesto configuration | Runs | Wall time (s) | Sampled peak RSS (MiB) | Retained annotations |
|---|---:|---:|---:|---:|
| Before optimization | 1 | 79.0 | 1753.7 | 1,632 |
| Two-search default | 1 | 143.4 | 1235.7 | 1,632 |

Sampled peak RSS decreased by **29.5%**, with a longer wall
time. All 4,900 HMM searches completed without errors. Six checked scientific
outputs were byte-identical between these two configurations: the long TSV,
summary, cluster summary, heatmap matrix, GFF3 and protein FASTA.
This checks equivalence for that optimization and snapshot, not every later
pipeline change. See [current resource behavior](performance.md).

## Measurement environment and limits

The runs used macOS 26.3.1 on an ARM64 Mac15,6 with 11 logical CPUs and 36 GiB
RAM, Python 3.13.13, HMMER 3.4, Prodigal 2.6.3 and BLAST+ 2.16.0.
The Efesto measurements refer to archived working-tree snapshots, not a clean
release identified by the current version number; the summary records hashes.

RSS is the sum for the timed wrapper and its live descendants, sampled using
`ps` approximately every 0.1 seconds plus sampling overhead. Shared pages can be
counted more than once; this is not PSS or exclusive physical memory. Brief peaks
can be missed. One MiB is 1,048,576 bytes. The monitor was outside the timed tree.

Runs were sequential with fresh output directories and no HMM result-cache reuse.
Filesystem caches were not flushed, CPU affinity and laptop load were not
controlled, and the matched-library and memory-optimization pilots had one run
per configuration. The ten reference genomes cover 40,197,383 bp; input identities,
ORF counts and SHA-256 checksums are in the downloadable summary.

No precision, recall or F1 estimate is supported by these timing experiments.
Training overlap and a limited reference panel also prevent treating additional
calls as independently validated discoveries. See [interpretation](interpretation.md).
