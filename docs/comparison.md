# Efesto and FeGenie: scope, model changes and comparisons

Efesto integrates multiple profile-HMM collections to annotate metal cycling,
acquisition, homeostasis and resistance. FeGenie contributes iron profiles and
the basis of the default context engine. MetHMMDB, MnOxGeneTool, Tabuteau et al.,
public protein-family databases and Efesto-built profiles contribute other
families and functional layers.

This page describes the functional coverage, model sources and workflow
differences between Efesto and FeGenie, followed by performance measurements.
Cite the software version using [CITATION.cff](../CITATION.cff), together with the
[original model and tool sources](references.md) used in the analysis.

## Current library at a glance

{{catalog_summary}}

Runtime status determines which profiles are searched. Source labels distinguish
import channels from the original family databases; model-specific scope is
listed in the [catalog](hmm_library.md).

## Functional scope beyond iron

The combined library includes manganese, copper, zinc, nickel, cobalt, cadmium,
chromium, mercury, silver, tungsten, lanthanide and metal/metalloid-associated
families. Acquisition, homeostasis, regulation, redox and resistance remain
separate functional categories. The following table gives examples of enabled profile families.

| Scope | Example profile families | Source contribution |
|---|---|---|
| Iron redox, acquisition, heme, storage and magnetosomes | Cyc2, Mtr/Mto, heme-handling and siderophore families | FeGenie, Tabuteau et al., public family databases and Efesto-built profiles |
| Manganese oxidation and acquisition | BoxA, MnxG, MopA; MntA and MntH | MnOxGeneTool, MetHMMDB and NCBI protein family models |
| Copper homeostasis and resistance | CopZ, copper-transporting ATPases, CueO | MetHMMDB |
| Other metal and metalloid functions | ZitB, cadmium-associated ATPases, mercuric reductase, ArsB/Acr3; Ni/Cr/Ag/Te-associated families | MetHMMDB, with role-specific categories |
| Lanthanide binding and handling | LanM, LanD | NCBI protein family models and Pfam |
| Tungsten and shared Mo/W acquisition | TupA/B/C, WtpA/B/C | NCBI protein family models, including TIGRFAM-derived profiles |
| Metallophore biosynthesis | CntK/L/M and MbnA/B/C | KOfam through Tabuteau et al.; NCBI protein family models |
| Fe–S cofactor biogenesis | SufA/B/C/D/E/S and IscS | NCBI protein family models and Pfam |

The eleven LanM/LanD, Wtp/Tup and MbnA/B/C profiles retain the official providers'
HMM parameters and GA cutoffs. Counts and runtime status are listed in the
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
absent from these counts. Exclusions include sparse training sets, overlapping
variants and functions outside the library scope. Per-profile status and category
records provide the detailed selection criteria.

These are not counts of models searched by the original FeGenie workflow.
Compare model stems/accessions and runtime status rather than treating a
difference in totals as newly trained models.

### Source contributions and modifications

| Source label | Contribution and changes in Efesto |
|---|---|
| `fegenie` | Iron-related profiles retain their imported scoring parameters. Selected sparse, overlapping and out-of-scope entries are excluded, including the mislabeled FpvI/FpvH duplicate. Efesto-built MtrA/MtoA/CymA profiles use the separate `curated` label. |
| `methmmdb` | Adds metal acquisition, homeostasis and resistance families beyond iron. Overlapping variants, low-seed entries and category mismatches were excluded; uptake and homeostasis families use categories distinct from resistance. Enabled profiles lack positive model-specific cutoffs and use the generic bitscore/E-value gates. |
| `tabuteau` | Adds acquisition and siderophore/metallophore biosynthetic profiles, including KOfam-derived models. Source cutoffs are retained; selected sparse or redundant entries are excluded. CntK/L/M are classified as metallophore synthesis. |
| `mnoxgenetool` | Adds 12 Mn-oxidation profile imports with published external cutoffs. Sparse source subclades include documented exceptions to seed-count exclusions. These were imported, not retrained by Efesto. |
| `ncbifam` | Adds official prokaryotic family profiles, including Fe–S machinery, regulation, stress, lanthanide handling, Mo/W uptake and methanobactin-related families. Source models and thresholds are retained, with model-specific cautions such as TupA coverage. |
| `pfam` | Adds direct domain/family imports, including SufE, LanD and a pilin-related profile. Context and interpretation guidance distinguish broad family detection from metal-specific function. |
| `interpro` | Records an import channel for member-database profiles. The MHB domain is assigned to heme handling, and MamH/MamJ to magnetosome formation. Three glyceraldehyde-dehydrogenase profiles were removed from the metal inventory; unrelated competence/virulence entries remain excluded. |
| `curated` | Six redox profiles have documented sequence inputs. Five are enabled; ImoA is excluded. Their functional scope is summarized below. |

Exclusion reasons and references are recorded per model in the
[source table](../hmm_library/model_licenses.tsv),
[category decisions](../hmm_library/category_curation.tsv) and
[curation guide](hmm_library_curation.md).

### Imported model versions

The original FeGenie `PF16525-MHB` file contains **Pfam PF16525.5**. Efesto uses
**PF16525.11**, imported through InterPro, under `Heme-Transport`. It has the
newer source model's parameters and GA of 27.4 bits, versus 27.3 in the old file.
Both accessions identify the same Pfam family; `.5` and `.11` denote model
revisions. Efesto imports revision 11 without retraining the HMM.

### Profiles built or rebuilt for Efesto

MtrA, MtoA and CymA replace the corresponding original FeGenie HMMs in Efesto's
ordinary library. CbcL and DA_402 add profile-based targets; ImoA is excluded
from ordinary searches. FeGenie also detects Geobacter/Thermincola reference
proteins through built-in BLAST searches, so these profile additions do not mean
that FeGenie lacked all detection for those organisms or electron-transfer systems.

In the FeGenie revision used for comparison, MtrA and MtoA each used a 140-bit full-sequence
cutoff, and CymA used 237 bits. The rebuilt Efesto profiles use the registry
cutoffs below. Their HMM parameters also changed, so the different score scales
must not be interpreted as a direct probability or accuracy gain. The matched
runtime experiment later on this page used the original profiles and thresholds,
not these rebuilt ones.

The following table is generated from registry rows labelled `curated`. NSEQ
counts training records, including any duplicate or closely related sequences;
it is not a count of independent functional experiments. Cutoffs are
full-sequence bitscores selected for Efesto, rather than inherited source GA
values.

{{curated_table}}

These profiles use separate, documented NCBI/UniProt sequence sets, aligned with
MAFFT (`--auto`) and converted into HMMs with HMMER `hmmbuild`. Threshold
selection is separate from profile construction; `hmmbuild` does not generate
these operating cutoffs. See [construction methods and cutoff provenance](hmm_library_curation.md#construction-of-efesto-built-profiles)
for the build recipe and training-set information.

| Profile | What changed and how to interpret it |
|---|---|
| MtrA | Rebuilt to separate the MtrA-side family from MtoA-related matches, with a manually selected cutoff. The family includes DmsE-related cytochromes; partner proteins and genomic context guide interpretation. |
| MtoA | Rebuilt with a distinct seed set and cutoff for the MtoA-side family. Six training records represent four unique sequences; electron-flow direction is not established by a homolog hit alone. |
| CymA | Rebuilt separately from ImoA. Retained as experimental because its characterized Fe(III)-respiration role coexists with other respiratory functions; it is not an Fe-specific terminal reductase. |
| ImoA | Excluded from ordinary runs. Its characterized quinol-oxidation role does not establish native Fe(II) oxidation. |
| CbcL | Adds an inner-membrane electron-transfer family associated with low-potential Fe(III) respiration in *G. sulfurreducens*. A homolog match alone does not establish iron reduction. |
| DA_402 | Adds a candidate electron-transfer family. The training set contains five records, including four same-genome paralogs; its functional scope remains provisional. |

Sequences, alignments, functional references and build instructions are available
in [`data/seeds/curated/`](https://github.com/l-gallucci/Efesto/tree/main/data/seeds/curated).

### Categories and thresholds

- Category names were standardized across sources, for example
  `iron_acquisition-iron_transport` → `IronAcquisition-Transport`.
  [Aliases](../hmm_library/category_aliases.tsv) support previous spellings.
- Functional reassignment is recorded separately from spelling changes.
  CntK/L/M describe `Metallophore-Synthesis`; the glyceraldehyde-dehydrogenase
  imports were removed rather than given another metal-function label.
- Official imported cutoffs and the Efesto-selected redox cutoffs retain their
  distinct provenance. The registry and runtime cutoff table preserve the
  same numeric precision.
- A missing or nonpositive model-specific cutoff uses the ordinary generic floor
  of 30 bits; the default E-value gate is 1e-5. These generic settings are not
  individually calibrated GA thresholds. The matched-library benchmark below
  uses FeGenie-compatible thresholds to compare the two workflows.

## Annotation and reporting changes

| Area | Current behavior |
|---|---|
| Genomic context | Coordinate-based grouping, optional strand awareness and system-specific rules; distinct systems in a mixed neighborhood are evaluated independently. |
| Retention | HMM gates, applicable sequence/context filters and restricted eggNOG rescue determine which candidates are retained. |
| Evidence | HMM statistics, context, UniOP pair probabilities, BGC overlap and optional structure/topology observations are reported as separate observations. |
| Traceability | `Efesto-filter-decisions.tsv` records actual decisions and final retention/rescue outcomes for selected HMM candidates, including exclusions. |
| Resource use | At most two HMM searches overlap by default; a bounded queue and disk-spooled decision records reduce avoidable memory use without setting a hard RAM ceiling. |
| Input handling | Genome-scoped identifiers, Bakta-aware contig resolution, external-tool error handling and input-bound caches protect multi-genome runs. |

See [pipeline](pipeline.md), [outputs](output_formats.md),
[evidence examples](evidence.md) and [performance](performance.md) for practical use.

## Performance comparisons

The **matched-library control** uses identical proteins and original FeGenie
HMMs. The **full-library comparison** uses each program's enabled library, with
different functional scope. Both comparisons measure runtime and memory, not
annotation accuracy.

**Benchmark configuration: 3–4 November 2026, before v1.0.0.** Results apply
to the software revisions and library configurations described below.
The [benchmark data](data/comparison_20261104.json) provide measurements, input
checksums and software identifiers.

## Matched-library benchmark

The control used the same **37,633 Prodigal-called proteins from 10 reference
genomes**, the same HMMER 3.4 executable and **191 distinct searched HMMs** from
unmodified FeGenie commit `8efc74cabfd09680050e35b4deafc5db33e381fb`.
Gene prediction was performed once and excluded from these timings.

| Program | Runs | Wall time (s) | Sampled process-tree peak RSS (MiB) | ORFs annotated through HMMs |
|---|---:|---:|---:|---:|
| Efesto | 1 | 38.0 | 562.8 | 1,304 |
| FeGenie | 1 | 182.7 | 886.4 | 1,304 |

In this configuration, Efesto completed the workflow **4.82 times faster**,
with **36.5% lower** sampled peak RSS. Each program was measured once; repeated
runs would be needed to estimate variability.

### Agreement of the search results

All **1,910 paired HMMER tables** (191 profiles × 10 genomes) contained identical
non-comment hit rows after sorting: **3,235 raw hit rows**, corresponding to
**1,304 ORFs annotated through HMMs** in each workflow. Downstream differences
concerned category assignment and the additional BLAST stage:

- Seven ORFs had different category labels: six involved the MtrB
  oxidation/reduction category alias, and one involved an equal-bitscore
  PchH/FpvD assignment.
- FeGenie's built-in reference BLAST searches added **12 ORFs**. These are
  results from the BLAST stage, separate from the shared HMM matches.

### Benchmark setup

- Both programs received the same FAA files; Efesto also received matching
  Prodigal GFF coordinates. Both used all-results mode to bypass context filters;
  their other sequence checks, best-hit selection and output logic remained.
- The control reproduced upstream FeGenie's integer-formatted search thresholds.
  Efesto used an E-value limit of `<=0.1`, versus upstream `<0.1`; this boundary
  difference did not change the observed results. These are control settings,
  not Efesto's ordinary default thresholds.
- Both had a nominal four-CPU budget. Efesto ran at most two concurrent searches,
  with HMMER `--cpu 1`; FeGenie ran searches serially with `--cpu 4`.
  Search scheduling and CPU allocation differed between workflows.
- The original library contained 216 category paths. The upstream metadata lookup
  enabled 192 paths representing 191 distinct searched names. A duplicate MtrB
  profile occurred in two categories; 24 paths lacked cutoff metadata. The
  comparison therefore covers the 191 models searched by the FeGenie workflow.

The measured wall time covers each complete invoked protein-input workflow.
It includes **FeGenie's built-in BLAST searches and deliberate waits** in the
upstream code. Optional external evidence and R plotting were disabled.
The ratio therefore describes these workflows; it does not isolate HMMER speed.

## Full libraries: different scope

The full-library benchmark used the same proteins with **490 Efesto profiles**
and **191 FeGenie profiles**, with three runs per program. The Efesto configuration
preceded the two-search concurrency limit. Median measurements were:

| Program | HMMs searched per genome | Runs | Median wall time (s) | Median sampled peak RSS (MiB) |
|---|---:|---:|---:|---:|
| Efesto | 490 | 3 | 73.3 | 1707.9 |
| FeGenie | 191 | 3 | 187.9 | 884.8 |

These libraries were unequal: Efesto searched 490 profiles, versus 191 in the
FeGenie workflow. Efesto also received GFF coordinates; FeGenie's protein
workflow lacked a contig-name map and used singleton neighborhoods. Consequently,
annotation-count differences combine library scope and context behavior. They
cannot be interpreted as sensitivity or specificity differences.

The [library catalog](hmm_library.md) lists the profiles enabled in the current
version of Efesto.

## Memory optimization within Efesto

The memory benchmark compared two search-scheduling configurations using the
same 490-profile library, FAA/GFF inputs and command arguments. The revised
configuration allowed at most two simultaneous HMM searches, bounded job
submission and loaded protein sequences after the searches.

| Efesto configuration | Runs | Wall time (s) | Sampled peak RSS (MiB) | Retained annotations |
|---|---:|---:|---:|---:|
| Before optimization | 1 | 79.0 | 1753.7 | 1,632 |
| Two-search default | 1 | 143.4 | 1235.7 | 1,632 |

Sampled peak RSS decreased by **29.5%**, while wall time increased from
79.0 to 143.4 seconds. All 4,900 HMM searches completed successfully. The long TSV,
summary, cluster summary, heatmap matrix, GFF3 and protein FASTA were byte-identical
between configurations. See [performance and memory](performance.md) for the
current scheduling behavior and configuration guidance.

## Measurement environment and limits

The runs used macOS 26.3.1 on an ARM64 Mac15,6 with 11 logical CPUs and 36 GiB
RAM, Python 3.13.13, HMMER 3.4, Prodigal 2.6.3 and BLAST+ 2.16.0.
The benchmark data identify the Efesto source revisions by content hashes
and FeGenie by its commit identifier.

RSS is the sum for the timed wrapper and its live descendants, sampled using
`ps` approximately every 0.1 seconds plus sampling overhead. Shared pages can be
counted more than once; this is not PSS or exclusive physical memory. Brief peaks
can be missed. One MiB is 1,048,576 bytes. The monitor was outside the timed tree.

Runs were sequential with fresh output directories and no HMM result-cache reuse.
Filesystem caches were not flushed, CPU affinity and laptop load were not
controlled, and the matched-library and memory benchmarks had one run
per configuration. The ten reference genomes cover 40,197,383 bp; input identities,
ORF counts and SHA-256 checksums are in the downloadable summary.

These timing experiments do not estimate precision, recall or F1. Additional
annotations from a larger library require their own functional interpretation;
see [interpreting results](interpretation.md).
