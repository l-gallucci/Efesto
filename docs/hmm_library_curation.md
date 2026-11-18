# Curation and provenance

Efesto organizes HMM profiles by biological role and retains their source,
references, search thresholds and runtime status. The [catalog](hmm_library.md)
lists the profiles available in each source revision; the
[source comparison](comparison.md) summarizes contributions from each collection.

## Library records

| File | Contents |
|---|---|
| `hmm_registry.tsv` | Model identifiers, sources, references, status and metadata |
| `category_curation.tsv` | Category assignments and supporting references |
| `category_aliases.tsv` | Alternative and previous category spellings |
| `HMM-bitcutoffs.txt` | Runtime full-sequence thresholds |
| `FeGenie-map.txt` | Readable labels indexed by model filename stem |
| `data/seeds/curated/` | Sequences, alignments and build information for Efesto-built profiles |
| `data/validation/` | Reference sequences and expected search behavior |

These files are available in the [source repository](https://github.com/l-gallucci/Efesto).

## Cutoff provenance

All runtime cutoffs are **full-sequence bitscores**. Their origin depends on
the model source:

| Profile source | Threshold used by Efesto |
|---|---|
| Official Pfam and NCBIfam/TIGRFAM profiles, including enabled member-database imports through InterPro | Full-sequence GA from the source HMM, retained without de novo estimation |
| Other published model collections | Model-specific cutoff from the source metadata, where available |
| Efesto-built profiles | Efesto-selected operating threshold recorded in the registry and runtime cutoff table |
| No positive model-specific cutoff | Generic bitscore floor, 30 bits by default |

An HMM can contain two GA values: one for the complete sequence and one for
an individual domain. Efesto's cutoff table uses the first. For example, SufD
TIGR01981.1 specifies 222.0 bits for the sequence and 221.5 bits for a domain;
Efesto uses **222.0 bits**. The E-value gate also applies to accepted matches.

## Construction of Efesto-built profiles

MtrA, MtoA, CymA, ImoA, CbcL and DA_402 use documented sequence sets from
NCBI/UniProt, combining literature-linked reference proteins with related homologs:

1. **Sequence selection:** each family has a separate FASTA set with accessions,
   functional references and sequence provenance.
2. **Alignment:** the sequences are aligned with MAFFT using `--auto`.
3. **Profile construction:** HMMER `hmmbuild` converts each alignment into a
   profile HMM under the corresponding model name.

The supplied build recipe uses **MAFFT 7.526 and HMMER 3.4**. FASTA files,
alignments, training-set composition, references and commands are available in
[`data/seeds/curated/`](https://github.com/l-gallucci/Efesto/tree/main/data/seeds/curated).
MtrA/MtoA and CymA/ImoA use separate family-specific training sets.
ImoA is excluded from ordinary searches.

The operating cutoffs for these profiles are selected separately from HMM
construction and stored in the registry and runtime cutoff table. They are
**Efesto-selected thresholds**, not GA values inherited from Pfam/NCBIfam or
automatically produced by `hmmbuild`. Training records include inferred homologs;
their number does not represent independent experimental validation.

## Model scope and status

| Model or family | Scope in Efesto |
|---|---|
| CbcL | Inner-membrane electron-transfer component associated with low-potential Fe(III) respiration in *Geobacter sulfurreducens*. A homolog match requires genomic context to support a physiological interpretation. |
| MtrA / MtoA | Related multiheme cytochrome families. Interpret assignments alongside partner proteins and genomic context; sequence similarity alone does not determine electron-flow direction. |
| CymA | Experimental profile for a multifunctional respiratory component involved in several electron-acceptor pathways, including Fe(III) respiration. |
| ImoA | Excluded from ordinary runs. Its characterized quinol-oxidation role does not establish native Fe(II) oxidation. |
| DA_402 | Candidate electron-transfer family with a provisional functional scope. |
| CntK/L/M | Opine-type metallophore biosynthesis. Normal context filtering requires CntL+CntM; CntK is optional across products. |
| LanM/LanD, Wtp/Tup and MbnA/B/C | Imported family profiles with source-specific scope and thresholds; see [lanthanides, tungsten and metallophores](metal_expansion.md). |

Profiles outside the metal-related scope and redundant entries are excluded from
ordinary searches. These include the glyceraldehyde-dehydrogenase imports
NF041018/19/20 and a mislabeled FpvI copy of FpvH. The registry records each
model's status and reason for exclusion.

## Related families and domains

Broad flavodoxin and pilin profiles describe protein families rather than a
specific metal function. Likewise, a SufE-related domain can occur in CsdE and
other sulfur-transfer proteins. Use partner genes, architecture and functional
references to distinguish related proteins.

Overlapping HMMs can represent different paralogs or biological roles. Model
accessions, training sequences, reference matches and genomic context help
determine whether profiles provide complementary coverage. Similar consensus
sequences alone do not establish functional equivalence.

See [Fe-S assembly](iron_sulfur_assembly.md) and
[interpreting results](interpretation.md) for practical guidance.
