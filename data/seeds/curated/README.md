# Sequence inputs for Efesto-built profiles

This directory contains FASTA sequences, alignments and source records for
MtrA, MtoA, CymA, ImoA, DA_402 and CbcL. Model status and search cutoffs are
listed in the [registry](../../../hmm_library/hmm_registry.tsv).

## Files

| File | Contents |
|---|---|
| `<model>_seeds.faa` | Ordered training sequences |
| `<model>_seeds.aln` | Multiple sequence alignment used for HMM construction |
| `curated_sequence_manifest.tsv` | Accessions, sequence hashes and provenance |
| `build_reproduction.tsv` | MAFFT/HMMER versions and scoring-parameter comparison |
| `recovered_history_manifest.json` | Source Git objects for available original input files |
| `CbcL_source_records.json` | UniProt records and NCBI cross-references for the CbcL inputs |
| `CbcL_reconstruction.json` | CbcL input provenance, build settings and hashes |

For CbcL, `CbcL_reconstruction.json` documents the derivation of the supplied
sequence set and alignment from public records.
`CbcL_verified_reference_only.faa` contains the single characterized reference
Q74GH2/GSU0274; the full training set is `CbcL_seeds.faa`.

## Building a profile

The documented environment uses MAFFT 7.526 and HMMER 3.4. Preserve sequence
order when regenerating an alignment:

```sh
mafft --thread 1 --auto --quiet MtrA_seeds.faa > /tmp/MtrA.aln
hmmbuild --cpu 1 -n MtrA /tmp/MtrA.hmm /tmp/MtrA.aln
```

To build directly from the supplied alignment:

```sh
hmmbuild --cpu 1 -n MtrA /tmp/MtrA.hmm MtrA_seeds.aln
```

Use the corresponding model name and files for other profiles.
`build_reproduction.tsv` compares scoring parameters, excluding build dates and
textual metadata. A new model or alignment requires its own threshold evaluation.

## Operating thresholds

The cutoffs for these six profiles are Efesto-selected full-sequence
operating thresholds. They are stored in `hmm_registry.tsv` and
`HMM-bitcutoffs.txt`, rather than inherited GA values from an external model.
`hmmbuild` constructs the profile; it does not select these biological acceptance
thresholds. Regenerating an existing HMM does not recalibrate its cutoff.

## Training-set composition and scope

| Model | Training records | Unique sequences | Scope |
|---|---:|---:|---|
| MtrA | 15 | 14 | MtrA/DmsE-related cytochromes; partner proteins and context guide orthology assignment |
| MtoA | 6 | 4 | MtoA-related cytochromes; individual homologs need not share a demonstrated electron-flow direction |
| CymA | 5 | 5 | Multifunctional respiratory electron transfer, including Fe(III) respiration |
| ImoA | 2 | 2 | Quinol oxidation; excluded from ordinary metal annotation |
| DA_402 | 5 | 5 | Candidate electron-transfer family; four records are same-genome paralogs |
| CbcL | 8 | 8 | Inner-membrane electron transfer associated with low-potential extracellular respiration |

Training records include experimentally characterized anchors and inferred
homologs. Their number does not represent independent functional experiments.
CbcL's characterized reference is associated with low-potential Fe(III)
respiration in *G. sulfurreducens*, rather than terminal Fe(III) reduction.
The manifest retains provider names such as `cbcY` for related sequences.

Repeated records are retained where they form part of the distributed model:
MtrA WP_011071900.1/AAN54830.1; MtoA ADE12722.1/WP_013030620.1 and
ADL56010.1/WP_013293942.1. Deduplication changes the training input and requires
rebuilding and threshold assessment.

## Functional references

- [MtrA](https://doi.org/10.1074/jbc.M302582200)
- [MtoA](https://doi.org/10.3389/fmicb.2012.00037)
- [CymA](https://doi.org/10.1128/jb.179.4.1143-1152.1997)
- [ImoA](https://doi.org/10.1128/mbio.02150-22)
- [DA_402](https://doi.org/10.1038/s41586-025-09467-0)
- [CbcL](https://doi.org/10.3389/fmicb.2022.898015)

## Source terms

The six explicitly listed Efesto model contributions use CC BY-SA 4.0; see
[HMM license scope](../../../hmm_library/LICENSES.md). Input sequence records
retain their NCBI/UniProt conditions and attributions. Preserve source accessions,
links and attribution metadata when reusing them.
