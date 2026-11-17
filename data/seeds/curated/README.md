# Curated HMM inputs: provenance audit, updated 2026-11-16

HMM licensing is described in [HMM license scope](../../../hmm_library/LICENSES.md).
The six explicitly listed Efesto model contributions use CC BY-SA 4.0;
this does not replace NCBI/UniProt conditions on the FASTA/MSA input records.
Preserve the source accessions, links and attributions in the manifest.

These inputs reproduce the current MtrA, MtoA, CymA, ImoA, DA_402 and CbcL HMM scoring parameters exactly with MAFFT 7.526 and HMMER 3.4. Dates and textual HMM metadata are excluded from the equality check. See `build_reproduction.tsv` and the accession/sequence SHA256 manifest. This establishes reproducibility, not independent classifier accuracy.

MtrA/MtoA FASTA and MSA files are original Git objects recovered at the commit recorded in `recovered_history_manifest.json`. All 21 sequences match official NCBI records. The other three MSAs were reconstructed with `mafft --auto --quiet` from the documented official sequences; their builds match the deployed matrices, checksums, lengths, effective counts and statistical parameters.

Reproduce a model into a separate output directory:

```sh
hmmbuild --cpu 1 -n MtrA /tmp/MtrA.reproduced.hmm MtrA_seeds.aln
```

Use the corresponding name and alignment for MtoA, CymA, ImoA, DA_402 or CbcL. Preserve FASTA order for reconstructing an alignment. Do not replace deployed models or choose cutoffs from this check alone.

| Model | Training records | Unique sequences | Provenance/function limit |
|---|---:|---:|---|
| MtrA | 15 | 14 | Six records have only a DmsE-family annotation; specific orthology needs genomic context/phylogeny |
| MtoA | 6 | 4 | Additional homologs have inferred, not individually assayed, functions |
| CymA | 5 | 5 | Respiratory hub supporting several acceptors, not exclusively iron |
| ImoA | 2 | 2 | Quinol oxidation demonstrated; pathway context does not prove direct Fe(II) oxidation |
| DA_402 | 5 | 5 | Four paralogs from one genome plus one unreviewed external homolog; candidate EET family |
| CbcL | 8 | 8 | Reconstructed from traceable public records; exact deployed-build reproduction, not independent classifier validation |

## CbcL reconstruction

`CbcL_seeds.faa` and `CbcL_seeds.aln` now contain all eight input sequences in
the order that reproduces the deployed model. `CbcL_source_records.json` stores
the official UniProt records and NCBI cross-references; every sequence was
independently verified byte-identical between those providers on 2026-11-16.
`CbcL_reconstruction.json` records commands, versions, file hashes and the result.

The historical input files were not found. The reconstruction used the old
genus/length description to retrieve public records, compared 792 eight-sequence
sets from 13 unique candidates against the existing build, and recovered input
order using the Easel alignment checksum. A fresh `mafft --thread 1 --auto --quiet`
alignment followed by HMMER 3.4 `hmmbuild --cpu 1 -n CbcL` reproduces the entire
current HMM apart from `DATE`: matrix, transitions, map, statistics, effective
sequence count and checksum **3842520616** all match. The deployed HMM and cutoff
**200 bits** have not been replaced or recalibrated.

The ordered UniProt accessions are Q74GH2, A0A0M4CZW0, A0ABM8EG83, B5E8E3,
A0A0C2DW34, A0A0B5FL41, Q39ZH5 and B9M704. All eight have nine CXXCH motifs.
Some provider records use the name `cbcY`; these annotations are preserved as
source metadata and are not treated as experimental evidence for every homolog.

`CbcL_verified_reference_only.faa` remains the original single-anchor convenience
file. It is not the training set. Q74GH2 is identical to AAR33608.1/GSU0274, linked
explicitly to CbcL in [Antunes et al. (2022)](https://doi.org/10.3389/fmicb.2022.898015).
It scores 1331.3 bits against the unchanged HMM. This anchor has experimental
support; the other seven are inferred homologs, not individually assayed in this
audit. Input reconstruction is not an independent sensitivity/specificity test.

CbcL is an inner-membrane electron-transfer component associated with respiration
of low-potential extracellular acceptors, including Fe(III) oxides, in
*G. sulfurreducens*. It is not a terminal Fe(III) reductase. A homolog hit alone
does not establish iron reduction, and absence does not exclude other EET routes.

Duplicate groups: MtrA WP_011071900.1/AAN54830.1; MtoA ADE12722.1/WP_013030620.1 and ADL56010.1/WP_013293942.1. They are retained to reproduce the deployed models; deduplication requires retraining and threshold reassessment. The manifest distinguishes experimentally linked anchors from inferred homologs. Source papers are family/context evidence, not per-accession experimental validation for every row.

Canonical references: [MtrA](https://doi.org/10.1074/jbc.M302582200), [MtoA](https://doi.org/10.3389/fmicb.2012.00037), [CymA](https://doi.org/10.1128/jb.179.4.1143-1152.1997), [ImoA](https://doi.org/10.1128/mbio.02150-22), [DA_402](https://doi.org/10.1038/s41586-025-09467-0). DA_402 is linked to WP_013162640.1 by the Nature Source Data coordinate table (Extended Figure 3, B16:H16), independently matching NC_014216.1:468410–469705, minus strand. Expression and structural evidence do not demonstrate the function of every homolog.
