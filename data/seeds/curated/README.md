# Curated HMM inputs: provenance audit, 2026-11-05

These inputs reproduce the current MtrA, MtoA, CymA, ImoA and DA_402 HMM scoring parameters exactly with MAFFT 7.526 and HMMER 3.4. Dates and textual HMM metadata are excluded from the equality check. See `build_reproduction.tsv` and the accession/sequence SHA256 manifest. This establishes reproducibility, not independent classifier accuracy.

MtrA/MtoA FASTA and MSA files are original Git objects recovered at the commit recorded in `recovered_history_manifest.json`. All 21 sequences match official NCBI records. The other three MSAs were reconstructed with `mafft --auto --quiet` from the documented official sequences; their builds match the deployed matrices, checksums, lengths, effective counts and statistical parameters.

Reproduce a model into a separate output directory:

```sh
hmmbuild --cpu 1 -n MtrA /tmp/MtrA.reproduced.hmm MtrA_seeds.aln
```

Use the corresponding name and alignment for MtoA, CymA, ImoA or DA_402. Preserve FASTA order for reconstructing an alignment. Do not replace deployed models or choose cutoffs from this check alone.

| Model | Training records | Unique sequences | Provenance/function limit |
|---|---:|---:|---|
| MtrA | 15 | 14 | Six records have only a DmsE-family annotation; specific orthology needs genomic context/phylogeny |
| MtoA | 6 | 4 | Additional homologs have inferred, not individually assayed, functions |
| CymA | 5 | 5 | Respiratory hub supporting several acceptors, not exclusively iron |
| ImoA | 2 | 2 | Quinol oxidation demonstrated; pathway context does not prove direct Fe(II) oxidation |
| DA_402 | 5 | 5 | Four paralogs from one genome plus one unreviewed external homolog; candidate EET family |
| CbcL | 8 reported in HMM | Unknown | One verified reference only; seven original accessions and original MSA missing |

`CbcL_verified_reference_only.faa` is deliberately not named as a complete seed set. Q74GH2 is identical to AAR33608.1/GSU0274. Antunes et al. explicitly link these identifiers in their cloning methods: [10.3389/fmicb.2022.898015](https://doi.org/10.3389/fmicb.2022.898015). The current HMM recognizes this reference at 1331.3 bits (cutoff 200).

Duplicate groups: MtrA WP_011071900.1/AAN54830.1; MtoA ADE12722.1/WP_013030620.1 and ADL56010.1/WP_013293942.1. They are retained to reproduce the deployed models; deduplication requires retraining and threshold reassessment. The manifest distinguishes experimentally linked anchors from inferred homologs. Source papers are family/context evidence, not per-accession experimental validation for every row.

Canonical references: [MtrA](https://doi.org/10.1074/jbc.M302582200), [MtoA](https://doi.org/10.3389/fmicb.2012.00037), [CymA](https://doi.org/10.1128/jb.179.4.1143-1152.1997), [ImoA](https://doi.org/10.1128/mbio.02150-22), [DA_402](https://doi.org/10.1038/s41586-025-09467-0). DA_402 is linked to WP_013162640.1 by the Nature Source Data coordinate table (Extended Figure 3, B16:H16), independently matching NC_014216.1:468410–469705, minus strand. Expression and structural evidence do not demonstrate the function of every homolog.
