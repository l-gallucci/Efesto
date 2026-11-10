# References and citations

Cite the sources applicable to the library and optional tools used in your run.
The runtime provenance table and registry identify the relevant model sources.


If you use Efesto, cite all of the following that apply:

**FeGenie** (iron cycling HMMs and operon logic):
> Garber AI et al. (2020) *FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods.* Front. Microbiol. 11:37. [doi:10.3389/fmicb.2020.00037](https://doi.org/10.3389/fmicb.2020.00037)

**MetHMMDB** (metal resistance gene HMMs):
> Ciuchcinski K, Dziurzynski M (2025) *Fast and accurate detection of metal resistance genes using MetHMMDb.* bioRxiv. [doi:10.1101/2024.12.26.629440](https://doi.org/10.1101/2024.12.26.629440)

**Tabuteau et al.** (additional iron acquisition HMMs):
> Tabuteau S et al. (2025) *Metagenomic profiling and genome-centric analysis reveal iron acquisition systems in cheese-associated bacteria and fungi.* Environmental Microbiology 27(12):e70218. [doi:10.1111/1462-2920.70218](https://doi.org/10.1111/1462-2920.70218)

**KOfam** (KEGG Ortholog HMMs used in Tabuteau et al.):
> Aramaki T et al. (2020) *KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold.* Bioinformatics 36:2251–2252. [doi:10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859)

**NCBI Protein Family Models** (NF* HMMs):
> Li W et al. (2021) *RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation.* Nucleic Acids Research 49:D1020–D1028. [doi:10.1093/nar/gkaa1105](https://doi.org/10.1093/nar/gkaa1105)

**MnOxGeneTool** (manganese oxidation HMMs):
> Wang Y et al. (2025) *MnOxGeneTool: a tool for identifying manganese oxidation genes in (meta)genomes.* Environmental Science & Technology. [doi:10.1021/acs.est.5c01235](https://doi.org/10.1021/acs.est.5c01235)

**UniOP** (operon prediction, if using `--operon_prediction`):
> Su H, Zhang R, Söding J (2024) *UniOP: a universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance.* bioRxiv. [doi:10.1101/2024.11.11.623000](https://doi.org/10.1101/2024.11.11.623000)

**eggNOG-mapper** (tier-2 confirmation for `needs_confirmation`-flagged HMM hits, if using `--eggnog_annotations` or `--run_eggnog`):
> Cantalapiedra CP et al. (2021) *eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale.* Mol. Biol. Evol. 38(12):5825–5829. [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293)

**Phobius** (topology/motif consistency check, if using `--topology_check`):
> Käll L, Krogh A, Sonnhammer ELL (2004) *A Combined Transmembrane Topology and Signal Peptide Prediction Method.* J. Mol. Biol. 338(5):1027–1036. [doi:10.1016/j.jmb.2004.03.016](https://doi.org/10.1016/j.jmb.2004.03.016)

**Foldseek + ProstT5** (tier-3 structural confirmation, if using `--run_struct_confirmation`):
> van Kempen M et al. (2023) *Fast and accurate protein structure search with Foldseek.* Nat. Biotechnol. [doi:10.1038/s41587-023-01773-0](https://doi.org/10.1038/s41587-023-01773-0)
> Heinzinger M et al. (2024) *Bilingual language model for protein sequence and structure.* NAR Genomics and Bioinformatics. [doi:10.1093/nargab/lqae150](https://doi.org/10.1093/nargab/lqae150)


## Development methods

- [Mankowski et al., supplied carbon-fixation preprint](https://doi.org/10.64898/2026.09.02.748891).
- [Tanabe and Dahl, HMS-S-S](https://doi.org/10.1111/1755-0998.13642).

An imported profile is not a newly trained Efesto model. Experimental model status
and registry source support do not imply an independent classifier benchmark.
