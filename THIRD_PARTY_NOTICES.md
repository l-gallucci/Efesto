# Third-party sources and conditions of use

Efesto is developed for academic, non-commercial research. Cite Efesto and the
original publications/databases relevant to the profiles and tools used in your
analysis; the [references](docs/references.md) and model registry identify these
sources.

The MIT license in `LICENSE` applies to Efesto's own code. Imported HMM profiles,
reference sequences, databases, and separately installed tools remain subject to
their original terms. Their inclusion and citation do not grant additional rights
or replace any permission required by the source.

Before non-academic or commercial use, consult **each original source**, check its
current license and conditions, and obtain permission where required. Academic
users must also comply with applicable source terms. In particular,
[MnOxGeneTool](https://github.com/wyh19990121/MnOxGeneTool) states that it is for
non-commercial academic use and directs commercial users to contact its authors.
This notice does not establish that every bundled third-party asset may be
redistributed under MIT.

The package metadata uses `LicenseRef-ThirdParty-Models` to refer to the
source-specific terms described here and in the original distributions; it is
not a new license or a grant of additional rights.

## Find the original sources

| Source | Original reference or provider |
|---|---|
| FeGenie | [Garber et al. (2020)](https://doi.org/10.3389/fmicb.2020.00037) |
| MetHMMDB | [Ciuchcinski and Dziurzynski](https://doi.org/10.1101/2024.12.26.629440) |
| Tabuteau et al. | [Iron acquisition profiles](https://doi.org/10.1111/1462-2920.70218) |
| KOfam / KEGG | [KOfamKOALA](https://doi.org/10.1093/bioinformatics/btz859) and [KEGG usage terms](https://www.kegg.jp/kegg/legal.html) |
| NCBI / RefSeq protein families | [NCBI Protein Family Models](https://www.ncbi.nlm.nih.gov/Structure/pgap/pgap.cgi) |
| InterPro / Pfam | [InterPro](https://www.ebi.ac.uk/interpro/) and the member database identified by each accession |
| MnOxGeneTool | [Repository and usage notice](https://github.com/wyh19990121/MnOxGeneTool) |

Check `source`, `acc`, `hmm_file`, and provenance fields in
`hmm_library/hmm_registry.tsv` to identify the applicable origin of an individual
profile. Consult the original distributions for their complete terms, including
when profiles were obtained through another compilation. Optional external tools
and their databases have separate licenses and citation requirements.
