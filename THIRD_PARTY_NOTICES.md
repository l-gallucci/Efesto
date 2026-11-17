# Third-party sources and conditions of use

All Efesto-owned code, including the Python application, standalone scripts,
setup utilities and tests, is distributed under **AGPLv3**. Project documentation
uses the same license. See
[license scope](LICENSES/README.md). These software licenses do not replace
the terms of imported HMMs, reference sequences, databases or external tools.
The logo is credited separately in the project README.

Cite Efesto and the original publications/databases relevant to an analysis.
Scientific references are listed in [References](docs/references.md).
Citation and a link to a source do not by themselves grant redistribution rights.
Academic use must also comply with the original conditions.

## Model sources

The registry's `source` records the import channel, which may differ from the
original database. [model_licenses.tsv](hmm_library/model_licenses.tsv) maps
**each registered HMM** to its source, applicable notice and any unresolved
conditions. It also covers excluded profiles: runtime status does not change
the terms for files included in a package. The table is a record of the source
notices checked on **2026-11-16**, not a new license or a permission grant.

| Source / original provider | Terms identified | Attribution and original source |
|---|---|---|
| FeGenie original/specialized profiles and context rules | AGPLv3 repository notice; no separate HMM exception located. The notice is preserved; model-specific scope/preferred source form is not independently certified. | FeGenie authors, [Garber et al. (2020)](https://doi.org/10.3389/fmicb.2020.00037); [repository and license](https://github.com/Arkadiy-Garber/FeGenie/blob/master/LICENSE) |
| Pfam, including original Pfam profiles imported through FeGenie/InterPro | CC0 on original Pfam data; preserve provenance and distinguish any collection-specific changes/metadata. A PF-like filename alone does not establish Pfam origin. | Pfam/EMBL-EBI; [source terms](https://pfam-docs.readthedocs.io/en/dev/pfam.html) |
| TIGRFAM, including models delivered through NCBI/InterPro | CC BY-SA 4.0. Preserve attribution, the license link and modification notices. FeGenie's MtrB/MtrC files require review of the historical model/version and modifications. | TIGR/JCVI and NCBI; [source and license](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/tigrfams/), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) |
| NCBIfam NF / NCBI-created family models | NCBI's public-domain/no-additional-restrictions policy for its own material, with the stated reservation for third-party rights. This is not a blanket public-domain claim for everything hosted by NCBI. | NCBI/NLM and original contributors; [data policy and disclaimer](https://www.ncbi.nlm.nih.gov/home/about/policies/) |
| Tabuteau original/other profiles | GPLv3 project notice. Imported KOfam/NCBI models retain their original source conditions. | Tabuteau, Herve, Irlinger and Monnet; [publication](https://doi.org/10.1111/1462-2920.70218), [repository and license](https://github.com/STabuteau/iron_cheese_metagenomes) |
| KOfam / KEGG profiles imported through Tabuteau | **Profile redistribution terms not established.** The KOfam README explicitly identifies KEGG rights/usage terms for KO definitions; it does not provide a clear general redistribution grant for the profiles. Do not infer one from KofamScan software licensing or Tabuteau's collection license. | Aramaki et al.; [KOfam README](https://www.genome.jp/ftp/db/kofam/README), [KEGG terms](https://www.kegg.jp/kegg/legal.html), [publication](https://doi.org/10.1093/bioinformatics/btz859) |
| MetHMMDB | **No project/data license specified in the checked sources; refer to the original repository/authors for applicable terms.** Efesto makes no new license grant for these HMMs and does not treat the absence of a notice as redistribution permission. | Ciuchcinski and Dziurzynski; [repository](https://github.com/kciuchcinski/MetHMMDB), [publication](https://doi.org/10.1101/2024.12.26.629440) |
| MnOxGeneTool | **Non-commercial academic use only**, with author consent required for commercial use. Redistribution/modification permission for this repackaging is not established by that use notice. | Wang et al.; [original notice](https://github.com/wyh19990121/MnOxGeneTool/blob/main/Readme.md), [publication](https://doi.org/10.1021/acs.est.5c01235) |
| Efesto rebuilt profiles | CC BY-SA 4.0 only for the six models explicitly listed in [HMM license scope](hmm_library/LICENSES.md), on copyrightable contributions held by Efesto contributors. Seed/source-sequence rights remain applicable. CbcL now has eight traceable input records and an exact reconstruction of its existing build; input rights remain applicable. | Efesto contributors; individual literature references in the model registry and source records in `data/seeds/curated/` |

Unspecified or unresolved terms are **not permission to redistribute**. This
notice records those gaps; it does not resolve them. Before publishing a release
that contains such assets, establish the required rights with the original
provider or obtain a separately documented authorization. Neither academic
intent nor adding this notice establishes that every bundled asset is cleared.
No blanket non-commercial restriction is added to Efesto's AGPLv3 code.

## Reference sequences and input data

NCBI records retain the [NCBI policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).
Copyrightable UniProt database material is provided under
[CC BY 4.0](https://www.uniprot.org/help/license); retain UniProt/source attribution
and identify changes such as alignment or sequence selection. Source accessions,
links and sequence hashes are recorded in the curated and validation manifests.
The license of an HMM does not replace the rights applicable to its seed inputs.

## Separately installed software

[Pyrodigal](https://github.com/althonos/pyrodigal) and
[Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) are GPLv3 Python libraries
imported for gene calling;
their own licenses and notices remain applicable.
[HMMER](https://github.com/EddyRivasLab/hmmer/blob/master/LICENSE) is a separate
executable under its upstream BSD three-clause terms. UniOP, eggNOG-mapper, Foldseek,
Phobius, their models/databases, and upstream annotation tools are not relicensed
by Efesto. Install and use them under their original conditions. In particular,
the checked UniOP repository does not specify a project license, and Phobius has
separate restrictive academic terms.

The wheel and source archive do not bundle the implementations of these
external tools or libraries. They are installed separately by the user. This
source-only workflow does not need to bundle every dependency's source merely
because it uses them. Distributing a combined executable, container or environment
would require satisfying the applicable dependency licenses and source obligations.

Efesto distributes its own standalone AGPLv3 R plotting script. It does not
redistribute the pheatmap or plotly packages or include their implementations.
Users install these packages separately; [pheatmap](https://cran.r-project.org/package=pheatmap)
retains GPL-2 and plotly retains MIT. Merely supplying this source workflow does
not create a blanket requirement to bundle every dependency's source. These
notices do not relicense the dependencies. A future distribution containing
combined software or bundled packages must meet their applicable licenses and
source obligations; see [license scope](LICENSES/README.md).

## Included license texts

- [GNU AGPLv3](LICENSE): all Efesto-owned code and documentation; also preserves the FeGenie notice.
- [GNU GPLv3](LICENSES/GPL-3.0.txt): preserved upstream GPL version 3 text.
- [MIT](LICENSES/MIT.txt): historical Efesto notice, retained for prior grants; not a separate current tooling license.
- [CC BY-SA 4.0 legal text](https://creativecommons.org/licenses/by-sa/4.0/legalcode): TIGRFAM and the explicitly identified Efesto-owned model contributions.

`LicenseRef-ThirdParty-Models` in the package metadata is a descriptive reference
to these source-specific conditions and `hmm_library/model_licenses.tsv`.
It is not a new license or an additional grant of rights.
