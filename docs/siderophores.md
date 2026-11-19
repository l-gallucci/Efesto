# Siderophores and metallophores

Efesto combines biosynthetic, acquisition and processing profiles from several
collections. A family hit identifies a candidate component. It does not by itself
identify the secreted molecule, a complete pathway, the transported metal or a
physiological phenotype. Siderophores are iron-binding members of the broader
metallophore group; some ligands can bind additional metals.

## Functional categories

| Category | Meaning |
|---|---|
| `Siderophore-Synthesis` | Biosynthetic enzymes, shared accessories and biosynthetic-pathway components, including precursor transport by PvdE |
| `Siderophore-Transport` | Uptake-associated components; broad families can have unresolved substrates or transport directions |
| `Siderophore-Export` | Outward secretion or ligand recycling; EntS, PvsC and PvdRT-related components |
| `Siderophore-Utilization` | Ligand processing and utilization-associated iron release; includes hydrolases and accessory proteins, not one shared reaction |
| `Siderophore-FunctionUnresolved` | Shared PF06276 domain: NIS biosynthesis versus FhuF-like reduction cannot be distinguished from the domain match alone |
| `Siderophore-TransportPotential` | Shared energy-coupling machinery such as TonB/Exb; not substrate-specific uptake evidence |
| `IronAcquisition-IronRelease` | Nutritional ferric-chelate reduction, including nonsiderophore chelates; distinct from respiratory `Iron-Reduction` |
| `Metallophore-Synthesis` | Cnt-type opine and Mbn-type peptide biosynthetic families; no single product or metal assignment from a generic family match |

`Fe`, `Fe-acquisition` and `Fe-metabolism` include the siderophore roles and
nutritional iron release. `metallophores` retains its synthesis-only meaning;
use the exact categories or `all` to include uptake, export and processing.

## Shared targets and domains

IroB from FeGenie and K23725 from KOfam (imported through Tabuteau et al.) both
target the C-glucosyltransferase involved in salmochelin biosynthesis. They use
different training sets and matrices. They are retained as overlapping models,
not as two independent biosynthetic steps. The reaction is characterized by
[Fischbach et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC545562/).

YqjH catalyses reductive iron release from several chelates, including ferric
dicitrate, using NADPH/FAD. Its scope includes siderophore utilization; the broad
category also covers iron release from other chelates. This is acquisition-associated iron assimilation,
not a membrane transporter or dissimilatory Fe(III) respiration
([Miethke et al.](https://doi.org/10.1021/bi201517h)).

PF06276 has a broader domain scope: it occurs in FhuF reductases and in
IucA/IucC-like NIS biosynthetic enzymes. The domain cannot by itself establish
either activity. PF06276 and the unresolved PHYLLO assignment are marked
`experimental`; both remain enabled in ordinary runs. Full protein architecture and specific family evidence are
needed ([Trindade et al.](https://doi.org/10.1007/s00775-021-01854-y),
[Bailey et al.](https://doi.org/10.1021/acs.biochem.6b00409),
[provider entry](https://www.ebi.ac.uk/interpro/entry/pfam/PF06276/)).

K02014 and K02016 identify broad receptor/binding components, not a named ligand.
Shared NRPS/NIS, phosphopantetheinyl-transferase and MbtH-like profiles also need
partner genes and genomic context. General precursor-supply enzymes and
locus-associated models remain less specific than experimentally characterized
pathway enzymes. For example, K14759/PHYLLO is a multifunctional quinone-biosynthesis
family and is not a specific siderophore marker.

VabS is a putative exporter: the vanchrobactin study supports a pathway-associated
role but does not directly measure specific transport. VabH has genetic support
for utilization; its precise reaction is not established by the IroD reaction
([Balado et al.](https://doi.org/10.1099/mic.0.29298-0)).

## Functional labels and neighborhood rules

Output categories describe biological roles. The bundled neighborhood rules use
an explicit, stable context grouping for reclassified profiles, preserving the
established candidate-retention policy while refining labels. For example,
BesA is reported as utilization but retains its original synthesis-neighborhood
membership. A minimum number of distinct HMMs is a contextual heuristic,
not a count of independent reactions or proof of pathway completeness.

The supplied JSON rules declare `category_basis: context` for the siderophore
count rules. Custom rules use annotation categories by default. See
[operon rules](operon_rules.md) and [separate evidence](evidence.md).

## Complete profile table

The table includes enabled and excluded profiles. Sources describe import
provenance; scope notes distinguish reference-family function from ligand or
pathway specificity. Excluded profiles are not searched in ordinary runs.
The machine-readable version is
[`siderophore_metadata.tsv`](../hmm_library/siderophore_metadata.tsv).

{{siderophore_table}}
