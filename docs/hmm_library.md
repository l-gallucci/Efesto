# HMM Library

The shipped library contains **652 registered HMMs**: **473 active**, **12 experimental** and **167 deprecated** (485 enabled). Category cleanup: 5 November 2026. Remaining profile bytes, accessions, seed counts and numeric search cutoffs were preserved.

## Category names and functional scope

Names follow `Process-Function`, for example `IronAcquisition-Transport`, `Siderophore-Synthesis`, `Siderophore-Transport`, `Heme-Transport`, `Iron-Reduction` and `MetalResistance-Copper`. Acquisition, homeostasis, regulation, stress response and resistance are distinct assignments. A model assigned to a process does not establish an organism's activity or resistance phenotype.

`IronAcquisition-ComplexTransport` contains families with alternative siderophore/heme substrates. `MetalAcquisition-MultiMetal` and `MetalAcquisition-CobaltNickel` preserve known substrate ambiguity. `Siderophore-Transport` does not assert that every family only imports: Pul3 also supports export of pulcherriminic acid. `Unclassified-MetalFunction` contains unresolved source assignments; these must not be counted as confirmed metal-resistance determinants.

The catalog below includes inactive historical records. A category with zero enabled profiles is retained only for provenance. Synthetic potential categories used by operon disambiguation and custom libraries may also appear at runtime.

| Category | Registered | Active | Experimental | Deprecated |
|---|---:|---:|---:|---:|
| `Arsenic-Oxidation` | 2 | 2 | 0 | 0 |
| `Competence-Pilus` | 1 | 0 | 0 | 1 |
| `Heme-Oxygenase` | 14 | 11 | 0 | 3 |
| `Heme-Transport` | 12 | 12 | 0 | 0 |
| `Iron-Oxidation` | 15 | 14 | 1 | 0 |
| `Iron-Reduction` | 34 | 30 | 1 | 3 |
| `Iron-Regulation` | 6 | 5 | 0 | 1 |
| `Iron-Resistance` | 3 | 0 | 0 | 3 |
| `Iron-Storage` | 3 | 3 | 0 | 0 |
| `Iron-StressResponse` | 2 | 2 | 0 | 0 |
| `IronAcquisition-ComplexTransport` | 3 | 3 | 0 | 0 |
| `IronAcquisition-Transport` | 38 | 37 | 0 | 1 |
| `IronSulfur-Assembly` | 7 | 7 | 0 | 0 |
| `IronSulfur-Regulation` | 1 | 1 | 0 | 0 |
| `Magnetosome-Formation` | 13 | 12 | 0 | 1 |
| `Manganese-Acquisition` | 2 | 2 | 0 | 0 |
| `Manganese-Oxidation` | 12 | 12 | 0 | 0 |
| `Manganese-Regulation` | 1 | 1 | 0 | 0 |
| `MetalAcquisition-CobaltNickel` | 2 | 2 | 0 | 0 |
| `MetalAcquisition-Copper` | 2 | 2 | 0 | 0 |
| `MetalAcquisition-IronManganese` | 2 | 2 | 0 | 0 |
| `MetalAcquisition-Magnesium` | 2 | 2 | 0 | 0 |
| `MetalAcquisition-Molybdenum` | 3 | 2 | 0 | 1 |
| `MetalAcquisition-MultiMetal` | 7 | 5 | 2 | 0 |
| `MetalAcquisition-Nickel` | 6 | 6 | 0 | 0 |
| `MetalAcquisition-Zinc` | 5 | 5 | 0 | 0 |
| `MetalHomeostasis-Copper` | 9 | 7 | 2 | 0 |
| `MetalHomeostasis-IronZinc` | 1 | 1 | 0 | 0 |
| `MetalRegulation-MultiMetal` | 3 | 3 | 0 | 0 |
| `MetalRegulation-Zinc` | 1 | 1 | 0 | 0 |
| `MetalResistance-Arsenic` | 31 | 10 | 0 | 21 |
| `MetalResistance-Chromium` | 6 | 1 | 1 | 4 |
| `MetalResistance-CobaltZincCadmium` | 16 | 3 | 0 | 13 |
| `MetalResistance-Copper` | 36 | 12 | 2 | 22 |
| `MetalResistance-Manganese` | 2 | 2 | 0 | 0 |
| `MetalResistance-Mercury` | 31 | 4 | 0 | 27 |
| `MetalResistance-MultiMetal` | 45 | 17 | 3 | 25 |
| `MetalResistance-Nickel` | 1 | 0 | 0 | 1 |
| `MetalResistance-NonSpecific` | 15 | 2 | 0 | 13 |
| `MetalResistance-Silver` | 4 | 2 | 0 | 2 |
| `MetalResistance-Tellurite` | 12 | 6 | 0 | 6 |
| `Siderophore-Synthesis` | 162 | 148 | 0 | 14 |
| `Siderophore-Transport` | 65 | 61 | 0 | 4 |
| `Siderophore-TransportPotential` | 6 | 6 | 0 | 0 |
| `StressResponse-Nitrosative` | 1 | 1 | 0 | 0 |
| `StressResponse-Oxidative` | 1 | 1 | 0 | 0 |
| `StressResponse-Redox` | 1 | 1 | 0 | 0 |
| `Unclassified-MetalFunction` | 4 | 4 | 0 | 0 |
| `Virulence-ActinNucleation` | 1 | 0 | 0 | 1 |

## Provenance and review status

`source` identifies the import channel, which can differ from the original model author/database. The preceding source audit verified upstream identity for 616 profiles, using exact bytes or equivalent scoring parameters against frozen snapshots. That is not biological validation. Exact-release identity for the remaining direct NCBIfam/Pfam imports is still incomplete. The five local builds MtrA, MtoA, CymA, ImoA and DA_402 have now been reproduced exactly; CbcL has its canonical reference verified but seven original seed accessions and its MSA remain missing. See `data/seeds/curated/`.

`category_curation.tsv` records explicit model decisions and evidence; the importer reads this table to avoid reintroducing known assignments. `category_aliases.tsv` maps historical spellings to canonical names. `category_catalog.tsv` records the inventory at cleanup time; regenerate it when the registry changes. Category aliases rename labels only: they do not apply biological recategorization to an external library or change original FeGenie profiles.

NF041018/NF041019/NF041020 (CutABC glyceraldehyde dehydrogenase) were removed from the shipped profiles, registry and runtime maps on 5 November 2026. Exclusion tombstones in `category_curation.tsv` prevent reimport. Metal cofactors alone do not place these carbon-metabolism enzymes in scope. NF040982 (competence pilin ComGD) remains archived as `deprecated_out_of_scope`. The FeGenie profile labeled FpvI is now `deprecated_misassigned_duplicate`: its numerical parameters are identical to FpvH and the reference-sequence check supports FpvH, which remains enabled.

Each registry row has `category_review_status`: `source_supported` means the assignment is supported by the cited source/header, `requires_review` means unresolved, and `not_individually_reviewed` means only naming has been normalized. None means an independently validated classifier. There are 211 source-supported records, 6 records requiring review (including the inactive mislabeled FpvI copy), and 435 records awaiting individual biological adjudication.

## Registry

`hmm_registry.tsv` is a tab-separated table with one row per model. `stem` matches the filename without `.hmm`; the internal HMM NAME can differ intentionally. `category` must match the parent directory in `hmm_file`, a path relative to the library. Other columns preserve model name, accession, source, NSEQ, cutoff, added date, status, bibliography, validation cohorts, confirmation settings and topology expectations. A blank `validated_in` is a metadata gap rather than proof of no historical validation.

The cleanup adds `category_original`, `category_review_status` and `category_note`. `category_original` preserves the pre-migration assignment. Review notes describe the functional evidence or the remaining uncertainty. Existing score/confirmation/topology columns retain their values.

## Search thresholds

`HMM-bitcutoffs.txt` and `FeGenie-map.txt` were unchanged by category renaming; the later physical removal of CutABC deleted only their three entries. In particular, 129 enabled models still lack a model-specific cutoff and use the existing generic filters. Threshold calibration requires an independent labeled evaluation and is separate from category naming. Historical rounding differences between registry and runtime thresholds remain documented in the source audit.

## Runtime and compatibility

`--annotate Fe`, `Cu`, `Ni`, `Mo`, etc. include the relevant new acquisition/homeostasis categories. Function-specific selectors such as `Ni-acquisition`, `Mo-acquisition`, `Cu-homeostasis` and `As-oxidation` are available. Resistance-specific selectors retain their resistance meaning. Exact canonical category names can also be passed to `--annotate`.

Legacy library-directory names and JSON rule categories are normalized on load. Unknown custom names are preserved. New output uses canonical names. Existing output files and frozen benchmark datasets are not rewritten. Comparison code normalizes spelling only; disagreements caused by biological reassignment remain visible. The private original-FeGenie adapter preserves original source categories and hashes in its manifest.

Default co-occurrence thresholds remain unchanged. Metal acquisition/homeostasis/regulation categories preserve the bypass used for their previous metal-resistance assignments; iron/siderophore assignments follow their corrected context group. Therefore a corrected category can change downstream retention and counts even with unchanged HMM scores. Old benchmark accuracy results must be rerun before attributing them to this library version.

## Maintenance

Verify paths and registry consistency:

```bash
python scripts/curate_hmm_library.py --verify hmm_library
```

Inspect a migration without changing files:

```bash
python scripts/migrate_categories.py --library hmm_library
```

The maintenance script uses `--apply` to execute its preflighted plan, checks file hashes during moves and rolls back moves if registry writing fails. Running it again after migration reports zero moves. To revise functional assignments, update the explicit curation table and registry together, rerun verification and regression tests, and record the supporting evidence.
