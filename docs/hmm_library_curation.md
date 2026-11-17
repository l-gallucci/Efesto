# Curation and provenance

The current inventory and runtime status come from the [registry-backed catalog](hmm_library.md).
This page describes the shipped library and its interpretation limits.
For a consolidated overview of source contributions, rebuilt profiles and
completed comparisons, read the [scope and model recap](comparison.md).

## Evidence records

| File | Purpose |
|---|---|
| `hmm_registry.tsv` | Sources, status, model metadata and references |
| `category_curation.tsv` | Explicit category/status decisions and evidence |
| `category_aliases.tsv` | Historical spelling normalization |
| `HMM-bitcutoffs.txt` | Runtime full-sequence thresholds |
| `FeGenie-map.txt` | Filename-stem to readable-label map |
| `data/seeds/curated/` | Recovered seed provenance and build-reproduction records |
| `data/validation/` | Targeted reference and regression panels |

These files are in the [source repository](https://github.com/l-gallucci/Efesto).
A source-supported category, reproducible build and independently tested classifier
are different kinds of evidence.

## Current decisions

- **ImoA:** deprecated and excluded from Mtr/Mto rules. The cited experiments
  support quinol oxidation/complementation, not direct native Fe(II) oxidation.
- **CymA:** retained as experimental with a multifunctional respiratory scope;
  direct Fe(III)-respiration evidence does not make it an Fe-specific terminal reductase.
- **CntK/L/M:** assigned to `Metallophore-Synthesis`. Normal context requires
  CntL+CntM; CntK is optional. A hit does not uniquely identify the final product.
- **CutA/B/C NF041018/19/20:** glyceraldehyde-dehydrogenase profiles excluded from
  the metal inventory; the old names did not justify an iron-transport assignment.
- **Mislabeled FpvI copy:** deprecated after sequence-parameter/reference checks
  showed it duplicated the FpvH model. No replacement true-FpvI model is claimed.
- **LanM/LanD, Wtp/Tup and MbnA/B/C:** imported official profiles preserve their
  source HMMs and GA cutoffs. See [the expansion page](metal_expansion.md).

## Redundancy

Duplicate files, overlapping families and paralogs require different decisions.
Consensus similarity or matching a small panel can identify candidates for review,
but does not establish functional interchangeability. Retain distinct paralogs
when their biological roles differ. Inspect accession, architecture, seed provenance,
real reference scores and context before deprecating a profile.

## Remaining limitations

Six local curated builds are now reproducible from documented inputs. CbcL's
eight public sequences were cross-verified against UniProt/NCBI, and a fresh
MAFFT/HMMER build reproduces its existing matrix, statistics and checksum exactly
(build date excluded). This is reconstruction of the inputs, not recovery of the
historical files or independent classifier validation. The HMM and cutoff are
unchanged. CbcL is an inner-membrane EET component linked to low-potential Fe(III)
respiration in the canonical organism; a homolog hit alone does not establish
iron reduction. See the source records in `data/seeds/curated/`.
Some MtrA seeds have only broad DmsE-family annotations; DA_402 remains a candidate
family with limited independent functional evidence.

Broad flavodoxin, pilin and SufE-related profiles require restrained interpretation.
SufE/CsdE and other near paralogs cannot be distinguished by a generic domain alone.
See [Fe-S assembly](iron_sulfur_assembly.md) and [interpretation](interpretation.md).
