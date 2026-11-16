# Installation

## Install from source

```bash
git clone https://github.com/l-gallucci/Efesto.git
cd Efesto
mamba env create -f environment.yml
conda activate efesto
efesto --help
efesto checkenv
```

The environment includes HMMER, Pyrodigal, Pyrodigal-GV, samtools and plotting
packages, and installs the local checkout with `pip install -e .`.
The package requires Python 3.9 or later; the supplied environment constrains
Python below 3.12 to accommodate optional eggNOG-mapper installations.

Editable installations read `hmm_library/` from the checkout, so keep that
checkout available. Regular wheel/source-package installations include the same
library inside the installed Python package and discover it automatically.
`--hmm_dir` can still select a custom library explicitly.

## Install into an existing environment

```bash
conda install -c conda-forge -c bioconda hmmer pyrodigal pyrodigal-gv
python -m pip install -e /path/to/Efesto
```

`pip install` alone does not install the external bioinformatics tools or the
optional gene-calling packages. For FAA-only runs, gene calling is unnecessary.
HMMER executables must be on `PATH`; Pyrodigal and Pyrodigal-GV are Python imports,
not a requirement to install a separate `prodigal` executable.

## Install a release wheel

Install HMMER in the target environment, then install the downloaded release file:

```bash
conda install -c conda-forge -c bioconda hmmer
python -m pip install /path/to/efesto-1.0.0-py3-none-any.whl
efesto --version
```

The wheel includes the registered HMM profiles, cutoff and gene maps, category
metadata, and source notices. It does not need the Git checkout at runtime.
For nucleotide inputs, install Pyrodigal/Pyrodigal-GV with conda as above or use
`python -m pip install 'efesto[gene-calling] @ file:///absolute/path/to/efesto-1.0.0-py3-none-any.whl'`.
The optional JSON context-rule example is packaged separately under
`efesto/data/examples/operon_rules.json`; it does not activate the JSON engine
unless copied into the selected HMM library. HMMER and optional external tools
remain separate installations.

## Optional tools

| Feature | Requirement | Efesto option |
|---|---|---|
| BAM coverage | samtools | `--bam` / `--bams` |
| Operon prediction | UniOP and its dependencies | `--operon_prediction` |
| Existing eggNOG annotations | Matching query IDs in an annotations file | `--eggnog_annotations` |
| Run eggNOG-mapper | eggNOG-mapper and a downloaded reference database | `--run_eggnog` |
| Topology checks | Phobius, installed under its license | `--topology_check` |
| Structural checks | Foldseek, ProstT5 weights and a curated reference database | `--run_struct_confirmation` |
| BGC context | Previously generated antiSMASH GFF3 | `--bgc_dir` |
| Heatmap figures | R packages in `environment.yml` | Plotting script |

Efesto does not download the large eggNOG or structure databases automatically.
Use `efesto checkenv` to inspect installed tools and optional dependencies.

## UniOP

Install [UniOP](https://github.com/hongsua/UniOP) following its own instructions,
then provide the script path:

```bash
efesto --fna_dir genomes/ --meta --out efesto_out \
  --operon_prediction --uniop_path /path/to/UniOP/src/UniOP
```

Efesto supplies its Prodigal-compatible protein output to UniOP with
`--faa_source prodigal_faa`. This is a file-format setting: Efesto's gene caller
remains Pyrodigal/Pyrodigal-GV. See [operon prediction](operon_rules.md#uniop).

## Optional confirmation tools

An existing annotation file is sufficient for the first example; eggNOG-mapper
does not need to be installed for that mode.

```bash
efesto --faa_dir proteins/ --gff_dir gff/ --out efesto_out \
  --eggnog_annotations annotations.emapper.annotations
```

For internal eggNOG-mapper execution, install it in a compatible environment and
supply `--eggnog_db_dir`. Structural checks additionally require a curated
Foldseek reference database and ProstT5 weights; topology checks require Phobius.
See the [command reference](command-reference.md) for paths and options and
[confidence scoring](cluster_confidence_scoring.md) for the scope of each check.

A missing optional result generally contributes a neutral score, not supporting
evidence. Inspect the run log for warnings and tool failures.
