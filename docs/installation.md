# Installation

---

## Requirements

| Dependency | Version | Role |
|---|---|---|
| HMMER | ≥ 3.3 | Profile HMM search |
| Prodigal / Pyrodigal | any | ORF calling (`--fna_dir`) |
| Python | ≥ 3.9 | Pipeline runtime |
| samtools | ≥ 1.10 | Coverage (`--bam` / `--bams`) |
| UniOP | any | Operon prediction (`--operon_prediction`) |
| eggNOG-mapper | 2.x | Tier-2 confirmation for `needs_confirmation`-flagged hits (`--run_eggnog`; not needed for `--eggnog_annotations`) — optional, see below |
| R + pheatmap + plotly | any | Heatmap visualisation (optional) |

---

## Option 1 — Conda / Mamba (recommended)

```bash
git clone https://github.com/l-gallucci/Efesto.git
cd Efesto
mamba env create -f environment.yml
conda activate efesto
efesto --help
```

`environment.yml` installs all Python and non-Python dependencies and registers
the `efesto` command via `pip install -e .`.

---

## Option 2 — Add to an existing conda environment

```bash
conda activate YOUR_ENV
pip install -e /path/to/Efesto
# Then install HMMER and samtools separately if not already present
conda install -c bioconda hmmer samtools
```

---

## Option 3 — Manual (no conda)

```bash
# Install HMMER system-wide or via your package manager
# Then:
pip install -e /path/to/Efesto
```

Ensure `hmmsearch`, `hmmpress`, and (if using `--fna_dir`) `prodigal` or
`pyrodigal-gv` are on `$PATH`.

---

## Installing UniOP (optional)

UniOP is required only for `--operon_prediction`. It has no conda package;
install from source:

```bash
git clone https://github.com/hongsua/UniOP.git
# No build needed — pure Python + dependencies (numpy, scipy)
pip install numpy scipy
```

Pass the path to the UniOP script with `--uniop_path /path/to/UniOP/src/UniOP`.

---

## Installing eggNOG-mapper (optional)

eggNOG-mapper is only needed for `--run_eggnog` (Efesto invoking it
internally on the small `needs_confirmation`-flagged ORF subset). If you
already have an `.emapper.annotations` file from a run you did elsewhere,
use `--eggnog_annotations /path/to/file` instead — that path does not
require eggNOG-mapper installed at all.

**Install it into the same `efesto` environment.** `environment.yml` caps
`python<3.12` specifically so this works — eggNOG-mapper's conda package
also pins `python<3.12`, and since the `efesto` environment never exceeds
that ceiling to begin with, adding eggNOG-mapper later never requires a
Python downgrade. Verified directly: a fresh solve of Efesto's full
dependency set together with `eggnog-mapper` resolves cleanly (conda picks
Python 3.11.x), and the test suite passes running under that exact Python
version.

```bash
conda activate efesto
conda install -c bioconda -c conda-forge eggnog-mapper
download_eggnog_data.py --data_dir /path/to/eggnog_db   # tens of GB, one-time
efesto ... --run_eggnog --eggnog_db_dir /path/to/eggnog_db
```

---

## Installing Phobius (optional)

Phobius is only needed for `--topology_check` (transmembrane topology /
signal-peptide prediction, compared against `expected_topology` set in the
registry — see `docs/cluster_confidence_scoring.md`, Component 7). Cheap
enough to run on every applicable hit, but off by default. No conda or pip
package (academic license, manual download):

```bash
# http://phobius.sbc.su.se/instructions.html — request/download the tarball,
# then unpack it somewhere on PATH or point --phobius_path at it directly
efesto ... --topology_check --phobius_path /path/to/phobius.pl
```

Missing Phobius doesn't break anything — `--topology_check` without it
prints a `[WARN]` and `topo_weight` falls back to whatever the heme-motif
signal alone provides (or neutral, for models with no `expected_heme_c`).

---

## Installing Foldseek (optional)

Foldseek is only needed for `--run_struct_confirmation` (tier-3 structural
confirmation, scoped identically to `--run_eggnog` — the small
`needs_confirmation`-flagged subset only, never the whole proteome — see
`docs/hmm_library_curation.md`, Tier 3). Foldseek translates amino-acid
sequences to the 3Di structural alphabet internally via ProstT5
(`--prostt5-model`), so **no separate `torch`/`transformers` Python
dependency is needed** — install Foldseek itself and download the ProstT5
weights once:

```bash
conda activate efesto
conda install -c bioconda -c conda-forge foldseek
foldseek databases ProstT5 /path/to/prostt5_weights tmp   # one-time, ~a few GB
```

You also need a small curated reference-structure DB (built with
`foldseek createdb`) — **not** the full AlphaFold DB/PDB, see
`docs/hmm_library_curation.md`, Tier 3, for why and for the recommended
first pilot target (MtrCAB):

```bash
efesto ... --run_struct_confirmation \
    --foldseek_db_dir /path/to/db --prostt5_weights /path/to/prostt5_weights
```

Missing Foldseek doesn't break anything — `--run_struct_confirmation`
without it prints a `[WARN]` and `struct_weight` stays neutral (1.0).

---

## Verifying the installation

```bash
# Check CLI is available
efesto --help

# Check every optional dependency at once (UniOP, eggNOG-mapper, Phobius,
# Foldseek) — reports [OK]/[WARN]/[FAIL] per tool, doesn't fail the check
# for optional ones you haven't installed yet
efesto checkenv

# Verify HMM library integrity
python scripts/curate_hmm_library.py --verify hmm_library/

# Run tests
pytest tests/ -v
```

---

## HMM library first run

On first use with a library that contains pre-HMMER3/f profiles, run:

```bash
efesto --normalize_hmms --faa_dir orfs/ --hmm_dir hmm_library/ --out results/
```

`--normalize_hmms` converts any legacy profiles in-place using `hmmconvert` before
the search. It is safe to run repeatedly (already-current files are skipped). After
the first conversion the flag is no longer needed.

---

## R packages (optional, for heatmaps)

```bash
conda install -c conda-forge r-pheatmap r-plotly r-htmlwidgets r-optparse r-rcolorbrewer r-scales
```
