# Installation

---

## Requirements

| Dependency | Version | Role |
|---|---|---|
| HMMER | ≥ 3.3 | Profile HMM search |
| Prodigal / Pyrodigal | any | ORF calling (`--fna_dir`) |
| Python | ≥ 3.10 | Pipeline runtime |
| samtools | ≥ 1.10 | Coverage (`--bam` / `--bams`) |
| UniOP | any | Operon prediction (`--operon_prediction`) |
| R + pheatmap + plotly | any | Heatmap visualisation (optional) |

---

## Option 1 — Conda / Mamba (recommended)

```bash
git clone https://github.com/YOUR_ORG/Efesto.git
cd Efesto
mamba env create -f environment.yml
conda activate efesto
Efesto --help
```

`environment.yml` installs all Python and non-Python dependencies and registers
the `Efesto` command via `pip install -e .`.

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

## Verifying the installation

```bash
# Check CLI is available
Efesto --help

# Verify HMM library integrity
python scripts/curate_hmm_library.py --verify hmm_library/

# Run tests
pytest tests/ -v
```

---

## HMM library first run

On first use with a library that contains pre-HMMER3/f profiles, run:

```bash
Efesto --normalize_hmms --faa_dir orfs/ --hmm_dir hmm_library/ --out results/
```

`--normalize_hmms` converts any legacy profiles in-place using `hmmconvert` before
the search. It is safe to run repeatedly (already-current files are skipped). After
the first conversion the flag is no longer needed.

---

## R packages (optional, for heatmaps)

```bash
conda install -c conda-forge r-pheatmap r-plotly r-htmlwidgets r-optparse r-rcolorbrewer r-scales
```
