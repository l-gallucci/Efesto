# Operon Rules

Efesto applies operon-context filters to remove spurious single-gene hits
and enforce co-occurrence requirements before reporting clusters. Two rule engines
are available: a default FeGenie-exact port, and a configurable JSON rule engine.

---

## Default filtering (no `operon_rules.json`)

Without `operon_rules.json` in `--hmm_dir`, Efesto uses an exact
reimplementation of FeGenie's per-category operon rules. Each cluster is routed
to one handler based on which special genes are present. For iron acquisition
categories, a per-ORF pass/break pattern is used: a gene that fails co-occurrence
is silently skipped, but other genes in the same cluster are still reported.

Use `--all_results` to bypass all filtering entirely.

---

## JSON rule engine

If `operon_rules.json` exists in `--hmm_dir`, the JSON engine replaces the default.
The file must be in the same directory as the HMM files.

### Top-level keys

```json
{
  "report_all_categories": ["metal_resistance-*", "iron_storage"],
  "rules": [ ... ]
}
```

**`report_all_categories`** — glob patterns. Clusters whose entire category set
matches bypass all rules and are reported unconditionally. MetHMMDB categories
(`metal_resistance-*`) and `iron_storage` bypass by default. `iron_sulfur_assembly`
and `iron_stress` are also in this list — single-gene hits are informative.

### Rule object schema

```json
{
  "name":            "RULE_NAME",
  "categories":      ["category_1", "category_2"],
  "genes":           ["GeneA", "GeneB"],
  "rule":            "require_n_of",
  "min_genes":       2,
  "on_fail":         "drop",
  "canonical_size":  5,
  "max_bp_gap":      2000,
  "neutral_energizer_stems": ["PF03544_TonB_C", "PF13103_TonB_2"]
}
```

| Field | Required | Description |
|-------|----------|-------------|
| `name` | yes | Rule identifier (shown in debug output) |
| `categories` | yes | Categories this rule governs |
| `genes` | rule-dependent | Gene names from `FeGenie-map.txt` |
| `rule` | yes | Rule type (see below) |
| `min_genes` | rule-dependent | Minimum count threshold |
| `on_fail` | yes | Action when rule fails |
| `canonical_size` | no | Expected gene count for completeness scoring |
| `max_bp_gap` | no | Override global `--max_bp_gap` for this rule's genes |
| `neutral_energizer_stems` | no | Stems that are substrate-agnostic (see below) |

---

## Rule types

| Rule | Behaviour |
|------|-----------|
| `require_n_of` | Cluster must contain ≥ `min_genes` unique members of `genes` list |
| `require_anchor` | A specific `anchor` gene must be present in the cluster |
| `require_n_cat` | Cluster must contain ≥ `min_genes` distinct HMMs from `categories` |
| `require_n_cat_or_lone_trusted` | As `require_n_cat`, or lone hit in `trusted_lone` list |
| `mtr_disambiguation` | Re-assigns `iron_oxidation` vs `iron_reduction` based on MtrA/MtoA co-presence |

---

## On-fail actions

| Action | Behaviour |
|--------|-----------|
| `passthrough_non_members` | Keep genes NOT in the rule's `genes` list; drop genes that are in it |
| `drop` | Remove the entire cluster from output |
| `keep_all` | Keep everything regardless of rule outcome |

---

## Per-rule bp-gap override (`max_bp_gap`)

When `max_bp_gap` is set in a rule, Efesto applies the more restrictive
of the two adjacent genes' rule-specific gaps during clustering. This makes tight
operons (e.g. MtrMto: 2 000 bp, FoxABC: 1 000 bp) stricter than the global
window without penalising loosely encoded systems (e.g. SIDERO_SYNTH: 5 000 bp).

---

## The TonB energizer guard (`neutral_energizer_stems`)

### Background — why TonB/ExbBD are substrate-agnostic

TonB, ExbB, and ExbD form the **Ton motor complex** — an inner membrane machinery
that harnesses the proton motive force (PMF) and delivers it to TonB-dependent
outer membrane transporters (TBDTs). The structural basis is well established:
ExbB forms a pentamer with an ExbD dimer inside its transmembrane pore; this
complex converts the electrochemical gradient into mechanical energy; TonB spans
the periplasm and cyclically contacts TBDTs via a conserved five-residue TonB box,
triggering plug movement and substrate release into the periplasm
([Ratliff, Celia & Buchanan 2022, *Front Microbiol*](https://doi.org/10.3389/fmicb.2022.852955);
[Celia et al. 2025, *Nat Commun*](https://doi.org/10.1038/s41467-025-61286-z)).

> **Note (2026-07-15):** the three citations previously here (Celia et al. 2016
> *Nature*; Silale & van den Berg 2023 *Annu Rev Microbiol*; Braun 2024
> *Mol Microbiol*) all had dead or mismatched DOIs — none resolved to a real
> paper matching the claimed author/year/journal. Replaced with verified papers
> covering the same structural claim (ExbB/ExbD/TonB architecture and mechanism).

**Critically, the Ton motor is substrate-agnostic.** The same TonB-ExbBD
complex energises ALL TBDTs regardless of their cargo. In *Bacteroidetes*, TBDTs
(SusC-like proteins) are the invariant component of every polysaccharide
utilisation locus (PUL), and TonB-ExbBD powers their carbohydrate import —
not iron uptake ([Koropatkin, Cameron & Martens 2012, *Nat Rev Microbiol*](https://doi.org/10.1038/nrmicro2746);
[Pollet et al. 2021, *Mol Microbiol*](https://doi.org/10.1111/mmi.14695)).
During phytoplankton blooms, TonB-dependent transporters were the most highly
expressed protein class in marine bacterioplankton (~16.7% of all detected proteins),
with the majority predicted to target polysaccharides
([Francis et al. 2021, *ISME J*](https://doi.org/10.1038/s41396-020-00858-x)).
In *Xanthomonas* plant pathogens, CAZyme clusters are organised in CUT systems
(**C**arbohydrate **U**tilisation with **T**onB-dependent transporters), directly
analogous to Bacteroidetes PULs
([Giuseppe et al. 2023, *Essays Biochem*](https://doi.org/10.1042/EBC20220128)).

### The problem

Efesto models the Ton motor with five HMMs:

| Stem | Pfam | Role |
|------|------|------|
| `PF03544_TonB_C` | PF03544 | TonB C-terminal periplasmic domain |
| `PF13103_TonB_2` | PF13103 | TonB domain 2 |
| `PF16031_TonB_N` | PF16031 | TonB N-terminal anchor |
| `PF01618-MotA_TolQ_ExbB_proton_channel_family` | PF01618 | ExbB / MotA / TolQ proton channel |
| `PF02472-ExbD` | PF02472 | ExbD periplasmic signalling |

All five are in `iron_acquisition-siderophore_transport_potential`. The `SIDERO_TRANSPORT`
rule requires `min_genes: 2` from `siderophore_transport_potential`. Therefore:

```
PF03544_TonB_C  +  PF01618_ExbB  →  2 hits  →  passes SIDERO_TRANSPORT
```

This is a **false positive** in Bacteroidetes-rich metagenomes: a Ton motor hit
without a siderophore receptor is evidence of carbohydrate transport, not iron
acquisition.

### The fix — `neutral_energizer_stems`

The `SIDERO_TRANSPORT` rule accepts an optional `neutral_energizer_stems` list.
If ALL hits in a cluster come exclusively from this list (no substrate-specific
receptor), the rule treats it as a fail → `drop`.

```json
{
  "name":       "SIDERO_TRANSPORT",
  "categories": ["iron_acquisition-siderophore_transport_potential",
                 "iron_acquisition-heme_transport",
                 "iron_acquisition-siderophore_transport"],
  "rule":       "require_n_cat_or_lone_trusted",
  "min_genes":  2,
  "neutral_energizer_stems": [
    "PF03544_TonB_C",
    "PF13103_TonB_2",
    "PF16031_TonB_N",
    "PF01618-MotA_TolQ_ExbB_proton_channel_family",
    "PF02472-ExbD",
    "HasB-TonB_paralog_heme_uptake-rep"
  ],
  "trusted_lone": [ ... ],
  "on_fail": "drop",
  "max_bp_gap": 2000
}
```

**Logic (proposed implementation in `operon.py`):**

```python
def passes_energizer_guard(cluster_hits, neutral_stems):
    """Return False (fail) if all hits are neutral energizers."""
    non_neutral = [h for h in cluster_hits if h["stem"] not in neutral_stems]
    return len(non_neutral) > 0
```

A cluster passes the guard when at least one hit is a substrate-specific gene
(receptor, periplasmic binding protein, ABC permease, etc.).

### Clusters the guard catches

| Composition | Guard outcome | Reason |
|---|---|---|
| `TonB_C` + `ExbB` | **Drop** | Ton motor only; substrate unknown |
| `TonB_C` + `ExbD` + `FutA1` | Pass | FutA1 is an iron-specific ABC binding protein |
| `TonB_C` + `PF00593-TonB_dependent_receptor` | Pass | TonB-dep receptor is substrate-specific |
| `ExbB` + `FepC-ATPase` | Pass | FepC is siderophore-specific |
| `HasB` + `HxuA` | Pass | HasB (heme-specific TonB paralog) + HxuA is heme-specific |

### HasB — the heme-specific TonB paralog

`HasB` (`HasB-TonB_paralog_heme_uptake-rep`) is a TonB paralog in *Serratia marcescens*
and related organisms. Unlike canonical TonB it is dedicated to the heme acquisition
system (HasASB), acting specifically with the HasR outer membrane heme receptor.
Its ExbB-interaction surface has a unique long periplasmic extension that discriminates
it from canonical TonB ([Biou et al. 2022, *Commun Biol*](https://doi.org/10.1038/s42003-022-03306-y)).
HasB is therefore substrate-specific (heme) — it is listed in `neutral_energizer_stems`
only when appearing WITHOUT HasR; when co-clustered with HasR or other heme acquisition
genes, it is a valid heme-transport signal.

> **Status:** The `neutral_energizer_stems` guard is the proposed implementation.
> The operon.py rule engine requires a ~20-line update to support it. See the
> open issue on the development tracker.

---

## Current rules (active `operon_rules.json`)

### FLEET

Iron oxidation via the electron shuttle complex in *Sideroxydans*-type organisms.
Requires ≥ 5 of 8 subunits (EetA/B, Ndh2, FmnB/A, DmkA/B, PplA).
`canonical_size: 8`, `max_bp_gap: 2000`. On fail: keep non-FLEET genes.

### MAM

Magnetosome island. Requires ≥ 5 of 10 MAM genes (MamA/B/E/K/P/M/Q/I/L/O).
`canonical_size: 10`, `max_bp_gap: 3000`. On fail: keep non-MAM genes.

### FOXABC

*Acidithiobacillus* FoxABC iron oxidation complex. Requires ≥ 2 of 3 subunits.
`canonical_size: 3`, `max_bp_gap: 1000`.

### FOXEYZ

FoxEYZ cytochrome *bc*1-type complex. FoxE is the mandatory anchor gene.
`canonical_size: 3`, `max_bp_gap: 1000`.

### DFE1 / DFE2

*Desulfosporosinus* iron reduction operons 1 and 2. Require ≥ 3 of 4 / 3 of 5
subunits respectively. `max_bp_gap: 1000`.

### MtrMto

Mtr/Mto disambiguation. Re-assigns `iron_oxidation` vs `iron_reduction` based
on co-presence of MtrA, MtrB, MtrC (reduction) vs MtoA (oxidation) and CymA.
`canonical_size: 3`, `max_bp_gap: 2000`. On fail: keep all.

### SIDERO_TRANSPORT

Siderophore/heme transport. Requires ≥ 2 distinct transport HMMs, or a lone
trusted receptor. `max_bp_gap: 2000`. On fail: drop.
Trusted lone genes: FutA1, FutA2, FutC, LbtU/LbtB variants, IroC.

### SIDERO_SYNTH

Siderophore biosynthesis. Requires ≥ 3 distinct synthesis HMMs.
`max_bp_gap: 5000` (biosynthesis operons can be large). On fail: drop.

### IRON_TRANSPORT

Iron and heme transport. Requires ≥ 2 distinct HMMs from `iron_acquisition-iron_transport`
and `iron_acquisition-heme_oxygenase`. `max_bp_gap: 2000`. On fail: drop.

---

## `--catalog_mode`

Bypasses co-occurrence count rules (FLEET ≥ 5, siderophore ≥ 2/3, iron_transport ≥ 2)
while keeping bitscore cutoffs, Cyc2 logic, and Mtr/Mto disambiguation. Use when
input is a deduplicated gene catalog (e.g. from CARD, KEGG) where genomic context
is unavailable and co-occurrence cannot be meaningfully tested.
