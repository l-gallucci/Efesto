# Cluster Confidence Scoring — Design Rationale

*Created: 2026-05-24*

---

## The formula

```
cluster_confidence = hmm_weight × co_occ_score × uniop_weight
                      × bgc_boost × annot_weight × struct_weight × topo_weight
```

Capped at 1.2. Recommended call threshold: ≥ 0.5.

Seven independent evidence sources multiplied together. Any single low score
brings the whole cluster down. This is intentional: a gene system with perfect
HMM hits but no spatial coherence is as suspicious as one with perfect spacing
but low-confidence HMMs.

`annot_weight` and `struct_weight` (added 2026-07-20) only ever affect
clusters containing a hit on a model flagged `needs_confirmation` in the
registry — everything else always sees 1.0 for both, identical to pre-2026-07-20
behaviour.

`topo_weight` (added 2026-07-23) is different in scope from the other two:
it checks every hit whose stem sets `expected_heme_c`/`expected_topology` in
the registry, not just `needs_confirmation`-flagged ones — see Component 7.
It is still gated behind `--topology_check`; without that flag it is always
1.0, identical to pre-2026-07-23 behaviour, even though `expected_heme_c`/
`expected_topology` and the underlying `heme_motifs` count are always present
in the registry/output.

---

## Component 1 — `hmm_weight`

```
hmm_weight = mean(w_i)  where w_i = 1.0 if calibrated, 0.5 if low_confidence
```

"calibrated" = hit passed the HMM's GA cutoff (equivalog-level models or
models with calibrated TC/GA/NC thresholds from the curation pipeline).

"low_confidence" = hit scored between NC and GA, or the model is isofunctional /
subfamily-level and expected to cross-hit.

A cluster of all calibrated hits → hmm_weight = 1.0.  
A cluster of all low_confidence hits → hmm_weight = 0.5.  
Mixed → proportional mean.

---

## Component 2 — `co_occ_score`

```
co_occ_score = completeness × distance_factor × edge_factor
```

### 2a. Completeness

```
completeness = min(n_observed / canonical_size, 1.0)
```

Only applies when the rule that matched the cluster defines a `canonical_size`.
If no canonical_size is known (siderophore biosynthesis, iron transport — gene
count is variable by biology), completeness = 1.0 (not penalized).

Single-gene clusters always return co_occ_score = 1.0 (co-occurrence is
undefined for one gene).

### 2b. Distance factor

```
distance_factor = mean( distance_decay(gap_i) )
                  over all consecutive gene pairs sorted by start coordinate
```

```
distance_decay(gap) = exp( -ln(2) × gap / half_life )
```

- gap=0 → decay=1.0
- gap=500 bp → decay=0.5
- gap=1000 bp → decay=0.25
- gap=5000 bp → decay=0.0045

**half_life = 500 bp** is the 90th-percentile intergenic distance within
bacterial operons (Zhang et al. 2006, doi:10.1016/j.compbiolchem.2006.03.002).
This means: a gap matching the typical upper bound of within-operon spacing
gets a 0.5 penalty; genes beyond that are increasingly unlikely to be in the
same operon.

The global `--max_bp_gap` (default 5000 bp) is the detection window — two genes
within this distance can be placed in the same cluster. The distance decay then
applies the scoring penalty within that window: genes 100 bp apart score
much higher than genes 4000 bp apart, even if both pass the detection threshold.
This design tolerates MAG fragmentation without rewarding it.

### 2c. Edge factor

```
edge_factor = 0.7 if any gene is within edge_margin (3000 bp) of a contig edge
              else 1.0
```

A gene near a contig edge means the rest of the operon may be on another contig
(truncation artifact of assembly). The 0.7 penalty is a soft discount, not
elimination, because truncation at a contig edge is expected in MAG work.

The 3000 bp edge margin was chosen to be larger than a typical single gene
(most bacterial genes are ≤ 3 kb) so any gene at the very edge of a contig
triggers the penalty regardless of gene length.

---

## Component 3 — `uniop_weight`

```
uniop_weight = min( uniop_prob(a, b) )  over all pairs (a, b) in cluster
```

**Weakest-link principle**: the pairwise UniOP operon probability is taken as
the minimum across all gene pairs in the cluster. A single low-probability pair
indicates the cluster spans an operon boundary — the system is likely a false
co-location. The minimum is more conservative than the mean but appropriate
because a wrong pair anywhere in the cluster contaminates the result.

Returns 1.0 (neutral) when:
- UniOP was not run (`--operon_prediction` not used)
- No pair probabilities are available for this cluster
- Single-gene cluster

UniOP probabilities come from parsing `uniop.pred` (raw pairwise scores, not
binarized). The threshold=0.5 used for operon membership assignment is separate
from the continuous probability stored here for scoring.

---

## Component 4 — `bgc_boost`

```
bgc_boost = 1.2 if cluster overlaps a siderophore BGC (antiSMASH)
            else 1.0
```

Optional. Requires `--bgc_dir` pointing to antiSMASH output; 1.0 (no boost)
when not provided. Implemented in `src/efesto/bgc.py`
(`bgc_boost_for_cluster`, `parse_antismash_gff`) — a siderophore synthesis
cluster inside an antiSMASH-predicted siderophore BGC is stronger evidence
than one without BGC support. `parse_antismash_gff` pre-filters to only
regions whose antiSMASH `product`/`type` matches siderophore keywords before
`bgc_boost_for_cluster` ever runs, so the boost can only trigger via genuine
overlap with a siderophore-annotated region — antiSMASH itself only detects
*secondary* metabolite biosynthesis (Blin et al. 2016, *Nucleic Acids Res*,
[doi:10.1093/nar/gkw960](https://doi.org/10.1093/nar/gkw960)), a categorically
different class of gene cluster from primary/housekeeping systems like SUF
Fe-S assembly or Tad pilus, so cross-category false boosts are not expected
in practice.

---

## Component 5 — `annot_weight`

```
annot_weight = 1.2 if every needs_confirmation-flagged hit in the cluster
                     is eggNOG-confirmed
             = 0.5 if any needs_confirmation-flagged hit is eggNOG-contradicted
             = 1.0 otherwise (neutral — no flagged hit, eggNOG not run,
                     or no informative eggNOG hit for the flagged ORF)
```

**Weakest-link, same principle as `uniop_weight`.** Only ever touches
clusters containing a hit on a model flagged `needs_confirmation` in the
registry (`sufA`/`sufD`/`sufE` currently — see
`docs/hmm_library_curation.md`, "SUF operon deployment"). Confirmation is via
`src/efesto/eggnog.py`, comparing the HMM's guessed gene identity against an
eggNOG-mapper annotation for the same ORF (KEGG Ortholog number where both
sides have one, falling back to a loose `Preferred_name` match otherwise).
Catches a false-positive mode co-occurrence alone can't: a gene that
coincidentally sits near real anchor genes without actually being what the
HMM guessed. Requires `--eggnog_annotations` (read an existing file) or
`--run_eggnog` (run `emapper.py` internally, scoped to just the flagged ORF
subset) — neither is on by default.

---

## Component 6 — `struct_weight`

Implemented 2026-07-24. Same confirm/contradict/neutral mechanism as
`annot_weight` and same scope (`needs_confirmation`-flagged hits dropped by
tier-1 operon filtering, `pending_confirmation`), but via **Baktfold**
(Foldseek+ProstT5 structural search, `src/efesto/baktfold.py`) instead of
eggNOG orthology — an independent evidence axis for the flagged hits eggNOG
can't resolve either (no KO assigned, uninformative annotation).

```
struct_weight = weakest_link( confirm_struct_hit(...) )
                over all needs_confirmation-flagged rows in the cluster
```

**No separate ProstT5 step.** Foldseek's `easy-search` accepts a plain
amino-acid FASTA and translates it to the 3Di structural alphabet internally
when given `--prostt5-model` (van Kempen et al. 2023, *Nat Biotechnol*,
doi:10.1038/s41587-023-01773-0; ProstT5 — Heinzinger et al. 2024,
doi:10.1093/nargab/lqae150), so `baktfold.py` never imports
torch/transformers itself — one subprocess call handles both translation and
search, same shape as `eggnog.run_eggnog`.

**Reference structure DB is small and curated, not the full AlphaFold
DB/PDB** — same narrow-scope philosophy as the rest of the escalation
pipeline (see the synteny-bootstrapping rejection in
`docs/hmm_library_curation.md`: heavy generic infrastructure was rejected
there for the same reason). Recommended first pilot target: **MtrCAB**,
since a real experimental structure exists (crystal structure, Edwards
et al. 2020, *Cell* 181(3):665-673, doi:10.1016/j.cell.2020.03.032) — a
stronger starting point than sufA/D/E, where structural data is sparser. Building that reference DB is separate curatorial follow-up work,
the same class of effort as curating HMM seeds — not part of the
`struct_weight` code path itself.

**Comparison is a loose substring match** against the reference DB's target
structure name (`baktfold.confirm_struct_hit`) — entries should be named so
the expected gene name appears in the target ID (e.g. `MtrA_ref_1zzz`), the
same naming convention used for HMM seed curation.

**Gated behind `--run_struct_confirmation`** (off by default, requires
`--foldseek_db_dir`). Mirrors `--run_eggnog`'s scope exactly: only the
`pending_confirmation` population (hits on a `needs_confirmation`-flagged
model, dropped by tier-1 filtering), never the whole proteome — the same
scoping argument that keeps tier-2 tractable at genome/metagenome-catalog
scale applies unchanged here.

---

## Component 7 — `topo_weight`

```
topo_weight = weakest_link( topology_verdict(...), heme_verdict(...) )
              over all cluster rows whose stem has an expectation set
```

**Runs on every applicable hit, not just `needs_confirmation`-flagged ones.**
Unlike `annot_weight`/`struct_weight`, this tier isn't scoped to narrow-margin
HMM models — it's scoped to whichever registry stems declare
`expected_heme_c`/`expected_topology` (see
`docs/hmm_library_curation.md`). That scoping is deliberately by *declared
expectation*, not by *total hit count*: as more stems get expectations
curated, the applicable set grows with curation effort, not with genome or
metagenome size — the same scaling argument that keeps tier-2/tier-3
tractable applies here for a different reason (cheap tool, not narrow
subset).

Pilot set (2026-07-24): all 12 `status=experimental` models in the registry
now have `expected_topology` (and, for the two cytochromes, `expected_heme_c`)
set — MtrA/MtoA (the original pilot) plus 10 metal-transport models added in
the same pass, each independently verified against literature before being
registered (see `docs/hmm_library_curation.md` for per-model citations and
confidence notes — one of the ten, `NiCoCd_resistance_nccB_1`, is flagged
lower-confidence there because its lipoprotein membrane anchor is a genuinely
harder call for Phobius than a classic Sec-cleaved signal peptide).

Two independent signals per applicable row, both from `src/efesto/topology.py`:

- **`topology_verdict`** — Phobius transmembrane-helix-count and
  signal-peptide prediction, run once per genome in a single batched call
  (`--topology_check`, requires Phobius on `PATH` or `--phobius_path`).
  Conservative by design: only contradicts on a clear architectural mismatch
  (e.g. `expected_topology=signal_peptide` but Phobius predicts ≥2 TM
  helices), stays neutral on borderline calls (a single predicted helix, or
  a lipoprotein signal Phobius' Sec-signal model may not catch cleanly).
- **`heme_verdict`** — compares `expected_heme_c` against `heme_motifs`,
  the CXXCH-style motif count `operon.count_heme` already computes for
  *every* hit today (previously purely descriptive output, now also an
  evidence signal). Confirmed if observed ≥ 0.7 × expected (tolerates minor
  regex undercount); contradicted only when the count is entirely absent
  (0) against a substantial expectation (≥4) — a decaheme cytochrome hit
  with zero heme motifs is a strong signal, a partial undercount is not.

**Weakest-link across both signals and all applicable rows**, same principle
as `annot_weight`: 1.0 neutral if no row has an applicable expectation or all
checks come back neutral, 1.2 if every applicable verdict confirms, 0.5 if
any contradicts.

**Why this tier exists**: a standalone hit with no operon neighbors, no
`--bam`, and no BGC/UniOP support falls back to just the raw HMM score under
every other component — no independent corroboration at all. `topo_weight`
is a genuinely orthogonal axis (biophysical structure prediction from
primary sequence, not another statistical summary of the same sequence data
an HMM profile already reads) that applies even to those standalone cases,
since it needs neither operon context nor coverage data. A sequence composition/ML-based classifier layer and a DIAMOND/BLAST
reciprocal-best-hit layer were both considered during design and rejected
for the same reason: both read the same primary-sequence signal the HMM
profile is already built from, so neither adds real independent evidence.

**Gated behind `--topology_check`** (off by default). Even though
`expected_heme_c`/`expected_topology` live in the registry unconditionally
and `heme_motifs` is always computed, `topo_weight` is only ever invoked
when the flag is passed — `cli.py` passes `topo_w=1.0` explicitly otherwise.
This keeps default-run output byte-identical to pre-2026-07-23 behaviour.

---

## Canonical sizes — biological basis

The `canonical_size` for each gene system defines how many genes a **complete**
instance is expected to have. An incomplete cluster (fewer observed genes) gets
a proportional completeness penalty.

| Rule | canonical_size | Biological basis |
|------|---------------|-----------------|
| **FLEET** | 8 | Complete FLEET electron shuttle cassette in *Sideroxydans lithotrophicus* ES-1: EetA, EetB, Ndh2, FmnB, FmnA, DmkA, DmkB, PplA. All 8 co-encoded in the *eet* locus. Barco et al. 2015 (*ISME J*, doi:10.1038/ismej.2014.212). |
| **MAM** | 10 | Conserved core magnetosome membrane proteins: MamA/B/E/K/P/M/Q/I/L/O. Minimal MAM set required for magnetosome biogenesis across Magnetococci and Deltaproteobacteria. Uebe & Schuler 2016 (*Nat Rev Microbiol*, doi:10.1038/nrmicro.2016.99). |
| **FOXABC** | 3 | Three-subunit cytochrome *bc*₁-like complex for acidophilic Fe(II) oxidation: FoxA (Rieske/cytochrome *b*), FoxB (cytochrome *c*₄), FoxC (high-potential Fe-S protein). Always co-encoded in *At. ferrooxidans*. Quatrini et al. 2009 (*BMC Genomics*, doi:10.1186/1471-2164-10-394). |
| **FOXEYZ** | 3 | Three-gene operon for photoferrotrophy: FoxE (reaction center cytochrome *c*), FoxY (flavoprotein), FoxZ (hypothetical). Always co-encoded in *Rhodobacter ferrooxidans* SW2. Croal et al. 2007 (*J Bacteriol*, doi:10.1128/JB.00929-06). |
| **DFE1** | 4 | Four-gene iron reduction operon in *Desulfosporosinus* sp. (DFE_0448–0451). Inherited from FeGenie (Garber et al. 2020). Canonical gene count = list length. |
| **DFE2** | 5 | Five-gene iron reduction operon in *Desulfosporosinus* sp. (DFE_0461–0465). Same rationale. |
| **MtrMto** | 3 | Minimal functional Mtr complex: MtrA–MtrB–MtrC for iron reduction (*Shewanella*; White et al. 2013, *PNAS*, doi:10.1073/pnas.1220074110); or the MtoA–MtoB–ImoA-associated locus in iron oxidizers (*Sideroxydans*; Liu et al. 2012, *Front. Microbiol.*, doi:10.3389/fmicb.2012.00037, characterized MtoA; ImoA renamed/characterized by a 2022 mBio paper, doi:10.1128/mbio.02150-22 — see below). The rule uses the MtrB-family model for the related outer-membrane scaffold. ImoA quinol oxidation is demonstrated; this rule does not establish direct Fe(II) oxidation by ImoA or physiological electron-flow direction. The full rule gene list has 5 entries (covering both systems); canonical_size=3 reflects the core trimer regardless of direction. **Not CymA**: Liu et al. 2012 originally called the oxidation-side gene "CymA-like" from loose NapC/NirT-superfamily sequence similarity; the 2022 paper showed it is only 33% identical to real Shewanella CymA (64%/54% identical to NirT/NapC instead) and renamed it ImoA (Slit_2495). Real CymA (Myers & Myers 1997, *J. Bacteriol.*, doi:10.1128/jb.179.4.1143-1152.1997) is a genuine, separate, well-characterized Shewanella gene — required for Fe(III)/nitrate/fumarate reduction, feeding electrons into MtrA — but it isn't part of this operon rule's oxidation side; see `docs/hmm_library_curation.md` for the full correction writeup and both new models. |
| **SIDERO_TRANSPORT** | — | Not set. Siderophore/heme transport systems are encoded by 2–6 genes depending on ABC transporter architecture and receptor type. Penalizing incomplete transport operons without a reference size would be arbitrary. |
| **SIDERO_SYNTH** | — | Not set. NRPS/PKS siderophore biosynthesis operons range from 3 to 20+ genes. Gene count reflects siderophore chemical complexity, not system incompleteness. |
| **IRON_TRANSPORT** | — | Not set. Same rationale as siderophore transport. |

---

## Per-rule `max_bp_gap` — biological basis

`max_bp_gap` in each rule is the **recommended maximum intergenic gap** for
genes in that system. Currently stored as metadata in the rule definition and
`operon_rules.json`. It is **not yet used for clustering** — the global
`--max_bp_gap` controls the detection window. Future work: per-category
clustering with these values.

The values below are informed estimates. The global `--max_bp_gap` default
of 5000 bp is a MAG-aware detection window. The per-rule values represent a
tighter biological expectation. The distance decay in co_occ_score already
penalizes large gaps within the detection window, so the per-rule max_bp_gap
serves mainly to document expected biology and guide future filtering.

| Rule | max_bp_gap | Rationale |
|------|-----------|-----------|
| **FLEET** | 2000 bp | The *eet* locus in *Sideroxydans* ES-1 spans ~16 kb for 8 genes (~2000 bp average per gene including ORFs). Intergenic gaps are typically 50–300 bp. 2000 bp adds MAG fragmentation tolerance while excluding random co-location across the chromosome. |
| **MAM** | 3000 bp | The mamAB operon in *Magnetospirillum magneticum* AMB-1 spans ~25 kb for ~10 core genes. Gaps are typically 100–500 bp. 3000 bp keeps genes within one operon sub-cluster of the magnetosome island without merging across the entire ~100 kb island. Murat et al. 2010 (*Mol Microbiol*, doi:10.1111/j.1365-2958.2010.07298.x). |
| **FOXABC** | 1000 bp | Compact three-gene cluster in *At. ferrooxidans*: foxA-foxB-foxC spans ~3 kb. Intergenic gaps are 50–200 bp. 1000 bp is generous for MAG work without being permissive. |
| **FOXEYZ** | 1000 bp | Same rationale as FOXABC — compact three-gene operon. |
| **DFE1** | 1000 bp | Four-gene operon expected to be tightly packed. No specific distance data; 1000 bp chosen by analogy with similarly tight iron reduction operons. |
| **DFE2** | 1000 bp | Same as DFE1. |
| **MtrMto** | 2000 bp | In *S. oneidensis* MR-1, mtrCAB spans ~5 kb for 3 genes. Intergenic gaps are ~100–300 bp. 2000 bp accommodates MAG fragmentation and variant architectures (e.g., mtrA and mtrB sometimes separated by mtrD). |
| **SIDERO_TRANSPORT** | 2000 bp | ABC transporter cassette (substrate-binding + permease + ATPase) is typically 3–5 kb. 2000 bp allows for variation in cassette architecture. |
| **SIDERO_SYNTH** | 5000 bp | NRPS/PKS genes are often >10 kb individually; biosynthesis operons can exceed 50 kb. Assembly fragmentation is common. Matches the global `--max_bp_gap` default deliberately — no tighter constraint is biologically justified here. |
| **IRON_TRANSPORT** | 2000 bp | Same as SIDERO_TRANSPORT (most iron transporters are ABC-type). |

---

## Interpretation guide

| cluster_confidence | Interpretation |
|--------------------|---------------|
| ≥ 0.8 | High confidence. Complete or near-complete gene set, tight spatial arrangement, UniOP-supported. |
| 0.5–0.8 | Moderate confidence. Incomplete gene set, or MAG-fragmented, or marginal UniOP support. Report but flag. |
| < 0.5 | Low confidence. Only with `--all_results`. Likely incomplete system or false co-location. |
| > 1.0 | BGC boost applied (siderophore cluster with antiSMASH support). Not yet active. |

The 0.5 threshold is a pragmatic choice. A cluster of 2/3 canonical genes
at typical spacing, all calibrated hits, without UniOP would score:
`hmm_weight=1.0 × co_occ=0.67×0.93×1.0 ≈ 0.62 × uniop=1.0 = 0.62` — above
threshold. Adding an edge penalty drops it to `0.62 × 0.7 ≈ 0.43` — below
threshold. This reflects the biological reality: a partial cluster at a contig
edge is not a confident call.

---

## Implementation

| Function | File | What it computes |
|----------|------|-----------------|
| `distance_decay(gap_bp)` | `scoring.py` | Exponential decay weight for one gap |
| `co_occurrence_score(...)` | `scoring.py` | completeness × distance_factor × edge_factor |
| `uniop_pair_score(...)` | `scoring.py` | min pairwise UniOP probability |
| `hmm_weight(...)` | `scoring.py` | mean calibrated/low_confidence weight |
| `cluster_confidence(...)` | `scoring.py` | final product, capped at 1.2 |
| `topo_weight(...)` | `scoring.py` | weakest-link over topology_verdict + heme_verdict |
| `run_phobius(...)` | `topology.py` | batched Phobius call, one per genome |
| `_parse_phobius_output(...)` | `topology.py` | returns {orf: {signal_peptide, tm_helices}} |
| `topology_verdict(...)` | `topology.py` | confirmed/contradicted/neutral vs expected_topology |
| `heme_verdict(...)` | `topology.py` | confirmed/contradicted/neutral vs expected_heme_c |
| `build_topology_expectation_map(registry)` | `io.py` | {stem → {expected_heme_c, expected_topology}} |
| `build_canonical_size_map(rules)` | `operon.py` | {gene_name → canonical_size} from rules |
| `_parse_uniop_pred(...)` | `uniop.py` | returns (orf_to_op, pair_probs) |
| `run_uniop(...)` | `uniop.py` | returns (genome_operon_map, genome_pair_probs) |
