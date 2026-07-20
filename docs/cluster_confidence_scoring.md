# Cluster Confidence Scoring — Design Rationale

*Created: 2026-05-24*

---

## The formula

```
cluster_confidence = hmm_weight × co_occ_score × uniop_weight
                      × bgc_boost × annot_weight × struct_weight
```

Capped at 1.2. Recommended call threshold: ≥ 0.5.

Six independent evidence sources multiplied together. Any single low score
brings the whole cluster down. This is intentional: a gene system with perfect
HMM hits but no spatial coherence is as suspicious as one with perfect spacing
but low-confidence HMMs.

`annot_weight` and `struct_weight` (added 2026-07-20) only ever affect
clusters containing a hit on a model flagged `needs_confirmation` in the
registry — everything else always sees 1.0 for both, identical to pre-2026-07-20
behaviour.

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

Designed, not yet implemented — placeholder parameter in `cluster_confidence`,
always 1.0. Intended as the same confirm/contradict/neutral mechanism as
`annot_weight`, but via Baktfold (ProstT5→Foldseek structural search) instead
of eggNOG, for the flagged hits that eggNOG can't resolve either.

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
| **MtrMto** | 3 | Minimal functional Mtr complex: MtrA–MtrB–MtrC for iron reduction (*Shewanella*; White et al. 2013, *PNAS*, doi:10.1073/pnas.1222358110); or MtoA–MtrB–CymA for iron oxidation (*Sideroxydans*; Liu et al. 2012, *Front. Microbiol.*, doi:10.3389/fmicb.2012.00037). Both share MtrB as the outer membrane scaffold. The full rule gene list has 5 entries (covering both systems); canonical_size=3 reflects the core trimer regardless of direction. |
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
| `build_canonical_size_map(rules)` | `operon.py` | {gene_name → canonical_size} from rules |
| `_parse_uniop_pred(...)` | `uniop.py` | returns (orf_to_op, pair_probs) |
| `run_uniop(...)` | `uniop.py` | returns (genome_operon_map, genome_pair_probs) |
