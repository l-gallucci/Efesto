# Publication Plan — MetalGenie-Evo

*Created: 2026-05-23*

---

## Target journal

**mSystems** (ASM) — first choice.
Rationale: accepts software + database as co-equal contributions; broad environmental
microbiology audience; no word limit; explicitly accepts metagenome case studies;
FeGenie was published in *Front. Microbiol.* — stepping up in venue is appropriate.

Alternative: *Nucleic Acids Research* database issue (if framed as a database resource).

---

## Framing — the ecological question that drives everything

> *"In a metagenome from a redox interface — where iron oxidation and reduction coexist,
> where iron limits primary production, where metal resistance co-selects with iron cycling
> — which organisms are doing what, and with what genes?"*

**Do NOT frame as "FeGenie but better."** FeGenie cannot answer this question.
MetalGenie-Evo can. Lead with the ecological question.

Suggested title direction:
> *MetalGenie-Evo: an expanded HMM library and metagenome-aware pipeline for
> dissimilatory metal cycling annotation*

---

## Paper structure

| Section | Content | Status |
|---------|---------|--------|
| **Introduction** | Gap: FeGenie cross-hits (MtrA/MtoA), no metagenome support, category errors, no metal resistance integration | Documented in curation MD |
| **Library curation** | Layer A/B methodology, category audit, 17 new models, MtrA/MtoA calibration | Fully documented |
| **Software** | Architecture: Pyrodigal-GV, operon detection (UniOP), TPM normalization, Anvi'o output | Code complete; needs architecture figure |
| **Benchmark — precision/recall** | MetalGenie-Evo vs FeGenie on reference genomes with ground truth | **MISSING — critical** |
| **Benchmark — MtrA/MtoA** | Cross-hit rate with FeGenie models; disambiguation with new models | Calibration report exists; needs tabulated FP comparison |
| **Case study** | ≥2 environmental metagenomes; show biological interpretation improvement | **MISSING — critical** |
| **Performance** | Runtime / memory across dataset sizes | Missing — simple to add |

---

## Priority tasks for publication readiness

### 1. Reference genome benchmark *(highest priority)*

Select 15–20 genomes with known metal cycling phenotype:

**Iron reducers:**
- *Shewanella oneidensis* MR-1
- *Geobacter sulfurreducens* PCA
- *Desulfovibrio vulgaris* Hildenborough

**Iron oxidizers:**
- *Acidithiobacillus ferrooxidans* ATCC 23270 (acidophilic — rusticyanin + Cyt579 test)
- *Sideroxydans lithotrophicus* ES-1 (MtoA canonical)
- *Gallionella capsiferriformans* ES-2 (MtoA canonical)
- *Leptospirillum ferrooxidans* (Cyt579 homolog, no rusticyanin)

**Metal resistance positive controls:**
- *Cupriavidus metallidurans* CH34 (multi-metal resistance)

**Negative controls:**
- *Escherichia coli* K-12
- *Bacillus subtilis* 168

Run FeGenie (original) and MetalGenie-Evo on same set. Compare:
- False positive rate — especially MtrA/MtoA cross-annotation
- Category correctness (iron transport miscalled as resistance)
- Recall on genes with known function

### 2. MtrA/MtoA false positive table

From existing `hmm_library/_calibration/mtr_mto/calibration_report.tsv`:
- N sequences misannotated by FeGenie model (scoring above 140)
- N correctly excluded by MetalGenie-Evo (below GA=580 / GA=520)
- Specific examples: MtrD family, betaproteobacteria DmsE

Already have the data — needs summarizing into a table (~1 day).

### 3. Environmental case study *(required for mSystems)*

Two contrasting environments:
- **Acid mine drainage** (SRA has many — e.g., Río Tinto): rusticyanin/Cyt579
  validation; expected high iron oxidation signal
- **Iron-rich groundwater or freshwater sediment**: MtrA vs MtoA separation
  matters most; iron stress markers (flavodoxin) expected if iron-limited fractions

**Key story to find:** one MAG/bin where FeGenie calls contradictory/ambiguous
MtrA+MtoA and MetalGenie-Evo correctly resolves them by score + operon context.
This is the concrete empirical validation of the cryptic cycling prediction
(Díaz-González et al. 2025 mSystems).

### 4. Runtime benchmark

Single script: run on 1, 10, 50, 100, 500 genomes. Report wall time + peak memory.
~0.5 day of work.

### 5. Figures needed

| Figure | Content |
|--------|---------|
| Fig 1 | Library composition — bar chart by category (active models, model sources) |
| Fig 2 | MtrA/MtoA score landscape — violin/dot plot across 3794-seq universe per class |
| Fig 3 | Software workflow diagram (input → gene calling → HMM search → operon → output) |
| Fig 4 | Benchmark precision/recall heatmap: MetalGenie-Evo vs FeGenie per category |
| Fig S1 | Category reclassification sankey or table |

### 6. Fix `validated_in` column

Fill in reference genomes where each model fires (from benchmark step 1).
Reviewers will ask; the column being empty is a liability.

---

## Suggested timeline

```
Week 1    MtrA/MtoA summary table from calibration_report.tsv
Week 1-2  Reference genome benchmark (15-20 genomes, FeGenie comparison)
Week 2    Fig 1 (library stats) + Fig 2 (score landscape) + runtime benchmark
Week 3-6  Environmental case study (pick SRA datasets, run, interpret)
Week 6-7  Fig 3 (workflow) + Fig 4 (benchmark) + fill validated_in
Week 8-10 Write
```

---

## Novelty claims (ranked by defensibility)

### 1 — MtrA/MtoA disambiguation (sharpest claim)

FeGenie's own developers flagged this as unresolved. Zero other tools address it.

Evidence in hand:
- Cross-hits quantified: 1.9× / 2.1× the calibrated cutoff
- Calibrated against 3794-sequence universe
- Score landscapes: distinct classes separated
- TC/GA/NC values with biological rationale per boundary

Iron reduction vs iron oxidation is not a subtle ambiguity — it is a 180° functional
difference. Any metagenome study using FeGenie models on mixed-redox environments
carries this error.

### 2 — First unified iron cycling + metal resistance pipeline

| Tool | Iron cycling | Metal resistance | Metagenome workflow | Operon context |
|------|-------------|-----------------|--------------------|----|
| FeGenie 2020 | ✅ | ❌ | ❌ genome-designed | basic |
| MetHMMDB 2025 | ❌ | ✅ | ❌ | ❌ |
| **MetalGenie-Evo** | ✅ | ✅ | ✅ | ✅ UniOP |

No published tool spans all four. The combination matters:
metal resistance without iron cycling context misses co-selection dynamics;
iron cycling without metal resistance misses organisms that do both.

### 3 — Flavodoxin:ferredoxin ratio as community iron stress index

LaRoche et al. 1996 established flavodoxin as in situ iron stress marker in
phytoplankton. No HMM pipeline has operationalized this for prokaryotic metagenomes.

MetalGenie-Evo has both models in `iron_stress`, correct category, outputs
per-sample counts. Enables a novel quantitative iron stress metric from
metagenomes. High value for oceanography, limnology, soil science.

### 4 — Complete acidophilic Fe(II)-oxidation marker set

FeGenie: Cyc2 only.
MetalGenie-Evo: Cyc2 + Rusticyanin + Cyt579 → full downhill electron chain.

In *At. ferrooxidans*, rusticyanin is the most abundant periplasmic protein
(~350 mg/mL). Annotating Cyc2 alone in AMD metagenomes underestimates acidophilic
iron oxidation capacity. First library covering the complete diagnostic set.

### 5 — Fe-S stress as ecological readout

SUF (sufB/sufC/sufS) vs ISC (IscS) ratio in metagenomes:
- SUF enrichment → organisms prepared for oxidative stress / iron limitation
- ISC dominance → housekeeping, iron-replete conditions

No existing tool separates these as distinct ecological signals. Application:
iron-limited ocean gyres, AMD oxidative gradients, soil redox transitions.

### 6 — Transparent, versioned, machine-readable library (FAIR)

FeGenie models: GitHub flat files, no provenance.
MetalGenie-Evo registry: per-model source, nseq, cutoff basis, added_date,
reference DOI, validated_in. FAIR data artifact citable independently of software.
Increasingly required by journals (NAR, mSystems data availability standards).

---

## Known flaws to address before submission

### Critical

| Flaw | Fix |
|------|-----|
| No benchmark — all claims are assertions | Run reference genome set (task 1 above) |
| MtoA nseq=6 → low sensitivity | State explicitly as limitation; model will improve as Gallionellaceae MAGs accumulate |
| Calibration universe is title-searched (circular) | Supplement with manually searched MtrD sequences (already done); discuss limitation in methods |
| MetHMMDB is preprint (115 models) | Acknowledge as limitation; track if published before submission |
| `iron_acquisition` typo (inherited from FeGenie) | Fix in category names, file paths, outputs — do before benchmark so paths don't change |

### Moderate

| Flaw | Fix |
|------|-----|
| Siderophore synthesis = 34% of active models (library imbalance) | Add narrative: siderophores are the numerically largest iron gene family in nature; imbalance reflects biology |
| No Archaea coverage | State as limitation; scope is bacteria |
| No Chlorobi / photoferrotrophic iron oxidation coverage | State as limitation; add FoxE/FoxY/FoxZ models in next version |
| Dissimilatory vs assimilatory iron reduction not distinguished | Add note in annotation guidance; future direction |
| `validated_in` empty for all new models | Fill from benchmark run |
| Source heterogeneity (4 sources, different training philosophies) | Document in methods; Layer A/B curation as harmonization step |

---

## The single experiment that would make the paper memorable

Find a MAG or enrichment metagenome bin where:
- FeGenie calls contradictory MtrA + MtoA (artifact of cross-hits)
- MetalGenie-Evo correctly assigns both by score margin + operon context
- The organism is demonstrably a dual-capacity iron cycler (cryptic cycling candidate)

This is what Díaz-González et al. 2025 (*mSystems*) predicted theoretically.
MetalGenie-Evo can make it empirical.
