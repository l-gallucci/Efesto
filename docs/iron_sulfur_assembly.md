# Iron-Sulfur Cluster Assembly — HMM Models and Biology

Efesto includes four HMMs targeting the bacterial iron-sulfur (Fe-S) cluster
assembly machinery, placed in the `iron_sulfur_assembly` category. This document covers
the biology of Fe-S assembly systems, the rationale for model selection, and guidance
for interpreting results.

---

## Biology

### What are Fe-S clusters?

Iron-sulfur ([Fe-S]) clusters — [2Fe-2S], [3Fe-4S], [4Fe-4S] — are ancient, ubiquitous
cofactors assembled from iron and inorganic sulfide. They are required by:

- **Electron transfer chains** (ferredoxins, respiratory complexes I and II)
- **Enzyme catalysis** (aconitase, endonuclease III, radical SAM enzymes)
- **Regulatory sensing** (IscR, FNR, NsrR, SoxR — all use Fe-S clusters as sensors)
- **DNA repair** (MutY, Nth, DNA primase)
- **Nitrogen fixation** (nitrogenase FeMoco assembly via the NIF system)

Because Fe-S cluster assembly requires both iron and sulfide, the machinery is
directly coupled to cellular iron availability. Under iron limitation or oxidative
stress, organisms switch assembly pathways — making Fe-S assembly genes reliable
biomarkers of iron status in environmental metagenomes.

### Three assembly systems in bacteria

| System | Operon | Distribution | Induction | Iron-stress signal? |
|--------|--------|-------------|-----------|---------------------|
| **ISC** | *iscRSUA-hscBA-fdx-iscX* | Primary housekeeping system; γ-/β-Proteobacteria, most gram-negatives | Repressed by IscR when Fe-S sufficient | Low — constitutive |
| **SUF** | *sufABCDSE* | Near-universal; sole system in many gram-positives and cyanobacteria | Induced by OxyR (H₂O₂), Fur (Fe starvation), SufR | **High — direct iron-stress marker** |
| **NIF** | *nifU, nifS, nifX, nifY...* | Restricted to diazotrophs | Co-regulated with nitrogenase | Nitrogen-fixation specific |

**Key ecological insight:** SUF is preferentially expressed under iron limitation
and oxidative stress. In iron-poor environments (oligotrophic oceans, iron-limited
groundwaters, host infection), SUF enrichment in metagenomes indicates that the
microbial community is experiencing Fe-stress and has activated the stress-tolerant
assembly pathway. ISC enrichment indicates housekeeping Fe-S metabolism without
iron stress. Detecting both allows distinction of these states.

---

## Model organism evidence

### *Escherichia coli* K-12 — ISC primary, SUF backup

The most studied dual-system organism. ISC is the primary housekeeping system.
SUF is the backup: it is induced by the Fur repressor (iron depletion) and OxyR
(H₂O₂), and its products are more resistant to oxidative damage than ISC
components. The *sufA* promoter contains a Fur box; iron chelation with
dipyridyl induces *suf* ~10-fold. Deletion of *isc* causes Fe-S protein
assembly defects; deletion of *suf* causes sensitivity to H₂O₂ and iron
chelation. Double *isc suf* deletion is lethal.

**References:**
- Outten FW, Djaman O, Storz G (2004). A *suf* operon requirement for Fe-S
  cluster assembly during iron starvation in *Escherichia coli*.
  *Mol Microbiol* 52:861–72.
  [doi:10.1111/j.1365-2958.2004.04025.x](https://doi.org/10.1111/j.1365-2958.2004.04025.x)
- Py B, Barras F (2010). Building Fe/S proteins: bacterial strategies.
  *Nat Rev Microbiol* 8:436–46.
  [doi:10.1038/nrmicro2356](https://doi.org/10.1038/nrmicro2356)

### *Bacillus subtilis* 168 — SUF only

*B. subtilis* has no ISC system. SUF (with SufU as the functional analog of
IscU) is the sole Fe-S assembly pathway and is essential for viability.
SufU is a zinc-containing scaffold with distinct mechanism from *E. coli* IscU.

**Reference:**
- Albrecht AG, Netz DJ, Miethke M, et al. (2010). SufU is an essential iron-sulfur
  cluster scaffold protein in *Bacillus subtilis*.
  *J Bacteriol* 192:1643–51.
  [doi:10.1128/JB.01536-09](https://doi.org/10.1128/JB.01536-09)

### *Synechocystis* sp. PCC 6803 — SUF essential for photosynthesis

Cyanobacteria encode SUF but not ISC. SUF is essential for assembly of the
multiple [4Fe-4S] clusters in photosystem I (PsaA/PsaB) and the [2Fe-2S]
cluster of the Rieske protein in the cytochrome *b*6*f* complex. Given that
photosynthesis generates reactive oxygen species, SUF's oxidative stress
tolerance is essential in oxygenic photosynthesizers.

**Reference:**
- Balasubramanian R, Shen G, Bryant DA, Golbeck JH (2006). Regulatory roles
  for IscA and SufA in iron homeostasis and redox stress responses in the
  cyanobacterium *Synechococcus* sp. strain PCC 7002.
  *J Bacteriol* 188:3182–91.
  [doi:10.1128/JB.188.9.3182-3191.2006](https://doi.org/10.1128/JB.188.9.3182-3191.2006)

### *Azotobacter vinelandii* — ISC + NIF

Encodes all three systems. The NIF system (NifU, NifS) is dedicated to
assembling the FeMoco and P-cluster of nitrogenase; ISC handles housekeeping
Fe-S assembly for the rest of the proteome. NIF genes are co-regulated with
the *nif* nitrogen-fixation operon and are expressed only under nitrogen-fixing
conditions.

### *Desulfovibrio vulgaris* Hildenborough — ISC + SUF, Fe-S proteome-dominant

Sulfate-reducing bacteria maintain both systems. Fe-S proteins dominate their
proteome (ferredoxins, hydrogenases, dissimilatory sulfite reductase). Both
ISC and SUF are expressed constitutively under standard anaerobic growth.

---

## SUF system — molecular mechanism

The SUF system operates through a dedicated sulfur-delivery and scaffold pathway
that is more resistant to oxidative inactivation than ISC:

```
cysteine  ──→  SufS (cysteine desulfurase)  ──→  persulfide on SufE
                                                        │
                                                        ▼
SufB-SufC-SufD complex (scaffold + ATPase)  ←── sulfur transfer
        │
        │  iron delivered (mechanism unclear; SufA may serve as iron chaperone)
        ▼
   [Fe-S] assembled on SufB scaffold
        │
        ▼
   transfer to target apo-protein via SufA
```

- **SufS** is a pyridoxal-phosphate (PLP)-dependent cysteine desulfurase. It
  is structurally related to IscS but forms a tighter complex with SufE rather
  than delivering sulfur directly to the scaffold.
- **SufE** is a sulfur-transfer intermediate that receives the persulfide from
  SufS and delivers it to SufB-SufC-SufD.
- **SufB-SufC-SufD** form an ATP-hydrolysing scaffold complex. SufC is the
  ATPase; SufB is the primary scaffold; SufD stabilises the complex.
- **SufA** (related to IscA) serves as an alternate scaffold or [Fe-S] carrier
  for delivery to client proteins.

**Reference:**
- Fontecave M, Ollagnier de Choudens S, Py B, Barras F (2005). Mechanisms of
  iron-sulfur cluster assembly: the SUF machinery.
  *J Biol Inorg Chem* 10:713–21.
  [doi:10.1007/s00775-005-0025-1](https://doi.org/10.1007/s00775-005-0025-1)
- Pérard J, Ollagnier de Choudens S (2017). Iron-sulfur clusters biogenesis by
  the SUF machinery: close to the molecular mechanism understanding.
  *J Biol Inorg Chem* 23:581–596.
  [doi:10.1007/s00775-017-1527-3](https://doi.org/10.1007/s00775-017-1527-3)

---

## ISC system — molecular mechanism

```
cysteine  ──→  IscS (cysteine desulfurase)  ──→  [Fe-S] assembled on IscU scaffold
                                                        │
                               iron via CyaY (bacterioferritin)?  ──→  IscU
                                                        │
                                                        ▼
                            HscA/HscB (Hsp70/Hsp40 chaperones) — trigger cluster release
                                                        │
                                                        ▼
                                          IscA / Fdx — transfer to target protein
```

- **IscS** delivers sulfur to IscU and is a key sulfur source for thio-modifications
  of tRNA as well as Fe-S cluster assembly.
- **IscU** is the primary scaffold; it binds three cysteine residues and assembles
  the nascent cluster.
- **HscA-HscB** are dedicated chaperones that bind holo-IscU and use ATP hydrolysis
  to release the cluster for downstream transfer.
- **IscR** (a [2Fe-2S]-containing regulator, not included in `iron_sulfur_assembly`)
  represses the *isc* operon when Fe-S clusters are sufficient; it autoregulates
  assembly based on cluster occupancy.

**Reference:**
- Roche B, Aussel L, Ezraty B, et al. (2013). Iron/sulfur proteins biogenesis in
  prokaryotes: formation, regulation and diversity.
  *Biochim Biophys Acta* 1827:455–69.
  [doi:10.1016/j.bbabio.2012.12.010](https://doi.org/10.1016/j.bbabio.2012.12.010)
- Bandyopadhyay S, Chandramouli K, Johnson MK (2008). Iron-sulfur cluster biosynthesis.
  *Biochem Soc Trans* 36:1112–9.
  [doi:10.1042/BST0361112](https://doi.org/10.1042/BST0361112)

---

## Active HMM models in Efesto

| Stem | Source | nseq | Cutoff (bits) | Target gene | System |
|------|--------|------|---------------|-------------|--------|
| `sufC_TIGR01978.1` | NCBIfam/TIGRFAM | ≥ 100 | calibrated | SufC — ATPase subunit of SufB-SufC-SufD scaffold complex | SUF |
| `sufB_TIGR01980.1` | NCBIfam/TIGRFAM | ≥ 100 | calibrated | SufB — scaffold protein, primary Fe-S assembly site | SUF |
| `sufS_TIGR01979.1` | NCBIfam/TIGRFAM | ≥ 100 | calibrated | SufS — cysteine desulfurase (sulfur donor to SufE) | SUF |
| `IscS_TIGR02006.1` | NCBIfam/TIGRFAM | ≥ 100 | calibrated | IscS — cysteine desulfurase (housekeeping, sulfur donor to IscU) | ISC |

All four models are from the NCBIfam/TIGRFAM database with curated bitscore cutoffs
(GA/TC/NC). They are searched with calibrated thresholds — no zero-cutoff fallback.

### Why these four and not others?

**Included:**
- `sufB`, `sufC`, `sufS` — the core SUF scaffold complex + sulfur donor. This triad
  is the minimal diagnostic set for the SUF system. SufB and SufC are almost always
  co-encoded; SufS is the sulfur entry point. Three models are sufficient to call
  a SUF operon with high confidence.
- `IscS` — the housekeeping cysteine desulfurase. It is the most conserved and
  diagnostic ISC component; it is also the sulfur source for tRNA thio-modifications
  even in organisms that lack ISC for Fe-S assembly, making it a reliable presence
  marker.

**Intentionally not included (yet):**
- `sufD`, `sufE`, `sufA` — diagnostic but highly redundant with `sufB`/`sufC` in
  operon context; adding them would inflate hit counts without improving specificity.
- `iscU`, `iscA`, `hscA`, `hscB`, `fdx` — ISC accessory components. IscU is the
  primary ISC scaffold but has no calibrated TIGRFAM model with sufficient coverage
  for metagenomics. Will be added when a suitable model is identified.
- NIF system (`nifU`, `nifS`) — nitrogen-fixation specific; overlap with ISC models
  is problematic. Will be added as a separate `nitrogen_fixation` category.

---

## Interpreting results

### Category: `iron_sulfur_assembly`

Models in this category are placed under `report_all_categories` — they pass
operon-context filtering unconditionally. Every cluster containing an Fe-S
assembly HMM hit is reported regardless of co-occurrence with other genes.
This is correct: SUF and ISC operons are well-defined systems where even a
single-gene hit is informative.

### Ecological interpretation

| Observation | Interpretation |
|---|---|
| Only `sufB`/`sufC`/`sufS` hits, no ISC | SUF-only organism (*B. subtilis*, cyanobacteria); or ISC-lacking anaerobe |
| `sufB/C/S` + `IscS` in same genome | Dual-system organism under potential Fe-stress (induced SUF) |
| High SUF / total Fe-S ratio across metagenome | Community experiencing iron limitation or oxidative stress |
| `IscS` only, no SUF hits | Housekeeping Fe-S assembly without iron stress; ISC-only γ-/β-Proteobacteria |
| SUF enriched in iron-poor biome | Consistent with Fe-S stress biomarker interpretation |

### Co-occurrence with regulatory models

Iron-sulfur assembly is tightly co-regulated with iron homeostasis:
- **IscR** (`iron_gene_regulation`) represses the ISC operon and induces SUF when
  clusters are limiting. IscR hits co-located with `IscS` indicate Fe-S-sensing
  feedback regulation.
- **Fur** (`iron_gene_regulation`) represses the SUF operon under iron sufficiency.
  Fur + SUF co-occurrence confirms iron-responsive regulation.
- **Flavodoxins** (`iron_stress`) substitute for ferredoxin (an Fe-S protein) under
  iron limitation. SUF + flavodoxin co-occurrence is strong evidence for active
  iron-stress response.

### Cluster confidence and `model_nseq`

All four Fe-S assembly models have `nseq ≥ 100` and calibrated bitscore cutoffs.
They contribute `hmm_weight = 1.0` (calibrated tier) to the cluster confidence score.
`model_nseq` in `results-long.tsv` reflects the training set size; values ≥ 100
indicate well-supported models.

---

## References

1. Outten FW, Djaman O, Storz G (2004). A *suf* operon requirement for Fe-S cluster
   assembly during iron starvation in *E. coli*. *Mol Microbiol* 52:861–72.
   [doi:10.1111/j.1365-2958.2004.04025.x](https://doi.org/10.1111/j.1365-2958.2004.04025.x)

2. Fontecave M, Ollagnier de Choudens S, Py B, Barras F (2005). Mechanisms of
   iron-sulfur cluster assembly: the SUF machinery. *J Biol Inorg Chem* 10:713–21.
   [doi:10.1007/s00775-005-0025-1](https://doi.org/10.1007/s00775-005-0025-1)

3. Pérard J, Ollagnier de Choudens S (2017). Iron-sulfur clusters biogenesis by
   the SUF machinery. *J Biol Inorg Chem* 23:581–596.
   [doi:10.1007/s00775-017-1527-3](https://doi.org/10.1007/s00775-017-1527-3)

4. Py B, Barras F (2010). Building Fe/S proteins: bacterial strategies.
   *Nat Rev Microbiol* 8:436–46.
   [doi:10.1038/nrmicro2356](https://doi.org/10.1038/nrmicro2356)

5. Roche B, Aussel L, Ezraty B, et al. (2013). Iron/sulfur proteins biogenesis
   in prokaryotes. *Biochim Biophys Acta* 1827:455–69.
   [doi:10.1016/j.bbabio.2012.12.010](https://doi.org/10.1016/j.bbabio.2012.12.010)

6. Bandyopadhyay S, Chandramouli K, Johnson MK (2008). Iron-sulfur cluster biosynthesis.
   *Biochem Soc Trans* 36:1112–9.
   [doi:10.1042/BST0361112](https://doi.org/10.1042/BST0361112)

7. Esquilin-Lebron K, Dubrac S, Barras F, Boyd JM (2021). Bacterial Approaches for
   Assembling Iron-Sulfur Proteins. *mBio* 12:e0242521.
   [doi:10.1128/mBio.02425-21](https://doi.org/10.1128/mBio.02425-21)
