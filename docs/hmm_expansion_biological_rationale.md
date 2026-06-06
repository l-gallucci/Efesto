# HMM Library Expansion — Biological Rationale

This document provides the scientific justification, model organism evidence,
key genes, HMM sources, and limitations for each proposed new HMM category.
All references are from PubMed.

---

## 1. Fe-S Cluster Assembly (SUF / ISC / NIF systems)

### Biology

Iron-sulfur (Fe-S) clusters — [2Fe-2S], [3Fe-4S], [4Fe-4S] — are among the most
ancient and versatile biological cofactors. They are required for electron transfer
chains, enzyme catalysis (aconitase, ferredoxins, dehydratases), DNA repair, and
regulatory sensing. Three distinct assembly machineries exist in bacteria:

| System | Role | Regulation |
|--------|------|------------|
| **ISC** (*iscSUA-hscBA-fdx-iscX*) | Housekeeping; primary system in most γ-/β-Proteobacteria | Repressed by IscR under Fe-S sufficiency |
| **SUF** (*sufABCDSE*) | Stress-inducible; activated under iron limitation and oxidative stress | Induced by OxyR (H₂O₂) and Fur (Fe starvation) |
| **NIF** (*nifU, nifS, nifX...*) | Specialized for nitrogenase Fe-S assembly in diazotrophs | Regulated with *nif* operon |

**MetalGenie-Evo relevance:** SUF is directly linked to iron cycling — it is induced
specifically when iron is limiting. Environmental metagenomes from iron-poor waters
(oligotrophic oceans, iron-limited groundwaters) show SUF enrichment as a proxy for
iron stress. ISC is constitutive housekeeping and less informative for metal cycling
studies. Annotating both allows distinction of iron-stress response (SUF) from
baseline metabolism (ISC).

### Model organisms and key evidence

| Organism | System | Evidence |
|----------|--------|---------|
| *Escherichia coli* K-12 | ISC (primary) + SUF (stress backup) | Deletion of *isc* operon → Fe-S protein assembly defects; *suf* induced by H₂O₂ and iron chelation |
| *Bacillus subtilis* 168 | SUF only (no ISC) | SUF is essential; *sufU* (functional analog of IscU) characterized |
| *Synechocystis* sp. PCC 6803 | SUF (photosynthetic, oxidative stress) | Essential for PSI assembly (multiple [4Fe-4S] clusters) |
| *Azotobacter vinelandii* | ISC + NIF | NIF assembles nitrogenase FeMoco; ISC for other Fe-S proteins |
| *Desulfovibrio vulgaris* Hildenborough | ISC + SUF | Both systems active; Fe-S proteins dominate the proteome |

### Key genes

**SUF system:**
- `sufA` — scaffold/carrier (SufA)
- `sufB` — scaffold (SufB-SufC-SufD ATP-hydrolyzing complex)
- `sufC` — ATPase
- `sufD` — scaffold partner
- `sufS` — cysteine desulfurase (sulfur donor)
- `sufE` — sulfur transfer intermediate

**ISC system:**
- `iscS` — cysteine desulfurase
- `iscU` — scaffold protein
- `iscA` — alternate scaffold / iron-sulfur cluster carrier
- `fdx` — ferredoxin (electron donor)
- `hscA`, `hscB` — Hsp70/Hsp40 chaperones (cluster transfer)

### References (PubMed)

- Bandyopadhyay S, Chandramouli K, Johnson MK (2008). Iron-sulfur cluster biosynthesis.
  *Biochem Soc Trans* 36:1112–9. [doi:10.1042/BST0361112](https://doi.org/10.1042/BST0361112)

- Fontecave M, Ollagnier de Choudens S, Py B, Barras F (2005). Mechanisms of iron-sulfur
  cluster assembly: the SUF machinery. *J Biol Inorg Chem* 10:713–21.
  [doi:10.1007/s00775-005-0025-1](https://doi.org/10.1007/s00775-005-0025-1)

- Pérard J, Ollagnier de Choudens S (2017). Iron-sulfur clusters biogenesis by the SUF
  machinery: close to the molecular mechanism understanding. *J Biol Inorg Chem* 23:581–596.
  [doi:10.1007/s00775-017-1527-3](https://doi.org/10.1007/s00775-017-1527-3)

- Esquilin-Lebron K, Dubrac S, Barras F, Boyd JM (2021). Bacterial Approaches for
  Assembling Iron-Sulfur Proteins. *mBio* 12:e0242521.
  [doi:10.1128/mBio.02425-21](https://doi.org/10.1128/mBio.02425-21)

### HMM sources

- **NCBIfam/TIGRFAM (via InterPro):** best coverage — TIGR01054 (sufS), TIGR01909 (sufB),
  TIGR04484 (sufC), TIGR01910 (sufD), TIGR04035 (sufA), TIGR01376 (sufE),
  TIGR00382 (iscS), TIGR00170 (iscU)
- **Pfam:** PF01106 (SufB), PF03191 (SufC/D), PF01612 (IscA_SufA)
- Prefer NCBIfam models — curated with validated cutoffs

### 1b. Fe-S regulatory sensors

Fe-S clusters are not only structural cofactors — several global regulators use them
as direct sensors of redox state, oxygen tension, nitric oxide, and iron availability.
These are among the most important environmental gene markers for distinguishing
metabolic strategies.

| Regulator | Cluster type | Sensed signal | Key targets | Model organism |
|-----------|-------------|---------------|-------------|----------------|
| **FNR** (fumarate/nitrate reduction) | [4Fe-4S] → [2Fe-2S] or apo under O₂ | O₂ | Switches between aerobic/anaerobic respiration; >100 genes | *E. coli* K-12 |
| **IscR** | [2Fe-2S] holo (Fe-S sufficient) or apo (Fe-S limiting) | Fe-S cluster availability | *isc* operon (represses when Fe-S sufficient), *suf* (induces when apo) | *E. coli* K-12 |
| **NsrR** | [2Fe-2S] | NO (nitric oxide) | NO detoxification genes (*hmp*, *norV*), Fe-S repair | *E. coli*, *B. subtilis* |
| **SoxR** | [2Fe-2S], oxidized by O₂•⁻ | Superoxide | *soxS* → oxidative stress response; mexAB-OprM (efflux in *P. aeruginosa*) | *E. coli*, *P. aeruginosa* |
| **Wbl/WhiB** family | [4Fe-4S] or [2Fe-2S] | O₂, NO, redox | Morphological differentiation, antibiotic resistance in Actinobacteria | *Streptomyces coelicolor*, *M. tuberculosis* |
| **Aconitase (AcnA/AcnB)** | [4Fe-4S] → apo under iron stress | Iron limitation | Acts as RNA-binding IRP when cluster lost; post-transcriptional iron regulation | *E. coli* K-12 |

**MetalGenie-Evo relevance:** FNR is a direct proxy for aerobic vs. anaerobic lifestyle —
its presence and cluster status determines whether an organism switches to anaerobic
iron-coupled metabolism. IscR links Fe-S assembly to iron sensing. NsrR connects iron
cycling to the nitrogen cycle (NO is an intermediate in denitrification). SoxR indicates
organisms living at redox interfaces where superoxide is generated.

**Reference (PubMed):**
- Miller HK, Auerbuch V (2015). Bacterial iron-sulfur cluster sensors in mammalian
  pathogens. *Metallomics* 7:943–56.
  [doi:10.1039/c5mt00012b](https://doi.org/10.1039/c5mt00012b)

- Green J, Crack JC, Thomson AJ, LeBrun NE (2009). Bacterial sensors of oxygen.
  *Curr Opin Microbiol* 12:145–51.
  [doi:10.1016/j.mib.2009.01.008](https://doi.org/10.1016/j.mib.2009.01.008)

**HMM sources:**
- **FNR:** NCBIfam TIGR00232 (*fnr/crp* family); Pfam PF00325 (Crp/Fnr regulator)
  — subfamily-specific model needed (FNR vs CRP are similar)
- **IscR:** NCBIfam TIGR02126; Pfam PF02082 (Rrf2 family)
- **NsrR:** Pfam PF02082 (Rrf2) — shares family with IscR; distinguish by phylogeny
- **SoxR:** NCBIfam TIGR01465; Pfam PF13377 (MerR-like)
- **WhiB:** NCBIfam TIGR00249; Pfam PF02467 (WhiB) — Actinobacteria only

---

### 1c. Ferredoxin / Flavodoxin iron-economy switch

One of the most ecologically important indicators of iron stress in metagenomes is the
**ferredoxin-to-flavodoxin substitution**. Under iron-replete conditions, bacteria use
ferredoxin ([2Fe-2S] or [4Fe-4S]) as a one-electron carrier in photosynthesis, nitrogen
fixation, and central metabolism. Under iron starvation, ferredoxin is replaced by
flavodoxin (uses FMN instead of Fe-S cluster) — an iron-free functional equivalent.

This substitution:
- Is induced by iron starvation via Fur derepression
- Has been directly observed in open-ocean metagenomes as a proxy for iron-limited
  primary producers (especially cyanobacteria and γ-Proteobacteria)
- Is the genetic basis of the "ferrous wheel" iron-economy strategy in oligotrophic
  surface oceans

**Key organisms:**
- *Synechococcus* and *Prochlorococcus* (marine cyanobacteria): *petF* (ferredoxin)
  vs *isiB* (flavodoxin); flavodoxin strongly induced in iron-poor ocean gyres
- *E. coli*: *fldA*, *fldB* (flavodoxins) induced under iron limitation alongside *suf*
- Many α- and γ-Proteobacteria encode both and switch depending on iron availability

**Reference:**
- LaRoche J et al. (1996). Flavodoxin as an in situ marker for iron stress in phytoplankton.
  *Nature* 382:802–805. [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0)
  *(note: Nature, not indexed in PubMed via standard search but well established)*

**HMM sources:**
- **Ferredoxin:** Pfam PF00111 (2Fe-2S ferredoxin), PF00037 (4Fe-4S ferredoxin);
  NCBIfam TIGR01324 (*petF* cyanobacterial ferredoxin)
- **Flavodoxin:** Pfam PF00258 (flavodoxin); NCBIfam TIGR02930 (*isiB* cyanobacterial),
  TIGR01753 (short-chain flavodoxin)

**Limitation:** Ferredoxins are among the most structurally diverse protein families;
a single model will not discriminate photosynthetic from respiratory from nitrogenase
ferredoxins. Use gene-neighborhood context (co-occurrence with *petA*, *nifH*, etc.)
for functional assignment.

---

### 1d. Rieske [2Fe-2S] proteins

Rieske proteins carry a [2Fe-2S] cluster coordinated by two His and two Cys residues
(distinct from plant-type ferredoxin coordination). They appear in two major contexts:

1. **Cytochrome bc₁ complex (Complex III)** and **cytochrome b₆f complex**:
   electron transfer in aerobic respiration and oxygenic photosynthesis. The Rieske
   subunit (PetC in cyanobacteria, QcrA in Actinobacteria) is the mobile head domain
   that shuttles electrons to cytochrome *c*₁/f.

2. **Rieske non-heme iron dioxygenases**: a large superfamily of oxygenases that
   use [2Fe-2S] + mononuclear Fe for aromatic ring hydroxylation. Relevant for:
   - Polycyclic aromatic hydrocarbon (PAH) degradation
   - Nitrobenzene degradation
   - Catechol biosynthesis (siderophore precursor via 2,3-dihydroxybenzoate pathway)

**MetalGenie-Evo relevance:** The bc₁ Rieske subunit is present in virtually all
aerobic heterotrophs — not discriminating. However, its absence (e.g., some anaerobes)
or replacement by alternative complexes indicates metabolic lifestyle. Rieske
dioxygenases are valuable markers of aromatic carbon cycling linked to iron cycling
(Fe²⁺ is the catalytic metal; reductase components use ferredoxin).

**HMM sources:**
- **bc₁ Rieske:** Pfam PF00355 (Rieske); NCBIfam TIGR01986 (*petC*, cyanobacterial)
- **Dioxygenase Rieske:** Pfam PF00355 + PF00910 (Ring_hydroxyl_B); subfamily models
  available in NCBIfam for specific dioxygenases

---

### Limitations (Fe-S category overall)

- ISC genes are near-universal in bacteria; alone not diagnostic of metal cycling
- NIF and ISC share gene names (*nifS/iscS*, *nifU/iscU*); subfamily-level models required
- *sufA*/*iscA* are paralogs with overlapping sequences — single model likely insufficient
- FNR and CRP (cAMP receptor) share the Crp/Fnr domain; require Fe-S-specific features
- IscR and NsrR share Rrf2 domain; distinguish by co-occurrence with *isc* vs *nsr* genes
- **Recommended category structure:**
  - `iron_stress-fes_cluster_suf` — SUF assembly system
  - `iron_metabolism-fes_cluster_isc` — ISC assembly system
  - `iron_metabolism-fes_sensors` — FNR, IscR, NsrR, SoxR, WhiB
  - `iron_stress-iron_economy` — ferredoxin / flavodoxin switch

---

## 2. Fur-Family Metalloregulators (Zur, Mur, Nur, PerR, Irr)

### Biology

The Fur (Ferric Uptake Regulator) superfamily is the primary metalloregulatory system
in bacteria. All members are homodimeric HTH transcription factors that use a metal
cofactor to sense intracellular metal availability. The currently covered **Fur** paralog
(already in MetalGenie-Evo: `iron_gene_regulation`) represses iron acquisition genes
when Fe²⁺ is sufficient.

The missing paralogs each sense a different metal or use a different mechanism:

| Paralog | Sensed signal | Key targets | First characterized in |
|---------|---------------|-------------|------------------------|
| **Zur** | Zn²⁺ | *znuABC* (zinc uptake), *zur* itself | *E. coli* K-12 |
| **Mur** | Mn²⁺ | *mntABC*, *sitABCD* (Mn uptake) | *Rhizobium leguminosarum* |
| **Nur** | Ni²⁺ | nickel uptake genes, *sodF* | *Streptomyces coelicolor* |
| **PerR** | H₂O₂ (via Fe or Mn oxidation) | *katA*, *mrgA*, *ahpCF* (peroxide defense) | *B. subtilis* 168 |
| **Irr** | Heme availability (not free Fe) | *hemB*, iron storage genes | *Bradyrhizobium japonicum* |
| **BosR** | Cu²⁺-like signal | RpoS regulation, oxidative stress | *Borrelia burgdorferi* |

**MetalGenie-Evo relevance:** Zur, Mur, and Nur directly link to Zn, Mn, and Ni
cycling. PerR is a key indicator of oxidative stress management coupled to Fe/Mn
homeostasis (relevant in environments with variable redox). Irr is critical for
understanding Fe regulation in α-Proteobacteria (rhizobia, Agrobacterium) which
lack canonical Fur.

### Model organisms and key evidence

| Organism | Paralog | Evidence |
|----------|---------|---------|
| *E. coli* K-12 | Zur, Fur | Zur deletion → constitutive *znuABC* expression; Zur crystal structure with Zn²⁺ cofactor |
| *B. subtilis* 168 | Fur, MntR (DtxR family), PerR | PerR senses Fe²⁺/Mn²⁺ ratio; H₂O₂ oxidizes metal → DNA release |
| *Streptomyces coelicolor* A3(2) | Nur | First Ni-sensing Fur paralog; controls *nikABCDE* |
| *Rhizobium leguminosarum* bv. viciae | Mur | Controls Mn homeostasis in nitrogen-fixing root symbiont |
| *Bradyrhizobium japonicum* USDA 110 | Irr | Heme-sensing; does NOT bind iron directly; unique mechanism |
| *Pseudomonas aeruginosa* PAO1 | Fur, Zur | Both characterized; *P. aeruginosa* Fur also controls virulence |

### References (PubMed)

- Lee J-W, Helmann JD (2007). Functional specialization within the Fur family of
  metalloregulators. *Biometals* 20:485–99.
  [doi:10.1007/s10534-006-9070-7](https://doi.org/10.1007/s10534-006-9070-7)

- Sevilla E, Bes MT, Peleato ML, Fillat MF (2021). Fur-like proteins: Beyond the ferric
  uptake regulator (Fur) paralog. *Arch Biochem Biophys* 701:108770.
  [doi:10.1016/j.abb.2021.108770](https://doi.org/10.1016/j.abb.2021.108770)

- Helmann JD (2014). Specificity of metal sensing: iron and manganese homeostasis in
  *Bacillus subtilis*. *J Biol Chem* 289:28112–20.
  [doi:10.1074/jbc.R114.587071](https://doi.org/10.1074/jbc.R114.587071)

### HMM sources

- **NCBIfam/TIGRFAM:** TIGR01781 (Fur), TIGR02086 (Zur), TIGR02783 (Mur/Fur family)
- **Pfam:** PF01475 (FUR) — covers whole family, insufficient for paralog discrimination
- **Critical:** Fur vs Zur vs PerR discrimination requires *subfamily-specific* models.
  PF01475 alone is not usable for paralog assignment.
- NCBIfam has paralog-discriminating models with calibrated cutoffs — prefer these.

### Limitations

- Fur/Zur/PerR are structurally very similar; HMM overlap is a real risk
- Metal selectivity is determined by subtle active-site differences not captured
  in global sequence-level models
- PerR in some organisms can use either Fe or Mn — functional annotation from
  sequence alone is ambiguous
- Irr is restricted to α-Proteobacteria; annotating it in other lineages = false positive
- **Recommendation:** separate models per paralog, annotated in categories
  `zinc_regulation-zur`, `manganese_regulation-mur`, `nickel_regulation-nur`,
  `oxidative_stress_regulation-perR`, `iron_regulation-irr`

---

## 3. Molybdenum Cofactor (Moco) Biosynthesis and Molybdate Transport

### Biology

Molybdenum (Mo) is an essential trace metal for most bacteria, required as the
molybdenum cofactor (Moco) in a wide range of enzymes: nitrate reductase (Nar),
sulfite oxidase, DMSO reductase, formate dehydrogenase, xanthine dehydrogenase,
and nitrogenase (Mo-Fe protein). Moco biosynthesis requires 6 conserved enzymatic
steps (MoaA–MoaE, MoeA, MogA) and subsequent assembly into enzyme-specific forms
(bis-MGD for DMSO reductase family, via MobA).

**Molybdate transport:**
- `modA` — periplasmic molybdate-binding protein (ABC transporter SBP)
- `modB` — permease
- `modC` — ATPase
- `modE` — transcriptional repressor, senses molybdate directly (binds ModE when
  Mo is sufficient → represses *modABC*)

**MetalGenie-Evo relevance:** ModABC is a direct molybdenum uptake system, and its
presence/absence indicates the capacity for Mo-dependent metabolism. Particularly
relevant for:
- N-cycling (nitrate reduction requires Mo-Nar)
- S-cycling (DMSO reductase in marine bacteria)
- Global biogeochemical surveys where Mo availability varies (ocean anoxic zones,
  euxinic basins where Mo is sequestered by sulfide)

### Model organisms and key evidence

| Organism | Key evidence |
|----------|-------------|
| *E. coli* K-12 | Complete Moco biosynthesis and *modABC* characterized; ModE crystal structure with molybdate |
| *Azotobacter vinelandii* | ModABC essential for Mo-nitrogenase; Mo-free nitrogenases (FeFe, VFe) induced when ModABC is repressed |
| *Rhodobacter capsulatus* | DMSO reductase (bis-MGD); *modABC* characterized |
| *Desulfovibrio* spp. | Formate dehydrogenase (Mo-dependent); relevant for sulfate-reducing bacteria |

### References (PubMed)

- Demtröder L, Narberhaus F, Masepohl B (2018). Coordinated regulation of nitrogen
  fixation and molybdate transport by molybdenum. *Mol Microbiol* 111:17–30.
  [doi:10.1111/mmi.14152](https://doi.org/10.1111/mmi.14152)

- Leimkühler S, Bühning M, Beilschmidt L (2017). Shared Sulfur Mobilization Routes
  for tRNA Thiolation and Molybdenum Cofactor Biosynthesis in Prokaryotes and
  Eukaryotes. *Biomolecules* 7:5.
  [doi:10.3390/biom7010005](https://doi.org/10.3390/biom7010005)

### HMM sources

- **NCBIfam/TIGRFAM:** TIGR01502 (moaA), TIGR01703 (moaC), TIGR01661 (moaD),
  TIGR01660 (moaE), TIGR02457 (moeA), TIGR02126 (mogA), TIGR00659 (modA),
  TIGR00680 (modB), TIGR00582 (modC)
- **Pfam:** PF00994 (MoeA), PF01077 (ModA)

### Limitations

- Moco biosynthesis genes (*moaA–moeA*) are near-universal in bacteria that use
  Mo-enzymes; alone not diagnostic of specific cycling function
- `modA` is most specific for molybdate uptake — prioritize this for cycling studies
- Mo-nitrogenase vs DMSO reductase vs nitrate reductase cannot be distinguished
  from Moco genes alone — require Nar/Nos/Dsr/NifHDK genes separately
- **Recommendation:** include `modABC` transport system in category
  `molybdenum_metabolism-transport`; include Moco biosynthesis as
  `molybdenum_metabolism-cofactor_biosynthesis`

---

## 4. Cobalamin (Vitamin B₁₂) Biosynthesis and Cobalt Trafficking

### Biology

Cobalamin is the only biological molecule that requires cobalt as an essential metal.
Its de novo biosynthesis requires ~25–30 enzymes and follows two distinct pathways
in bacteria:

| Pathway | Model organism | Key difference |
|---------|---------------|----------------|
| **Aerobic** (*cob* genes) | *Pseudomonas denitrificans*, *Rhodobacter capsulatus* | Cobalt inserted late (after ring contraction) |
| **Anaerobic** (*cbi* genes) | *Salmonella enterica* LT2, *Bacillus megaterium* | Cobalt inserted early; oxygen-sensitive intermediates |

Only ~30% of sequenced bacteria can synthesize cobalamin de novo. The majority
are auxotrophs that rely on import via BtuB (outer membrane, TonB-dependent) and
BtuCD (ABC transporter). Cobalt homeostasis is directly tied to B12 biosynthesis.

**Cobalt transport:**
- `cbtA`/`corA`/`nikABCDE` (some NikABCDE systems transport Co²⁺ as well as Ni²⁺)
- `cbiN`, `cbiO`, `cbiQ` — cobalt-specific ABC transporter (anaerobic pathway)

**MetalGenie-Evo relevance:**
- Cobalt is already in the library as resistance (`metal_resistance-cobalt_zinc_cadmium`)
  but cobalamin biosynthesis is the main biological *use* of cobalt
- Distinguishing cobalt resistance from cobalt utilization (B12 biosynthesis) is
  ecologically critical
- In cobalt-rich environments (ophiolites, ultramafic soils), presence of the
  complete *cbi/cob* pathway indicates cobalt-utilizing organisms vs cobalt-resistant organisms

### Model organisms and key evidence

| Organism | Pathway | Evidence |
|----------|---------|---------|
| *Salmonella enterica* LT2 | Anaerobic | Roth group: complete anaerobic pathway mapped; *cbiA–cbiZ* |
| *Bacillus megaterium* DSM319 | Anaerobic | Used for commercial B12 production; complete pathway |
| *Pseudomonas denitrificans* | Aerobic | First aerobic pathway organism; *cobA–cobX* |
| *Rhodobacter capsulatus* SB1003 | Aerobic | Photosynthetic organism; B12 required for light-harvesting proteins |
| *Sinorhizobium meliloti* 1021 | Aerobic | B12 required for symbiotic nitrogen fixation |

### References (PubMed)

- Moore SJ, Warren MJ (2012). The anaerobic biosynthesis of vitamin B12.
  *Biochem Soc Trans* 40:581–6.
  [doi:10.1042/BST20120066](https://doi.org/10.1042/BST20120066)

- Scott AI, Roessner CA (2002). Biosynthesis of cobalamin (vitamin B12).
  *Biochem Soc Trans* 30:613–20.
  [doi:10.1042/bst0300613](https://doi.org/10.1042/bst0300613)

### HMM sources

- **NCBIfam/TIGRFAM:** extensive — TIGR00738 (*cobI*/*cbiG*), TIGR02463 (*cbiA*),
  TIGR02472 (*cbiB*), TIGR02467 (*cbiC*), TIGR02464 (*cbiD*)... and ~20 more
- **Pfam:** PF02614 (CbiA), PF01923 (CbiX), PF02315 (CobL)
- Coverage is good in TIGRFAM; many models have validated cutoffs
- Key differentiator: aerobic (*cob*) vs anaerobic (*cbi*) models are distinct —
  annotate separately

### Limitations

- ~30 genes per complete pathway → high model count
- Partial pathways are common (import of intermediates); incomplete gene sets
  do not necessarily mean absence of B12
- Some *cbi/cob* genes have paralogs in other tetrapyrrole pathways (siroheme, heme);
  careful with model selection
- BtuB (cobalamin outer membrane receptor) is already partially in FeGenie
  (heme_transport section) — check for overlap before adding
- **Recommendation:** add as `cobalt_metabolism-cobalamin_biosynthesis_aerobic`,
  `cobalt_metabolism-cobalamin_biosynthesis_anaerobic`, `cobalt_metabolism-cobalamin_transport`

---

## 5. Zinc and Manganese Uptake Systems

### Biology

Zinc and manganese are essential micronutrients but toxic in excess. Bacteria
maintain homeostasis through tightly regulated high-affinity uptake systems:

**Zinc uptake:**
- `znuA` — periplasmic Zn²⁺-binding protein (cluster A-I SBP; binds Zn via His residues)
- `znuB` — permease
- `znuC` — ATPase
- Regulated by **Zur** (Zn²⁺-sensing Fur paralog, see category 2 above)
- `zupT` — ZupT (ZIP-family Zn²⁺ importer; lower affinity, constitutive)
- `znuD` / `zur`-regulated outer membrane TonB-dependent Zn receptor — in some bacteria

**Manganese uptake:**
- `mntA`/`psaA`/`sitA` — periplasmic Mn²⁺-binding protein (cluster A-I SBP)
- `mntB`/`psaB`/`sitB` — permease
- `mntC`/`psaC`/`sitC` — ATPase
- `mntH` / `mntP` — MntH (NRAMP family H⁺-coupled Mn²⁺ importer; single component)
- Regulated by **MntR** (DtxR family) in Firmicutes, **Mur** (Fur family) in Proteobacteria

**Note on naming confusion:** Multiple ABC transporter systems transport Mn²⁺:
- MntABC (*B. subtilis*) — primary Mn uptake
- PsaBCA (*Streptococcus pneumoniae*) — Mn/Zn
- SitABCD (*E. coli/Salmonella*) — Mn/Fe (broad specificity)
- All share the cluster A-I SBP fold; subfamily-specific models needed

**MetalGenie-Evo relevance:**
- Zinc and manganese are the most common enzyme cofactors after iron
- ZnuABC and MntABC presence indicates need for high-affinity uptake (limiting conditions)
- MntH (NRAMP) is ecologically relevant: NRAMP transporters are key players in
  Mn bioavailability in soil and ocean environments
- These systems connect directly to oxidative stress defense (Mn-SOD,
  Mn-catalase in peroxide stress)

### Model organisms and key evidence

| Organism | System | Metal | Evidence |
|----------|--------|-------|---------|
| *E. coli* K-12 | ZnuABC | Zn | Patzer & Hantke 1998, 2000; ZnuA crystal with Zn²⁺ |
| *B. subtilis* 168 | MntABC | Mn | Que & Helmann 2000; MntR regulator characterized |
| *Salmonella* Typhimurium | SitABCD | Mn/Fe | Required for virulence; mutants attenuated |
| *S. pneumoniae* D39 | PsaBCA | Mn (Zn) | Required for virulence; Mn critical for Mn-SOD |
| *E. coli* K-12 | MntH | Mn | Kehres et al. 2000; induced under Mn limitation and H₂O₂ |

### References (PubMed)

- Porcheron G, Garenaux A, Proulx J, Sabri M, Dozois CM (2013). Iron, copper, zinc,
  and manganese transport and regulation in pathogenic Enterobacteria. *Front Cell Infect
  Microbiol* 3:90.
  [doi:10.3389/fcimb.2013.00090](https://doi.org/10.3389/fcimb.2013.00090)

- Charbonnier M et al. (2022). Battle for Metals: Regulatory RNAs at the Front Line.
  *Front Cell Infect Microbiol* 12:952948.
  [doi:10.3389/fcimb.2022.952948](https://doi.org/10.3389/fcimb.2022.952948)

- Helmann JD (2014). Specificity of metal sensing: iron and manganese homeostasis
  in *Bacillus subtilis*. *J Biol Chem* 289:28112–20.
  [doi:10.1074/jbc.R114.587071](https://doi.org/10.1074/jbc.R114.587071)

### HMM sources

- **ZnuABC:** NCBIfam TIGR00801 (*znuA*), TIGR00800 (*znuC*); Pfam PF00496 (SBP_bac_1)
- **MntH:** NCBIfam TIGR00361 (*mntH*); Pfam PF01566 (Nramp)
- **MntABC/SitABCD:** TIGRFAM has subfamily models — important to use *specific*
  Mn-SBP models, not generic cluster A-I models (too broad)
- Pfam PF01048 (SBP_bac_3) covers Mn-SBP but non-specifically

### Limitations

- ZnuA / MntA / SitA all belong to cluster A-I ABC transporter SBP superfamily;
  metal selectivity determined by binding site residues, not overall fold
- Generic SBP models will produce massive false positives; use only
  subfamily-specific NCBIfam models
- SitABCD in *E. coli* / *Salmonella* transports both Mn²⁺ and Fe²⁺ — functional
  ambiguity at annotation level
- MntH (NRAMP) is ubiquitous in Firmicutes and some Proteobacteria; use with caution
  in broad surveys as it will hit in almost every genome in those phyla
- **Recommendation:** add as `zinc_metabolism-zinc_uptake` (ZnuABC + ZupT),
  `manganese_metabolism-manganese_uptake` (MntABC + MntH + PsaBCA)

---

## 6. Cytochrome 579 and Rusticyanin — Acidophilic Iron Oxidation Chain

### Biology

Acidophilic iron-oxidizing bacteria (principally *Acidithiobacillus ferrooxidans*)
use a unique downhill electron transfer chain to oxidize Fe²⁺ at low pH where chemical
reoxidation is slow. This chain is composed of distinct proteins not found in
neutrophilic iron oxidizers:

| Component | Type | Function |
|-----------|------|----------|
| **Cyc2** | Outer membrane cytochrome c | Accepts electrons from periplasmic Fe²⁺; already in FeGenie library |
| **Rusticyanin** | Blue copper protein (cupredoxin) | Most abundant periplasmic protein in *At. ferrooxidans* (~350 mg/mL); shuttles electrons from Cyc2 to terminal oxidase; redox potential +680 mV — among the highest known |
| **Cytochrome 579** | High-potential c-type cytochrome | Periplasmic; accepts electrons from rusticyanin; part of downhill chain to *aa₃* oxidase |
| **CycA1** | Diheme cytochrome c₄ | Bridges rusticyanin and the terminal *aa₃* copper oxidase |

Based on PubMed-indexed work by Blake & White (2020), at least **six distinct
electron carrier arrangements** are used across 8+ phyla of iron-oxidizing microorganisms,
suggesting that rusticyanin/Cyt579 is only one solution to the problem of Fe²⁺ oxidation
at low pH. The FeGenie library captured Cyc2 but lacked rusticyanin and Cyt579.
**Both are now added (2026-05-23).**

**MetalGenie-Evo relevance:**
- Rusticyanin + Cyt579 are diagnostic of *acidophilic* iron oxidation (acid mine drainage,
  sulfidic mine waste, bioleaching systems)
- Their presence in a metagenome indicates Fe-S mineral weathering capacity
- Cyc1 and Cyc2 are already in the library (FeGenie); rusticyanin and Cyt579 complete
  the acidophilic Fe(II)-oxidation marker set

### Model organisms

| Organism | Components | pH optimum |
|----------|-----------|------------|
| *Acidithiobacillus ferrooxidans* ATCC 23270 | Cyc2–rusticyanin–Cyt579–CycA1–*aa₃* | 1.5–3.0 |
| *Leptospirillum ferrooxidans* | Cyc2 + Cyt579 homolog (no rusticyanin) | 1.5–3.5 |
| *Ferroplasma acidiphilum* (archaeon) | Different cupredoxin; no rusticyanin | 1.0–2.5 |

### Reference (PubMed)

- Blake RC, White RA (2020). In situ absorbance measurements: a new means to study
  respiratory electron transfer in chemolithotrophic microorganisms.
  *Adv Microb Physiol* 76:81–127.
  [doi:10.1016/bs.ampbs.2020.01.003](https://doi.org/10.1016/bs.ampbs.2020.01.003)

### HMM sources — DEPLOYED

| Protein | HMM | Source | GA | LENG | NSEQ | File | Key reference |
|---------|-----|--------|----|------|------|------|---------------|
| Rusticyanin | TIGR03095.1 | NCBIfam/JCVI | 146.7 | 148 | 4 | `hmm_library/iron_oxidation/rusticyanin_TIGR03095.hmm` | [doi:10.1099/mic.0.26966-0](https://doi.org/10.1099/mic.0.26966-0) |
| Cytochrome 579 | NF033864.1 | NCBIfam | 280.0 | 179 | 4 | `hmm_library/iron_oxidation/cytochrome579_NF033864.hmm` | [doi:10.1128/AEM.02799-07](https://doi.org/10.1128/AEM.02799-07) |

Note: NF033156 is *not* rusticyanin (it is a bleomycin-binding protein). Correct model
is TIGR03095 (`rusti_cyanin` equivalog, JCVI, NSEQ=4, GA=TC=146.7, NC=68.5).

Both models have NCBIfam-calibrated GA=TC cutoffs (no borderline sequences in NCBI
training universe). Added to `iron_oxidation` category alongside Cyc1/Cyc2 (consistent
placement — no separate subcategory created).

---

## 7. MtrA vs MtoA — Iron Reduction vs Iron Oxidation Discrimination

### Biology and the problem

MtrA (Metal-reducing) and MtoA (Metal-oxidizing) are both **periplasmic decaheme
cytochrome c proteins** involved in extracellular electron transfer — but they
catalyze **opposite reactions**:

| Protein | System | Direction | Reference organism |
|---------|--------|-----------|-------------------|
| **MtrA** | Mtr pathway (MtrCAB + OmcA) | **Fe³⁺ → Fe²⁺** (iron REDUCTION) | *Shewanella oneidensis* MR-1 |
| **MtoA** | MtoAB complex | **Fe²⁺ → Fe³⁺** (iron OXIDATION) | *Sideroxydans lithotrophicus* ES-1 |

The problem: MtrA and MtoA share significant sequence similarity (both are decaheme
cytochromes in the same MHC superfamily). A single HMM trained on one can
cross-hit the other. This is a known unresolved issue in FeGenie flagged by its
developers.

**Quantified cross-hit (measured in MetalGenie-Evo v0.x):**

```
MtrA HMM (nseq=9, cutoff=140) vs MtoA consensus: 265.8 bits  — 1.9× cutoff
MtoA HMM (nseq=10, cutoff=140) vs MtrA consensus: 289.8 bits  — 2.1× cutoff
```

At the calibrated cutoff of 140 bits, both models would annotate any decaheme
MHC of this family as BOTH MtrA AND MtoA. The "best-hit wins" assignment in the
pipeline reduces error rate but does not eliminate it.

**Ecological consequence:** Misannotating MtoA as MtrA (or vice versa) in a metagenome
would flip the functional prediction from iron oxidation to iron reduction — a
fundamental error in biogeochemical interpretation.

### Evidence for MtoAB in diverse environments

From PubMed (He et al. 2016), MtoAB was identified in metagenomic analysis of the
nitrate-dependent Fe(II)-oxidizing culture KS: MtoAB homologs found in *Gallionellaceae*
alongside Cyc2 and OmpB, confirming that MtoAB is widespread in neutrophilic iron oxidizers
and not limited to *Sideroxydans*.

### Model organisms

| Organism | Protein | Direction | Notes |
|----------|---------|-----------|-------|
| *Shewanella oneidensis* MR-1 | MtrA | Reduction | Best-characterized DIR system; MtrCAB + OmcA |
| *Shewanella loihica* PV-4 | MtrA homolog | Reduction | Deep-sea iron reducer |
| *Sideroxydans lithotrophicus* ES-1 | MtoA | Oxidation | First MtoA characterized; microaerobic |
| *Gallionellaceae* sp. (enrichment KS) | MtoAB homolog | Oxidation | Nitrate-dependent Fe(II) oxidation |
| *Geobacter sulfurreducens* PCA | OmcB/OmcZ | Reduction | Different MHC topology; Geobacter uses different architecture |

### Reference (PubMed)

- He S, Tominski C, Kappler A, Behrens S, Roden EE (2016). Metagenomic Analyses of
  the Autotrophic Fe(II)-Oxidizing, Nitrate-Reducing Enrichment Culture KS.
  *Appl Environ Microbiol* 82:2656–2668.
  [doi:10.1128/AEM.03493-15](https://doi.org/10.1128/AEM.03493-15)

### Recommended solution

1. **Build separate subfamily HMMs** from curated seed alignments using
   `scripts/build_mtr_mto_subfamily_hmms.py` (seeds in `data/seeds/mtr_mto_seeds.tsv`).
   Anchor on the heme-binding motifs and the distinct N-terminal signal sequences.
   Seed accessions curated from primary literature:
   - MtrA seeds: *Shewanella oneidensis* MR-1 WP_011071012.1, *S. loihica* PV-4
     WP_012254745.1, *S. baltica* OS185 WP_006080780.1, *S.* sp. ANA-3 YP_875003.1
   - MtoA seeds: *Sideroxydans lithotrophicus* ES-1 ADE57625.1,
     *Gallionella capsiferriformans* ES-2 AEX59023.1, *Gallionellaceae* sp. KS WP_044467591.1

2. **Use gene neighborhood** as an interim tie-breaker:
   - `MtrC_TIGR03507` (outer-membrane cytochrome, iron_reduction-only) adjacent
     to the hit → iron reduction assignment.
   - `Cyc2`, `FoxA`, `FoxE` (iron_oxidation) adjacent → oxidation assignment.

3. **Set high bitscore cutoffs** after cross-validation — from `build_mtr_mto_subfamily_hmms.py`
   output, set cutoff = correct-hit minimum score − 10 bits.

4. **Operon architecture confirms assignment:**
   - Reduction: MtrC → MtrA → MtrB (or MtrA → MtrB → MtrC depending on strain)
   - Oxidation: MtoA → MtoB (two-gene cluster in *Sideroxydans/Gallionella*)

5. **Recommended categories remain:** `iron_reduction` (MtrA/B/C, OmcA, OmcZ)
   vs `iron_oxidation` (MtoA/B, Cyc2, FoxA-Z); existing FeGenie category split
   is correct — only the model boundaries need improvement.

---

## 8. Iron Cycling in Fully Oxic, Iron-Rich Environments

### The question

In environments with abundant iron but fully oxic conditions (e.g., iron-rich soils,
mining drainage settling ponds, hydrothermal vent plumes in oxygenated water), how
do both iron oxidation AND iron reduction coexist?

### The answer: cryptic iron cycling

Based on Consensus-retrieved literature [Berg et al. 2016; Peng et al. 2019;
Emerson et al. 2012; Díaz-González et al. 2025], the dominant mechanism is
**cryptic iron cycling** — rapid oxidation and reduction occurring simultaneously
in tight spatial coupling, with no net measurable change in bulk Fe²⁺/Fe³⁺ but
high turnover rates.

**Mechanisms enabling both oxidation and reduction in bulk-oxic systems:**

| Mechanism | Scale | Key genes |
|-----------|-------|-----------|
| **Biofilm micro-anoxic gradients** | Biofilm interior (µm scale) | Iron reducers (MtrCAB, OmcA) in anoxic core; iron oxidizers (MtoAB, Cyc2) at surface |
| **Fe-OM complexes + photoreduction** | Surface water (photic zone) | Light → Fe³⁺-OM → Fe²⁺-OM; microbes reoxidize via MtoAB or multicopper oxidases |
| **OmcA "capacitor" effect** | Mineral surface | OmcA stores photogenerated electrons on hematite surface, releases for Fe reduction in dark |
| **Dual-capacity organisms** | Single cell | *Shewanella* spp. can switch between oxidation/reduction; *Gallionellaceae* tolerate narrow O₂ gradients |
| **Siderophore-mediated Fe mobilization** | Cell vicinity | Siderophores dissolve Fe(OH)₃ → bioavailable Fe³⁺ even in fully oxic conditions |
| **Dps/ferritin iron sequestration** | Intracellular | Store Fe²⁺ → prevent Fenton; release as needed → local cycling |

### Key literature (Consensus)

- [The microbial ferrous wheel: iron cycling in terrestrial, freshwater, and marine
  environments](https://consensus.app/papers/details/2c40d32352b857e1b85385cf1d3d7e98/?utm_source=claude_desktop)
  — Emerson D et al. 2012, *Front Microbiol*. Conceptual framework for the "ferrous wheel" —
  iron cycling as continuous rapid turnover at oxic-anoxic interfaces and within biofilms.

- [Cryptic Cycling of Complexes Containing Fe(III) and Organic Matter by Phototrophic
  Fe(II)-Oxidizing Bacteria](https://consensus.app/papers/details/b9881804775e5669b5e34abe6c85439f/?utm_source=claude_desktop)
  — Peng C et al. 2019, *Appl Environ Microbiol*. Light-driven Fe(III)-OM → Fe(II)-OM
  → microbial Fe(II) reoxidation = cryptic cycle in photic zones.

- [Iron Cycle Tuned by Outer-Membrane Cytochromes of DMRB](https://consensus.app/papers/details/6b28eb0e2f075451a8c60413b1f7a56d/?utm_source=claude_desktop)
  — Yu S-S et al. 2021, *Environ Sci Technol*. OmcA as a photon-electron capacitor on
  hematite — stores and releases electrons for iron cycling independent of bulk O₂.

- [Trait-based meta-analysis of microbial guilds in the iron redox cycle](https://consensus.app/papers/details/2e8a2870992a5cc8a21bf065d32100e4/?utm_source=claude_desktop)
  — Díaz-González F et al. 2025, *mSystems*. Dual-capacity Fe oxidizers/reducers as
  overlooked mediators of cryptic cycling; guild framework for environmental metagenomes.

- [Iron is not everything: unexpected complex metabolic responses between iron-cycling
  microorganisms](https://consensus.app/papers/details/fac85fcca74d5d6c848b1bb52b6395c5/?utm_source=claude_desktop)
  — Cooper R et al. 2020, *ISME J*. *Sideroxydans* + *Shewanella* co-culture: iron cycling
  genes NOT the most differentially expressed — metabolic crosstalk (biofilm formation,
  amino acid metabolism) drives the interaction.

### Gene annotation targets for oxic iron-rich environments

**Already in library (relevant):**
- Siderophore synthesis + transport — mobilize Fe(OH)₃ even in oxic conditions
- Iron storage (Bfr, FtnA) — buffer free iron, prevent Fenton
- Heme oxygenase — iron recycling from heme
- OmcA/MtrCAB — iron reduction capacity even in organisms co-occurring with O₂

**Missing (to add):**
- **Dps** (DNA-binding protein from starved cells) — ferritin superfamily; sequesters
  Fe²⁺ using H₂O₂ as oxidant (kills two toxins: free Fe²⁺ + H₂O₂). Key Fenton prevention.
  Category: `iron_storage` (extend existing)
- **Multicopper oxidases (MCO)** — CueO, CopA-linked MCO oxidize Fe²⁺ in periplasm.
  Relevant in oxic iron-rich soils where these function as iron oxidases incidentally.
  Category: new `iron_oxidation-multicopper_oxidase`
- **Ferrireductase (Fre/FRE homologs)** — NADH-dependent Fe³⁺ → Fe²⁺ reduction
  inside cell for iron assimilation from transferrin-like chelates. Not dissimilatory.
  Category: `iron_metabolism-assimilatory_reduction`

### Limitation

**Most cryptic cycling is not directly gene-annotatable from static sequence data.**
The balance between oxidation and reduction in oxic environments is determined by:
- Local oxygen gradients (not encoded in genes)
- Biofilm architecture (not sequence-interpretable)
- Light availability (drives photo-Fenton; purely abiotic)

MetalGenie-Evo can identify the *capacity* for both iron oxidation and iron reduction
genes co-occurring in a single genome or metagenome bin — the co-presence of MtrCAB
(reduction) + Cyc2/MtoAB (oxidation) in one organism is a strong indicator of
dual-capacity cryptic cycling potential. This is the most powerful annotation strategy
for this ecological context.

---

## Priority and implementation order

| Priority | Category | Rationale | Key references (PubMed) |
|----------|----------|-----------|-------------------------|
| 1 | **MtrA / MtoA disambiguation** | ✅ DONE (2026-05-23) — curated seeds, TC/GA/NC calibrated | [doi:10.1128/AEM.03493-15](https://doi.org/10.1128/AEM.03493-15) |
| 2 | **Rusticyanin + Cytochrome 579** | ✅ DONE (2026-05-23) — TIGR03095 + NF033864 deployed to `iron_oxidation/` | [doi:10.1099/mic.0.26966-0](https://doi.org/10.1099/mic.0.26966-0); [doi:10.1128/AEM.02799-07](https://doi.org/10.1128/AEM.02799-07) |
| 3 | **Fe-S sensors (IscR, SoxR, NsrR)** | ✅ DONE (2026-05-23) — TIGR02010.1/TIGR01950.1/NF008240.0 deployed. FNR skipped: no clean NCBIfam equivalog (CRP-FNR domain non-specific). | IscR: [doi:10.1016/s1369-5274(03)00039-0](https://doi.org/10.1016/s1369-5274(03)00039-0); SoxR: [doi:10.1016/j.jinorgbio.2013.11.008](https://doi.org/10.1016/j.jinorgbio.2013.11.008); NsrR: [doi:10.1371/journal.pone.0003623](https://doi.org/10.1371/journal.pone.0003623) |
| 4 | **Ferredoxin / Flavodoxin switch** | ✅ DONE (2026-05-23) — TIGR01752.1 (flav_long) + TIGR01753.1 (flav_short) added to `iron_stress`. No clean petF equivalog in NCBIfam. | [doi:10.1038/382802a0](https://doi.org/10.1038/382802a0) |
| 5 | **Fur-family (Zur, PerR)** | ✅ DONE (2026-05-23) — NF008646.0 (Zur) + NF052545.1 (PerR). Mur/Nur: no clean NCBIfam equivalog found. | Zur: [doi:10.3389/fcimb.2013.00059](https://doi.org/10.3389/fcimb.2013.00059); PerR: [doi:10.1111/j.1365-2958.2008.06192.x](https://doi.org/10.1111/j.1365-2958.2008.06192.x) |
| 6 | **ZnuABC / MntH** | ✅ ALREADY IN LIBRARY — ZnuABC in `metal_resistance-cobalt_zinc_cadmium`, MntH in `metal_resistance-non-specific`. No action needed. | — |
| 7 | **SUF / ISC** | ✅ DONE (2026-05-23) — sufC/sufB/sufS (TIGR01978/01980/01979) + IscS (TIGR02006.1) in new `iron_sulfur_assembly` category. sufD skipped (7-bit GA-NC gap). | SUF: [doi:10.1074/jbc.M308004200](https://doi.org/10.1074/jbc.M308004200); IscS: [doi:10.1074/jbc.M401261200](https://doi.org/10.1074/jbc.M401261200) |
| 8 | **Dps iron storage + Multicopper oxidase** | ✅ DONE (2026-05-23) — Dps (NF009990.0) added to `iron_storage`; multicopper oxidase already in `metal_resistance-copper` (CueO/Mco/MmcO). Dps overlaps with PF00210 (ferritin domain) — intentional, Dps is functionally distinct. | [doi:10.1128/MMBR.66.4.630-670.2002](https://doi.org/10.1128/MMBR.66.4.630-670.2002) |
| 9 | **ModABC (molybdate transport)** | ✅ DONE (2026-05-23) — ModA already in library; added ModC (TIGR02142.1). ModB skipped (GA-NC gap = 0.75 bits). | [doi:10.1111/mmi.14152](https://doi.org/10.1111/mmi.14152) |
| 10 | **Cobalamin biosynthesis** | Complex (30 genes), two pathways; highest implementation cost | [doi:10.1042/BST20120066](https://doi.org/10.1042/BST20120066) |

For each category, the workflow is:
1. Download target NCBIfam/TIGRFAM models by accession from NCBI FTP
2. Check version tag (HMMER3/f) — convert if needed via `normalize_hmm_versions.py`
3. Add to `hmm_library/<new_category>/`
4. Add entries to `hmm_registry.tsv` with source, reference DOI, nseq, cutoff
5. Add directory names to `_ANNOTATE_MAP` in `src/metalgenie_evo/io.py`
6. Run Layer B clustering to check for cross-model redundancy
