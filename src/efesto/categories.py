# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Efesto contributors.

"""Canonical category names and aliases for existing libraries/configuration."""

LEGACY_CATEGORY_NAMES = {'iron_acquisition-iron_transport': 'IronAcquisition-Transport',
 'iron_acquisition-heme_transport': 'Heme-Transport',
 'iron_acquisition-heme_oxygenase': 'Heme-Oxygenase',
 'iron_acquisition-siderophore_synthesis': 'Siderophore-Synthesis',
 'iron_acquisition-siderophore_transport': 'Siderophore-Transport',
 'iron_acquisition-siderophore_transport_potential': 'Siderophore-TransportPotential',
 'iron_reduction': 'Iron-Reduction',
 'iron_oxidation': 'Iron-Oxidation',
 'iron_storage': 'Iron-Storage',
 'iron_gene_regulation': 'Iron-Regulation',
 'iron_resistance': 'Iron-Resistance',
 'iron_stress': 'Iron-StressResponse',
 'iron_sulfur_assembly': 'IronSulfur-Assembly',
 'magnetosome_formation': 'Magnetosome-Formation',
 'possible_iron_oxidation_and_possible_iron_reduction': 'Iron-RedoxPotential',
 'probable_iron_reduction': 'Iron-ReductionPotential',
 'manganese_acquisition': 'Manganese-Acquisition',
 'manganese_oxidation': 'Manganese-Oxidation',
 'mn_gene_regulation': 'Manganese-Regulation',
 'metal_resistance-*': 'MetalResistance-*',
 'metal_resistance-arsenic': 'MetalResistance-Arsenic',
 'metal_resistance-copper': 'MetalResistance-Copper',
 'metal_resistance-mercury': 'MetalResistance-Mercury',
 'metal_resistance-cobalt_zinc_cadmium': 'MetalResistance-CobaltZincCadmium',
 'metal_resistance-chromium': 'MetalResistance-Chromium',
 'metal_resistance-lead': 'MetalResistance-Lead',
 'metal_resistance-nickel': 'MetalResistance-Nickel',
 'metal_resistance-silver': 'MetalResistance-Silver',
 'metal_resistance-gold': 'MetalResistance-Gold',
 'metal_resistance-tellurite': 'MetalResistance-Tellurite',
 'metal_resistance-tellurium': 'MetalResistance-Tellurium',
 'metal_resistance-antimony': 'MetalResistance-Antimony',
 'metal_resistance-bismuth': 'MetalResistance-Bismuth',
 'metal_resistance-unknown': 'MetalResistance-Unknown',
 'metal_resistance-non-specific': 'MetalResistance-NonSpecific',
 'metal_resistance-multimetal': 'MetalResistance-MultiMetal',
 'metal_resistance-molybdenum': 'MetalResistance-Molybdenum',
 'metal_resistance-magnesium': 'MetalResistance-Magnesium',
 'metal_resistance-manganese': 'MetalResistance-Manganese'}

def canonical_category(name):
    """Rename known categories; preserve custom category names and wildcard patterns."""
    name = name.replace("iron_aquisition", "iron_acquisition")
    return LEGACY_CATEGORY_NAMES.get(name, name)


CANONICAL_CATEGORIES = frozenset(['Metallophore-Synthesis',
 'MetalHomeostasis-Lanthanides',
 'MetalAcquisition-Tungsten',
 'MetalAcquisition-MolybdenumTungsten',
 'Arsenic-Oxidation',
 'CarbonMetabolism-Dehydrogenase',
 'Competence-Pilus',
 'Heme-Oxygenase',
 'Heme-Transport',
 'Iron-Oxidation',
 'Iron-RedoxPotential',
 'Iron-Reduction',
 'Iron-ReductionPotential',
 'Iron-Regulation',
 'Iron-Resistance',
 'Iron-Storage',
 'Iron-StressResponse',
 'IronAcquisition-ComplexTransport',
 'IronAcquisition-Transport',
 'IronSulfur-Assembly',
 'IronSulfur-Regulation',
 'Magnetosome-Formation',
 'Manganese-Acquisition',
 'Manganese-Oxidation',
 'Manganese-Regulation',
 'MetalAcquisition-CobaltNickel',
 'MetalAcquisition-Copper',
 'MetalAcquisition-IronManganese',
 'MetalAcquisition-Magnesium',
 'MetalAcquisition-Molybdenum',
 'MetalAcquisition-MultiMetal',
 'MetalAcquisition-Nickel',
 'MetalAcquisition-Zinc',
 'MetalHomeostasis-Copper',
 'MetalHomeostasis-IronZinc',
 'MetalRegulation-MultiMetal',
 'MetalRegulation-Zinc',
 'MetalResistance-Antimony',
 'MetalResistance-Arsenic',
 'MetalResistance-Bismuth',
 'MetalResistance-Chromium',
 'MetalResistance-CobaltZincCadmium',
 'MetalResistance-Copper',
 'MetalResistance-Gold',
 'MetalResistance-Lead',
 'MetalResistance-Magnesium',
 'MetalResistance-Manganese',
 'MetalResistance-Mercury',
 'MetalResistance-Molybdenum',
 'MetalResistance-MultiMetal',
 'MetalResistance-Nickel',
 'MetalResistance-NonSpecific',
 'MetalResistance-Silver',
 'MetalResistance-Tellurite',
 'MetalResistance-Tellurium',
 'MetalResistance-Unknown',
 'Siderophore-Export',
 'Siderophore-Utilization',
 'Siderophore-FunctionUnresolved',
 'IronAcquisition-IronRelease',
 'Siderophore-Synthesis',
 'Siderophore-Transport',
 'Siderophore-TransportPotential',
 'StressResponse-Nitrosative',
 'StressResponse-Oxidative',
 'StressResponse-Redox',
 'Unclassified-MetalFunction',
 'Virulence-ActinNucleation'])


# Preserve the established neighborhood membership when functional labels change.
# Keys are stable bundled model identifiers; custom categories are not remapped.
_CONTEXT_CATEGORY_OVERRIDES = {
    'Sid_IroB_Shigella_dysenteriae_AAR13266.1_180606': ('Siderophore-Synthesis', 'Siderophore-Transport'),
    'Sid_BesA_Bacillibactin_hydrolyze_Bsubtilis_YP_003867459.1_180603': ('Siderophore-Utilization', 'Siderophore-Synthesis'),
    'IroD-VabH-enterochelin_esterase-rep': ('Siderophore-Utilization', 'Siderophore-Synthesis'),
    'PvsC_vibrioferrin_biosynthesis-rep': ('Siderophore-Export', 'Siderophore-Synthesis'),
    'Sid_PchH_pyochelin_biosynthesis_Pseudomonas_aeruginosa_PA4223_180623': ('Siderophore-Transport', 'Siderophore-Synthesis'),
    'Sid_PchI_pyochelin_biosynthesis_Pseudomonas_aeruginosa_PA4222_180623': ('Siderophore-Transport', 'Siderophore-Synthesis'),
    'Sid_IroE_Shigella_dysenteriae_AAR13269.1_180606': ('Siderophore-Utilization', 'Siderophore-Transport'),
    'Sid_ViuB_Vibriobactin_siderophore_transport_utilization_protein_Vibrio_cholerae_A5F660.txt-nr': ('Siderophore-Utilization', 'Siderophore-Transport'),
    'Sid_FpvG_iron_reductase_PA2403_Paeruginosa_PAO1_180602': ('Siderophore-Utilization', 'Siderophore-Transport'),
    'Sid_FpvH_iron_release_PA2404_Paeruginosa_PAO1_180602': ('Siderophore-Utilization', 'Siderophore-Transport'),
    'Sid_YqjH_NADPH-dep_ferric_chelate_reductase_enterobactin_Q46871_Escherichia_coli_180606': ('IronAcquisition-IronRelease', 'Siderophore-Transport'),
    'PF06276-RhbCF-AcsD-IucA-IucC-PvsB-Ferric_iron_reductase_FhuF_like_transporter': ('Siderophore-FunctionUnresolved', 'Siderophore-Transport'),
    'Sid_EntS_enterobactin_exporter_Ecoli_P24077_180603': ('Siderophore-Export', 'Siderophore-Transport'),
    'Sid_PvdR_Pyoverdine_export_Paeruginosa_PA2389_180623': ('Siderophore-Export', 'Siderophore-Transport'),
    'Sid_PvdT_Pyoverdine_export_Paeruginosa_PA2390_180623': ('Siderophore-Export', 'Siderophore-Transport'),
}


def context_category(row):
    """Return the neighborhood group independently of the output role label."""
    category = canonical_category(row["cat"])
    override = _CONTEXT_CATEGORY_OVERRIDES.get(row["hmm_stem"])
    return override[1] if override and category == override[0] else category


def rule_category(row, rule):
    """Custom JSON rules use output categories unless context basis is explicit."""
    if rule.get("category_basis", "annotation") == "context":
        return context_category(row)
    return canonical_category(row["cat"])
