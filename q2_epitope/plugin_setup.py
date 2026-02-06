# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import q2_epitope
from qiime2.plugin import Plugin, Str, Choices, Range, Collection, Bool, Float
from q2_types.feature_data import FeatureData
from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import (Zscore, Epitope, MappedEpitope, GMT,
                                     Enriched, PSEAScores)


plugin = Plugin(
    name='epitope',
    version=q2_epitope.__version__,
    website='https://github.com/qiime2/q2-epitope',
    package='q2_epitope',
    short_description=('Used to manipulate pepsirf epitope data in QIIME 2'),
    description=('Used to manipulate pepsirf epitope data in QIIME 2')
)

plugin.methods.register_function(
    function=q2_epitope.create_epitope_map,
    inputs={'epitope': FeatureData[Epitope]},
    parameters={
        'collapse': Str % Choices(['Bacterial', 'Viral', 'Both'])
    },
    outputs=[
        ('epitope_map', FeatureData[MappedEpitope])
    ],
    input_descriptions={'epitope': 'FeatureTable containing at least '
                        'CodeName, SpeciesID, ClusterID, EpitopeWindow, '
                        'Species, and Subtype columns'},
    parameter_descriptions={},
    output_descriptions={'epitope_map': 'FeatureTable containing columns '
                         'described in action descriptions.'},
    name='create epitope map',
    description='Creates the fully defined epitope name '
                'species_clusterID_EpitopeWindow mapped to peptide code names '
                'and species/subtypes the epitope is associated with.'
)

plugin.methods.register_function(
    function=q2_epitope.epitope_zscore,
    inputs={
        'zscores': FeatureTable[Zscore],
        'epitope_map': FeatureData[MappedEpitope]
    },
    parameters={},
    outputs=[
        ('epitope_zscore', FeatureTable[Zscore]),
    ],
    input_descriptions={
        'zscores': 'FeatureTable containing the code names of peptides and '
                   'their per sample z scores',
        'epitope_map': 'FeatureTable containing epitopes and their associated '
                       'peptides and subtypes'
    },
    parameter_descriptions={},
    output_descriptions={
        'epitope_zscore': 'FeatureTable containing the epitopes and their per '
                          'sample z scores.'
    },
    name='zscore',
    description='Creates a map of epitopes to their max z-score within each '
                'sample. The maxes are taken by finding the per sample maxes '
                'among z scores of peptides associated with a given epitope',
)

plugin.methods.register_function(
    function=q2_epitope.taxa_to_epitope,
    inputs={
        'epitope': FeatureData[Epitope],
    },
    parameters={
        'collapse': Str % Choices(['Bacterial', 'Viral', 'Both'])
    },
    outputs=[
        ('epitope_gmt', GMT),
    ],
    input_descriptions={
        'epitope': 'Feature table containing at least SpeciesID, ClusterID, '
                   'and EpitopeWindow columns'
    },
    parameter_descriptions={},
    output_descriptions={
        'epitope_gmt': 'GMT mapping SpeciesIDs to associated epitopes.'
    },
    name='taxa to epitope',
    description='Creates a GMT file mapping SpeciesIDs to their associated '
                'epitopes',
)

plugin.methods.register_function(
    function=q2_epitope.enriched_subtypes,
    inputs={
        'scores': Collection[FeatureData[PSEAScores]],
        'subtypes': FeatureData[MappedEpitope],
    },
    parameters={
        'p_value': Float % Range(0, None),
        'enrichment_score': Float % Range(0, None),
        'include_negative_enrichment': Bool,
        'split_column': Str,
        'peptide_library': Str
    },
    outputs=[
        ('enriched', Collection[FeatureData[Enriched]]),
    ],
    input_descriptions={
        'scores': 'PSEAScores of peptides/epitopes.',
        'subtypes': 'subtypes'
    },
    parameter_descriptions={},
    output_descriptions={
        'enriched': 'Enriched subtypes.'
    },
    name='enriched subtypes',
    description='Counts which subtypes have been enriched',
)
