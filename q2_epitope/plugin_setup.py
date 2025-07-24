# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import q2_epitope
from qiime2.plugin import Plugin
from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import Zscore, Epitope, MappedEpitope, GMT


plugin = Plugin(
    name='epitope',
    version=q2_epitope.__version__,
    website='https://github.com/qiime2/q2-epitope',
    package='q2_epitope',
    short_description=(''),
    description=('')
)

plugin.methods.register_function(
    function=q2_epitope.create_epitope_map,
    inputs={'epitope': FeatureTable[Epitope]},
    parameters={},
    outputs=[
        ('epitope_map', FeatureTable[MappedEpitope])
    ],
    input_descriptions={'epitope': ''},
    parameter_descriptions={},
    output_descriptions={'epitope_map': ''},
    name='create epitope map',
    description='create epitope map'
)

plugin.methods.register_function(
    function=q2_epitope.zscore,
    inputs={
        'scores': FeatureTable[Zscore],
        'epitope': FeatureTable[MappedEpitope]
    },
    parameters={},
    outputs=[
        ('zscore_map', FeatureTable[Zscore]),
    ],
    input_descriptions={'scores': '', 'epitope': ''},
    parameter_descriptions={},
    output_descriptions={'zscore_map': ''},
    name='zscore',
    description='zscore',
)

plugin.methods.register_function(
    function=q2_epitope.taxa_to_epitope,
    inputs={
        'epitope': FeatureTable[Epitope],
    },
    parameters={},
    outputs=[
        ('epitope_map', GMT),
    ],
    input_descriptions={'epitope': ''},
    parameter_descriptions={},
    name='epitope_map',
    description='epitop_map',
)
