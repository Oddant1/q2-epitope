# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import q2_epitope
from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch, TypeMap)
from q2_types.feature_table import FeatureTable
from q2_pepsirf.format_types import Zscore


plugin = Plugin(
    name='epitope',
    version=q2_epitope.__version__,
    website='https://github.com/qiime2/q2-epitope',
    package='q2_epitope',
    short_description=(''),
    description=('')
)

plugin.methods.register_function(
    function=q2_epitope.zscore,
    inputs={
        "scores": FeatureTable[Zscore],
    },
    parameters={
        "metadata": qiime2.plugin.Metadata,
    },
    outputs=[
        ("zscore_map", FeatureTable[Zscore]),
    ],
    input_descriptions={"scores": ""},
    parameter_descriptions={"metadata": ""},
    output_descriptions={"zscore_map": ""},
    name="zscore",
    description="zscore",
)
