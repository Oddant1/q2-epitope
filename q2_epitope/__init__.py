# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

__all__ = ['epitope_zscore', 'create_epitope_map', 'taxa_to_epitope',
           'enriched_subtypes']

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'

from .epitope import (epitope_zscore, create_epitope_map, taxa_to_epitope,
                      enriched_subtypes)
