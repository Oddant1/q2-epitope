# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import qiime2
import numpy as np
import pandas as pd

from biom.table import Table

from q2_pepsirf.format_types import (PepsirfContingencyTSVFormat,
                                     EpitopeFormat, MappedEpitopeFormat)
from q2_types.feature_table import BIOMV210Format


def create_epitope_map(epitope: pd.DataFrame) -> pd.DataFrame:
    epitope['SpeciesID'] = epitope['SpeciesID'].str.split(';')
    epitope['ClusterID'] = epitope['ClusterID'].fillna('clusterNA')
    epitope['ClusterID'] = epitope['ClusterID'].str.split(';')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].fillna('Peptide_NA')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].str.split(';')

    epitope = epitope.explode(['SpeciesID', 'ClusterID', 'EpitopeWindow'])

    epitope.drop_duplicates(inplace=True)

    def combine(row):
        return f"{row['SpeciesID']}_{row['ClusterID']}_{row['EpitopeWindow']}"

    epitope['Epitope'] = epitope.apply(combine, axis=1)

    mapped = epitope[['Epitope', 'CodeName']]

    mapped = epitope.groupby('Epitope')['CodeName'].agg(list).reset_index()
    mapped['CodeName'] = mapped['CodeName'].transform(lambda x: ';'.join(x))

    return mapped


# We need a good format for the cleaned up metadata
def zscore(
        scores: pd.DataFrame,
        epitope: pd.DataFrame
    ) -> Table:
    observations = list(epitope.index)
    samples = list(scores.index)

    data = []
    for _, row in epitope.iterrows():
        max_z_scores_per_sample = []
        split_row = row['CodeName'].split(';')
        z_scores = scores.columns[scores.columns.isin(split_row)]
        for _, row in scores[z_scores.values].iterrows():
            max_z_scores_per_sample.append(max(row.values, key=abs))

        data.append(max_z_scores_per_sample)

    data = np.array(data)
    table = Table(data, observations, samples)

    return table
