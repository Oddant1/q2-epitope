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

from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from q2_types.feature_table import BIOMV210Format


def create_epitope_map() -> pd.DataFrame:
    df = pd.read_csv("/home/anthony/Downloads/IN2_metadata.tsv", sep="\t", low_memory=False)

    df['SpeciesID'] = df['SpeciesID'].str.split(';')
    df['ClusterID'] = df['ClusterID'].fillna('clusterNA')
    df['ClusterID'] = df['ClusterID'].str.split(';')
    df['EpitopeWindow'] = df['EpitopeWindow'].fillna('Peptide_NA')
    df['EpitopeWindow'] = df['EpitopeWindow'].str.split(';')

    df = df.explode(['SpeciesID', 'ClusterID', 'EpitopeWindow'])

    df.drop_duplicates(inplace=True)

    def combine(row):
        return f"{row['SpeciesID']}_{row['ClusterID']}_{row['EpitopeWindow']}"

    df['Epitope'] = df.apply(combine, axis=1)

    df2 = df[['Epitope', 'CodeName']]

    df2 = df.groupby('Epitope')['CodeName'].agg(list).reset_index()
    df2['CodeName'] = df2['CodeName'].transform(lambda x: ';'.join(x))

    df2.to_csv('~/df2.tsv', sep='\t')


# We need a good format for the cleaned up metadata
def zscore(
        scores: PepsirfContingencyTSVFormat,
        metadata: qiime2.Metadata
    ) -> Table:
    scores = "%s" % os.path.abspath(str(scores))
    df_z = pd.read_csv(scores, sep='\t', index_col=0)
    df_epitope = metadata.to_dataframe()

    observations = list(df_epitope.index)
    samples = list(df_z.columns[1:])

    data = []

    for _, row in df_epitope.iterrows():
        max_z_scores_per_sample = []
        z_scores = df_z[df_z.index.isin(row['CodeName'].split(';'))]

        for column in z_scores.columns[1:]:
            max_z_scores_per_sample.append(max(z_scores[column].values, key=abs))

        data.append(max_z_scores_per_sample)

    data = np.array(data)
    table = Table(data, observations, samples)

    with open('./delete.tsv', 'w') as fh:
        table.to_tsv(direct_io=fh, observation_column_name='Sequence name')

    return table
