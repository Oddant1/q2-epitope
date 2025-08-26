# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from biom import Table

from q2_types.feature_table import BIOMV210Format


def create_epitope_map(epitope: pd.DataFrame) -> pd.DataFrame:
    epitope = _create_EpitopeID_row(epitope)
    epitope = epitope.reset_index()
    epitope = _create_SpeciesSubtype_row(epitope)

    mapped = epitope[['EpitopeID', 'CodeName', 'SpeciesSubtype']]
    mapped = mapped.groupby('EpitopeID')[['CodeName', 'SpeciesSubtype']].agg(list).reset_index()
    mapped['CodeName'] = mapped['CodeName'].transform(lambda x: ';'.join(x))
    mapped['SpeciesSubtype'] = mapped['SpeciesSubtype'].transform(lambda x: ';'.join(x))
    mapped.set_index('EpitopeID', inplace=True)

    return mapped


def zscore(
            scores: pd.DataFrame,
            epitope: pd.DataFrame
        ) -> BIOMV210Format:
    scores.fillna(value=0, axis=1, inplace=True)
    samples = list(scores.index)
    observations = list(epitope.index)

    data = []
    for _, row in epitope.iterrows():
        max_z_scores_per_sample = []
        split_row = row['CodeName'].split(';')

        # Filter the scores dataframe to only include columns corresponding to
        # the peptides we're looking at
        z_scores = scores.columns[scores.columns.isin(split_row)]

        for _, row in scores[z_scores.values].iterrows():
            max_z_scores_per_sample.append(max(row.values, key=abs))

        data.append(max_z_scores_per_sample)

    data = np.array(data)
    table = Table(data, observations, samples)

    result = BIOMV210Format()
    with result.open() as fh:
        table.to_hdf5(fh, generated_by="q2-pepsirf for pepsirf")

    return result


def taxa_to_epitope(epitope: pd.DataFrame) -> pd.DataFrame:
    epitope = _create_EpitopeID_row(epitope)

    mapped = epitope[['EpitopeID', 'SpeciesID']]
    mapped = mapped.groupby(['SpeciesID'])
    mapped = mapped['EpitopeID'].unique()
    mapped = mapped.reset_index()
    mapped.set_index('SpeciesID', inplace=True)

    return mapped


def _create_EpitopeID_row(epitope):
    epitope['SpeciesID'] = epitope['SpeciesID'].str.split(';')
    epitope['ClusterID'] = epitope['ClusterID'].fillna('clusterNA')
    epitope['ClusterID'] = epitope['ClusterID'].str.split(';')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].fillna('Peptide_NA')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].str.split(';')

    epitope = epitope.explode(['SpeciesID', 'ClusterID', 'EpitopeWindow'])
    epitope.drop_duplicates(inplace=True)

    def combine(row):
        return f"{row['SpeciesID']}_{row['ClusterID']}_{row['EpitopeWindow']}"

    epitope['EpitopeID'] = epitope.apply(combine, axis=1)

    return epitope


def _create_SpeciesSubtype_row(epitope):
    epitope['Species'] = epitope['Species'].fillna('speciesNA')
    epitope['Subtype'] = epitope['Subtype'].fillna('subtypeNA')

    def combine(row):
        combined_row = ''
        species = row['Species'].split(';')
        subtype = row['Subtype'].split(';')
        zipped = zip(species, subtype)

        for species, subtype in zipped:
            combined_row += f'{species}:{subtype};'

        return combined_row[:-1]

    epitope['SpeciesSubtype'] = epitope.apply(combine, axis=1)

    return epitope
