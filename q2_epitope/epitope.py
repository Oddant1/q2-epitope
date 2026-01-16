# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from biom import Table

from q2_types.feature_table import BIOMV210Format


def create_epitope_map(epitope: pd.DataFrame) -> pd.DataFrame:
    epitope = _create_EpitopeID_row(epitope)
    epitope = epitope.reset_index()
    epitope = _create_SpeciesSubtype_row(epitope)

    epitope_map = epitope[['EpitopeID', 'CodeName', 'SpeciesSubtype']]
    epitope_map = \
        epitope_map.groupby(
            'EpitopeID')[
                ['CodeName', 'SpeciesSubtype']].agg(list).reset_index()
    epitope_map['CodeName'] = \
        epitope_map['CodeName'].transform(lambda x: ';'.join(x))
    epitope_map['SpeciesSubtype'] = \
        epitope_map['SpeciesSubtype'].transform(lambda x: ';'.join(x))
    epitope_map.set_index('EpitopeID', inplace=True)

    return epitope_map


def epitope_zscore(
            zscores: pd.DataFrame,
            epitope_map: pd.DataFrame
        ) -> BIOMV210Format:
    zscores.fillna(value=0, axis=1, inplace=True)
    samples = list(zscores.index)
    observations = list(epitope_map.index)

    data = []
    for _, row in epitope_map.iterrows():
        max_z_scores_per_sample = []
        split_row = row['CodeName'].split(';')

        # Filter the scores dataframe to only include columns corresponding to
        # the peptides we're looking at
        z_scores = zscores.columns[zscores.columns.isin(split_row)]

        for _, row in zscores[z_scores.values].iterrows():
            max_z_scores_per_sample.append(max(row.values, key=abs))

        data.append(max_z_scores_per_sample)

    data = np.array(data)
    table = Table(data, observations, samples)

    result = BIOMV210Format()
    with result.open() as fh:
        table.to_hdf5(fh, generated_by="q2-pepsirf for pepsirf")

    return result


def taxa_to_epitope(
        epitope: pd.DataFrame, collapse: str='Viral') -> pd.DataFrame:
    epitope = _create_EpitopeID_row(epitope, collapse)

    mapped = epitope[['EpitopeID', 'SpeciesID']]
    mapped = mapped.groupby(['SpeciesID'])
    mapped = mapped['EpitopeID'].unique()
    mapped = mapped.reset_index()
    mapped.set_index('SpeciesID', inplace=True)

    return mapped


def _create_EpitopeID_row(epitope, collapse):
    epitope['SpeciesID'] = epitope['SpeciesID'].str.split(';')
    epitope['ClusterID'] = epitope['ClusterID'].fillna('clusterNA')
    epitope['ClusterID'] = epitope['ClusterID'].str.split(';')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].fillna('Peptide_NA')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].str.split(';')

    epitope = epitope.explode(['SpeciesID', 'ClusterID', 'EpitopeWindow'])
    epitope.drop_duplicates(inplace=True)

    def combine(row):
        if collapse == 'Both' or row['Category'] == collapse:
            return \
                f"{row['SpeciesID']}_{row['ClusterID']}_{row['EpitopeWindow']}"

        return row.name

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


# Run PSEA and get the significant species based on p.adjust < .05 or some specified threshold
# Look at core_enrichment column of that output for / separated peptide values for all output tables which are per sample
def filter_psea_outputs(in_dir: str, out_path: str, threshold: float=.05):
    filtered_df = pd.DataFrame()
    for fp in os.listdir(in_dir):
        fp = os.path.join(in_dir, fp)
        df = pd.read_csv(fp, sep='\t')
        df = df.loc[(df['p.adjust'] < threshold)]
        filtered_df = pd.concat([filtered_df, df])

    filtered_df.set_index('ID', inplace=True)
    filtered_df.to_csv(out_path, sep='\t')

    global_enriched = []
    for row in filtered_df['core_enrichment']:
        enriched = row.split('/')
        global_enriched.extend(enriched)

    global_enriched = set(global_enriched)


# Take those peptides and only collapse those to epitope and get subtypes (should be same workflow but on already filtered data)

# NOTE: At the end we are going to run tens out thousands of samples so we will want to make things efficient
