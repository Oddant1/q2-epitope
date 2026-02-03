# ----------------------------------------------------------------------------
# Copyright (c) 2025-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import ast
import numpy as np
import pandas as pd
from biom import Table

from q2_types.feature_table import BIOMV210Format


def create_epitope_map(epitope: pd.DataFrame, collapse:str='Viral') -> pd.DataFrame:
    epitope = _create_EpitopeID_row(epitope, collapse)
    epitope = epitope.reset_index()

    epitope = \
        epitope.groupby(
            'EpitopeID').agg(list).reset_index()
    epitope.set_index('EpitopeID', inplace=True)

    def validate_categories(row):
        '''
        Ensure the categories column is actually valid. After aggregating some
        rows will have a list of multiple categories, these should just be the
        same category multiple times.

        1. Ensure that is the case
        2. Make it just a single value not a list
        '''
        category_set = set(row['Category'])
        if len(category_set) > 1:
            raise ValueError(
                'Collapsed epitope mapped some subtypes to one category and '
                f'some to another. Offending row is: {row}'
            )

        return row['Category'][0]

    epitope['Category'] = epitope.apply(validate_categories, axis=1)

    return epitope


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
    epitope['Species'] = epitope['Species'].str.split(';')
    epitope['Subtype'] = epitope['Subtype'].str.split(';')
    epitope['SpeciesID'] = epitope['SpeciesID'].str.split(';')
    epitope['ClusterID'] = epitope['ClusterID'].str.split(';')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].str.split(';')

    epitope = epitope.explode(
        ['Species', 'Subtype', 'SpeciesID', 'ClusterID', 'EpitopeWindow']
    )

    epitope['Subtype'] = epitope['Subtype'].fillna('subtypeNA')
    # Happens because some subtype rows have an empty value along with valid
    # values indicating one or more missing
    epitope['Subtype'] = epitope['Subtype'].replace('', 'subtypeNA')
    epitope['ClusterID'] = epitope['ClusterID'].fillna('clusterNA')
    epitope['EpitopeWindow'] = epitope['EpitopeWindow'].fillna('Peptide_NA')

    epitope.drop_duplicates(inplace=True)

    def combine(row):
        if collapse == 'Both' or row['Category'] == collapse:
            return \
                f"{row['SpeciesID']}_{row['ClusterID']}_{row['EpitopeWindow']}"

        return row.name

    epitope['EpitopeID'] = epitope.apply(combine, axis=1)

    return epitope


# Take those peptides and only collapse those to epitope and get subtypes (should be same workflow but on already filtered data)

# NOTE: At the end we are going to run tens out thousands of samples so we will want to make things efficient
def enriched_subtypes(
        scores: pd.DataFrame, subtypes: pd.DataFrame, p_value: float = .05,
        enrichment_score: float = 1, include_negative_enrichment: bool = True,
        split_column: str = None) -> pd.DataFrame:

    split_values = []
    keys = ['species-peptide', 'species-epitope', 'subspecies-peptide']
    subtypes['Subtype'] = subtypes['Subtype'].apply(ast.literal_eval)
    subtypes['CodeName'] = subtypes['CodeName'].apply(ast.literal_eval)

    if split_column is not None:
        if split_column not in subtypes:
            raise KeyError(
                f"The requested split_column: '{split_column}' does not "
                f"exist. Valid columns are: {list(subtypes.columns)}"
            )

        split_values = subtypes[split_column].unique()

        new_keys = []
        for value in split_values:
            for key in keys:
                new_keys.append(f'{value}-{key}')

        keys = new_keys

    scores = pd.concat(list(scores.values()))
    scores = scores.loc[scores['p.adjust'] <= p_value]

    if include_negative_enrichment:
        scores = scores.loc[abs(scores['enrichmentScore'] >= enrichment_score)]
    else:
        scores = scores.loc[scores['enrichmentScore'] >= enrichment_score]

    counts = {key: {} for key in keys}
    for _, row in scores.iterrows():
        peptides = row['core_enrichment'].split('/')
        species = row['species_name']

        for peptide in peptides:
            possibles = \
                subtypes.loc[
                    subtypes['CodeName'].apply(
                        lambda peptides: peptide in peptides)]
            for _, possible in possibles.iterrows():
                index = possible['CodeName'].index(peptide)
                subtype = possible['Subtype'][index]
                species_subtype = f'{species}:{subtype}'
                epitope = possible.name
                if split_column is not None:
                    found_split_value = possible[split_column]
                    if isinstance(found_split_value, list):
                        found_split_value = found_split_value[index]

                    if not f'{found_split_value}-{species}-{peptide}' in counts[f'{found_split_value}-species-peptide']:
                        counts[f'{found_split_value}-species-peptide'][f'{found_split_value}-{species}-{peptide}'] = 0
                    counts[f'{found_split_value}-species-peptide'][f'{found_split_value}-{species}-{peptide}'] += 1

                    if not f'{found_split_value}-{species}-{epitope}' in counts[f'{found_split_value}-species-epitope']:
                        counts[f'{found_split_value}-species-epitope'][f'{found_split_value}-{species}-{epitope}'] = 0
                    counts[f'{found_split_value}-species-epitope'][f'{found_split_value}-{species}-{epitope}'] +=1

                    if not f'{found_split_value}-{species_subtype}-{peptide}' in counts[f'{found_split_value}-subspecies-peptide']:
                        counts[f'{found_split_value}-subspecies-peptide'][f'{found_split_value}-{species_subtype}-{peptide}'] = 0
                    counts[f'{found_split_value}-subspecies-peptide'][f'{found_split_value}-{species_subtype}-{peptide}'] += 1
                else:
                    if not f'{species}-{peptide}' in counts['species-peptide']:
                        counts['species-peptide'][f'{species}-{peptide}'] = 0
                    counts['species-peptide'][f'{species}-{peptide}'] += 1

                    if not f'{species}-{epitope}' in counts['species-epitope']:
                        counts['species-epitope'][f'{species}-{epitope}'] = 0
                    counts['species-epitope'][f'{species}-{epitope}'] +=1

                    if not f'{species_subtype}-{peptide}' in counts['subspecies-peptide']:
                        counts['subspecies-peptide'][f'{species_subtype}-{peptide}'] = 0
                    counts['subspecies-peptide'][f'{species_subtype}-{peptide}'] += 1

    # TODO: Sort this by values
    for key, value in counts.items():
        counts[key] = pd.DataFrame({'Counts': value.values()}, index=value.keys())

    return counts
