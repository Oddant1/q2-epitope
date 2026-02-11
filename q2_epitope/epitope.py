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


def create_epitope_map(
            epitope: pd.DataFrame, collapse:str='Viral'
        ) -> pd.DataFrame:
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

        # Filter the scores dataframe to only include columns corresponding to
        # the peptides we're looking at
        z_scores = zscores.columns[zscores.columns.isin(row['CodeName'])]

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
            epitope: pd.DataFrame, collapse: str='Viral'
        ) -> pd.DataFrame:
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


def enriched_subtypes(
            scores: pd.DataFrame, subtypes: pd.DataFrame, p_value: float = .05,
            enrichment_score: float = 1,
            include_negative_enrichment: bool = True,
            split_column: str = None,
            peptide_library: str = 'IN2'
        ) -> pd.DataFrame:
    split_values = []
    keys = ['species-peptide', 'species-epitope', 'subspecies-peptide']

    # NOTE: This will definitely work if Category is chosen as split column
    # which was the intended usage. If a column that is formatted wildly
    # different from that is chosen, things may get hairy.
    #
    # Here we assume that a split column will have cells that are either a
    # single value or a list of values of the same length as the list of
    # peptides for its row
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
        # Get core_enrichement from psea_out
        enriched_elements = row['core_enrichment'].split('/')
        species = row['species_name']

        for enriched in enriched_elements:
            # Determine if we are looking at something that has been collapsed
            # to epitope or not
            if enriched.startswith(peptide_library):
                # uncollapsed
                peptide = enriched
                hits = \
                    subtypes.loc[
                        subtypes['CodeName'].apply(
                            lambda peptides: peptide in peptides)]

                for _, hit in hits.iterrows():
                    index = hit['CodeName'].index(peptide)
                    subtype = hit['Subtype'][index]
                    species_subtype = f'{species}:{subtype}'
                    # NOTE: If this is uncollapsed the epitope will just be the
                    # peptide (if this is a bacterium and bacteria was not
                    # collapsed for instance)
                    epitope = hit.name

                    found_split_value = \
                        _find_split_value(hit, split_column, index)
                    _count_enriched(counts, species, species_subtype, epitope,
                                    peptide, found_split_value)
            else:
                # collapsed
                epitope = enriched
                hit = subtypes.loc[epitope]
                for index, (peptide, subtype) in enumerate(
                        zip(hit['CodeName'], hit['Subtype'])):
                    species_subtype = f'{species}:{subtype}'

                    found_split_value = \
                        _find_split_value(hit, split_column, index)
                    _count_enriched(counts, species, species_subtype, epitope,
                                    peptide, found_split_value)

    for key, value in counts.items():
        _sorted = \
            dict(sorted(value.items(), key=lambda item: item[1], reverse=True))
        counts[key] = \
            pd.DataFrame({'Counts': _sorted.values()}, index=_sorted.keys())

    return counts


def _find_split_value(hit, split_column, index):
    found_split_value = ''
    if split_column is not None:
        found_split_value = hit[split_column]
        if isinstance(found_split_value, list):
            found_split_value = found_split_value[index]
        found_split_value += '-'

    return found_split_value


def _count_enriched(counts, species, species_subtype, epitope, peptide,
                    found_split_value):
    # Track species and peptide including split value if relevant
    split_species_peptide = f'{found_split_value}species-peptide'
    found_split_species_peptide = \
        f'{found_split_value}{species}-{peptide}'

    if not found_split_species_peptide in \
            counts[split_species_peptide]:
        counts[split_species_peptide]\
            [found_split_species_peptide] = 0
    counts[split_species_peptide][found_split_species_peptide] += 1

    # Track species and epitope including split value if relevant
    split_species_epitope = f'{found_split_value}species-epitope'
    found_split_species_epitope = \
        f'{found_split_value}{species}-{epitope}'

    if not found_split_species_epitope in \
            counts[split_species_epitope]:
        counts[split_species_epitope]\
            [found_split_species_epitope] = 0
    counts[split_species_epitope][found_split_species_epitope] +=1

    # Track subspecies and peptide including split value if
    # relevant
    split_sub_peptide = f'{found_split_value}subspecies-peptide'
    found_split_sub_peptide = \
        f'{found_split_value}{species_subtype}-{peptide}'

    if not found_split_species_peptide in \
            counts[split_sub_peptide]:
        counts[split_sub_peptide][found_split_sub_peptide] = 0
    counts[split_sub_peptide][found_split_sub_peptide] += 1
