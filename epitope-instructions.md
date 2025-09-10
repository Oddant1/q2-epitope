# Install using conda, activate the environment

These install instructions assume that you have installed conda, e.g., through [Miniforge](https://github.com/conda-forge/miniforge).

Run the following command to create a new environment and install q2-epitope into it.

```
conda env create -n q2-epitope -f https://raw.githubusercontent.com/Oddant1/q2-epitope/refs/heads/main/env.yml
```

Activate the new environment:

```
conda activate q2-epitope
```

# Prepping epitope data with q2-epitope

The following steps use provided example data. We haven't yet put this into a GitHub repo/etc, as we didn't know if these were ok to share.

1. Import peptipe/epitope metadata (I think this artifact only really needs to be created once per peptide/epitope library?):

`qiime tools import --type FeatureData[Epitope] --input-path IN2_metadata.tsv --output-path epitope.qza`

2. Import peptide z-scores file:

`qiime tools import --type FeatureTable[Zscore] --input-path IM0206_IN2_raw_3mm_Z-HDI75.tsv --output-path zscores.qza`

3. Run `qiime epitope create-epitope-map`. This will put together the full epitope identifiers species_clusterID_EpitopeWindow and map them to the peptide codenames and species/subtypes those epitopes are associated with:

`qiime epitope create-epitope-map --i-epitope epitope.qza --o-epitope-map epitope-map.qza`

4. Run `qiime epitope epitope-zscore`. This maps the epitopes to their max z-scores within each sample. The maxes are taken by finding the per sample maxes among z-scores of peptides associated with a given epitope:

`qiime epitope epitope-zscore --i-zscores zscores.qza --i-epitope-map epitope-map.qza --o-epitope-zscore epitope-zscores.qza`

5. Run `qiime epitope taxa-to-epitope`. Creates `gmt` file (as a QIIME 2 Artifact) mapping the speciesIDs to associated epitopes (not associated peptides):

`qiime epitope taxa-to-epitope --i-epitope epitope.qza --o-epitope-gmt epitope-gmt.qza`

# Using epitope data in q2-PSEA

For the time being, you will unfortunately have to unzip your artifacts to use them with q2-PSEA. You can do this using `qiime tools export`.

When we got to it, q2-PSEA did not take QIIME 2 Artifacts as inputs, but rather used raw filepaths and we did not refactor that as part of this work. You will need to export the three artifacts made with the Actions above.

`qiime tools export --input-path epitope-map.qza --output-path epitope-map`

`qiime tools export --input-path epitope-zscores.qza --output-path epitope-zscores`

`qiime tools export --input-path epitope-gmt.qza --output-path epitope-gmt`

Then run q2-psea make-psea-table:

`qiime psea make-psea-table --p-scores-file ./epitope-zscores/pepsirf-table.tsv --p-peptide-sets-file ./epitope-gmt/map.gmt --p-mapped-epitope-file ./epitope-map/mapped-epitope.tsv --p-species-taxa-file PV2species.tsv --p-pairs-file pairs.tsv --p-threshold .05 --p-nes-thresh 0 --output-dir make-psea-table-out`

All other parameters set as desired/as would be necessary previously. The outputs will now be generated using per epitope data. The new enriched subtype outputs will be saved to the folder specified by the `--p-enriched-subtypes-dir` (defaults to `./psea_enriched_subtypes_tables`).
