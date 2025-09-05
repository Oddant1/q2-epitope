# Install and activate environment

Run the following two commands:

```
conda env create -n q2-epitope -f https://raw.githubusercontent.com/Oddant1/q2-epitope/refs/heads/main/env.yml
conda activate q2-epitope
```

# Prepping epitope data with q2-epitope

The file names in the import are the file names I was sent. I wasn't sure if it was ok to put the files in a public GitHub repo or not, so I didn't.

1. Import peptipe/epitope metadata (I think this artifact only really needs to be created once per peptide/epitope library?):

`qiime tools import --type FeatureTable[Epitope] --input-path IN2_metadata.tsv --output-path epitope.qza`

2. Import peptide z-scores file:

`qiime tools import --type FeatureTable[Zscore] --input-path IM0206_IN2_raw_3mm_Z-HDI75.tsv --output-path zscores.qza`

3. Run qiime epitope create-epitope-map:

`qiime epitope create-epitope-map --i-epitope epitope.qza --o-epitope-map epitope-map.qza`

4. Run qiime epitope epitope-zscore:

`qiime epitope epitope-zscore --i-zscores zscores.qza --i-epitope-map epitope-map.qza --o-epitope-zscore epitope-zscores.qza`

5. Run qiime epitope taxa-to-epitope:

`qiime epitope taxa-to-epitope --i-epitope epitope.qza --o-epitope-gmt epitope-gmt.qza`

# Using epitope data in q2-PSEA

For the time being, you will unfortunately have to unzip your artifacts to use them with q2-PSEA. When I got to it, q2-PSEA did not take any Artifacts as inputs, it used raw filepaths. I did not refactor that as part of this work. You can do this usings `qiime tools export`. You will need to export the three artifacts made with the Actions above.

`qiime tools export --input-path epitope-map.qza --output-path epitope-map`

`qiime tools export --input-path epitope-zscores.qza --output-path epitope-zscores`

`qiime tools export --input-path epitope-gmt.qza --output-path epitope-gmt`

Then run q2-psea make-psea-table:

`qiime psea make-psea-table --p-scores-file ./epitope-zscores/pepsirf-table.tsv --p-peptide-sets-file ./epitope-gmt/map.gmt --p-mapped-epitope-file ./epitope-map/mapped-epitope.tsv --p-species-taxa-file ./PV2species.tsv`

all other parameters set as desired/as would be necessary previously. The outputs will now be generated using per epitope data. The new enriched subtype outputs will be saved to the folder specified by the `--p-enriched-subtypes-dir` defaults to `./psea_enriched_subtypes_tables`.
