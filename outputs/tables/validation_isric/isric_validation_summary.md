# ISRIC External Salinity Validation Summary

## Inputs
- Years: 2000, 2002, 2005, 2009, 2016
- Salinity directory: `D:\Codex\260724 小论文\revise\external_validation_inputs\isric_global_soil_salinity_remote_crops`
- Eco-components directory: `D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024`
- HARSEI directory: `D:\Codex\260724 小论文\revise\annual_harsei_outputs\rasters\HARSEI`

## RSEI reconstruction
- RSEI pooled PC1 explained variance ratio: 0.7429
- RSEI PC1 loadings: NDVI=0.6603, WET=0.3985, NDBSI=0.5169, LST=0.3715

## Output tables
- `isric_spearman_correlations.csv`
- `isric_salinity_class_index_statistics.csv`
- `isric_saline_vs_nonsaline_harsei_minus_rsei.csv`
- `isric_harsei_vs_rsei_salinity_gradient_test.csv`

## Manuscript interpretation guide
- SRSI should show a positive Spearman correlation with external salinity.
- HARSEI should show a negative Spearman correlation with external salinity.
- HARSEI - RSEI is a diagnostic difference; its sign depends on how salinity stress and human-activity adjustment interact.
- A clearer HARSEI salinity gradient is supported when HARSEI has a more negative median slope across salinity classes than RSEI.
