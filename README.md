# HARSEI Altai Mountains 2000-2024

This repository archives the Google Earth Engine (GEE) and Python scripts, key
parameters, and statistical outputs used for the revised HARSEI analysis in the
Altai Mountains.

## Contents

- `gee/`: GEE scripts for annual remote-sensing input export, HAI inputs, water
  and snow sensitivity, ISRIC salinity export, and irrigation-mask export.
- `src/01_preprocessing_harsei/`: local preprocessing, AWRSEI/HAI/HARSEI
  reconstruction, fixed normalization, pooled PCA, entropy fusion, and
  diagnostic validation scripts.
- `src/02_spatiotemporal_analysis/`: trend, CV, Hurst, transition-matrix, Moran
  and LISA scripts.
- `src/03_geodetector/`: GeoDetector factor and interaction detection scripts.
- `src/04_rf_shap/`: random forest and SHAP scripts.
- `src/05_validation/`: ISRIC salinity, pressure-zone overlay, classic RSEI, and
  Manas field-salinity validation scripts.
- `outputs/tables/`: key statistical outputs used in the revised manuscript and
  reviewer response.
- `outputs/figures/`: non-raster generated figures and figure notes.
- `config/`: run notes, parameter settings, and data-source descriptions.
- `data/data_manifest.csv`: public data sources, local data requirements, and
  excluded large-file notes.

## Main Reproducibility Parameters

- Study area GEE asset: `projects/jlu2024621038/assets/yjq`
- Period: 2000-2024 annual series
- Main compositing window: May-September
- Working resolution: 1000 m
- Water mask: large-water pixels removed using `NDWI > 0.2`
- AWRSEI inputs: NDVI, WET, NDBSI, LST, and SRSI
- HAI formula: `(LUCC + POP + LIGHT) / 3`
- HARSEI fusion: annual entropy-weight fusion of AWRSEI and reverse HAI
- GeoDetector: continuous drivers discretized into seven quantile classes; WRB
  retained as a categorical factor
- RF-SHAP: up to 100,000 valid pixels per target year, 25% test set, 2,000 SHAP
  samples, random forest with 160 trees

Full parameter details are provided in `config/parameters.md`.

## Large Files Not Included

Large GeoTIFF rasters and large intermediate pixel-sample tables are not stored
directly in this GitHub package. They can be regenerated from the GEE scripts
and public data sources listed in `config/data_sources.md` and
`data/data_manifest.csv`. This keeps the GitHub repository lightweight while
preserving the code and statistical outputs needed to audit the analysis.

The pooled PCA sample exported from GEE,
`YJQ_pooled_pca_sample_2000_2024.csv`, is approximately 36 MB in the local run
and was therefore not copied into the upload package. Re-export it using
`gee/gee_01_export_annual_inputs.js` if full rerunning is required.

## Notes For Reuse

Some Python scripts retain the local directory constants used in the manuscript
revision run. Before rerunning on another computer, update the `ROOT`,
`DRIVER_ROOT`, input, and output paths at the top of the relevant scripts, or
adapt them to your own folder structure. The GEE scripts are intended to be
copied directly into the GEE Code Editor after confirming access to the ROI
asset.

## Suggested Citation Statement In The Manuscript

The Google Earth Engine and Python scripts used for data preprocessing, HARSEI
construction, spatiotemporal analysis, GeoDetector analysis, RF-SHAP modelling,
key parameter settings, and main statistical outputs have been archived in a
public GitHub repository. Large raster data are not stored directly in the
repository because of file-size limitations; instead, their public sources, GEE
asset identifiers, processing parameters, and derived statistical tables are
provided.
