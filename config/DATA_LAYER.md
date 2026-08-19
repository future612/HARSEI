# Data Layer Notes

## ROI

- GEE asset: `projects/jlu2024621038/assets/yjq`
- Use as exact export region and analysis boundary.

## Annual Ecological Components

| Component | Dataset | Bands | Unit / scale | Notes |
|---|---|---|---|---|
| NDVI | `MODIS/061/MOD13A1` | `NDVI`, `SummaryQA` | scale factor 0.0001 | Growing-season mean after `SummaryQA <= 1` mask |
| Surface reflectance | `MODIS/061/MOD09A1` | `sur_refl_b01`-`sur_refl_b07`, `QA`, `StateQA` | scale factor 0.0001 | Used for WET, NDBSI, SRSI and NDWI water mask; masked for cloud, cloud shadow, cirrus and snow/ice |
| WET | derived from MOD09A1 | red, NIR, blue, green, SWIR1, SWIR2, SWIR3 | unitless | Uses MODIS tasseled-cap wetness coefficients from Lobser and Cohen (2007): `0.1147, 0.2489, 0.2408, 0.3132, -0.3122, -0.6416, -0.5087` |
| NDBSI | derived from MOD09A1 | red, NIR, blue, green, SWIR1 | unitless | Average of SI and IBI |
| LST | `MODIS/061/MOD11A2` | `LST_Day_1km`, `QC_Day` | K * 0.02, then Celsius | Uses product LST, not Landsat thermal retrieval |
| SRSI | derived from MOD09A1 | green, red, NIR | unitless | Candidate retained from submitted manuscript: `sqrt((NDVI_SR - 1)^2 + SI1^2)`, where `SI1 = sqrt(green * red)`; final formula/citation must be confirmed |
| Salinity diagnostics | derived from MOD09A1 | green, red, blue, NIR, SWIR1 | unitless | Extra bands exported for formula checking: `SI1`, `NDSI_SAL`, `SI3_GRB`, `SI_SWIR_BLUE`, `NDVI_SR` |
| Water mask | derived from MOD09A1 | green, NIR | unitless | `NDWI <= 0.2` retained as land in the main run; `gee_03` tests thresholds 0.0, 0.1, 0.2, 0.3 |
| Snow diagnostic | `MODIS/061/MOD13A1` and `MODIS/061/MOD09A1` | `SummaryQA`, `StateQA` | fraction | `snow_frac` from MOD13A1 and `mod09_snow_frac` from MOD09A1 StateQA are exported to evaluate high-elevation snow contamination |

## Required HAI Inputs

| Component | Dataset | Bands | Notes |
|---|---|---|---|
| Population | `WorldPop/GP/100m/pop` | `population` | Current downloaded stack uses annual data through 2020; 2021-2024 carry 2020 forward unless replaced |
| Night lights, DMSP | `NOAA/DMSP-OLS/NIGHTTIME_LIGHTS` | `stable_lights` | Used through 2013 |
| Night lights, VIIRS | `NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG` | `avg_rad`, `cf_cvg` | Annual VIIRS composite is coverage-weighted by `cf_cvg`; used from 2014 onward after local overlap harmonization |
| Land-use disturbance, submitted-code mode | `projects/sat-io/open-datasets/GLC-FCS30D/annual` | annual bands `b1`-`b23` | Matches the user's pasted LUCC scoring code and Zhang et al. (2024); b1-b23 correspond to 2000-2022, so 2023-2024 require carry-forward or replacement |
| Land-use disturbance, optional CLCD mode | user-uploaded CLCD asset | class band | Use only if the manuscript retains CLCD; the latest public CLCD record includes 2024, but CLCD must be uploaded because it is not a built-in GEE catalog layer |
| Land-use disturbance, fallback mode | `MODIS/061/MCD12Q1` | `LC_Type1` | Official fallback product available through 2024; lower spatial resolution and different land-cover taxonomy |

## HAI Land-Use Scoring

The LUCC component follows the Chinese manuscript's rule before fixed normalization:

| Land-cover type | GLC-FCS30D code | CLCD code if using CLCD | HAI LUCC score |
|---|---:|---:|---:|
| Built-up / impervious land | 190 | 8 | 10 |
| Cropland / farmland | 10, 11, 12 | 1 | 7 |
| Grassland | 130 | 4 | 4 |
| Water | 210 | 5 | 1 |
| Remaining classes | all other classes | 2, 3, 6, 7, 9, other | 0 |

The 0-10 LUCC score is then normalized together with population and harmonized night-time lights in the local HAI workflow. The NDWI water mask is used to define the final common valid-analysis pixels for ecological components and HARSEI fusion; the water score of 1 documents the low HAI pressure of water in the raw LUCC component and does not mean that large open-water pixels are retained in the final ecological PCA.

## Local Download Folder

Put every file exported from GEE into:

```text
D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024
```

Required files:

- `YJQ_ecocomponents_2000.tif` through `YJQ_ecocomponents_2024.tif`
- `YJQ_hai_inputs_2000.tif` through `YJQ_hai_inputs_2024.tif`
- `YJQ_pooled_pca_sample_2000_2024.csv`
- `YJQ_water_snow_sensitivity_2000_2024.csv`

## Items Requiring Manuscript Confirmation

1. SRSI: cite Alhammadi and Glenn (2008) and describe SRSI as a vegetation-salinity response indicator; do not describe it as identical to the Khan-family standard salinity indices.
2. LUCC source: the user's pasted code uses GLC-FCS30D, not CLCD. If the manuscript retains this code, revise the data-source statement from CLCD to GLC-FCS30D and disclose the 2024 handling. If the manuscript must retain CLCD, upload annual CLCD assets and set `LANDCOVER_MODE = 'CLCD_ASSET'`.
3. WorldPop 2024: the downloaded stack uses 2020 as the latest available WorldPop layer; 2021-2024 carry 2020 forward. This must be disclosed or replaced by a current population product.
4. Compositing window: the main code uses May-Sep. Use `gee_03` to decide whether the manuscript should fully switch from Apr-Sep to May-Sep or report Apr-Sep only as a sensitivity comparison.

## Fixed Comparability Rules

For the revised manuscript, do these after download:

1. Pool all years before normalization.
2. Use fixed percentile ranges across 2000-2024.
3. Direction-standardize ecological variables:
   - positive: NDVI, WET
   - negative: NDBSI, LST, SRSI
4. Compute one pooled PCA model for all years.
5. Use one fixed PC1 direction and report loadings.
6. Use fixed entropy weights for AWRSEI and reverse-HAI fusion.
