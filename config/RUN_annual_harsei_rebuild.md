# Annual HARSEI Rebuild Runbook

## Research Design

| Decision | Choice | Reason | Needs confirmation |
|---|---|---|---|
| Study area | `projects/jlu2024621038/assets/yjq` | User-provided GEE asset for the study boundary | Confirm this is the final manuscript ROI |
| Period | 2000-2024 annual sequence | Directly addresses reviewer criticism that six time points are insufficient | None |
| Main growing-season window | May-September (`START_MONTH = 5`, `END_MONTH_EXCLUSIVE = 10`) | Avoids April high-elevation snow contamination; can be compared with Apr-Sep sensitivity | Confirm whether manuscript will use May-Sep or Apr-Sep sensitivity |
| Output target | Google Drive GeoTIFF exports from GEE | User will download from GEE independently | Confirm Drive folder if you want a different name |
| Analysis scale | 1000 m | Matches original HARSEI grid and MOD11A2 native scale | None |
| Core method | Export annual physical/ecological component stacks, then compute fixed normalization and pooled PCA locally | Avoids year-by-year min-max/PCA and supports reviewer response on temporal comparability | None |

## Files

- `gee_01_export_annual_inputs.js`  
  Copy into the GEE Code Editor. It exports:
  - `YJQ_ecocomponents_YYYY.tif`, 2000-2024
  - `YJQ_pooled_pca_sample_2000_2024.csv`
  - Optional diagnostic bands inside each GeoTIFF: `NDWI`, `land_mask`, `snow_frac`, `mod09_snow_frac`, `SI1`, `NDSI_SAL`, `SI3_GRB`, `SI_SWIR_BLUE`, `NDVI_SR`

- `gee_02_export_annual_hai_inputs.js`  
  Copy into the GEE Code Editor after or alongside the ecological export. It exports:
  - `YJQ_hai_inputs_YYYY.tif`, 2000-2024

- `gee_03_water_snow_sensitivity.js`  
  Copy into the GEE Code Editor for the reviewer-requested sensitivity check. It exports:
  - `YJQ_water_snow_sensitivity_2000_2024.csv`

- `local_02_build_harsei_from_downloads.py`  
  Run after downloading the GEE outputs. It builds fixed-normalized AWRSEI and HARSEI.

- `local_03_summarize_water_snow_sensitivity.py`  
  Run after downloading `YJQ_water_snow_sensitivity_2000_2024.csv`. It writes a compact summary table and report-ready sentences for the reviewer response.

## GEE Steps

If an older `gee_01_export_annual_inputs.js` task is already running, use the revised script for the final manuscript outputs. If the old tasks only just started, cancel and rerun; if they are almost finished, keep them as a backup but do not use them as the main revision result because they lack the strengthened StateQA snow mask and salinity diagnostic bands.

1. Open Google Earth Engine Code Editor.
2. Paste `gee_01_export_annual_inputs.js`.
3. Keep these defaults for the main revision:
   - `ROI_ASSET = 'projects/jlu2024621038/assets/yjq'`
   - `START_YEAR = 2000`
   - `END_YEAR = 2024`
   - `START_MONTH = 5`
   - `END_MONTH_EXCLUSIVE = 10`
   - `SCALE = 1000`
   - `EXPORT_ECO_COMPONENTS = true`
   - `EXPORT_POOLED_SAMPLE = true`
   - `EXPORT_HAI_INPUTS = false` initially
4. Run the script.
5. In the Tasks panel, start the 25 image exports and the sample CSV export.
6. Download outputs from Google Drive folder `YJQ_HARSEI_annual_inputs_2000_2024`.
7. Put all downloaded files here:

```text
D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024
```

## Required HAI Export

HAI is required for the final annual HARSEI. Use the second GEE script:

1. Open `gee_02_export_annual_hai_inputs.js`.
2. Keep the default:
   - `LANDCOVER_MODE = 'MCD12Q1'`
3. Run it and start the 25 `YJQ_hai_inputs_YYYY.tif` tasks.
4. Download them into the same local folder as the ecological components:

```text
D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024
```

The default annual HAI uses:

- WorldPop population; 2021-2024 use 2020 as the nearest available population layer in the current downloaded stack.
- DMSP-OLS and VIIRS night-time lights; VIIRS monthly composites are weighted by `cf_cvg`, and local harmonization uses the 2012-2013 overlap.
- GLC-FCS30D annual land cover, matching the user's submitted LUCC code; annual bands currently cover 2000-2022, so 2023-2024 use 2022 unless replaced by CLCD/MCD12Q1.
- LUCC score follows the manuscript's 0-10 HAI rule: built-up = 10, cropland = 7, grassland = 4, water = 1, all other classes = 0.

## Land-Cover Mode Choice

The pasted original LUCC code uses:

```js
var annual = ee.ImageCollection("projects/sat-io/open-datasets/GLC-FCS30D/annual");
```

Therefore, if you keep that code path, revise the manuscript data source from CLCD to GLC-FCS30D and report that 2024 uses the nearest available 2022 land-cover layer unless you provide a true 2024 LUCC layer.

## Optional CLCD Mode

If you upload annual CLCD images to GEE and want to retain the original manuscript's LUCC source, change:

```js
var LANDCOVER_MODE = 'CLCD_ASSET';
var CLCD_ASSET_PREFIX = 'projects/jlu2024621038/assets/CLCD/CLCD_v01_';
var CLCD_ASSET_SUFFIX = '';
var CLCD_LAST_YEAR = 2024;
```

The latest public CLCD Zenodo record includes `CLCD_v01_2024_albert.tif`, so retaining CLCD is defensible if the revised manuscript cites the updated record and documents that 2024 was used. If your uploaded asset names include `_albert`, set:

```js
var CLCD_ASSET_SUFFIX = '_albert';
```

If the manuscript must keep CLCD rather than switching to MCD12Q1, provide the actual annual CLCD asset prefix/suffix before running HAI exports.

For the final manuscript route:

- Use `LANDCOVER_MODE = 'GLC_FCS30D'` if the paper should match your pasted code.
- Use `LANDCOVER_MODE = 'CLCD_ASSET'` if the paper should retain the CLCD wording and you upload annual CLCD assets including 2024.
- Use `LANDCOVER_MODE = 'MCD12Q1'` only as a reproducible fallback.

## Water Mask And Snow Sensitivity

Run `gee_03_water_snow_sensitivity.js` once. It exports a CSV comparing:

- Apr-Sep (`start_month = 4`) versus May-Sep (`start_month = 5`)
- NDWI water thresholds `0.0`, `0.1`, `0.2`, and `0.3`
- Total water/land area and snow frequency in the whole ROI and high-elevation pixels (`SRTM elevation >= 2500 m`)

Use this table in the response to justify either:

- retaining `NDWI > 0.2` as the water threshold if land/water area and HARSEI inputs are stable across thresholds; or
- revising the threshold if the sensitivity table shows large threshold dependence.

Use the Apr-Sep versus May-Sep comparison to justify the revised May-Sep main window, or to document that April snow contamination is negligible if you decide to keep Apr-Sep.

After downloading the sensitivity CSV:

```powershell
python .\runs\20260724-annual-harsei-rebuild\local_03_summarize_water_snow_sensitivity.py `
  --csv ".\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024\YJQ_water_snow_sensitivity_2000_2024.csv" `
  --out-dir ".\revise\annual_harsei_outputs\tables"
```

Outputs:

- `water_snow_sensitivity_summary.csv`
- `water_snow_sensitivity_summary.md`

## Local Rebuild Steps

After all `YJQ_ecocomponents_YYYY.tif` files are downloaded into one folder:

```powershell
python .\runs\20260724-annual-harsei-rebuild\local_02_build_harsei_from_downloads.py `
  --eco-dir ".\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024" `
  --hai-dir ".\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024" `
  --out-dir ".\revise\annual_harsei_outputs" `
  --hai-component-weights entropy
```

Expected outputs:

- `rasters/YJQ_AWRSEI_YYYY.tif`
- `rasters/YJQ_HAI_reverse_YYYY.tif` if HAI inputs are available
- `rasters/YJQ_HARSEI_YYYY.tif` if HAI inputs are available
- `tables/fixed_normalization_ranges.csv`
- `tables/pooled_pca_pc1.csv`
- `tables/entropy_weights_and_light_harmonization.csv`
- `tables/annual_index_summary.csv`

## Manuscript Use

Use the output tables to update:

- Methods: fixed-reference normalization, pooled PCA, fixed entropy weights.
- Methods: MODIS-specific WET coefficients; MOD11A2 LST product conversion; DMSP/VIIRS harmonization; WorldPop/land-cover temporal availability.
- Results 3.1: annual mean HARSEI/AWRSEI trajectory.
- Results 3.3/3.5: annual CV and Theil-Sen/Mann-Kendall trend.
- Results 3.6: Hurst only if annual HARSEI is completed; otherwise keep Hurst removed.
- Supplementary: PC1 loadings, explained variance, normalization ranges, entropy weights, night-light harmonization parameters, and NDWI/snow sensitivity table.
