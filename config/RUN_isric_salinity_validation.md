# ISRIC Global Soil Salinity External Validation

## Research Design

| Decision | Choice | Reason | Needs confirmation |
|---|---|---|---|
| Study area | `projects/jlu2024621038/assets/yjq` | Same study-area boundary as the HARSEI reconstruction. | Confirm in GEE if the asset loads. |
| External product | `projects/sat-io/open-datasets/global_soil_salinity` | ISRIC/Ivushkin et al. Global Soil Salinity Maps, 250 m, CC BY 4.0. | Print metadata in GEE before export. |
| Validation years | 2000, 2002, 2005, 2009, 2016 | Years overlapping external product and annual HARSEI/SRSI. | None. |
| Export target | Google Drive folder `YJQ_ISRIC_salinity_validation` | User can download GeoTIFFs manually and place them locally. | Download to the exact local folder below. |
| Local input folder | `D:\Codex\260724 小论文\revise\external_validation_inputs\isric_global_soil_salinity` | Local script searches this folder by year. | Required before running local validation. |
| Analysis grid | 1000 m HARSEI grid | SRSI, RSEI and HARSEI are compared on the same annual index grid. | None. |

## Workflow

1. Open `gee_04_export_isric_salinity.js` in the GEE Code Editor.
2. Run the script and inspect the Console output. It prints collection size, `system:index` values and selected image metadata.
3. Start the five export tasks for 2000, 2002, 2005, 2009 and 2016.
4. Download the exported GeoTIFFs from Google Drive.
5. Put them in:

```text
D:\Codex\260724 小论文\revise\external_validation_inputs\isric_global_soil_salinity
```

Expected file names can be:

```text
ISRIC_global_soil_salinity_2000_YJQ.tif
ISRIC_global_soil_salinity_2002_YJQ.tif
ISRIC_global_soil_salinity_2005_YJQ.tif
ISRIC_global_soil_salinity_2009_YJQ.tif
ISRIC_global_soil_salinity_2016_YJQ.tif
```

The local script only requires that the filename contains the year.

6. Run:

```powershell
python runs\20260730-isric-salinity-validation\local_14_isric_external_validation.py
```

## Outputs

The script writes results to:

```text
D:\Codex\260724 小论文\revise\isric_salinity_validation
```

Main tables:

- `isric_spearman_correlations.csv`
- `isric_salinity_class_index_statistics.csv`
- `isric_saline_vs_nonsaline_harsei_minus_rsei.csv`
- `isric_harsei_vs_rsei_salinity_gradient_test.csv`
- `rsei_pooled_pca_pc1.csv`

If matplotlib is available, it also creates:

- `fig_isric_salinity_gradient_medians.png`

