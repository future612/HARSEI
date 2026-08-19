# Data Layer Notes

## ISRIC / Ivushkin Global Soil Salinity Maps

| Item | Value |
|---|---|
| GEE collection | `projects/sat-io/open-datasets/global_soil_salinity` |
| Provider/listing | Community Catalog, derived from Ivushkin et al. / ISRIC data |
| Type | Image collection |
| Years used | 2000, 2002, 2005, 2009, 2016 |
| Nominal resolution | 250 m |
| License | CC BY 4.0 |
| Band used | First band, renamed to `isric_salinity` during export |
| Units/classes | Treated as an ordinal external salinity value/class. The local script detects whether values are categorical. |
| Resampling to HARSEI grid | Auto: categorical maps use mode/nearest; continuous maps use average. |
| Saline/non-saline rule | Auto: lowest ordinal salinity class is non-saline, higher classes are saline. |

## HARSEI/SRSI/RSEI inputs

| Layer | Source |
|---|---|
| SRSI | Band 5 of `YJQ_ecocomponents_YYYY.tif`; higher values indicate stronger salinity/surface-structure stress. |
| HARSEI | `revise\annual_harsei_outputs\rasters\HARSEI\YJQ_HARSEI_YYYY.tif`; higher values indicate better pressure-adjusted ecological quality. |
| RSEI | Reconstructed locally from NDVI, WET, NDBSI and LST using fixed-reference normalization and pooled PCA. |

## Validation logic

- SRSI should be positively correlated with external salinity.
- RSEI and HARSEI should be negatively correlated with external salinity.
- HARSEI is considered more salinity-sensitive than RSEI if it shows a clearer negative median gradient across external salinity classes and `HARSEI - RSEI` is more negative in saline than non-saline pixels.

