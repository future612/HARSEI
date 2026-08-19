# Key Parameter Settings

## Study Area And Time Range

| Item | Setting |
|---|---|
| Study area | `projects/jlu2024621038/assets/yjq` |
| Years | 2000-2024 |
| Target years for selected maps and driver analysis | 2000, 2005, 2010, 2015, 2020, 2024 |
| Spatial resolution | 1000 m |
| Main seasonal window | May-September |

## Data Preprocessing

| Item | Setting |
|---|---|
| MODIS surface reflectance | MOD09A1, scale factor 0.0001 |
| MODIS NDVI | MOD13A1, scale factor 0.0001 |
| MODIS LST | MOD11A2, `LST = 0.02 * LST_Day_1km - 273.15` |
| Cloud/snow screening | MODIS QA bands |
| Water mask | `NDWI > 0.2` treated as large water and removed from ecological PCA |
| Sensitivity test | NDWI thresholds 0.0, 0.1, 0.2 and 0.3; April-September versus May-September |

## AWRSEI And HARSEI Construction

| Component | Setting |
|---|---|
| AWRSEI variables | NDVI, WET, NDBSI, LST, SRSI |
| WET coefficients | MODIS tasseled-cap wetness coefficients of Lobser and Cohen (2007): `0.1147, 0.2489, 0.2408, 0.3132, -0.3122, -0.6416, -0.5087` |
| SRSI | Salinity/structure response index following Alhammadi and Glenn (2008) |
| Normalization | Fixed cross-year reference ranges for 2000-2024 |
| Dimensionality reduction | Pooled PCA model for annual comparability |
| HAI formula | `(LUCC + POP + LIGHT) / 3` |
| LUCC score | built-up = 10, cropland = 7, grassland = 4, water = 1, others = 0 |
| Night-light harmonization | DMSP-OLS and VIIRS harmonized using the 2012-2013 overlap |
| HARSEI fusion | Annual entropy-weight fusion of AWRSEI and reverse HAI |

## Spatiotemporal Analysis

| Method | Setting |
|---|---|
| Trend | Theil-Sen slope and Mann-Kendall significance test based on 25 annual observations |
| Variation | Coefficient of variation |
| Persistence | Hurst exponent based on annual HARSEI |
| Grade transition | Five-class HARSEI grade transition and three-group improved/stable/worsened summary |
| Spatial autocorrelation | Global Moran's I and LISA for six target years |

## GeoDetector

| Item | Setting |
|---|---|
| Dependent variable | Revised annual HARSEI |
| Explanatory variables | DEM, SLOPE, ASPECT, WRB, PRE, TEMMAX, TEMMIN, TEM, SOIL, AET, SWE |
| Excluded variables | HARSEI component variables and direct HAI constituents |
| Continuous-variable discretization | Seven quantile classes |
| Categorical variable | WRB retained as a categorical soil-class factor |
| Interaction detection | Standard GeoDetector interaction rules based on `q(X1)`, `q(X2)`, `q(X1 and X2)`, and `q(X1) + q(X2)` |

## Random Forest And SHAP

| Item | Setting |
|---|---|
| Target years | 2000, 2005, 2010, 2015, 2020, 2024 |
| Maximum model samples | 100,000 valid pixels per year |
| Train/test split | 75% / 25% |
| SHAP samples | 2,000 validation samples |
| Random seed | 20260811 + year |
| Model | `RandomForestRegressor` |
| Trees | 160 |
| Max depth | 20 |
| Max features | `sqrt` |
| Min samples leaf | 6 |
| Bootstrap | true |
| Metrics | RMSE and R2 |

