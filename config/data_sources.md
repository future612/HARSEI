# Data Sources

## GEE And Public Datasets

| Dataset | Use |
|---|---|
| `MODIS/061/MOD13A1` | NDVI |
| `MODIS/061/MOD09A1` | WET, NDBSI, SRSI, NDWI mask, QA screening |
| `MODIS/061/MOD11A2` | Land surface temperature |
| `WorldPop/GP/100m/pop` | Population density for HAI |
| `NOAA/DMSP-OLS/NIGHTTIME_LIGHTS` | Night-time lights before VIIRS period |
| `NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG` | VIIRS night-time lights |
| `projects/sat-io/open-datasets/GLC-FCS30D/annual` | Land-use disturbance scoring |
| `IDAHO_EPSCOR/TERRACLIMATE` | PRE, TEM, TEMMAX, TEMMIN, SOIL, AET, SWE |
| `USGS/SRTMGL1_003` | DEM and high-elevation snow diagnostic |
| ISRIC Global Soil Salinity Maps | External salinity-product validation |
| GFSAD/LGRIP30 | Irrigated-area pressure-zone overlay |
| Global Mining Areas v2 | Mining-pressure-zone overlay |

## Key Literature Sources

- Lobser, S.E., Cohen, W.B., 2007. MODIS tasseled cap: land cover characteristics expressed through transformed MODIS data. International Journal of Remote Sensing 28, 5079-5101.
- Alhammadi, M.S., Glenn, E.P., 2008. Detecting date palm trees health and vegetation greenness change on the eastern coast of the United Arab Emirates using SAVI. International Journal of Remote Sensing 29, 1745-1765.
- Wang, J.F., Xu, C.D., 2017. Geodetector: Principle and prospective. Acta Geographica Sinica 72, 116-134.
- Breiman, L., 2001. Random forests. Machine Learning 45, 5-32.
- Lundberg, S.M., Lee, S.I., 2017. A unified approach to interpreting model predictions.

## Data Availability Notes

Large derived GeoTIFF rasters are not included in the GitHub package. They can be
regenerated from the GEE scripts and local Python scripts. If a journal or
reviewer requests the actual raster layers, archive the rasters separately using
Zenodo, Figshare, an institutional repository, or GitHub Releases/Git LFS.

