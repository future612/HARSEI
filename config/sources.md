# Sources To Cite Or Check

Primary data/catalog pages:

- MOD13A1 NDVI: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13A1
- MOD09A1 surface reflectance: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09A1
- MOD11A2 land surface temperature: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A2
- WorldPop population: https://developers.google.com/earth-engine/datasets/catalog/WorldPop_GP_100m_pop
- DMSP-OLS night-time lights: https://developers.google.com/earth-engine/datasets/catalog/NOAA_DMSP-OLS_NIGHTTIME_LIGHTS
- VIIRS monthly night-time lights: https://developers.google.com/earth-engine/datasets/catalog/NOAA_VIIRS_DNB_MONTHLY_V1_VCMCFG
- MODIS MCD12Q1 annual land cover: https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1
- SRTM elevation: https://developers.google.com/earth-engine/datasets/catalog/USGS_SRTMGL1_003
- CLCD data paper: https://essd.copernicus.org/articles/13/3907/2021/
- CLCD latest Zenodo record: https://zenodo.org/records/4417810/latest
- GLC-FCS30D SAT-IO GEE asset: `projects/sat-io/open-datasets/GLC-FCS30D/annual`
- GLC_FCS30D data paper: https://essd.copernicus.org/articles/16/1353/2024/
- GLC_FCS30D Zenodo record: https://zenodo.org/records/8239305

Methods:

- Lobser, S.E., Cohen, W.B., 2007. MODIS tasseled cap: land cover characteristics expressed through transformed MODIS data. International Journal of Remote Sensing 28, 5079-5101. DOI: https://doi.org/10.1080/01431160701253303
- Xu, H., 2013. Remote sensing ecological index / RSEI reference for the general RSEI framework.
- Alhammadi, M.S., Glenn, E.P., 2008. Detecting date palm trees health and vegetation greenness change on the eastern coast of the United Arab Emirates using SAVI. International Journal of Remote Sensing 29, 1745-1765. DOI: https://doi.org/10.1080/01431160701281089
- Khan, N.M., Rastoskuev, V.V., Shalina, E.V., Sato, Y., 2001. Mapping salt-affected soils using remote sensing indicators. Check as the original NDSI / salinity-index family source before final citation.
- Yang, J., Huang, X., 2021. The 30 m annual land cover dataset and its dynamics in China from 1990 to 2019. Earth System Science Data 13, 3907-3925. DOI: https://doi.org/10.5194/essd-13-3907-2021
- Zhang, X., Zhao, T., Xu, H., Liu, W., Wang, J., Chen, X., Liu, L., 2024. GLC_FCS30D: the first global 30 m land-cover dynamics monitoring product with a fine classification system for the period from 1985 to 2022 generated using dense-time-series Landsat imagery and the continuous change-detection method. Earth System Science Data 16, 1353-1381. DOI: https://doi.org/10.5194/essd-16-1353-2024
- Mann, H.B., 1945; Sen, P.K., 1968; Kendall, M.G., 1975 for trend tests.

Reviewer-facing notes:

- Cite MOD11A2 as product LST with scale factor, not a Landsat radiance/K1/K2 retrieval.
- Cite MODIS tasseled-cap wetness as sensor-specific; do not cite Landsat wetness coefficients for MOD09A1.
- Report that VIIRS annual composites are weighted by `cf_cvg`; then report the local DMSP = a + b * log1p(VIIRS) overlap fit using 2012-2013.
- State WorldPop 2024 is not available in the current downloaded stack; carry 2020 forward for 2021-2024 with disclosure and, if time allows, add a POP sensitivity test.
- State whether LUCC is GLC-FCS30D, user-uploaded CLCD, or MCD12Q1. The user's pasted code uses GLC-FCS30D; do not describe that workflow as CLCD.
- State that HARSEI annual reconstruction uses fixed-reference normalization and pooled PCA.
