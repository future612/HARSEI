GeoDetector interaction detection for revised annual HARSEI

Dependent variable: revised HARSEI rasters from D:\Codex\260724 小论文\revise\annual_harsei_outputs\rasters\HARSEI.
Driver variables: DEM, SLOPE, ASPECT, WRB, PRE, TEMMAX, TEMMIN, TEM, SOIL, AET, and SWE from E:\wl24\hys to wl\数据清单\1数据下载.
Excluded from the main driver-interaction figure: HARSEI component variables and HAI/LUCC-related variables, to avoid circularity in driver interpretation.
Continuous drivers were discretized into 7 quantile classes for each year; WRB was retained as a categorical factor.
Interaction type follows standard GeoDetector rules based on q(X1), q(X2), q(X1∩X2), and q(X1)+q(X2).

Key outputs:
- figures/Fig_HARSEI_GeoDetector_interaction_2000_2024.png/.tif/.pdf
- tables/GD_factor_revised_HARSEI_YYYY.csv
- tables/GD_interaction_q_revised_HARSEI_YYYY.csv
- tables/GD_interaction_type_revised_HARSEI_YYYY.csv
- tables/GD_interaction_long_revised_HARSEI_2000_2024.csv
- tables/HARSEI_GeoDetector_interaction_data.xlsx