RF-SHAP variable importance and model explanation for revised HARSEI

Dependent variable: revised annual HARSEI rasters.
Drivers: DEM, SLOPE, ASPECT, WRB, PRE, TEMMAX, TEMMIN, TEM, SOIL, AET, SWE.
Excluded: HARSEI component variables and HAI/LUCC-related variables to avoid circular interpretation.
Sampling: up to 100,000 valid pixels per year; test size=0.25; SHAP samples=2,000.

Key outputs:
- figures/Fig_HARSEI_RF_SHAP_beeswarm_2000_2024.png/.tif/.pdf
- figures/Fig_HARSEI_RF_SHAP_importance_heatmap_2000_2024.png/.tif/.pdf
- tables/RF_SHAP_model_metrics.csv
- tables/RF_SHAP_variable_importance.csv
- tables/RF_SHAP_beeswarm_plot_data.csv
- tables/RF_SHAP_sample_wide_YYYY.csv
- tables/HARSEI_RF_SHAP_variable_importance_data.xlsx