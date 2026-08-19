# Manuscript/Response Text: Pressure-Zone Overlay Validation

## 中文建议表述

为回应审稿人关于“已知盐渍区、灌溉区及矿区区域的 HARSEI 与 RSEI 差异图定量叠加分析”的意见，我们补充构建了三类外部压力区掩膜，并在 2000-2024 年年度栅格上统计 `HARSEI - RSEI` 的区域差异。盐渍区由 ISRIC Global Soil Salinity Maps 的研究区裁剪结果定义，包括任意盐渍等级区（salinity class > 0）和中高盐等级区（salinity class >= 2）；灌溉区由 GFSAD LGRIP30 的 irrigated croplands 类别聚合至 1 km 后，以灌溉农田比例 >= 25% 定义；矿区采用 Global Mining Areas v2 栅格，并重采样到 HARSEI 参考网格。

叠加分析显示，不同压力区的 `HARSEI - RSEI` 差异具有明确的区域分异。盐渍区的差异最明显：2000-2024 年，任意盐渍等级区的 `Median(HARSEI - RSEI)` 相对背景区平均高 0.216，中高盐等级区平均高 0.185；同时，盐渍区内 HARSEI 和 RSEI 均低于背景区，但 RSEI 下降幅度更大。矿区也表现出稳定差异，矿区 `Median(HARSEI - RSEI)` 相对背景区平均高 0.048，且 HARSEI 与 RSEI 均低于非矿区背景。相比之下，灌溉区 `Median(HARSEI - RSEI)` 相对背景区平均为 -0.013，接近零且略为负值；灌溉区内 HARSEI 与 RSEI 均略高于背景区，说明该类区域在本研究区更多代表绿洲农田生态状态，而非单纯负向扰动区。

因此，我们不再将 HARSEI 与 RSEI 的高相关性表述为“验证”，而将其作为模型一致性证据；真正的外部支撑由三部分组成：一是 ISRIC 盐度产品对 SRSI/RSEI/HARSEI 盐度梯度方向的验证，二是玛纳斯河流域实测土壤盐分对 SRSI 盐度响应的独立佐证，三是盐渍区、灌溉区和矿区掩膜与 `HARSEI - RSEI` 差异图的定量叠加分析。

## Suggested English Text

To address the reviewer's request for a quantitative overlay analysis of the HARSEI-RSEI difference map with known saline, irrigated and mining areas, we constructed three independent pressure-zone masks and summarized `HARSEI - RSEI` for annual rasters from 2000 to 2024. Saline areas were derived from the clipped ISRIC Global Soil Salinity Maps, including any saline class (salinity class > 0) and moderate-to-high salinity classes (salinity class >= 2). Irrigated areas were derived from the irrigated-cropland class of GFSAD LGRIP30, aggregated to 1 km and thresholded at an irrigated-cropland fraction of >= 25%. Mining areas were derived from the Global Mining Areas v2 raster and resampled to the HARSEI reference grid.

The overlay analysis revealed distinct patterns among pressure-zone types. The strongest difference occurred in saline areas. Over 2000-2024, the median `HARSEI - RSEI` in any saline class was on average 0.216 higher than in the background area, and the corresponding value for moderate-to-high salinity classes was 0.185. Both HARSEI and RSEI were lower in saline areas than in the background, but RSEI showed a stronger decline. Mining areas also showed a stable difference: the median `HARSEI - RSEI` in mining areas was on average 0.048 higher than in the background, while both HARSEI and RSEI were lower than in non-mining areas. By contrast, irrigated areas showed a near-zero and slightly negative median difference relative to the background (-0.013). HARSEI and RSEI were both slightly higher in irrigated areas than in the background, suggesting that irrigated areas in this study region represent oasis agricultural ecological conditions rather than purely negative disturbance zones.

Accordingly, we no longer present the high correlation between HARSEI and RSEI as validation. Instead, it is treated as model-consistency evidence. External support is now based on three complementary analyses: validation against ISRIC salinity gradients for SRSI/RSEI/HARSEI, field-salinity support for SRSI from the Manas River Basin dataset, and quantitative overlays of the `HARSEI - RSEI` difference maps with saline, irrigated and mining-area masks.

## Key Output Files

- `pressure_zone_mask_manifest.csv`
- `pressure_zone_harsei_minus_rsei_overlay_summary.csv`
- `pressure_zone_overlay_cross_year_summary.csv`
- `fig_pressure_zone_harsei_minus_rsei_overlay.png`
