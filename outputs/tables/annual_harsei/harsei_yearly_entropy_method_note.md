# HARSEI Yearly Entropy Fusion Note

HARSEI was recalculated year by year following the provided R workflow.

For each year, AWRSEI and the HAI-derived layer were stacked, complete-case pixels were extracted, entropy weights were calculated from that year's valid pixels, and the final raster was calculated as:

`HARSEI = AWRSEI * w1 + HAI_layer * w2`

The active HAI layer is `HAI_reverse`. In inverse mode, `HAI_reverse = 1 - HAI`, so higher HARSEI remains ecological-positive.

The annual fusion uses the original AWRSEI and HAI-layer values after weights are obtained, matching `data_valid %*% weights` in the R code.

Mean annual weights over 2000-2024:

| Variable | Mean weight |
| --- | ---: |
| AWRSEI | 0.56641085 |
| HAI_reverse | 0.43358915 |

Detailed annual weights are saved in `harsei_yearly_entropy_weights.csv`.

Backup of replaced HARSEI outputs:

revise\annual_harsei_outputs\backups\before_yearly_entropy_harsei_20260727_163644
