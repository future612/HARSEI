# Revised HAI/HARSEI Rebuild Note

## Revised HAI

HAI was rebuilt as the equal-weight mean of three dimensionless human-activity components:

`HAI = (POP_norm + LIGHT_norm + LUCC_SCORE_norm) / 3`

Each component was normalized using the pooled 2000-2024 valid-pixel min-max range and `1e-6` in the denominator.

## Night-Light Harmonization

`LIGHT` uses DMSP for 2000-2013 and DMSP-equivalent VIIRS for 2014-2024.
The overlap-period fit is `DMSP = -2.74682210 + 18.23736726 * log1p(VIIRS)`.

## HARSEI Fusion

HARSEI was fused from `AWRSEI` and `HAI_reverse`. For the default inverse mode, `HAI_reverse = 1 - HAI`, so higher HARSEI consistently indicates better ecological conditions.
The entropy weights follow the user-provided R function exactly, including row-wise proportions.

| Variable | Weight | Fusion min | Fusion max |
| --- | ---: | ---: | ---: |
| AWRSEI | 0.56819105 | 0.00000000 | 1.00000000 |
| HAI_reverse | 0.43180895 | 0.04596901 | 1.00000000 |

## HAI Component Ranges

| Component | Min | Max |
| --- | ---: | ---: |
| POP | 0.00000000 | 61.02000046 |
| LIGHT | 0.00000000 | 90.68939209 |
| LUCC_SCORE | 0.00000000 | 10.00000000 |

## Backup

revise\annual_harsei_outputs\backups\before_equal_hai_20260727_161245
