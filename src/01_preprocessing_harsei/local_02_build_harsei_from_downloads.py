#!/usr/bin/env python
"""
Build annual AWRSEI/HARSEI from downloaded GEE annual input stacks.

Expected ecological component files:
  YJQ_ecocomponents_2000.tif ... YJQ_ecocomponents_2024.tif

Expected first five bands:
  1 NDVI
  2 WET
  3 NDBSI
  4 LST
  5 SRSI

Optional HAI input files:
  YJQ_hai_inputs_2000.tif ... YJQ_hai_inputs_2024.tif

Expected first five bands:
  1 POP
  2 LIGHT_RAW
  3 LIGHT_DMSP
  4 LIGHT_VIIRS
  5 LUCC_SCORE

The script uses:
  - fixed pooled percentile normalization
  - pooled PCA for AWRSEI
  - optional DMSP/VIIRS overlap harmonization
  - entropy weights for AWRSEI and reverse-HAI fusion

Install requirements in a normal Python environment if needed:
  pip install rasterio numpy pandas
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


YEARS = list(range(2000, 2025))
ECO_BANDS = ["NDVI", "WET", "NDBSI", "LST", "SRSI"]
ECO_DIRECTIONS = {
    "NDVI": 1,
    "WET": 1,
    "NDBSI": -1,
    "LST": -1,
    "SRSI": -1,
}
HAI_BANDS = ["POP", "LIGHT_RAW", "LIGHT_DMSP", "LIGHT_VIIRS", "LUCC_SCORE"]
NODATA = -9999.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--eco-dir", required=True, help="Folder with YJQ_ecocomponents_YYYY.tif files.")
    p.add_argument("--hai-dir", default=None, help="Optional folder with YJQ_hai_inputs_YYYY.tif files.")
    p.add_argument("--out-dir", required=True, help="Output folder for AWRSEI/HARSEI rasters and tables.")
    p.add_argument("--prefix", default="YJQ")
    p.add_argument("--range-low", type=float, default=1.0, help="Lower pooled percentile for robust normalization.")
    p.add_argument("--range-high", type=float, default=99.0, help="Upper pooled percentile for robust normalization.")
    p.add_argument("--sample-per-year", type=int, default=50000, help="Samples per year for percentiles/PCA/weights.")
    p.add_argument("--hai-component-weights", choices=["entropy", "equal"], default="entropy")
    p.add_argument("--seed", type=int, default=20260724)
    return p.parse_args()


def import_rasterio():
    try:
        import rasterio
    except Exception as exc:  # pragma: no cover
        raise SystemExit(
            "Missing dependency rasterio. Install with: pip install rasterio numpy pandas"
        ) from exc
    return rasterio


def finite_mask(arr: np.ndarray) -> np.ndarray:
    return np.all(np.isfinite(arr), axis=0)


def sample_rows(arr: np.ndarray, mask: np.ndarray, max_n: int, rng: np.random.Generator) -> np.ndarray:
    # arr is bands x rows x cols. Return n x bands.
    flat = arr.reshape(arr.shape[0], -1).T
    valid = mask.reshape(-1)
    idx = np.flatnonzero(valid)
    if idx.size == 0:
        return np.empty((0, arr.shape[0]), dtype=np.float32)
    if idx.size > max_n:
        idx = rng.choice(idx, size=max_n, replace=False)
    return flat[idx, :].astype(np.float32)


def read_multiband(path: Path, count: int):
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(list(range(1, count + 1))).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def load_eco(eco_dir: Path):
    data = {}
    profile = None
    for year in YEARS:
        path = eco_dir / f"YJQ_ecocomponents_{year}.tif"
        if not path.exists():
            raise FileNotFoundError(path)
        arr, prof = read_multiband(path, 5)
        data[year] = arr
        if profile is None:
            profile = prof
    return data, profile


def load_hai(hai_dir: Path | None):
    if hai_dir is None:
        return None
    data = {}
    for year in YEARS:
        path = hai_dir / f"YJQ_hai_inputs_{year}.tif"
        if not path.exists():
            raise FileNotFoundError(path)
        arr, _ = read_multiband(path, 5)
        data[year] = arr
    return data


def pooled_ranges(raw_by_year: dict[int, np.ndarray], names: list[str], low: float, high: float,
                  sample_per_year: int, rng: np.random.Generator):
    samples = []
    for year, arr in raw_by_year.items():
        mask = finite_mask(arr)
        samples.append(sample_rows(arr, mask, sample_per_year, rng))
    x = np.vstack(samples)
    out = {}
    for i, name in enumerate(names):
        vals = x[:, i]
        vals = vals[np.isfinite(vals)]
        out[name] = (float(np.nanpercentile(vals, low)), float(np.nanpercentile(vals, high)))
    return out


def normalize_band(values: np.ndarray, lo: float, hi: float, direction: int = 1) -> np.ndarray:
    den = hi - lo
    if abs(den) < 1e-12:
        z = np.zeros_like(values, dtype=np.float32)
    else:
        z = (values - lo) / den
    z = np.clip(z, 0.0, 1.0).astype(np.float32)
    if direction < 0:
        z = 1.0 - z
    return z


def normalize_eco(raw_by_year: dict[int, np.ndarray], ranges: dict[str, tuple[float, float]]):
    norm = {}
    for year, arr in raw_by_year.items():
        out = np.empty_like(arr, dtype=np.float32)
        for i, name in enumerate(ECO_BANDS):
            lo, hi = ranges[name]
            out[i] = normalize_band(arr[i], lo, hi, ECO_DIRECTIONS[name])
        out[:, ~finite_mask(arr)] = np.nan
        norm[year] = out
    return norm


def fit_pooled_pca(norm_by_year: dict[int, np.ndarray], sample_per_year: int, rng: np.random.Generator):
    samples = []
    for arr in norm_by_year.values():
        samples.append(sample_rows(arr, finite_mask(arr), sample_per_year, rng))
    x = np.vstack(samples)
    x = x[np.all(np.isfinite(x), axis=1)]
    mean = x.mean(axis=0)
    cov = np.cov(x - mean, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pc1 = eigvecs[:, 0].astype(np.float64)
    if pc1.sum() < 0:
        pc1 *= -1.0
    explained = eigvals / eigvals.sum()
    scores = (x - mean) @ pc1
    score_lo = float(np.nanpercentile(scores, 1.0))
    score_hi = float(np.nanpercentile(scores, 99.0))
    return mean, pc1, eigvals, explained, score_lo, score_hi


def compute_awrsei_for_year(arr: np.ndarray, mean: np.ndarray, pc1: np.ndarray,
                            score_lo: float, score_hi: float) -> np.ndarray:
    mask = finite_mask(arr)
    flat = arr.reshape(arr.shape[0], -1).T
    scores = (flat - mean) @ pc1
    aw = normalize_band(scores.reshape(arr.shape[1], arr.shape[2]), score_lo, score_hi, 1)
    aw[~mask] = np.nan
    return aw.astype(np.float32)


def fit_light_harmonization(hai_raw: dict[int, np.ndarray], sample_per_year: int,
                            rng: np.random.Generator):
    xs = []
    ys = []
    for year in [2012, 2013]:
        if year not in hai_raw:
            continue
        arr = hai_raw[year]
        dmsp = arr[2]
        viirs = arr[3]
        mask = np.isfinite(dmsp) & np.isfinite(viirs) & (viirs > 0)
        idx = np.flatnonzero(mask.reshape(-1))
        if idx.size == 0:
            continue
        if idx.size > sample_per_year:
            idx = rng.choice(idx, size=sample_per_year, replace=False)
        xs.append(np.log1p(viirs.reshape(-1)[idx]))
        ys.append(dmsp.reshape(-1)[idx])
    if not xs:
        return None
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    if x.size < 100:
        return None
    b, a = np.polyfit(x, y, 1)
    return float(a), float(b)


def build_hai_pressure(hai_raw: dict[int, np.ndarray], sample_per_year: int,
                       weight_mode: str, rng: np.random.Generator):
    light_fit = fit_light_harmonization(hai_raw, sample_per_year, rng)
    components = {}
    for year, arr in hai_raw.items():
        pop = arr[0]
        dmsp = arr[2]
        viirs = arr[3]
        lucc = arr[4]
        if year <= 2013:
            light = dmsp
        elif light_fit is not None:
            a, b = light_fit
            light = a + b * np.log1p(viirs)
            light = np.clip(light, 0.0, None)
        else:
            light = np.log1p(arr[1])
        comp = np.stack([pop, light, lucc]).astype(np.float32)
        components[year] = comp

    ranges = pooled_ranges(components, ["POP", "LIGHT", "LUCC_SCORE"], 1.0, 99.0,
                           sample_per_year, rng)
    norm = {}
    for year, comp in components.items():
        out = np.empty_like(comp, dtype=np.float32)
        for i, name in enumerate(["POP", "LIGHT", "LUCC_SCORE"]):
            lo, hi = ranges[name]
            out[i] = normalize_band(comp[i], lo, hi, 1)
        out[:, ~finite_mask(comp)] = np.nan
        norm[year] = out

    samples = []
    for arr in norm.values():
        samples.append(sample_rows(arr, finite_mask(arr), sample_per_year, rng))
    x = np.vstack(samples)
    x = x[np.all(np.isfinite(x), axis=1)]
    if weight_mode == "equal":
        weights = np.array([1 / 3, 1 / 3, 1 / 3], dtype=np.float64)
    else:
        weights = entropy_weights(x)

    pressure = {}
    for year, arr in norm.items():
        p = np.tensordot(weights, arr, axes=(0, 0)).astype(np.float32)
        p[~finite_mask(arr)] = np.nan
        pressure[year] = p
    return pressure, ranges, weights, light_fit


def entropy_weights(x: np.ndarray) -> np.ndarray:
    x = x.astype(np.float64)
    x = np.clip(x, 1e-12, None)
    x = x[np.all(np.isfinite(x), axis=1)]
    col_sum = x.sum(axis=0)
    valid = col_sum > 0
    if not np.all(valid):
        return np.ones(x.shape[1], dtype=np.float64) / x.shape[1]
    p = x / col_sum
    n = x.shape[0]
    e = -np.sum(p * np.log(p), axis=0) / np.log(n)
    d = 1.0 - e
    if d.sum() <= 0:
        return np.ones(x.shape[1], dtype=np.float64) / x.shape[1]
    return d / d.sum()


def write_single(path: Path, arr: np.ndarray, profile: dict):
    rasterio = import_rasterio()
    out_profile = profile.copy()
    out_profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
    data = arr.astype(np.float32).copy()
    data[~np.isfinite(data)] = NODATA
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(data, 1)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    eco_dir = Path(args.eco_dir)
    hai_dir = Path(args.hai_dir) if args.hai_dir else None
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raster_dir = out_dir / "rasters"
    table_dir = out_dir / "tables"
    raster_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)

    eco_raw, profile = load_eco(eco_dir)
    eco_ranges = pooled_ranges(eco_raw, ECO_BANDS, args.range_low, args.range_high,
                               args.sample_per_year, rng)
    eco_norm = normalize_eco(eco_raw, eco_ranges)
    mean, pc1, eigvals, explained, score_lo, score_hi = fit_pooled_pca(
        eco_norm, args.sample_per_year, rng
    )

    awrsei = {}
    for year, arr in eco_norm.items():
        aw = compute_awrsei_for_year(arr, mean, pc1, score_lo, score_hi)
        awrsei[year] = aw
        write_single(raster_dir / f"{args.prefix}_AWRSEI_{year}.tif", aw, profile)

    hai_raw = load_hai(hai_dir)
    hai_reverse = None
    harsei = None
    fusion_weights = None
    hai_ranges = None
    hai_component_weights = None
    light_fit = None
    if hai_raw is not None:
        hai_pressure, hai_ranges, hai_component_weights, light_fit = build_hai_pressure(
            hai_raw, args.sample_per_year, args.hai_component_weights, rng
        )
        hai_reverse = {year: (1.0 - arr).astype(np.float32) for year, arr in hai_pressure.items()}

        fusion_samples = []
        for year in YEARS:
            a = awrsei[year]
            h = hai_reverse[year]
            mask = np.isfinite(a) & np.isfinite(h)
            pair = np.stack([a, h])
            fusion_samples.append(sample_rows(pair, mask, args.sample_per_year, rng))
        fx = np.vstack(fusion_samples)
        fusion_weights = entropy_weights(fx)
        harsei = {}
        for year in YEARS:
            out = fusion_weights[0] * awrsei[year] + fusion_weights[1] * hai_reverse[year]
            out[~np.isfinite(out)] = np.nan
            harsei[year] = out.astype(np.float32)
            write_single(raster_dir / f"{args.prefix}_HAI_reverse_{year}.tif", hai_reverse[year], profile)
            write_single(raster_dir / f"{args.prefix}_HARSEI_{year}.tif", harsei[year], profile)

    range_rows = []
    for name in ECO_BANDS:
        lo, hi = eco_ranges[name]
        range_rows.append({"variable": name, "p_low": lo, "p_high": hi, "direction": ECO_DIRECTIONS[name]})
    if hai_ranges:
        for name, (lo, hi) in hai_ranges.items():
            range_rows.append({"variable": name, "p_low": lo, "p_high": hi, "direction": 1})
    write_csv(table_dir / "fixed_normalization_ranges.csv", range_rows,
              ["variable", "p_low", "p_high", "direction"])

    pca_rows = []
    for i, name in enumerate(ECO_BANDS):
        pca_rows.append({
            "variable": name,
            "pc1_loading": pc1[i],
            "pooled_mean_after_direction_norm": mean[i],
            "eigenvalue_pc1": eigvals[0],
            "pc1_explained_variance_ratio": explained[0],
            "pc1_score_p1": score_lo,
            "pc1_score_p99": score_hi,
        })
    write_csv(table_dir / "pooled_pca_pc1.csv", pca_rows,
              ["variable", "pc1_loading", "pooled_mean_after_direction_norm",
               "eigenvalue_pc1", "pc1_explained_variance_ratio", "pc1_score_p1", "pc1_score_p99"])

    weight_rows = []
    if hai_component_weights is not None:
        for name, weight in zip(["POP", "LIGHT", "LUCC_SCORE"], hai_component_weights):
            weight_rows.append({"level": "HAI_components", "variable": name, "weight": weight})
    if fusion_weights is not None:
        weight_rows.append({"level": "HARSEI_fusion", "variable": "AWRSEI", "weight": fusion_weights[0]})
        weight_rows.append({"level": "HARSEI_fusion", "variable": "HAI_reverse", "weight": fusion_weights[1]})
    if light_fit is not None:
        weight_rows.append({"level": "night_light_harmonization", "variable": "DMSP=a+b*log1p(VIIRS):a", "weight": light_fit[0]})
        weight_rows.append({"level": "night_light_harmonization", "variable": "DMSP=a+b*log1p(VIIRS):b", "weight": light_fit[1]})
    if weight_rows:
        write_csv(table_dir / "entropy_weights_and_light_harmonization.csv", weight_rows,
                  ["level", "variable", "weight"])

    stat_rows = []
    for year in YEARS:
        row = {
            "year": year,
            "AWRSEI_mean": float(np.nanmean(awrsei[year])),
            "AWRSEI_std": float(np.nanstd(awrsei[year])),
        }
        if harsei is not None:
            row.update({
                "HAI_reverse_mean": float(np.nanmean(hai_reverse[year])),
                "HARSEI_mean": float(np.nanmean(harsei[year])),
                "HARSEI_std": float(np.nanstd(harsei[year])),
            })
        stat_rows.append(row)
    fieldnames = list(stat_rows[0].keys())
    write_csv(table_dir / "annual_index_summary.csv", stat_rows, fieldnames)

    print("Done.")
    print("Rasters:", raster_dir)
    print("Tables:", table_dir)


if __name__ == "__main__":
    main()
