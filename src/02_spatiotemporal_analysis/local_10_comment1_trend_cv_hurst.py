#!/usr/bin/env python
"""
Reviewer comment 1 reanalysis: annual HARSEI trend/CV/Hurst for 2000-2024.

Inputs are the final annual HARSEI rasters generated after the revised HAI and
yearly entropy fusion steps.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


YEARS = np.arange(2000, 2025, dtype=np.int16)
NODATA = -9999.0

TREND_LABELS = {
    1: "Significant degradation",
    2: "Non-significant degradation",
    3: "Stable/no trend",
    4: "Non-significant improvement",
    5: "Significant improvement",
}
CV_LABELS = {
    1: "Low fluctuation (<0.05)",
    2: "Relatively low fluctuation (0.05-0.10)",
    3: "Medium fluctuation (0.10-0.15)",
    4: "High fluctuation (0.15-0.20)",
    5: "Very high fluctuation (>0.20)",
}
HURST_LABELS = {
    1: "Anti-persistent (H < 0.5)",
    2: "Weak persistence (0.5-0.6)",
    3: "Moderate persistence (0.6-0.7)",
    4: "Strong persistence (H > 0.7)",
}
PERSIST_LABELS = {
    1: "Persistent degradation",
    2: "Possible reversal to degradation",
    3: "Weak/no persistence",
    4: "Possible reversal to improvement",
    5: "Persistent improvement",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harsei-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--prefix", default="YJQ")
    parser.add_argument("--chunk-size", type=int, default=50000)
    return parser.parse_args()


def import_rasterio():
    import rasterio

    return rasterio


def read_stack(harsei_dir: Path, prefix: str) -> tuple[np.ndarray, dict]:
    rasterio = import_rasterio()
    arrays = []
    profile = None
    for year in YEARS:
        path = harsei_dir / f"{prefix}_HARSEI_{int(year)}.tif"
        if not path.exists():
            raise FileNotFoundError(path)
        with rasterio.open(path) as src:
            arr = src.read(1).astype(np.float32)
            nodata = src.nodata
            if profile is None:
                profile = src.profile.copy()
            else:
                same = (
                    profile["height"] == src.height
                    and profile["width"] == src.width
                    and profile["crs"] == src.crs
                    and profile["transform"] == src.transform
                )
                if not same:
                    raise ValueError(f"Raster profile differs in {path}")
        if nodata is not None:
            arr[arr == nodata] = np.nan
        arr[arr == NODATA] = np.nan
        arr[~np.isfinite(arr)] = np.nan
        arrays.append(arr)
    assert profile is not None
    return np.stack(arrays, axis=0), profile


def write_single(path: Path, arr: np.ndarray, profile: dict, dtype: str, nodata, description: str) -> None:
    rasterio = import_rasterio()
    out_profile = profile.copy()
    out_profile.update(count=1, dtype=dtype, nodata=nodata, compress="deflate")
    data = arr.copy()
    if np.issubdtype(np.dtype(dtype), np.floating):
        data = data.astype(dtype)
        data[~np.isfinite(data)] = nodata
    else:
        data = data.astype(dtype)
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)


def mk_z_from_s(s: np.ndarray, n: int) -> np.ndarray:
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    z = np.zeros_like(s, dtype=np.float32)
    pos = s > 0
    neg = s < 0
    z[pos] = (s[pos] - 1.0) / math.sqrt(var_s)
    z[neg] = (s[neg] + 1.0) / math.sqrt(var_s)
    return z


def theilsen_mk(stack: np.ndarray, valid: np.ndarray, chunk_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n, h, w = stack.shape
    flat = stack.reshape(n, -1)
    valid_flat = valid.reshape(-1)
    slope_flat = np.full(flat.shape[1], np.nan, dtype=np.float32)
    z_flat = np.full(flat.shape[1], np.nan, dtype=np.float32)
    s_flat = np.full(flat.shape[1], np.nan, dtype=np.float32)
    valid_idx = np.flatnonzero(valid_flat)
    year_vals = YEARS.astype(np.float32)
    pair_count = n * (n - 1) // 2
    for start in range(0, valid_idx.size, chunk_size):
        idx = valid_idx[start : start + chunk_size]
        x = flat[:, idx]
        slopes = np.empty((pair_count, idx.size), dtype=np.float32)
        s = np.zeros(idx.size, dtype=np.float32)
        k = 0
        for i in range(n - 1):
            for j in range(i + 1, n):
                diff = x[j] - x[i]
                slopes[k] = diff / (year_vals[j] - year_vals[i])
                s += np.sign(diff).astype(np.float32)
                k += 1
        slope_flat[idx] = np.median(slopes, axis=0).astype(np.float32)
        s_flat[idx] = s
        z_flat[idx] = mk_z_from_s(s, n)
    return slope_flat.reshape(h, w), z_flat.reshape(h, w), s_flat.reshape(h, w)


def hurst_rs_chunked(stack: np.ndarray, valid: np.ndarray, chunk_size: int) -> np.ndarray:
    n, h, w = stack.shape
    flat = stack.reshape(n, -1).astype(np.float32)
    valid_idx = np.flatnonzero(valid.reshape(-1))
    hurst_flat = np.full(flat.shape[1], np.nan, dtype=np.float32)
    scales = np.array([5, 6, 8, 10, 12, 15, 20, 25], dtype=np.int16)
    log_scales = np.log(scales.astype(np.float64))
    for start in range(0, valid_idx.size, chunk_size):
        idx = valid_idx[start : start + chunk_size]
        x = flat[:, idx].astype(np.float64)
        rs_values = np.full((scales.size, idx.size), np.nan, dtype=np.float64)
        for si, scale in enumerate(scales):
            seg_count = n // int(scale)
            if seg_count < 1:
                continue
            trimmed = x[: seg_count * int(scale), :]
            seg = trimmed.reshape(seg_count, int(scale), idx.size)
            mean = np.mean(seg, axis=1, keepdims=True)
            centered = seg - mean
            cumdev = np.cumsum(centered, axis=1)
            r = np.max(cumdev, axis=1) - np.min(cumdev, axis=1)
            s = np.std(seg, axis=1, ddof=1)
            rs = np.where(s > 1e-12, r / s, np.nan)
            rs_values[si] = np.nanmean(rs, axis=0)
        good = np.isfinite(rs_values) & (rs_values > 0)
        log_rs = np.where(good, np.log(rs_values), np.nan)
        count = np.sum(good, axis=0).astype(np.float64)
        sx = np.nansum(np.where(good, log_scales[:, None], np.nan), axis=0)
        sy = np.nansum(log_rs, axis=0)
        sxx = np.nansum(np.where(good, log_scales[:, None] ** 2, np.nan), axis=0)
        sxy = np.nansum(np.where(good, log_scales[:, None] * log_rs, np.nan), axis=0)
        den = count * sxx - sx * sx
        slope = np.where((count >= 2) & (np.abs(den) > 1e-12), (count * sxy - sx * sy) / den, np.nan)
        hurst_flat[idx] = slope.astype(np.float32)
    return hurst_flat.reshape(h, w)


def classify_trend(slope: np.ndarray, z: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.zeros(slope.shape, dtype=np.uint8)
    sig = np.abs(z) >= 1.96
    eps = 1e-6
    out[valid & (slope < -eps) & sig] = 1
    out[valid & (slope < -eps) & ~sig] = 2
    out[valid & (np.abs(slope) <= eps)] = 3
    out[valid & (slope > eps) & ~sig] = 4
    out[valid & (slope > eps) & sig] = 5
    return out


def classify_cv(cv: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.zeros(cv.shape, dtype=np.uint8)
    out[valid & (cv < 0.05)] = 1
    out[valid & (cv >= 0.05) & (cv < 0.10)] = 2
    out[valid & (cv >= 0.10) & (cv < 0.15)] = 3
    out[valid & (cv >= 0.15) & (cv < 0.20)] = 4
    out[valid & (cv >= 0.20)] = 5
    return out


def classify_hurst(hurst: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.zeros(hurst.shape, dtype=np.uint8)
    out[valid & (hurst < 0.5)] = 1
    out[valid & (hurst >= 0.5) & (hurst < 0.6)] = 2
    out[valid & (hurst >= 0.6) & (hurst < 0.7)] = 3
    out[valid & (hurst >= 0.7)] = 4
    return out


def classify_persistence(slope: np.ndarray, hurst: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.zeros(slope.shape, dtype=np.uint8)
    eps = 1e-6
    weak = np.abs(hurst - 0.5) <= 0.02
    out[valid & weak] = 3
    out[valid & (hurst > 0.52) & (slope < -eps)] = 1
    out[valid & (hurst < 0.48) & (slope > eps)] = 2
    out[valid & (hurst < 0.48) & (slope < -eps)] = 4
    out[valid & (hurst > 0.52) & (slope > eps)] = 5
    return out


def class_summary(arr: np.ndarray, labels: dict[int, str], valid: np.ndarray) -> list[dict]:
    total = int(valid.sum())
    rows = []
    for code, label in labels.items():
        count = int(np.sum(arr[valid] == code))
        rows.append({
            "code": code,
            "class": label,
            "pixels": count,
            "percent": 100.0 * count / total if total else np.nan,
        })
    return rows


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def theilsen_1d(y: np.ndarray, years: np.ndarray) -> float:
    slopes = []
    for i in range(len(y) - 1):
        for j in range(i + 1, len(y)):
            slopes.append((float(y[j]) - float(y[i])) / (float(years[j]) - float(years[i])))
    return float(np.median(slopes))


def mk_1d(y: np.ndarray) -> tuple[float, float, bool]:
    n = len(y)
    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += int(np.sign(y[j] - y[i]))
    z = float(mk_z_from_s(np.array([s], dtype=np.float32), n)[0])
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return z, p, p < 0.05


def save_timeseries_figure(path: Path, annual_rows: list[dict], slope: float, mk_p: float) -> None:
    years = np.array([row["year"] for row in annual_rows], dtype=np.float64)
    means = np.array([row["mean"] for row in annual_rows], dtype=np.float64)
    intercept = float(np.median(means - slope * years))
    trend = intercept + slope * years
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    ax.plot(years, means, color="#377eb8", marker="o", lw=2.2, ms=4.2, label="Annual mean HARSEI")
    ax.plot(years, trend, color="#d95f02", lw=2.0, ls="--", label="Theil-Sen trend")
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean HARSEI")
    ax.set_title("Annual HARSEI trajectory (2000-2024)")
    ax.grid(True, color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False)
    ax.text(
        0.02,
        0.04,
        f"Slope = {slope:.5f} yr$^{{-1}}$; MK p = {mk_p:.3f}",
        transform=ax.transAxes,
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.92, "pad": 5},
    )
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_class_map(path: Path, arr: np.ndarray, labels: dict[int, str], colors: list[str], title: str) -> None:
    cmap = ListedColormap(["#f7f7f7"] + colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 1.5, 1), cmap.N)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    ax.imshow(arr, cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    legend = [Patch(facecolor=colors[code - 1], edgecolor="none", label=label) for code, label in labels.items()]
    ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=8)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_continuous_map(path: Path, arr: np.ndarray, title: str, cmap: str, vmin=None, vmax=None) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.8), constrained_layout=True)
    masked = np.ma.masked_invalid(arr)
    im = ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, shrink=0.78)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    harsei_dir = Path(args.harsei_dir)
    out_dir = Path(args.out_dir)
    raster_dir = out_dir / "rasters"
    table_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"
    for folder in [raster_dir, table_dir, fig_dir]:
        folder.mkdir(parents=True, exist_ok=True)

    stack, profile = read_stack(harsei_dir, args.prefix)
    valid = np.all(np.isfinite(stack), axis=0)
    mean_ts = np.nanmean(stack, axis=(1, 2))
    std_ts = np.nanstd(stack, axis=(1, 2))
    annual_rows = [
        {"year": int(year), "mean": float(mean), "std": float(std)}
        for year, mean, std in zip(YEARS, mean_ts, std_ts)
    ]
    write_csv(table_dir / "annual_harsei_mean_timeseries.csv", annual_rows, ["year", "mean", "std"])

    slope, z, s = theilsen_mk(stack, valid, args.chunk_size)
    cv = (np.nanstd(stack, axis=0) / (np.nanmean(stack, axis=0) + 1e-6)).astype(np.float32)
    cv[~valid] = np.nan
    hurst_raw = hurst_rs_chunked(stack, valid, args.chunk_size)
    hurst_raw[~valid] = np.nan
    hurst = np.clip(hurst_raw, 0.0, 1.0).astype(np.float32)
    hurst[~valid] = np.nan

    trend_class = classify_trend(slope, z, valid)
    cv_class = classify_cv(cv, valid)
    hurst_class = classify_hurst(hurst, valid & np.isfinite(hurst))
    persist_class = classify_persistence(slope, hurst, valid & np.isfinite(hurst))

    write_single(raster_dir / "YJQ_HARSEI_theilsen_slope_2000_2024.tif", slope, profile, "float32", NODATA, "Theil-Sen slope")
    write_single(raster_dir / "YJQ_HARSEI_mk_z_2000_2024.tif", z, profile, "float32", NODATA, "Mann-Kendall Z")
    write_single(raster_dir / "YJQ_HARSEI_cv_2000_2024.tif", cv, profile, "float32", NODATA, "Coefficient of variation")
    write_single(raster_dir / "YJQ_HARSEI_hurst_rs_raw_2000_2024.tif", hurst_raw, profile, "float32", NODATA, "Raw Hurst exponent by R/S")
    write_single(raster_dir / "YJQ_HARSEI_hurst_rs_2000_2024.tif", hurst, profile, "float32", NODATA, "Hurst exponent by R/S clipped to 0-1")
    write_single(raster_dir / "YJQ_HARSEI_trend_class_2000_2024.tif", trend_class, profile, "uint8", 0, "Theil-Sen/MK trend class")
    write_single(raster_dir / "YJQ_HARSEI_cv_class_2000_2024.tif", cv_class, profile, "uint8", 0, "CV class")
    write_single(raster_dir / "YJQ_HARSEI_hurst_class_2000_2024.tif", hurst_class, profile, "uint8", 0, "Hurst class")
    write_single(raster_dir / "YJQ_HARSEI_persistence_class_2000_2024.tif", persist_class, profile, "uint8", 0, "Trend-Hurst persistence class")

    trend_rows = class_summary(trend_class, TREND_LABELS, valid)
    cv_rows = class_summary(cv_class, CV_LABELS, valid)
    hurst_rows = class_summary(hurst_class, HURST_LABELS, valid & np.isfinite(hurst))
    persist_rows = class_summary(persist_class, PERSIST_LABELS, valid & np.isfinite(hurst))
    write_csv(table_dir / "trend_class_summary.csv", trend_rows, ["code", "class", "pixels", "percent"])
    write_csv(table_dir / "cv_class_summary.csv", cv_rows, ["code", "class", "pixels", "percent"])
    write_csv(table_dir / "hurst_class_summary.csv", hurst_rows, ["code", "class", "pixels", "percent"])
    write_csv(table_dir / "persistence_class_summary.csv", persist_rows, ["code", "class", "pixels", "percent"])

    regional_slope = theilsen_1d(mean_ts, YEARS.astype(np.float32))
    regional_z, regional_p, regional_sig = mk_1d(mean_ts)
    summary_rows = [
        {"metric": "n_years", "value": int(len(YEARS))},
        {"metric": "mean_2000", "value": float(mean_ts[0])},
        {"metric": "mean_2024", "value": float(mean_ts[-1])},
        {"metric": "net_change", "value": float(mean_ts[-1] - mean_ts[0])},
        {"metric": "net_change_percent", "value": float((mean_ts[-1] - mean_ts[0]) / mean_ts[0] * 100.0)},
        {"metric": "minimum_year", "value": int(YEARS[int(np.argmin(mean_ts))])},
        {"metric": "minimum_mean", "value": float(np.min(mean_ts))},
        {"metric": "maximum_year", "value": int(YEARS[int(np.argmax(mean_ts))])},
        {"metric": "maximum_mean", "value": float(np.max(mean_ts))},
        {"metric": "regional_theilsen_slope", "value": regional_slope},
        {"metric": "regional_mk_z", "value": regional_z},
        {"metric": "regional_mk_p", "value": regional_p},
        {"metric": "regional_mk_significant_0.05", "value": str(regional_sig)},
        {"metric": "pixel_slope_mean", "value": float(np.nanmean(slope[valid]))},
        {"metric": "pixel_slope_median", "value": float(np.nanmedian(slope[valid]))},
        {"metric": "cv_min", "value": float(np.nanmin(cv[valid]))},
        {"metric": "cv_mean", "value": float(np.nanmean(cv[valid]))},
        {"metric": "cv_median", "value": float(np.nanmedian(cv[valid]))},
        {"metric": "cv_max", "value": float(np.nanmax(cv[valid]))},
        {"metric": "hurst_raw_min", "value": float(np.nanmin(hurst_raw[valid]))},
        {"metric": "hurst_raw_mean", "value": float(np.nanmean(hurst_raw[valid]))},
        {"metric": "hurst_raw_median", "value": float(np.nanmedian(hurst_raw[valid]))},
        {"metric": "hurst_raw_max", "value": float(np.nanmax(hurst_raw[valid]))},
        {"metric": "hurst_clipped_min", "value": float(np.nanmin(hurst[valid]))},
        {"metric": "hurst_clipped_mean", "value": float(np.nanmean(hurst[valid]))},
        {"metric": "hurst_clipped_median", "value": float(np.nanmedian(hurst[valid]))},
        {"metric": "hurst_clipped_max", "value": float(np.nanmax(hurst[valid]))},
    ]
    write_csv(table_dir / "regional_summary_metrics.csv", summary_rows, ["metric", "value"])

    save_timeseries_figure(fig_dir / "fig_comment1_annual_harsei_timeseries.png", annual_rows, regional_slope, regional_p)
    save_class_map(
        fig_dir / "fig_comment1_theilsen_mk_trend_classes.png",
        trend_class,
        TREND_LABELS,
        ["#b2182b", "#ef8a62", "#cfcfcf", "#67a9cf", "#2166ac"],
        "Theil-Sen/Mann-Kendall trend classes (n = 25)",
    )
    save_class_map(
        fig_dir / "fig_comment1_cv_classes.png",
        cv_class,
        CV_LABELS,
        ["#1a9850", "#91cf60", "#fee08b", "#fc8d59", "#d73027"],
        "Temporal fluctuation classes based on CV (n = 25)",
    )
    save_continuous_map(
        fig_dir / "fig_comment1_hurst_continuous.png",
        hurst,
        "Hurst exponent by rescaled range analysis (n = 25; clipped to 0-1)",
        "viridis",
        vmin=0.0,
        vmax=1.0,
    )
    save_class_map(
        fig_dir / "fig_comment1_trend_hurst_persistence.png",
        persist_class,
        PERSIST_LABELS,
        ["#a50026", "#f46d43", "#f7f7f7", "#74add1", "#313695"],
        "Trend-Hurst persistence classes (n = 25)",
    )

    print("Done comment 1 annual trend/CV/Hurst reanalysis.")
    print(f"Outputs: {out_dir}")
    print(f"Regional slope: {regional_slope:.6f}; MK p={regional_p:.4f}; significant={regional_sig}")
    print(f"Mean 2000={mean_ts[0]:.4f}; mean 2024={mean_ts[-1]:.4f}")


if __name__ == "__main__":
    main()
