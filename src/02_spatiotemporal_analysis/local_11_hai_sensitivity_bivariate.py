#!/usr/bin/env python
"""HAI construct sensitivity and AWRSEI-HAI bivariate zoning."""

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


YEARS = list(range(2000, 2025))
NODATA = -9999.0


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--awrsei-dir", required=True)
    p.add_argument("--hai-main-dir", required=True)
    p.add_argument("--hai-input-dir", required=True)
    p.add_argument("--harsei-main-dir", required=True)
    p.add_argument("--tables-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--prefix", default="YJQ")
    return p.parse_args()


def import_rasterio():
    import rasterio

    return rasterio


def read_single(path: Path):
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def read_hai_input(path: Path):
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read([1, 2, 3, 4, 5]).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def write_single(path: Path, arr: np.ndarray, profile: dict, description: str):
    rasterio = import_rasterio()
    out_profile = profile.copy()
    out_profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
    data = arr.astype(np.float32).copy()
    data[~np.isfinite(data)] = NODATA
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)


def read_ranges(path: Path):
    out = {}
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            out[row["component"]] = (
                float(row["raw_min_used_for_norm"]),
                float(row["raw_max_used_for_norm"]),
            )
    return out


def read_light_fit(path: Path):
    a = b = None
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row["variable"] == "DMSP=a+b*log1p(VIIRS):a":
                a = float(row["weight"])
            elif row["variable"] == "DMSP=a+b*log1p(VIIRS):b":
                b = float(row["weight"])
    if a is None or b is None:
        raise ValueError("Night-light harmonization coefficients not found.")
    return a, b


def normalize(x, lo, hi):
    out = (x - lo) / (hi - lo + 1e-6)
    return np.clip(out, 0.0, 1.0).astype(np.float32)


def entropy_weights_user(data_matrix: np.ndarray):
    x = data_matrix.astype(np.float64)
    col_min = np.nanmin(x, axis=0)
    col_max = np.nanmax(x, axis=0)
    norm_data = (x - col_min) / (col_max - col_min + 1e-6)
    p = norm_data / np.sum(norm_data + 1e-12, axis=1)[:, None]
    k = 1.0 / math.log(p.shape[0])
    e = -k * np.sum(p * np.log(p + 1e-12), axis=0)
    d = 1.0 - e
    if abs(float(np.sum(d))) < 1e-12:
        return np.ones(x.shape[1]) / x.shape[1]
    return d / np.sum(d)


def grade(x):
    out = np.zeros(x.shape, dtype=np.uint8)
    out[np.isfinite(x) & (x < 0.2)] = 1
    out[np.isfinite(x) & (x >= 0.2) & (x < 0.4)] = 2
    out[np.isfinite(x) & (x >= 0.4) & (x < 0.6)] = 3
    out[np.isfinite(x) & (x >= 0.6) & (x < 0.8)] = 4
    out[np.isfinite(x) & (x >= 0.8)] = 5
    return out


def save_bivariate(path: Path, aw: np.ndarray, hai: np.ndarray):
    mask = np.isfinite(aw) & np.isfinite(hai)
    aw_thr = float(np.nanmedian(aw[mask]))
    hai_thr = float(np.nanmedian(hai[mask]))
    cls = np.zeros(aw.shape, dtype=np.uint8)
    cls[mask & (aw < aw_thr) & (hai < hai_thr)] = 1
    cls[mask & (aw < aw_thr) & (hai >= hai_thr)] = 2
    cls[mask & (aw >= aw_thr) & (hai < hai_thr)] = 3
    cls[mask & (aw >= aw_thr) & (hai >= hai_thr)] = 4
    labels = {
        1: "Low state / low pressure",
        2: "Low state / high pressure",
        3: "High state / low pressure",
        4: "High state / high pressure",
    }
    colors = ["#fddbc7", "#b2182b", "#2166ac", "#92c5de"]
    cmap = ListedColormap(["#f7f7f7"] + colors)
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    ax.imshow(cls, cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_title("AWRSEI-HAI state-pressure zoning (2024)")
    ax.set_xticks([])
    ax.set_yticks([])
    legend = [Patch(facecolor=colors[k - 1], edgecolor="none", label=v) for k, v in labels.items()]
    ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=8)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    rows = []
    total = int(mask.sum())
    for k, label in labels.items():
        count = int(np.sum(cls == k))
        rows.append({"class": label, "pixels": count, "percent": 100.0 * count / total, "awrsei_median_threshold": aw_thr, "hai_median_threshold": hai_thr})
    return rows


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    aw_dir = Path(args.awrsei_dir)
    hai_main_dir = Path(args.hai_main_dir)
    hai_input_dir = Path(args.hai_input_dir)
    harsei_main_dir = Path(args.harsei_main_dir)
    tables_dir = Path(args.tables_dir)
    out_dir = Path(args.out_dir)
    raster_dir = out_dir / "rasters"
    table_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"
    for d in [raster_dir, table_dir, fig_dir]:
        d.mkdir(parents=True, exist_ok=True)

    ranges = read_ranges(tables_dir / "hai_equal_component_ranges.csv")
    light_fit = read_light_fit(tables_dir / "entropy_weights_and_light_harmonization.csv")
    annual_rows = []
    weight_rows = []
    aw_2024 = hai_2024 = None
    for year in YEARS:
        aw, profile = read_single(aw_dir / f"{args.prefix}_AWRSEI_{year}.tif")
        hai_main, _ = read_single(hai_main_dir / f"{args.prefix}_HAI_{year}.tif")
        harsei_main, _ = read_single(harsei_main_dir / f"{args.prefix}_HARSEI_{year}.tif")
        inp, _ = read_hai_input(hai_input_dir / f"{args.prefix}_hai_inputs_{year}.tif")
        pop = normalize(inp[0], *ranges["POP"])
        if year <= 2013:
            light_raw = inp[2]
        else:
            a, b = light_fit
            light_raw = np.clip(a + b * np.log1p(inp[3]), 0.0, None).astype(np.float32)
        light = normalize(light_raw, *ranges["LIGHT"])
        mask = np.isfinite(aw) & np.isfinite(pop) & np.isfinite(light)
        hai_nolucc = ((pop + light) / 2.0).astype(np.float32)
        hai_nolucc[~mask] = np.nan
        hai_reverse = (1.0 - hai_nolucc).astype(np.float32)
        data = np.column_stack([aw[mask], hai_reverse[mask]])
        weights = entropy_weights_user(data)
        harsei_nolucc = np.full(aw.shape, np.nan, dtype=np.float32)
        harsei_nolucc[mask] = (data @ weights).astype(np.float32)
        write_single(raster_dir / f"{args.prefix}_HAI_noLUCC_{year}.tif", hai_nolucc, profile, "HAI without LUCC")
        write_single(raster_dir / f"{args.prefix}_HARSEI_noLUCC_{year}.tif", harsei_nolucc, profile, "HARSEI sensitivity without LUCC in HAI")
        common = np.isfinite(harsei_main) & np.isfinite(harsei_nolucc)
        corr = float(np.corrcoef(harsei_main[common].ravel(), harsei_nolucc[common].ravel())[0, 1])
        g_main = grade(harsei_main)
        g_alt = grade(harsei_nolucc)
        changed = int(np.sum(common & (g_main != g_alt)))
        annual_rows.append({
            "year": year,
            "main_HAI_mean": float(np.nanmean(hai_main[common])),
            "noLUCC_HAI_mean": float(np.nanmean(hai_nolucc[common])),
            "main_HARSEI_mean": float(np.nanmean(harsei_main[common])),
            "noLUCC_HARSEI_mean": float(np.nanmean(harsei_nolucc[common])),
            "HARSEI_mean_diff_noLUCC_minus_main": float(np.nanmean(harsei_nolucc[common] - harsei_main[common])),
            "HARSEI_pearson_r": corr,
            "grade_changed_pixels": changed,
            "grade_changed_percent": 100.0 * changed / int(common.sum()),
        })
        weight_rows.extend([
            {"year": year, "variable": "AWRSEI", "weight": float(weights[0])},
            {"year": year, "variable": "HAI_reverse_noLUCC", "weight": float(weights[1])},
        ])
        if year == 2024:
            aw_2024 = aw
            hai_2024 = hai_main

    write_csv(table_dir / "hai_no_lucc_sensitivity_annual.csv", annual_rows, list(annual_rows[0].keys()))
    write_csv(table_dir / "hai_no_lucc_sensitivity_weights.csv", weight_rows, ["year", "variable", "weight"])
    summary = {
        "mean_abs_harsei_diff": float(np.mean([abs(r["HARSEI_mean_diff_noLUCC_minus_main"]) for r in annual_rows])),
        "mean_pearson_r": float(np.mean([r["HARSEI_pearson_r"] for r in annual_rows])),
        "max_grade_changed_percent": float(max(r["grade_changed_percent"] for r in annual_rows)),
        "mean_grade_changed_percent": float(np.mean([r["grade_changed_percent"] for r in annual_rows])),
        "mean_main_hai": float(np.mean([r["main_HAI_mean"] for r in annual_rows])),
        "mean_nolucc_hai": float(np.mean([r["noLUCC_HAI_mean"] for r in annual_rows])),
    }
    write_csv(table_dir / "hai_no_lucc_sensitivity_summary.csv", [summary], list(summary.keys()))
    bivar_rows = save_bivariate(fig_dir / "fig_awrsei_hai_bivariate_2024.png", aw_2024, hai_2024)
    write_csv(table_dir / "awrsei_hai_bivariate_2024_summary.csv", bivar_rows, list(bivar_rows[0].keys()))
    print("HAI sensitivity done")
    print(summary)


if __name__ == "__main__":
    main()
