#!/usr/bin/env python
"""
External validation against ISRIC / Ivushkin et al. Global Soil Salinity Maps.

Required inputs:
  - ISRIC salinity GeoTIFFs exported from GEE:
      ISRIC_global_soil_salinity_2000_YJQ.tif
      ISRIC_global_soil_salinity_2002_YJQ.tif
      ISRIC_global_soil_salinity_2005_YJQ.tif
      ISRIC_global_soil_salinity_2009_YJQ.tif
      ISRIC_global_soil_salinity_2016_YJQ.tif
  - Annual ecological component stacks:
      YJQ_ecocomponents_YYYY.tif
      band 1 NDVI, band 2 WET, band 3 NDBSI, band 4 LST, band 5 SRSI
  - Annual HARSEI rasters:
      YJQ_HARSEI_YYYY.tif

Outputs:
  - Spearman correlations between salinity and SRSI/RSEI/HARSEI.
  - Means/medians of SRSI, RSEI and HARSEI by external salinity class.
  - HARSEI - RSEI difference in saline vs. non-saline areas.
  - Gradient test: whether HARSEI declines more clearly than RSEI with salinity.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np


YEARS_DEFAULT = [2000, 2002, 2005, 2009, 2016]
ECO_BANDS = ["NDVI", "WET", "NDBSI", "LST", "SRSI"]
RSEI_BANDS = ["NDVI", "WET", "NDBSI", "LST"]
NODATA = -9999.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--salinity-dir", default=r"D:\Codex\260724 小论文\revise\external_validation_inputs\isric_global_soil_salinity")
    p.add_argument("--eco-dir", default=r"D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024")
    p.add_argument("--harsei-dir", default=r"D:\Codex\260724 小论文\revise\annual_harsei_outputs\rasters\HARSEI")
    p.add_argument("--fixed-ranges", default=r"D:\Codex\260724 小论文\revise\annual_harsei_outputs\tables\fixed_normalization_ranges.csv")
    p.add_argument("--pca-sample", default=r"D:\Codex\260724 小论文\revise\gee_downloads\YJQ_HARSEI_annual_inputs_2000_2024\YJQ_pooled_pca_sample_2000_2024.csv")
    p.add_argument("--out-dir", default=r"D:\Codex\260724 小论文\revise\isric_salinity_validation")
    p.add_argument("--years", default=",".join(str(y) for y in YEARS_DEFAULT))
    p.add_argument("--salinity-resampling", choices=["auto", "nearest", "mode", "average", "bilinear"], default="auto")
    p.add_argument("--max-unique-classes", type=int, default=15)
    p.add_argument("--continuous-class-count", type=int, default=5)
    p.add_argument("--saline-threshold", default="auto", help="auto, or numeric threshold. Auto uses the lowest ordinal class as non-saline.")
    return p.parse_args()


def import_rasterio():
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.merge import merge
    from rasterio.warp import reproject

    return rasterio, Resampling, merge, reproject


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def find_year_tifs(folder: Path, year: int) -> list[Path]:
    patterns = [
        f"*{year}*.tif",
        f"*{year}*.tiff",
        f"*{year}*.TIF",
        f"*{year}*.TIFF",
    ]
    files: list[Path] = []
    for pat in patterns:
        files.extend(folder.glob(pat))
    files = sorted(set(files))
    return [p for p in files if p.is_file()]


def read_single(path: Path) -> tuple[np.ndarray, dict]:
    rasterio, _, _, _ = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def read_multiband(path: Path, count: int) -> tuple[np.ndarray, dict]:
    rasterio, _, _, _ = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(list(range(1, count + 1))).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def read_salinity(files: list[Path]) -> tuple[np.ndarray, dict]:
    rasterio, _, merge, _ = import_rasterio()
    if len(files) == 1:
        return read_single(files[0])
    datasets = [rasterio.open(str(p)) for p in files]
    try:
        merged, transform = merge(datasets)
        profile = datasets[0].profile.copy()
        profile.update(height=merged.shape[1], width=merged.shape[2], transform=transform, count=1)
        arr = merged[0].astype(np.float32)
        nodata = profile.get("nodata")
        if nodata is not None:
            arr[arr == nodata] = np.nan
        arr[arr == NODATA] = np.nan
        arr[~np.isfinite(arr)] = np.nan
        return arr, profile
    finally:
        for ds in datasets:
            ds.close()


def profile_matches(a: dict, b: dict) -> bool:
    return all(a.get(key) == b.get(key) for key in ["height", "width", "crs", "transform"])


def is_ordinal_class(arr: np.ndarray, max_unique: int = 15) -> tuple[bool, np.ndarray]:
    vals = arr[np.isfinite(arr)]
    if vals.size == 0:
        return False, np.array([])
    sample = vals if vals.size <= 300000 else vals[np.linspace(0, vals.size - 1, 300000).astype(int)]
    uniques = np.unique(sample)
    uniques = uniques[np.isfinite(uniques)]
    if uniques.size <= max_unique and np.all(np.isclose(uniques, np.round(uniques), atol=1e-6)):
        return True, np.sort(uniques)
    return False, np.array([])


def resampling_enum(name: str, categorical: bool):
    _, Resampling, _, _ = import_rasterio()
    if name == "auto":
        if categorical and hasattr(Resampling, "mode"):
            return Resampling.mode
        return Resampling.nearest if categorical else Resampling.average
    return {
        "nearest": Resampling.nearest,
        "mode": Resampling.mode if hasattr(Resampling, "mode") else Resampling.nearest,
        "average": Resampling.average,
        "bilinear": Resampling.bilinear,
    }[name]


def reproject_to_reference(arr: np.ndarray, src_profile: dict, ref_profile: dict, resampling) -> np.ndarray:
    if profile_matches(src_profile, ref_profile):
        return arr.astype(np.float32)
    _, _, _, reproject = import_rasterio()
    dst = np.full((ref_profile["height"], ref_profile["width"]), np.nan, dtype=np.float32)
    src_nodata = src_profile.get("nodata")
    if src_nodata is None or not np.isfinite(src_nodata):
        src_nodata = NODATA
    src_data = arr.astype(np.float32).copy()
    src_data[~np.isfinite(src_data)] = src_nodata
    reproject(
        source=src_data,
        destination=dst,
        src_transform=src_profile["transform"],
        src_crs=src_profile["crs"],
        src_nodata=src_nodata,
        dst_transform=ref_profile["transform"],
        dst_crs=ref_profile["crs"],
        dst_nodata=np.nan,
        resampling=resampling,
    )
    dst[~np.isfinite(dst)] = np.nan
    return dst.astype(np.float32)


def normalize_band(values: np.ndarray, lo: float, hi: float, direction: int) -> np.ndarray:
    den = hi - lo
    if abs(den) < 1e-12:
        z = np.zeros_like(values, dtype=np.float32)
    else:
        z = (values - lo) / den
    z = np.clip(z, 0.0, 1.0).astype(np.float32)
    if direction < 0:
        z = 1.0 - z
    return z


def load_ranges(path: Path) -> dict[str, tuple[float, float, int]]:
    rows = read_csv_rows(path)
    out = {}
    for row in rows:
        var = row["variable"]
        if var in ECO_BANDS:
            out[var] = (float(row["p_low"]), float(row["p_high"]), int(float(row["direction"])))
    return out


def pca_from_matrix(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = x[np.all(np.isfinite(x), axis=1)]
    mean = x.mean(axis=0)
    xc = x - mean
    cov = np.cov(xc, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    explained = eigvals / eigvals.sum()
    pc1 = eigvecs[:, 0].astype(np.float64)
    if pc1.sum() < 0:
        pc1 *= -1.0
    return mean.astype(np.float64), pc1, eigvals, explained


def fit_rsei_pca(sample_csv: Path, ranges: dict[str, tuple[float, float, int]]) -> dict:
    rows = []
    with sample_csv.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            vals = []
            ok = True
            for var in RSEI_BANDS:
                try:
                    raw = float(row[var])
                except Exception:
                    ok = False
                    break
                lo, hi, direction = ranges[var]
                vals.append(float(normalize_band(np.array([raw], dtype=np.float32), lo, hi, direction)[0]))
            if ok:
                rows.append(vals)
    x = np.asarray(rows, dtype=np.float64)
    mean, pc1, eigvals, explained = pca_from_matrix(x)
    scores = (x - mean) @ pc1
    score_p1 = float(np.nanpercentile(scores, 1))
    score_p99 = float(np.nanpercentile(scores, 99))
    return {
        "mean": mean,
        "pc1": pc1,
        "eigenvalue_pc1": float(eigvals[0]),
        "explained_pc1": float(explained[0]),
        "score_p1": score_p1,
        "score_p99": score_p99,
        "sample_size": int(x.shape[0]),
    }


def compute_rsei_from_ecostack(eco: np.ndarray, ranges: dict[str, tuple[float, float, int]], pca: dict) -> np.ndarray:
    norm = []
    for band_i, var in enumerate(RSEI_BANDS):
        lo, hi, direction = ranges[var]
        norm.append(normalize_band(eco[band_i], lo, hi, direction))
    x = np.stack(norm)
    mask = np.all(np.isfinite(x), axis=0)
    flat = x.reshape(x.shape[0], -1).T
    scores = (flat - pca["mean"]) @ pca["pc1"]
    rsei = normalize_band(scores.reshape(x.shape[1], x.shape[2]), pca["score_p1"], pca["score_p99"], 1)
    rsei[~mask] = np.nan
    return rsei.astype(np.float32)


def rankdata_average(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    sorted_vals = values[order]
    i = 0
    while i < values.size:
        j = i + 1
        while j < values.size and sorted_vals[j] == sorted_vals[i]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return float("nan")
    x = x - x.mean()
    y = y - y.mean()
    den = math.sqrt(float(np.sum(x * x) * np.sum(y * y)))
    if den == 0:
        return float("nan")
    return float(np.sum(x * y) / den)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[mask], dtype=np.float64)
    y = np.asarray(y[mask], dtype=np.float64)
    if x.size < 3:
        return float("nan")
    return pearson(rankdata_average(x), rankdata_average(y))


def median(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.nanmedian(values)) if values.size else float("nan")


def mean(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.nanmean(values)) if values.size else float("nan")


def std(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.nanstd(values)) if values.size else float("nan")


def class_labels(sal: np.ndarray, categorical: bool, unique_classes: np.ndarray, class_count: int):
    if categorical:
        classes = np.full(sal.shape, -9999, dtype=np.int16)
        label_map = {}
        for rank, val in enumerate(unique_classes, start=1):
            classes[np.isclose(sal, val, atol=1e-6)] = rank
            label_map[rank] = f"class_{int(val) if float(val).is_integer() else val:g}"
        return classes, label_map
    vals = sal[np.isfinite(sal)]
    qs = np.unique(np.nanquantile(vals, np.linspace(0, 1, class_count + 1)))
    classes = np.full(sal.shape, -9999, dtype=np.int16)
    label_map = {}
    for i in range(len(qs) - 1):
        lo, hi = qs[i], qs[i + 1]
        if i == len(qs) - 2:
            m = (sal >= lo) & (sal <= hi)
        else:
            m = (sal >= lo) & (sal < hi)
        classes[m] = i + 1
        label_map[i + 1] = f"Q{i+1}_{lo:.3g}_{hi:.3g}"
    return classes, label_map


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask].astype(np.float64)
    y = y[mask].astype(np.float64)
    if x.size < 2:
        return float("nan")
    x0 = x - x.mean()
    den = np.sum(x0 * x0)
    if den == 0:
        return float("nan")
    return float(np.sum(x0 * (y - y.mean())) / den)


def monotonic_decrease_steps(y: Iterable[float]) -> int:
    vals = [v for v in y if np.isfinite(v)]
    return int(sum(vals[i + 1] <= vals[i] for i in range(len(vals) - 1)))


def try_make_figures(class_stats_path: Path, out_dir: Path):
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except Exception:
        return
    df = pd.read_csv(class_stats_path)
    if df.empty:
        return
    pooled = df.groupby("class_rank")[["SRSI_median", "RSEI_median", "HARSEI_median"]].mean().reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4), dpi=300)
    for ax, col, color, ylabel in [
        (axes[0], "SRSI_median", "#d95f02", "SRSI median"),
        (axes[1], "RSEI_median", "#1b9e77", "RSEI median"),
        (axes[2], "HARSEI_median", "#7570b3", "HARSEI median"),
    ]:
        ax.plot(pooled["class_rank"], pooled[col], marker="o", color=color, lw=1.8)
        ax.set_xlabel("External salinity class rank")
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#dddddd", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_isric_salinity_gradient_medians.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    years = [int(y.strip()) for y in args.years.split(",") if y.strip()]
    sal_dir = Path(args.salinity_dir)
    eco_dir = Path(args.eco_dir)
    harsei_dir = Path(args.harsei_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    missing = []
    for year in years:
        if not find_year_tifs(sal_dir, year):
            missing.append(year)
    if missing:
        raise SystemExit(
            "Missing ISRIC salinity GeoTIFF(s) for years "
            + ", ".join(str(y) for y in missing)
            + f". Put files in: {sal_dir}"
        )

    ranges = load_ranges(Path(args.fixed_ranges))
    rsei_pca = fit_rsei_pca(Path(args.pca_sample), ranges)
    pca_rows = []
    for var, loading, mean_val in zip(RSEI_BANDS, rsei_pca["pc1"], rsei_pca["mean"]):
        pca_rows.append({
            "variable": var,
            "pc1_loading": loading,
            "mean_after_direction_norm": mean_val,
            "eigenvalue_pc1": rsei_pca["eigenvalue_pc1"],
            "pc1_explained_variance_ratio": rsei_pca["explained_pc1"],
            "score_p1": rsei_pca["score_p1"],
            "score_p99": rsei_pca["score_p99"],
            "sample_size": rsei_pca["sample_size"],
        })
    write_csv(out_dir / "rsei_pooled_pca_pc1.csv", pca_rows, list(pca_rows[0].keys()))

    correlation_rows: list[dict] = []
    class_rows: list[dict] = []
    diff_rows: list[dict] = []
    gradient_rows: list[dict] = []
    pooled_arrays: dict[str, list[np.ndarray]] = {"salinity": [], "SRSI": [], "RSEI": [], "HARSEI": [], "DIFF": []}

    for year in years:
        eco_path = eco_dir / f"YJQ_ecocomponents_{year}.tif"
        harsei_path = harsei_dir / f"YJQ_HARSEI_{year}.tif"
        if not eco_path.exists():
            raise FileNotFoundError(eco_path)
        if not harsei_path.exists():
            raise FileNotFoundError(harsei_path)

        eco, ref_profile = read_multiband(eco_path, 5)
        harsei, harsei_profile = read_single(harsei_path)
        if not profile_matches(ref_profile, harsei_profile):
            resampling = resampling_enum("bilinear", False)
            harsei = reproject_to_reference(harsei, harsei_profile, ref_profile, resampling)

        sal_files = find_year_tifs(sal_dir, year)
        sal_raw, sal_profile = read_salinity(sal_files)
        categorical, unique_vals = is_ordinal_class(sal_raw, args.max_unique_classes)
        sal = reproject_to_reference(
            sal_raw,
            sal_profile,
            ref_profile,
            resampling_enum(args.salinity_resampling, categorical),
        )
        if categorical:
            sal = np.rint(sal).astype(np.float32)

        srsi = eco[4].astype(np.float32)
        rsei = compute_rsei_from_ecostack(eco, ranges, rsei_pca)
        diff = harsei - rsei
        valid = np.isfinite(sal) & np.isfinite(srsi) & np.isfinite(rsei) & np.isfinite(harsei)

        sal_v = sal[valid]
        srsi_v = srsi[valid]
        rsei_v = rsei[valid]
        harsei_v = harsei[valid]
        diff_v = diff[valid]

        pooled_arrays["salinity"].append(sal_v)
        pooled_arrays["SRSI"].append(srsi_v)
        pooled_arrays["RSEI"].append(rsei_v)
        pooled_arrays["HARSEI"].append(harsei_v)
        pooled_arrays["DIFF"].append(diff_v)

        correlation_rows.extend([
            {"year": year, "variable": "SRSI", "n": int(sal_v.size), "spearman_rho_with_salinity": spearman(srsi_v, sal_v), "expected_direction": "positive"},
            {"year": year, "variable": "RSEI", "n": int(sal_v.size), "spearman_rho_with_salinity": spearman(rsei_v, sal_v), "expected_direction": "negative"},
            {"year": year, "variable": "HARSEI", "n": int(sal_v.size), "spearman_rho_with_salinity": spearman(harsei_v, sal_v), "expected_direction": "negative"},
            {"year": year, "variable": "HARSEI_minus_RSEI", "n": int(sal_v.size), "spearman_rho_with_salinity": spearman(diff_v, sal_v), "expected_direction": "diagnostic; sign depends on HAI adjustment"},
        ])

        sal_classes, labels = class_labels(sal, categorical, unique_vals, args.continuous_class_count)
        valid_class = valid & (sal_classes > 0)
        med_harsei_by_class = []
        med_rsei_by_class = []
        class_rank_values = []
        for rank in sorted(labels):
            m = valid_class & (sal_classes == rank)
            if not np.any(m):
                continue
            class_rank_values.append(rank)
            med_harsei_by_class.append(median(harsei[m]))
            med_rsei_by_class.append(median(rsei[m]))
            class_rows.append({
                "year": year,
                "class_rank": rank,
                "external_salinity_class": labels[rank],
                "n": int(m.sum()),
                "salinity_mean": mean(sal[m]),
                "salinity_median": median(sal[m]),
                "SRSI_mean": mean(srsi[m]),
                "SRSI_median": median(srsi[m]),
                "RSEI_mean": mean(rsei[m]),
                "RSEI_median": median(rsei[m]),
                "HARSEI_mean": mean(harsei[m]),
                "HARSEI_median": median(harsei[m]),
                "HARSEI_minus_RSEI_mean": mean(diff[m]),
                "HARSEI_minus_RSEI_median": median(diff[m]),
            })

        ranks_arr = np.asarray(class_rank_values, dtype=np.float64)
        harsei_med_arr = np.asarray(med_harsei_by_class, dtype=np.float64)
        rsei_med_arr = np.asarray(med_rsei_by_class, dtype=np.float64)
        h_slope = linear_slope(ranks_arr, harsei_med_arr)
        r_slope = linear_slope(ranks_arr, rsei_med_arr)
        h_spearman = spearman(harsei_med_arr, ranks_arr)
        r_spearman = spearman(rsei_med_arr, ranks_arr)
        gradient_rows.append({
            "year": year,
            "class_count": int(len(class_rank_values)),
            "HARSEI_median_slope_per_class": h_slope,
            "RSEI_median_slope_per_class": r_slope,
            "HARSEI_spearman_across_class_medians": h_spearman,
            "RSEI_spearman_across_class_medians": r_spearman,
            "HARSEI_monotonic_decrease_steps": monotonic_decrease_steps(harsei_med_arr),
            "RSEI_monotonic_decrease_steps": monotonic_decrease_steps(rsei_med_arr),
            "HARSEI_clearer_decline_than_RSEI": (
                "Yes" if np.isfinite(h_slope) and np.isfinite(r_slope)
                and h_slope < 0 and abs(h_slope) > abs(r_slope) else "No"
            ),
        })

        if args.saline_threshold == "auto":
            if categorical and unique_vals.size > 0:
                threshold = float(np.min(unique_vals))
                saline_mask = valid & (sal > threshold)
                nonsaline_mask = valid & np.isclose(sal, threshold, atol=1e-6)
                threshold_note = f"auto categorical: saline > lowest class ({threshold:g})"
            else:
                threshold = float(np.nanpercentile(sal_v, 75))
                saline_mask = valid & (sal >= threshold)
                nonsaline_mask = valid & (sal < threshold)
                threshold_note = f"auto continuous: saline >= P75 ({threshold:.6g})"
        else:
            threshold = float(args.saline_threshold)
            saline_mask = valid & (sal > threshold)
            nonsaline_mask = valid & (sal <= threshold)
            threshold_note = f"user threshold: saline > {threshold:g}"

        for group_name, m in [("non_saline", nonsaline_mask), ("saline", saline_mask)]:
            diff_rows.append({
                "year": year,
                "group": group_name,
                "threshold_rule": threshold_note,
                "n": int(m.sum()),
                "SRSI_mean": mean(srsi[m]),
                "SRSI_median": median(srsi[m]),
                "RSEI_mean": mean(rsei[m]),
                "RSEI_median": median(rsei[m]),
                "HARSEI_mean": mean(harsei[m]),
                "HARSEI_median": median(harsei[m]),
                "HARSEI_minus_RSEI_mean": mean(diff[m]),
                "HARSEI_minus_RSEI_median": median(diff[m]),
            })
        if np.any(saline_mask) and np.any(nonsaline_mask):
            diff_rows.append({
                "year": year,
                "group": "saline_minus_non_saline",
                "threshold_rule": threshold_note,
                "n": int(saline_mask.sum() + nonsaline_mask.sum()),
                "SRSI_mean": mean(srsi[saline_mask]) - mean(srsi[nonsaline_mask]),
                "SRSI_median": median(srsi[saline_mask]) - median(srsi[nonsaline_mask]),
                "RSEI_mean": mean(rsei[saline_mask]) - mean(rsei[nonsaline_mask]),
                "RSEI_median": median(rsei[saline_mask]) - median(rsei[nonsaline_mask]),
                "HARSEI_mean": mean(harsei[saline_mask]) - mean(harsei[nonsaline_mask]),
                "HARSEI_median": median(harsei[saline_mask]) - median(harsei[nonsaline_mask]),
                "HARSEI_minus_RSEI_mean": mean(diff[saline_mask]) - mean(diff[nonsaline_mask]),
                "HARSEI_minus_RSEI_median": median(diff[saline_mask]) - median(diff[nonsaline_mask]),
            })

    sal_pool = np.concatenate(pooled_arrays["salinity"])
    pooled_rows = [
        {"year": "pooled", "variable": "SRSI", "n": int(sal_pool.size), "spearman_rho_with_salinity": spearman(np.concatenate(pooled_arrays["SRSI"]), sal_pool), "expected_direction": "positive"},
        {"year": "pooled", "variable": "RSEI", "n": int(sal_pool.size), "spearman_rho_with_salinity": spearman(np.concatenate(pooled_arrays["RSEI"]), sal_pool), "expected_direction": "negative"},
        {"year": "pooled", "variable": "HARSEI", "n": int(sal_pool.size), "spearman_rho_with_salinity": spearman(np.concatenate(pooled_arrays["HARSEI"]), sal_pool), "expected_direction": "negative"},
        {"year": "pooled", "variable": "HARSEI_minus_RSEI", "n": int(sal_pool.size), "spearman_rho_with_salinity": spearman(np.concatenate(pooled_arrays["DIFF"]), sal_pool), "expected_direction": "diagnostic; sign depends on HAI adjustment"},
    ]
    correlation_rows.extend(pooled_rows)

    write_csv(out_dir / "isric_spearman_correlations.csv", correlation_rows, ["year", "variable", "n", "spearman_rho_with_salinity", "expected_direction"])
    write_csv(out_dir / "isric_salinity_class_index_statistics.csv", class_rows, [
        "year", "class_rank", "external_salinity_class", "n", "salinity_mean", "salinity_median",
        "SRSI_mean", "SRSI_median", "RSEI_mean", "RSEI_median", "HARSEI_mean", "HARSEI_median",
        "HARSEI_minus_RSEI_mean", "HARSEI_minus_RSEI_median",
    ])
    write_csv(out_dir / "isric_saline_vs_nonsaline_harsei_minus_rsei.csv", diff_rows, [
        "year", "group", "threshold_rule", "n",
        "SRSI_mean", "SRSI_median", "RSEI_mean", "RSEI_median", "HARSEI_mean", "HARSEI_median",
        "HARSEI_minus_RSEI_mean", "HARSEI_minus_RSEI_median",
    ])
    write_csv(out_dir / "isric_harsei_vs_rsei_salinity_gradient_test.csv", gradient_rows, [
        "year", "class_count", "HARSEI_median_slope_per_class", "RSEI_median_slope_per_class",
        "HARSEI_spearman_across_class_medians", "RSEI_spearman_across_class_medians",
        "HARSEI_monotonic_decrease_steps", "RSEI_monotonic_decrease_steps", "HARSEI_clearer_decline_than_RSEI",
    ])

    try_make_figures(out_dir / "isric_salinity_class_index_statistics.csv", out_dir)

    summary_lines = [
        "# ISRIC External Salinity Validation Summary",
        "",
        "## Inputs",
        f"- Years: {', '.join(str(y) for y in years)}",
        f"- Salinity directory: `{sal_dir}`",
        f"- Eco-components directory: `{eco_dir}`",
        f"- HARSEI directory: `{harsei_dir}`",
        "",
        "## RSEI reconstruction",
        f"- RSEI pooled PC1 explained variance ratio: {rsei_pca['explained_pc1']:.4f}",
        "- RSEI PC1 loadings: " + ", ".join(f"{v}={rsei_pca['pc1'][i]:.4f}" for i, v in enumerate(RSEI_BANDS)),
        "",
        "## Output tables",
        "- `isric_spearman_correlations.csv`",
        "- `isric_salinity_class_index_statistics.csv`",
        "- `isric_saline_vs_nonsaline_harsei_minus_rsei.csv`",
        "- `isric_harsei_vs_rsei_salinity_gradient_test.csv`",
        "",
        "## Manuscript interpretation guide",
        "- SRSI should show a positive Spearman correlation with external salinity.",
        "- HARSEI should show a negative Spearman correlation with external salinity.",
        "- HARSEI - RSEI is a diagnostic difference; its sign depends on how salinity stress and human-activity adjustment interact.",
        "- A clearer HARSEI salinity gradient is supported when HARSEI has a more negative median slope across salinity classes than RSEI.",
        "",
    ]
    (out_dir / "isric_validation_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")

    print("External salinity validation complete.")
    print("Output:", out_dir)


if __name__ == "__main__":
    main()
