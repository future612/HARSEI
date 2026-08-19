#!/usr/bin/env python
"""
Rebuild HAI and HARSEI while preserving the existing AWRSEI rasters.

Revised method:
  HAI = (POP_norm + LIGHT_norm + LUCC_SCORE_norm) / 3

HARSEI is fused from AWRSEI and the ecological-positive human activity term
(1 - HAI) using the user's entropy-weight function shown in the revision notes:

  norm_data = apply(data_matrix, 2, function(x) (x-min(x))/(max(x)-min(x)+1e-6))
  P = norm_data / rowSums(norm_data + 1e-12)
  k = 1 / log(nrow(P))
  E = -k * colSums(P * log(P + 1e-12))
  d = 1 - E
  w = d / sum(d)

The script moves old HAI/HARSEI GeoTIFFs into a timestamped backup folder before
writing revised outputs, and copies replaced CSV tables into the same backup.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np


YEARS = list(range(2000, 2025))
NODATA = -9999.0
COMPONENT_NAMES = ["POP", "LIGHT", "LUCC_SCORE"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--awrsei-dir", required=True, help="Folder containing YJQ_AWRSEI_YYYY.tif.")
    parser.add_argument("--hai-input-dir", required=True, help="Folder containing YJQ_hai_inputs_YYYY.tif.")
    parser.add_argument("--out-rasters-dir", required=True, help="Parent raster folder with AWRSEI/HAI/HARSEI subfolders.")
    parser.add_argument("--out-tables-dir", required=True, help="Output table folder.")
    parser.add_argument("--prefix", default="YJQ")
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--light-fit-sample-per-year", type=int, default=250000)
    parser.add_argument(
        "--fusion-hai-orientation",
        choices=["inverse", "direct"],
        default="inverse",
        help="Use inverse so higher HARSEI always means better ecological condition.",
    )
    parser.add_argument("--no-backup", action="store_true", help="Overwrite outputs without moving old files first.")
    return parser.parse_args()


def import_rasterio():
    try:
        import rasterio
    except Exception as exc:  # pragma: no cover
        raise SystemExit("rasterio is required for GeoTIFF I/O in this script.") from exc
    return rasterio


def finite(arr: np.ndarray) -> np.ndarray:
    return np.isfinite(arr)


def safe_resolve(path: Path, allowed_root: Path) -> Path:
    resolved = path.resolve()
    root = allowed_root.resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError(f"Refusing to modify path outside allowed root: {resolved}")
    return resolved


def read_single(path: Path) -> tuple[np.ndarray, dict]:
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def read_hai_input(path: Path) -> tuple[np.ndarray, dict, tuple[str | None, ...]]:
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        if src.count < 5:
            raise ValueError(f"{path} has {src.count} band(s), expected at least 5.")
        arr = src.read([1, 2, 3, 4, 5]).astype(np.float32)
        profile = src.profile.copy()
        descriptions = src.descriptions
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile, descriptions


def profiles_match(reference: dict, other: dict) -> bool:
    keys = ["height", "width", "crs", "transform"]
    return all(reference.get(k) == other.get(k) for k in keys)


def write_single(path: Path, arr: np.ndarray, profile: dict, description: str) -> None:
    rasterio = import_rasterio()
    out_profile = profile.copy()
    out_profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
    data = arr.astype(np.float32).copy()
    data[~np.isfinite(data)] = NODATA
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)


def minmax(values: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    valid = values[mask & np.isfinite(values)]
    if valid.size == 0:
        raise ValueError("No valid values available for min-max calculation.")
    return float(np.nanmin(valid)), float(np.nanmax(valid))


def normalize(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
    out = (values - lo) / (hi - lo + 1e-6)
    return np.clip(out, 0.0, 1.0).astype(np.float32)


def fit_light_harmonization(
    hai_inputs: dict[int, np.ndarray],
    masks: dict[int, np.ndarray],
    sample_per_year: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for year in [2012, 2013]:
        arr = hai_inputs[year]
        dmsp = arr[2]
        viirs = arr[3]
        mask = masks[year] & np.isfinite(dmsp) & np.isfinite(viirs) & (viirs > 0)
        idx = np.flatnonzero(mask.reshape(-1))
        if idx.size > sample_per_year:
            idx = rng.choice(idx, size=sample_per_year, replace=False)
        if idx.size:
            xs.append(np.log1p(viirs.reshape(-1)[idx]).astype(np.float64))
            ys.append(dmsp.reshape(-1)[idx].astype(np.float64))
    if not xs:
        raise ValueError("Cannot fit night-light harmonization: no valid 2012-2013 overlap samples.")
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    if x.size < 100:
        raise ValueError(f"Cannot fit night-light harmonization robustly: only {x.size} samples.")
    b, a = np.polyfit(x, y, 1)
    return float(a), float(b)


def build_light(year: int, hai_arr: np.ndarray, light_fit: tuple[float, float]) -> np.ndarray:
    dmsp = hai_arr[2]
    viirs = hai_arr[3]
    if year <= 2013:
        return dmsp.astype(np.float32)
    a, b = light_fit
    out = a + b * np.log1p(viirs)
    out = np.clip(out, 0.0, None)
    return out.astype(np.float32)


def backup_replaced_outputs(
    out_rasters_dir: Path,
    out_tables_dir: Path,
    prefix: str,
    allowed_root: Path,
) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_root = out_rasters_dir.parent / "backups" / f"before_equal_hai_{timestamp}"
    safe_resolve(backup_root, allowed_root)
    backup_root.mkdir(parents=True, exist_ok=False)

    for subdir, patterns in {
        "HAI": [f"{prefix}_HAI*.tif"],
        "HARSEI": [f"{prefix}_HARSEI*.tif"],
    }.items():
        src_dir = safe_resolve(out_rasters_dir / subdir, allowed_root)
        dst_dir = backup_root / "rasters" / subdir
        dst_dir.mkdir(parents=True, exist_ok=True)
        for pattern in patterns:
            for src in src_dir.glob(pattern):
                shutil.move(str(src), str(dst_dir / src.name))

    table_backup = backup_root / "tables"
    table_backup.mkdir(parents=True, exist_ok=True)
    for name in [
        "annual_index_summary.csv",
        "entropy_weights_and_light_harmonization.csv",
        "hai_equal_component_ranges.csv",
        "hai_harsei_rebuild_validation.csv",
        "hai_harsei_rebuild_method_note.md",
    ]:
        src = out_tables_dir / name
        if src.exists():
            shutil.copy2(src, table_backup / name)
    return backup_root


def user_entropy_weights_stream(
    first: dict[int, np.ndarray],
    second: dict[int, np.ndarray],
    masks: dict[int, np.ndarray],
    first_range: tuple[float, float],
    second_range: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    entropy_sum = np.zeros(2, dtype=np.float64)
    n = 0
    for year in YEARS:
        mask = masks[year] & np.isfinite(first[year]) & np.isfinite(second[year])
        a = normalize(first[year], *first_range)[mask].astype(np.float64)
        b = normalize(second[year], *second_range)[mask].astype(np.float64)
        denom = a + b + 2e-12
        p_a = a / denom
        p_b = b / denom
        entropy_sum[0] += np.sum(p_a * np.log(p_a + 1e-12))
        entropy_sum[1] += np.sum(p_b * np.log(p_b + 1e-12))
        n += int(a.size)
    if n <= 1:
        raise ValueError("Not enough valid rows to compute entropy weights.")
    k = 1.0 / math.log(n)
    entropy = -k * entropy_sum
    diversity = 1.0 - entropy
    if not np.isfinite(diversity).all() or abs(float(diversity.sum())) < 1e-12:
        weights = np.array([0.5, 0.5], dtype=np.float64)
    else:
        weights = diversity / diversity.sum()
    return weights.astype(np.float64), entropy, diversity, n


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_method_note(
    path: Path,
    light_fit: tuple[float, float],
    component_ranges: dict[str, tuple[float, float]],
    fusion_label: str,
    fusion_ranges: dict[str, tuple[float, float]],
    fusion_weights: np.ndarray,
    backup_root: Path | None,
) -> None:
    backup_text = str(backup_root) if backup_root else "No backup requested."
    lines = [
        "# Revised HAI/HARSEI Rebuild Note",
        "",
        "## Revised HAI",
        "",
        "HAI was rebuilt as the equal-weight mean of three dimensionless human-activity components:",
        "",
        "`HAI = (POP_norm + LIGHT_norm + LUCC_SCORE_norm) / 3`",
        "",
        "Each component was normalized using the pooled 2000-2024 valid-pixel min-max range and `1e-6` in the denominator.",
        "",
        "## Night-Light Harmonization",
        "",
        "`LIGHT` uses DMSP for 2000-2013 and DMSP-equivalent VIIRS for 2014-2024.",
        f"The overlap-period fit is `DMSP = {light_fit[0]:.8f} + {light_fit[1]:.8f} * log1p(VIIRS)`.",
        "",
        "## HARSEI Fusion",
        "",
        f"HARSEI was fused from `AWRSEI` and `{fusion_label}`. For the default inverse mode, `{fusion_label} = 1 - HAI`, so higher HARSEI consistently indicates better ecological conditions.",
        "The entropy weights follow the user-provided R function exactly, including row-wise proportions.",
        "",
        "| Variable | Weight | Fusion min | Fusion max |",
        "| --- | ---: | ---: | ---: |",
        f"| AWRSEI | {fusion_weights[0]:.8f} | {fusion_ranges['AWRSEI'][0]:.8f} | {fusion_ranges['AWRSEI'][1]:.8f} |",
        f"| {fusion_label} | {fusion_weights[1]:.8f} | {fusion_ranges[fusion_label][0]:.8f} | {fusion_ranges[fusion_label][1]:.8f} |",
        "",
        "## HAI Component Ranges",
        "",
        "| Component | Min | Max |",
        "| --- | ---: | ---: |",
    ]
    for name in COMPONENT_NAMES:
        lo, hi = component_ranges[name]
        lines.append(f"| {name} | {lo:.8f} | {hi:.8f} |")
    lines.extend(["", "## Backup", "", backup_text, ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    awrsei_dir = Path(args.awrsei_dir)
    hai_input_dir = Path(args.hai_input_dir)
    out_rasters_dir = Path(args.out_rasters_dir)
    out_tables_dir = Path(args.out_tables_dir)
    allowed_root = Path.cwd()

    safe_resolve(awrsei_dir, allowed_root)
    safe_resolve(hai_input_dir, allowed_root)
    safe_resolve(out_rasters_dir, allowed_root)
    safe_resolve(out_tables_dir, allowed_root)

    hai_out_dir = out_rasters_dir / "HAI"
    harsei_out_dir = out_rasters_dir / "HARSEI"
    hai_out_dir.mkdir(parents=True, exist_ok=True)
    harsei_out_dir.mkdir(parents=True, exist_ok=True)
    out_tables_dir.mkdir(parents=True, exist_ok=True)

    awrsei: dict[int, np.ndarray] = {}
    hai_inputs: dict[int, np.ndarray] = {}
    masks: dict[int, np.ndarray] = {}
    awrsei_profile: dict | None = None
    validation_rows: list[dict] = []

    for year in YEARS:
        aw_path = awrsei_dir / f"{args.prefix}_AWRSEI_{year}.tif"
        hai_path = hai_input_dir / f"{args.prefix}_hai_inputs_{year}.tif"
        if not aw_path.exists():
            raise FileNotFoundError(aw_path)
        if not hai_path.exists():
            raise FileNotFoundError(hai_path)

        aw, aw_profile = read_single(aw_path)
        hai_arr, hai_profile, descriptions = read_hai_input(hai_path)
        if awrsei_profile is None:
            awrsei_profile = aw_profile
        if not profiles_match(awrsei_profile, aw_profile):
            raise ValueError(f"AWRSEI profile differs from first year in {year}.")
        if not profiles_match(awrsei_profile, hai_profile):
            raise ValueError(f"HAI input profile does not match AWRSEI profile in {year}.")

        light_valid = np.isfinite(hai_arr[2]) if year <= 2013 else np.isfinite(hai_arr[3])
        raw_valid = (
            np.isfinite(aw)
            & np.isfinite(hai_arr[0])
            & light_valid
            & np.isfinite(hai_arr[4])
        )
        awrsei[year] = aw
        hai_inputs[year] = hai_arr
        masks[year] = raw_valid
        validation_rows.append({
            "year": year,
            "valid_pixels": int(raw_valid.sum()),
            "awrsei_min": float(np.nanmin(aw[raw_valid])),
            "awrsei_mean": float(np.nanmean(aw[raw_valid])),
            "awrsei_max": float(np.nanmax(aw[raw_valid])),
            "hai_input_band_count": hai_arr.shape[0],
            "hai_input_band_descriptions": ";".join("" if d is None else d for d in descriptions[:5]),
        })

    assert awrsei_profile is not None
    light_fit = fit_light_harmonization(hai_inputs, masks, args.light_fit_sample_per_year, rng)

    raw_components: dict[int, dict[str, np.ndarray]] = {}
    for year in YEARS:
        arr = hai_inputs[year]
        raw_components[year] = {
            "POP": arr[0].astype(np.float32),
            "LIGHT": build_light(year, arr, light_fit),
            "LUCC_SCORE": arr[4].astype(np.float32),
        }

    component_ranges: dict[str, tuple[float, float]] = {}
    for name in COMPONENT_NAMES:
        mins = []
        maxs = []
        for year in YEARS:
            lo, hi = minmax(raw_components[year][name], masks[year])
            mins.append(lo)
            maxs.append(hi)
        component_ranges[name] = (float(min(mins)), float(max(maxs)))

    hai_pressure: dict[int, np.ndarray] = {}
    hai_reverse: dict[int, np.ndarray] = {}
    for year in YEARS:
        norm_components = []
        for name in COMPONENT_NAMES:
            norm_components.append(normalize(raw_components[year][name], *component_ranges[name]))
        stack = np.stack(norm_components)
        out = np.nanmean(stack, axis=0).astype(np.float32)
        out[~masks[year]] = np.nan
        hai_pressure[year] = out
        hai_reverse[year] = (1.0 - out).astype(np.float32)
        hai_reverse[year][~masks[year]] = np.nan

    fusion_label = "HAI_reverse" if args.fusion_hai_orientation == "inverse" else "HAI"
    fusion_second = hai_reverse if args.fusion_hai_orientation == "inverse" else hai_pressure
    aw_range = (
        float(min(np.nanmin(awrsei[y][masks[y]]) for y in YEARS)),
        float(max(np.nanmax(awrsei[y][masks[y]]) for y in YEARS)),
    )
    second_range = (
        float(min(np.nanmin(fusion_second[y][masks[y]]) for y in YEARS)),
        float(max(np.nanmax(fusion_second[y][masks[y]]) for y in YEARS)),
    )
    fusion_ranges = {"AWRSEI": aw_range, fusion_label: second_range}
    fusion_weights, entropy, diversity, entropy_n = user_entropy_weights_stream(
        awrsei, fusion_second, masks, aw_range, second_range
    )

    backup_root = None
    if not args.no_backup:
        backup_root = backup_replaced_outputs(out_rasters_dir, out_tables_dir, args.prefix, allowed_root)

    harsei: dict[int, np.ndarray] = {}
    for year in YEARS:
        norm_aw = normalize(awrsei[year], *aw_range)
        norm_second = normalize(fusion_second[year], *second_range)
        out = (fusion_weights[0] * norm_aw + fusion_weights[1] * norm_second).astype(np.float32)
        out[~masks[year]] = np.nan
        harsei[year] = out
        write_single(hai_out_dir / f"{args.prefix}_HAI_{year}.tif", hai_pressure[year], awrsei_profile, "HAI_pressure_equal_POP_LIGHT_LUCC")
        write_single(harsei_out_dir / f"{args.prefix}_HARSEI_{year}.tif", harsei[year], awrsei_profile, f"HARSEI_entropy_AWRSEI_{fusion_label}")

    component_rows = []
    for name in COMPONENT_NAMES:
        lo, hi = component_ranges[name]
        component_rows.append({
            "component": name,
            "raw_min_used_for_norm": lo,
            "raw_max_used_for_norm": hi,
            "hai_component_weight": 1.0 / 3.0,
        })
    write_csv(
        out_tables_dir / "hai_equal_component_ranges.csv",
        component_rows,
        ["component", "raw_min_used_for_norm", "raw_max_used_for_norm", "hai_component_weight"],
    )

    annual_rows = []
    for year in YEARS:
        mask = masks[year] & np.isfinite(harsei[year])
        annual_rows.append({
            "year": year,
            "valid_pixels": int(mask.sum()),
            "AWRSEI_mean": float(np.nanmean(awrsei[year][mask])),
            "AWRSEI_std": float(np.nanstd(awrsei[year][mask])),
            "HAI_mean": float(np.nanmean(hai_pressure[year][mask])),
            "HAI_std": float(np.nanstd(hai_pressure[year][mask])),
            "HAI_reverse_mean": float(np.nanmean(hai_reverse[year][mask])),
            "HAI_reverse_std": float(np.nanstd(hai_reverse[year][mask])),
            "HARSEI_mean": float(np.nanmean(harsei[year][mask])),
            "HARSEI_std": float(np.nanstd(harsei[year][mask])),
        })
    write_csv(
        out_tables_dir / "annual_index_summary.csv",
        annual_rows,
        [
            "year",
            "valid_pixels",
            "AWRSEI_mean",
            "AWRSEI_std",
            "HAI_mean",
            "HAI_std",
            "HAI_reverse_mean",
            "HAI_reverse_std",
            "HARSEI_mean",
            "HARSEI_std",
        ],
    )

    weight_rows = [
        {"level": "HAI_components_equal", "variable": "POP", "weight": 1.0 / 3.0},
        {"level": "HAI_components_equal", "variable": "LIGHT", "weight": 1.0 / 3.0},
        {"level": "HAI_components_equal", "variable": "LUCC_SCORE", "weight": 1.0 / 3.0},
        {"level": "HARSEI_fusion_user_entropy", "variable": "AWRSEI", "weight": float(fusion_weights[0])},
        {"level": "HARSEI_fusion_user_entropy", "variable": fusion_label, "weight": float(fusion_weights[1])},
        {"level": "HARSEI_fusion_user_entropy", "variable": "entropy_AWRSEI", "weight": float(entropy[0])},
        {"level": "HARSEI_fusion_user_entropy", "variable": f"entropy_{fusion_label}", "weight": float(entropy[1])},
        {"level": "HARSEI_fusion_user_entropy", "variable": "diversity_AWRSEI", "weight": float(diversity[0])},
        {"level": "HARSEI_fusion_user_entropy", "variable": f"diversity_{fusion_label}", "weight": float(diversity[1])},
        {"level": "HARSEI_fusion_user_entropy", "variable": "n_rows", "weight": int(entropy_n)},
        {"level": "night_light_harmonization", "variable": "DMSP=a+b*log1p(VIIRS):a", "weight": light_fit[0]},
        {"level": "night_light_harmonization", "variable": "DMSP=a+b*log1p(VIIRS):b", "weight": light_fit[1]},
    ]
    write_csv(out_tables_dir / "entropy_weights_and_light_harmonization.csv", weight_rows, ["level", "variable", "weight"])

    for row in validation_rows:
        year = int(row["year"])
        mask = masks[year]
        row.update({
            "hai_min": float(np.nanmin(hai_pressure[year][mask])),
            "hai_mean": float(np.nanmean(hai_pressure[year][mask])),
            "hai_max": float(np.nanmax(hai_pressure[year][mask])),
            "harsei_min": float(np.nanmin(harsei[year][mask])),
            "harsei_mean": float(np.nanmean(harsei[year][mask])),
            "harsei_max": float(np.nanmax(harsei[year][mask])),
        })
    write_csv(
        out_tables_dir / "hai_harsei_rebuild_validation.csv",
        validation_rows,
        [
            "year",
            "valid_pixels",
            "awrsei_min",
            "awrsei_mean",
            "awrsei_max",
            "hai_min",
            "hai_mean",
            "hai_max",
            "harsei_min",
            "harsei_mean",
            "harsei_max",
            "hai_input_band_count",
            "hai_input_band_descriptions",
        ],
    )

    write_method_note(
        out_tables_dir / "hai_harsei_rebuild_method_note.md",
        light_fit,
        component_ranges,
        fusion_label,
        fusion_ranges,
        fusion_weights,
        backup_root,
    )

    print("Rebuilt HAI/HARSEI with equal HAI components.")
    print(f"HAI rasters: {hai_out_dir}")
    print(f"HARSEI rasters: {harsei_out_dir}")
    print(f"Tables: {out_tables_dir}")
    print(f"Fusion label: {fusion_label}")
    print(f"Fusion weights: AWRSEI={fusion_weights[0]:.8f}, {fusion_label}={fusion_weights[1]:.8f}")
    if backup_root:
        print(f"Backup: {backup_root}")


if __name__ == "__main__":
    main()
