#!/usr/bin/env python
"""
Rebuild HARSEI/OURS with year-by-year entropy weights.

This follows the user's R workflow:

  r2_resampled <- resample(r2, r1, method = "bilinear")
  stack_r <- stack(r1, r2_resampled)
  values <- getValues(stack_r)
  valid_rows <- complete.cases(values)
  data_valid <- values[valid_rows, ]
  weights <- entropy_weights(data_valid)
  fusion[valid_rows] <- data_valid %*% weights

For ecological direction consistency the active default uses HAI_reverse = 1 - HAI
as the second layer. Use --hai-orientation direct only if the HAI raster is
already ecological-positive rather than human-pressure-positive.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--awrsei-dir", required=True)
    parser.add_argument("--hai-dir", required=True)
    parser.add_argument("--out-harsei-dir", required=True)
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--prefix", default="YJQ")
    parser.add_argument("--hai-orientation", choices=["inverse", "direct"], default="inverse")
    parser.add_argument("--no-backup", action="store_true")
    return parser.parse_args()


def import_rasterio():
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import reproject

    return rasterio, Resampling, reproject


def read_single(path: Path) -> tuple[np.ndarray, dict]:
    rasterio, _, _ = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def profile_matches(a: dict, b: dict) -> bool:
    return all(a.get(key) == b.get(key) for key in ["height", "width", "crs", "transform"])


def resample_to_reference(arr: np.ndarray, src_profile: dict, ref_profile: dict) -> np.ndarray:
    if profile_matches(src_profile, ref_profile):
        return arr
    _, Resampling, reproject = import_rasterio()
    dst = np.full((ref_profile["height"], ref_profile["width"]), np.nan, dtype=np.float32)
    src_nodata = src_profile.get("nodata")
    if src_nodata is None:
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
        resampling=Resampling.bilinear,
    )
    dst[~np.isfinite(dst)] = np.nan
    return dst.astype(np.float32)


def write_single(path: Path, arr: np.ndarray, profile: dict, description: str) -> None:
    rasterio, _, _ = import_rasterio()
    out_profile = profile.copy()
    out_profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
    data = arr.astype(np.float32).copy()
    data[~np.isfinite(data)] = NODATA
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)


def entropy_weights_user(data_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Equivalent to the user's R function, including row-wise P."""
    x = data_matrix.astype(np.float64)
    col_min = np.nanmin(x, axis=0)
    col_max = np.nanmax(x, axis=0)
    norm_data = (x - col_min) / (col_max - col_min + 1e-6)
    row_sum = np.sum(norm_data + 1e-12, axis=1)
    p = norm_data / row_sum[:, None]
    k = 1.0 / math.log(p.shape[0])
    entropy = -k * np.sum(p * np.log(p + 1e-12), axis=0)
    diversity = 1.0 - entropy
    if not np.isfinite(diversity).all() or abs(float(np.sum(diversity))) < 1e-12:
        weights = np.ones(x.shape[1], dtype=np.float64) / x.shape[1]
    else:
        weights = diversity / np.sum(diversity)
    return weights.astype(np.float64), entropy.astype(np.float64), diversity.astype(np.float64)


def backup_existing_harsei(out_harsei_dir: Path, tables_dir: Path, prefix: str) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_root = out_harsei_dir.parent.parent / "backups" / f"before_yearly_entropy_harsei_{timestamp}"
    backup_harsei = backup_root / "rasters" / "HARSEI"
    backup_tables = backup_root / "tables"
    backup_harsei.mkdir(parents=True, exist_ok=False)
    backup_tables.mkdir(parents=True, exist_ok=True)
    for src in out_harsei_dir.glob(f"{prefix}_HARSEI_*.tif"):
        shutil.move(str(src), str(backup_harsei / src.name))
    for name in [
        "annual_index_summary.csv",
        "entropy_weights_and_light_harmonization.csv",
        "harsei_yearly_entropy_weights.csv",
        "harsei_yearly_entropy_formula_validation.csv",
        "harsei_yearly_entropy_method_note.md",
    ]:
        src = tables_dir / name
        if src.exists():
            shutil.copy2(src, backup_tables / name)
    return backup_root


def read_light_rows(tables_dir: Path) -> list[dict[str, str]]:
    path = tables_dir / "entropy_weights_and_light_harmonization.csv"
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return [row for row in rows if row["level"] in {"HAI_components_equal", "night_light_harmonization"}]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_method_note(
    path: Path,
    second_label: str,
    mean_weights: tuple[float, float],
    backup_root: Path | None,
) -> None:
    lines = [
        "# HARSEI Yearly Entropy Fusion Note",
        "",
        "HARSEI was recalculated year by year following the provided R workflow.",
        "",
        "For each year, AWRSEI and the HAI-derived layer were stacked, complete-case pixels were extracted, entropy weights were calculated from that year's valid pixels, and the final raster was calculated as:",
        "",
        "`HARSEI = AWRSEI * w1 + HAI_layer * w2`",
        "",
        f"The active HAI layer is `{second_label}`. In inverse mode, `HAI_reverse = 1 - HAI`, so higher HARSEI remains ecological-positive.",
        "",
        "The annual fusion uses the original AWRSEI and HAI-layer values after weights are obtained, matching `data_valid %*% weights` in the R code.",
        "",
        "Mean annual weights over 2000-2024:",
        "",
        "| Variable | Mean weight |",
        "| --- | ---: |",
        f"| AWRSEI | {mean_weights[0]:.8f} |",
        f"| {second_label} | {mean_weights[1]:.8f} |",
        "",
        "Detailed annual weights are saved in `harsei_yearly_entropy_weights.csv`.",
        "",
        "Backup of replaced HARSEI outputs:",
        "",
        str(backup_root) if backup_root else "No backup requested.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    awrsei_dir = Path(args.awrsei_dir)
    hai_dir = Path(args.hai_dir)
    out_harsei_dir = Path(args.out_harsei_dir)
    tables_dir = Path(args.tables_dir)
    out_harsei_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    backup_root = None
    if not args.no_backup:
        backup_root = backup_existing_harsei(out_harsei_dir, tables_dir, args.prefix)

    annual_rows: list[dict] = []
    weight_rows: list[dict] = []
    validation_rows: list[dict] = []
    harsei_by_year: dict[int, np.ndarray] = {}
    second_label = "HAI_reverse" if args.hai_orientation == "inverse" else "HAI"

    for year in YEARS:
        aw_path = awrsei_dir / f"{args.prefix}_AWRSEI_{year}.tif"
        hai_path = hai_dir / f"{args.prefix}_HAI_{year}.tif"
        if not aw_path.exists():
            raise FileNotFoundError(aw_path)
        if not hai_path.exists():
            raise FileNotFoundError(hai_path)

        aw, aw_profile = read_single(aw_path)
        hai, hai_profile = read_single(hai_path)
        hai_aligned = resample_to_reference(hai, hai_profile, aw_profile)
        second = (1.0 - hai_aligned).astype(np.float32) if args.hai_orientation == "inverse" else hai_aligned
        second[~np.isfinite(hai_aligned)] = np.nan

        mask = np.isfinite(aw) & np.isfinite(second)
        data_valid = np.column_stack([aw[mask], second[mask]])
        if data_valid.shape[0] < 2:
            raise ValueError(f"Not enough valid pixels for yearly entropy fusion in {year}.")
        weights, entropy, diversity = entropy_weights_user(data_valid)
        fusion = np.full(aw.shape, np.nan, dtype=np.float32)
        fusion[mask] = (data_valid @ weights).astype(np.float32)
        harsei_by_year[year] = fusion
        write_single(out_harsei_dir / f"{args.prefix}_HARSEI_{year}.tif", fusion, aw_profile, f"HARSEI_yearly_entropy_AWRSEI_{second_label}")

        expected = np.full(aw.shape, np.nan, dtype=np.float32)
        expected[mask] = (data_valid @ weights).astype(np.float32)
        validation_rows.append({
            "year": year,
            "valid_pixels": int(mask.sum()),
            "profiles_originally_matched": profile_matches(aw_profile, hai_profile),
            "formula_max_abs_diff": float(np.nanmax(np.abs(expected[mask] - fusion[mask]))),
            "harsei_min": float(np.nanmin(fusion[mask])),
            "harsei_mean": float(np.nanmean(fusion[mask])),
            "harsei_max": float(np.nanmax(fusion[mask])),
        })
        weight_rows.append({
            "year": year,
            "variable": "AWRSEI",
            "weight": float(weights[0]),
            "entropy": float(entropy[0]),
            "diversity": float(diversity[0]),
            "valid_pixels": int(mask.sum()),
        })
        weight_rows.append({
            "year": year,
            "variable": second_label,
            "weight": float(weights[1]),
            "entropy": float(entropy[1]),
            "diversity": float(diversity[1]),
            "valid_pixels": int(mask.sum()),
        })
        annual_rows.append({
            "year": year,
            "valid_pixels": int(mask.sum()),
            "AWRSEI_mean": float(np.nanmean(aw[mask])),
            "AWRSEI_std": float(np.nanstd(aw[mask])),
            "HAI_mean": float(np.nanmean(hai_aligned[mask])),
            "HAI_std": float(np.nanstd(hai_aligned[mask])),
            "HAI_reverse_mean": float(np.nanmean((1.0 - hai_aligned)[mask])),
            "HAI_reverse_std": float(np.nanstd((1.0 - hai_aligned)[mask])),
            "HARSEI_mean": float(np.nanmean(fusion[mask])),
            "HARSEI_std": float(np.nanstd(fusion[mask])),
        })

    write_csv(
        tables_dir / "harsei_yearly_entropy_weights.csv",
        weight_rows,
        ["year", "variable", "weight", "entropy", "diversity", "valid_pixels"],
    )
    write_csv(
        tables_dir / "harsei_yearly_entropy_formula_validation.csv",
        validation_rows,
        ["year", "valid_pixels", "profiles_originally_matched", "formula_max_abs_diff", "harsei_min", "harsei_mean", "harsei_max"],
    )
    write_csv(
        tables_dir / "annual_index_summary.csv",
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

    mean_awrsei_weight = float(np.mean([row["weight"] for row in weight_rows if row["variable"] == "AWRSEI"]))
    mean_second_weight = float(np.mean([row["weight"] for row in weight_rows if row["variable"] == second_label]))
    combined_weight_rows = read_light_rows(tables_dir)
    combined_weight_rows.extend([
        {"level": "HARSEI_fusion_yearly_entropy_mean", "variable": "AWRSEI", "weight": mean_awrsei_weight},
        {"level": "HARSEI_fusion_yearly_entropy_mean", "variable": second_label, "weight": mean_second_weight},
    ])
    for row in weight_rows:
        combined_weight_rows.append({
            "level": "HARSEI_fusion_yearly_entropy",
            "variable": f"{row['year']}_{row['variable']}",
            "weight": row["weight"],
        })
    write_csv(tables_dir / "entropy_weights_and_light_harmonization.csv", combined_weight_rows, ["level", "variable", "weight"])

    write_method_note(
        tables_dir / "harsei_yearly_entropy_method_note.md",
        second_label,
        (mean_awrsei_weight, mean_second_weight),
        backup_root,
    )

    print("Rebuilt HARSEI with year-by-year entropy fusion.")
    print(f"HAI orientation: {args.hai_orientation} ({second_label})")
    print(f"Mean weights: AWRSEI={mean_awrsei_weight:.8f}, {second_label}={mean_second_weight:.8f}")
    print(f"HARSEI rasters: {out_harsei_dir}")
    print(f"Tables: {tables_dir}")
    if backup_root:
        print(f"Backup: {backup_root}")


if __name__ == "__main__":
    main()
