#!/usr/bin/env python
"""Validate revised HAI and HARSEI GeoTIFFs against the documented formulas."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


NODATA = -9999.0
COMPONENTS = ["POP", "LIGHT", "LUCC_SCORE"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--awrsei-dir", required=True)
    parser.add_argument("--hai-dir", required=True)
    parser.add_argument("--harsei-dir", required=True)
    parser.add_argument("--hai-input-dir", required=True)
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--years", default="2000,2013,2014,2024")
    parser.add_argument("--prefix", default="YJQ")
    return parser.parse_args()


def import_rasterio():
    import rasterio

    return rasterio


def read_single(path: Path) -> tuple[np.ndarray, dict]:
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


def read_hai_input(path: Path) -> tuple[np.ndarray, dict]:
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read([1, 2, 3, 4, 5]).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def read_component_ranges(path: Path) -> dict[str, tuple[float, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return {
        row["component"]: (float(row["raw_min_used_for_norm"]), float(row["raw_max_used_for_norm"]))
        for row in rows
    }


def read_weights(path: Path) -> tuple[tuple[float, float], tuple[float, float]]:
    light_a = light_b = None
    w_a = w_h = None
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            variable = row["variable"]
            value = float(row["weight"])
            if variable == "DMSP=a+b*log1p(VIIRS):a":
                light_a = value
            elif variable == "DMSP=a+b*log1p(VIIRS):b":
                light_b = value
            elif row["level"] == "HARSEI_fusion_user_entropy" and variable == "AWRSEI":
                w_a = value
            elif row["level"] == "HARSEI_fusion_user_entropy" and variable == "HAI_reverse":
                w_h = value
    missing = [name for name, value in {
        "light_a": light_a,
        "light_b": light_b,
        "w_AWRSEI": w_a,
        "w_HAI_reverse": w_h,
    }.items() if value is None]
    if missing:
        raise ValueError(f"Missing weight values: {missing}")
    return (float(light_a), float(light_b)), (float(w_a), float(w_h))


def normalize(arr: np.ndarray, lo: float, hi: float) -> np.ndarray:
    out = (arr - lo) / (hi - lo + 1e-6)
    return np.clip(out, 0.0, 1.0).astype(np.float32)


def profile_key(profile: dict) -> tuple:
    return profile["height"], profile["width"], profile["crs"], profile["transform"]


def main() -> None:
    args = parse_args()
    years = [int(part.strip()) for part in args.years.split(",") if part.strip()]
    awrsei_dir = Path(args.awrsei_dir)
    hai_dir = Path(args.hai_dir)
    harsei_dir = Path(args.harsei_dir)
    hai_input_dir = Path(args.hai_input_dir)
    tables_dir = Path(args.tables_dir)

    component_ranges = read_component_ranges(tables_dir / "hai_equal_component_ranges.csv")
    light_fit, fusion_weights = read_weights(tables_dir / "entropy_weights_and_light_harmonization.csv")

    all_years = list(range(2000, 2025))
    aw_min, aw_max = np.inf, -np.inf
    hr_min, hr_max = np.inf, -np.inf
    for year in all_years:
        aw, _ = read_single(awrsei_dir / f"{args.prefix}_AWRSEI_{year}.tif")
        hai, _ = read_single(hai_dir / f"{args.prefix}_HAI_{year}.tif")
        mask = np.isfinite(aw) & np.isfinite(hai)
        aw_min = min(aw_min, float(np.nanmin(aw[mask])))
        aw_max = max(aw_max, float(np.nanmax(aw[mask])))
        hr = 1.0 - hai
        hr_min = min(hr_min, float(np.nanmin(hr[mask])))
        hr_max = max(hr_max, float(np.nanmax(hr[mask])))

    rows = []
    reference_profile = None
    for year in years:
        aw, aw_profile = read_single(awrsei_dir / f"{args.prefix}_AWRSEI_{year}.tif")
        hai, hai_profile = read_single(hai_dir / f"{args.prefix}_HAI_{year}.tif")
        harsei, harsei_profile = read_single(harsei_dir / f"{args.prefix}_HARSEI_{year}.tif")
        hai_input, input_profile = read_hai_input(hai_input_dir / f"{args.prefix}_hai_inputs_{year}.tif")
        if reference_profile is None:
            reference_profile = aw_profile
        profiles_ok = (
            profile_key(reference_profile)
            == profile_key(aw_profile)
            == profile_key(hai_profile)
            == profile_key(harsei_profile)
            == profile_key(input_profile)
        )

        light = hai_input[2] if year <= 2013 else light_fit[0] + light_fit[1] * np.log1p(hai_input[3])
        light = np.clip(light, 0.0, None).astype(np.float32)
        pop_norm = normalize(hai_input[0], *component_ranges["POP"])
        light_norm = normalize(light, *component_ranges["LIGHT"])
        lucc_norm = normalize(hai_input[4], *component_ranges["LUCC_SCORE"])
        expected_hai = ((pop_norm + light_norm + lucc_norm) / 3.0).astype(np.float32)
        source_valid = (
            np.isfinite(aw)
            & np.isfinite(hai_input[0])
            & np.isfinite(light)
            & np.isfinite(hai_input[4])
        )
        expected_hai[~source_valid] = np.nan

        norm_aw = normalize(aw, aw_min, aw_max)
        norm_hai_reverse = normalize(1.0 - expected_hai, hr_min, hr_max)
        expected_harsei = (fusion_weights[0] * norm_aw + fusion_weights[1] * norm_hai_reverse).astype(np.float32)
        expected_harsei[~source_valid] = np.nan

        compare_hai = np.isfinite(expected_hai) & np.isfinite(hai)
        compare_harsei = np.isfinite(expected_harsei) & np.isfinite(harsei)
        rows.append({
            "year": year,
            "profiles_match": profiles_ok,
            "valid_pixels": int(compare_hai.sum()),
            "hai_min": float(np.nanmin(hai[compare_hai])),
            "hai_mean": float(np.nanmean(hai[compare_hai])),
            "hai_max": float(np.nanmax(hai[compare_hai])),
            "hai_formula_max_abs_diff": float(np.nanmax(np.abs(expected_hai[compare_hai] - hai[compare_hai]))),
            "harsei_min": float(np.nanmin(harsei[compare_harsei])),
            "harsei_mean": float(np.nanmean(harsei[compare_harsei])),
            "harsei_max": float(np.nanmax(harsei[compare_harsei])),
            "harsei_formula_max_abs_diff": float(np.nanmax(np.abs(expected_harsei[compare_harsei] - harsei[compare_harsei]))),
        })

    out_path = tables_dir / "revised_hai_harsei_formula_validation.csv"
    with out_path.open("w", encoding="utf-8", newline="") as f:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out_path}")
    for row in rows:
        print(
            f"{row['year']}: profiles={row['profiles_match']}, "
            f"HAI diff={row['hai_formula_max_abs_diff']:.8g}, "
            f"HARSEI diff={row['harsei_formula_max_abs_diff']:.8g}"
        )


if __name__ == "__main__":
    main()
