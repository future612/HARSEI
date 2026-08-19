#!/usr/bin/env python
"""
Validate downloaded HARSEI inputs and summarize index sensitivity for stricter
NDWI thresholds that are available within the current NDWI <= 0.2 export mask.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import rasterio


YEARS = list(range(2000, 2025))
THRESHOLDS = [0.0, 0.1, 0.2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-dir", required=True)
    parser.add_argument("--raster-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def read_band(path: Path, band: int) -> np.ndarray:
    with rasterio.open(path) as src:
        arr = src.read(band).astype(np.float32)
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr


def finite_min_max(arr: np.ndarray) -> tuple[float, float, float]:
    vals = arr[np.isfinite(arr)]
    if vals.size == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.nanmin(vals)), float(np.nanmean(vals)), float(np.nanmax(vals))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    if value != value:
        return "NA"
    return f"{value:.6f}"


def main() -> None:
    args = parse_args()
    download_dir = Path(args.download_dir)
    raster_dir = Path(args.raster_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    validation_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []

    for year in YEARS:
        eco_path = download_dir / f"YJQ_ecocomponents_{year}.tif"
        hai_path = download_dir / f"YJQ_hai_inputs_{year}.tif"
        awrsei_path = raster_dir / f"YJQ_AWRSEI_{year}.tif"
        hai_rev_path = raster_dir / f"YJQ_HAI_reverse_{year}.tif"
        harsei_path = raster_dir / f"YJQ_HARSEI_{year}.tif"

        ndwi = read_band(eco_path, 6)
        lucc_score = read_band(hai_path, 5)
        lucc_class = read_band(hai_path, 6)
        pop_year = read_band(hai_path, 8)
        lucc_year = read_band(hai_path, 9)
        awrsei = read_band(awrsei_path, 1)
        hai_rev = read_band(hai_rev_path, 1)
        harsei = read_band(harsei_path, 1)

        lucc_min, lucc_mean, lucc_max = finite_min_max(lucc_score)
        class_min, class_mean, class_max = finite_min_max(lucc_class)
        pop_vals = pop_year[np.isfinite(pop_year)]
        lucc_year_vals = lucc_year[np.isfinite(lucc_year)]
        validation_rows.append(
            {
                "year": year,
                "pop_year_used": int(np.nanmedian(pop_vals)) if pop_vals.size else "",
                "lucc_year_used": int(np.nanmedian(lucc_year_vals)) if lucc_year_vals.size else "",
                "lucc_score_min": lucc_min,
                "lucc_score_mean": lucc_mean,
                "lucc_score_max": lucc_max,
                "lucc_class_min": class_min,
                "lucc_class_mean": class_mean,
                "lucc_class_max": class_max,
            }
        )

        for threshold in THRESHOLDS:
            mask = np.isfinite(ndwi) & (ndwi <= threshold)
            row: dict[str, object] = {
                "year": year,
                "ndwi_threshold": threshold,
                "valid_pixels": int(mask.sum()),
            }
            for name, arr in [
                ("AWRSEI", awrsei),
                ("HAI_reverse", hai_rev),
                ("HARSEI", harsei),
            ]:
                vals = arr[mask & np.isfinite(arr)]
                row[f"{name}_mean"] = float(np.nanmean(vals)) if vals.size else float("nan")
            sensitivity_rows.append(row)

    write_csv(
        out_dir / "input_stack_validation_summary.csv",
        validation_rows,
        [
            "year",
            "pop_year_used",
            "lucc_year_used",
            "lucc_score_min",
            "lucc_score_mean",
            "lucc_score_max",
            "lucc_class_min",
            "lucc_class_mean",
            "lucc_class_max",
        ],
    )
    write_csv(
        out_dir / "water_threshold_index_sensitivity_current_mask.csv",
        sensitivity_rows,
        [
            "year",
            "ndwi_threshold",
            "valid_pixels",
            "AWRSEI_mean",
            "HAI_reverse_mean",
            "HARSEI_mean",
        ],
    )

    by_threshold: dict[float, list[dict[str, object]]] = {t: [] for t in THRESHOLDS}
    for row in sensitivity_rows:
        by_threshold[float(row["ndwi_threshold"])].append(row)

    lines = ["# Current-Stack NDWI Threshold Sensitivity", ""]
    lines.append(
        "This summary uses only pixels available in the current main export (`NDWI <= 0.2`). "
        "It therefore evaluates stricter thresholds 0.0 and 0.1 against the main 0.2 threshold."
    )
    lines.append("")
    lines.append("| NDWI threshold | Mean valid pixels | Mean AWRSEI | Mean HAI_reverse | Mean HARSEI |")
    lines.append("|---:|---:|---:|---:|---:|")
    summary = {}
    for threshold in THRESHOLDS:
        rows = by_threshold[threshold]
        valid = np.array([float(r["valid_pixels"]) for r in rows], dtype=float)
        aw = np.array([float(r["AWRSEI_mean"]) for r in rows], dtype=float)
        hr = np.array([float(r["HAI_reverse_mean"]) for r in rows], dtype=float)
        hs = np.array([float(r["HARSEI_mean"]) for r in rows], dtype=float)
        summary[threshold] = {
            "valid": float(np.nanmean(valid)),
            "awrsei": float(np.nanmean(aw)),
            "hai_reverse": float(np.nanmean(hr)),
            "harsei": float(np.nanmean(hs)),
        }
        lines.append(
            f"| {threshold:.1f} | {fmt(summary[threshold]['valid'])} | "
            f"{fmt(summary[threshold]['awrsei'])} | {fmt(summary[threshold]['hai_reverse'])} | "
            f"{fmt(summary[threshold]['harsei'])} |"
        )

    base = summary[0.2]["harsei"]
    lines.extend(["", "## Reportable Sentence", ""])
    for threshold in [0.0, 0.1]:
        diff = summary[threshold]["harsei"] - base
        lines.append(
            f"Compared with the main NDWI > 0.2 water mask, using NDWI > {threshold:.1f} "
            f"changed the multi-year mean HARSEI by {diff:.6f}."
        )
    (out_dir / "water_threshold_index_sensitivity_current_mask.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
