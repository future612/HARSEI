#!/usr/bin/env python
"""Compute equal-interval grade areas and transition matrices from final HARSEI."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


YEARS = [2000, 2005, 2010, 2015, 2020, 2024]
LABELS = {
    1: "Poor",
    2: "Relatively poor",
    3: "Fair",
    4: "Good",
    5: "Excellent",
}
NODATA = -9999.0


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--harsei-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--prefix", default="YJQ")
    return p.parse_args()


def import_rasterio():
    import rasterio

    return rasterio


def read(path: Path):
    rasterio = import_rasterio()
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr


def grade(arr):
    g = np.zeros(arr.shape, dtype=np.uint8)
    m = np.isfinite(arr)
    g[m & (arr < 0.2)] = 1
    g[m & (arr >= 0.2) & (arr < 0.4)] = 2
    g[m & (arr >= 0.4) & (arr < 0.6)] = 3
    g[m & (arr >= 0.6) & (arr < 0.8)] = 4
    g[m & (arr >= 0.8)] = 5
    return g


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    harsei_dir = Path(args.harsei_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    grades = {}
    area_rows = []
    for year in YEARS:
        arr = read(harsei_dir / f"{args.prefix}_HARSEI_{year}.tif")
        g = grade(arr)
        grades[year] = g
        total = int(np.sum(g > 0))
        for code, label in LABELS.items():
            count = int(np.sum(g == code))
            area_rows.append({"year": year, "grade": code, "class": label, "pixels": count, "percent": 100.0 * count / total})
    write_csv(out_dir / "final_harsei_grade_area_summary.csv", area_rows, ["year", "grade", "class", "pixels", "percent"])

    trans_rows = []
    three_rows = []
    pairs = list(zip(YEARS[:-1], YEARS[1:])) + [(2000, 2024)]
    for y1, y2 in pairs:
        g1, g2 = grades[y1], grades[y2]
        mask = (g1 > 0) & (g2 > 0)
        total = int(mask.sum())
        worse = int(np.sum(mask & (g2 < g1)))
        same = int(np.sum(mask & (g2 == g1)))
        better = int(np.sum(mask & (g2 > g1)))
        three_rows.append({
            "period": f"{y1}-{y2}",
            "worse_pixels": worse,
            "same_pixels": same,
            "better_pixels": better,
            "worse_percent": 100.0 * worse / total,
            "same_percent": 100.0 * same / total,
            "better_percent": 100.0 * better / total,
        })
        for from_code in LABELS:
            for to_code in LABELS:
                count = int(np.sum(mask & (g1 == from_code) & (g2 == to_code)))
                trans_rows.append({
                    "period": f"{y1}-{y2}",
                    "from_grade": from_code,
                    "from_class": LABELS[from_code],
                    "to_grade": to_code,
                    "to_class": LABELS[to_code],
                    "pixels": count,
                    "percent_of_period_valid": 100.0 * count / total,
                })
    write_csv(out_dir / "final_harsei_transition_5class_long.csv", trans_rows, ["period", "from_grade", "from_class", "to_grade", "to_class", "pixels", "percent_of_period_valid"])
    write_csv(out_dir / "final_harsei_transition_3group_summary.csv", three_rows, ["period", "worse_pixels", "same_pixels", "better_pixels", "worse_percent", "same_percent", "better_percent"])
    print("Transition done")


if __name__ == "__main__":
    main()
