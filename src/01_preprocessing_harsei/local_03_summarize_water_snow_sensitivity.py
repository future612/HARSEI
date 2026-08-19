#!/usr/bin/env python
"""
Summarize the GEE water-mask and snow sensitivity CSV.

Input:
  YJQ_water_snow_sensitivity_2000_2024.csv

Outputs:
  water_snow_sensitivity_summary.csv
  water_snow_sensitivity_summary.md
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean


NUMERIC_FIELDS = [
    "ndwi_threshold",
    "land_area_km2",
    "water_area_km2",
    "high_elevation_area_km2",
    "snow_fraction_all",
    "snow_fraction_high_elevation",
    "snow_area_fraction_gt_0_1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        required=True,
        help="Path to YJQ_water_snow_sensitivity_2000_2024.csv downloaded from GEE.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Folder for summary outputs.",
    )
    parser.add_argument("--main-threshold", type=float, default=0.2)
    parser.add_argument("--main-start-month", type=int, default=5)
    return parser.parse_args()


def to_float(value: str) -> float:
    if value is None or value == "":
        return float("nan")
    return float(value)


def read_rows(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            parsed: dict[str, object] = dict(row)
            parsed["year"] = int(float(str(row["year"])))
            parsed["start_month"] = int(float(str(row["start_month"])))
            for field in NUMERIC_FIELDS:
                if field in row:
                    parsed[field] = to_float(row[field])
            rows.append(parsed)
    return rows


def avg(values: list[float]) -> float:
    clean = [v for v in values if v == v]
    return mean(clean) if clean else float("nan")


def pct_change(new: float, old: float) -> float:
    if old != old or abs(old) < 1e-12:
        return float("nan")
    return (new - old) / old * 100.0


def summarize(rows: list[dict[str, object]], main_threshold: float) -> list[dict[str, float | int]]:
    grouped: dict[tuple[int, float], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["start_month"]), float(row["ndwi_threshold"]))].append(row)

    out = []
    for (start_month, threshold), group in sorted(grouped.items()):
        out.append(
            {
                "start_month": start_month,
                "ndwi_threshold": threshold,
                "mean_land_area_km2": avg([float(r["land_area_km2"]) for r in group]),
                "mean_water_area_km2": avg([float(r["water_area_km2"]) for r in group]),
                "mean_snow_fraction_all": avg([float(r["snow_fraction_all"]) for r in group]),
                "mean_snow_fraction_high_elevation": avg(
                    [float(r["snow_fraction_high_elevation"]) for r in group]
                ),
                "mean_snow_area_fraction_gt_0_1": avg(
                    [float(r["snow_area_fraction_gt_0_1"]) for r in group]
                ),
            }
        )

    baseline_by_month = {
        int(r["start_month"]): r
        for r in out
        if abs(float(r["ndwi_threshold"]) - main_threshold) < 1e-9
    }
    for row in out:
        base = baseline_by_month.get(int(row["start_month"]))
        if base:
            row["land_area_pct_vs_threshold_0_2"] = pct_change(
                float(row["mean_land_area_km2"]), float(base["mean_land_area_km2"])
            )
            row["water_area_pct_vs_threshold_0_2"] = pct_change(
                float(row["mean_water_area_km2"]), float(base["mean_water_area_km2"])
            )
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "start_month",
        "ndwi_threshold",
        "mean_land_area_km2",
        "mean_water_area_km2",
        "mean_snow_fraction_all",
        "mean_snow_fraction_high_elevation",
        "mean_snow_area_fraction_gt_0_1",
        "land_area_pct_vs_threshold_0_2",
        "water_area_pct_vs_threshold_0_2",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fmt(value: float) -> str:
    if value != value:
        return "NA"
    return f"{value:.4f}"


def write_md(path: Path, rows: list[dict[str, object]], main_start_month: int, main_threshold: float) -> None:
    by_key = {
        (int(r["start_month"]), float(r["ndwi_threshold"])): r
        for r in rows
    }
    main = by_key.get((main_start_month, main_threshold))
    apr = by_key.get((4, main_threshold))
    may = by_key.get((5, main_threshold))

    lines = ["# Water Mask And Snow Sensitivity Summary", ""]
    lines.append("| Start month | NDWI threshold | Mean land area km2 | Mean water area km2 | Snow fraction all | Snow fraction high elevation |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            "| {start_month} | {threshold:.1f} | {land} | {water} | {snow_all} | {snow_high} |".format(
                start_month=int(row["start_month"]),
                threshold=float(row["ndwi_threshold"]),
                land=fmt(float(row["mean_land_area_km2"])),
                water=fmt(float(row["mean_water_area_km2"])),
                snow_all=fmt(float(row["mean_snow_fraction_all"])),
                snow_high=fmt(float(row["mean_snow_fraction_high_elevation"])),
            )
        )

    lines.extend(["", "## Reportable Sentences", ""])
    if apr and may:
        snow_drop = pct_change(
            float(may["mean_snow_fraction_high_elevation"]),
            float(apr["mean_snow_fraction_high_elevation"]),
        )
        lines.append(
            "At NDWI > {thr:.1f}, mean high-elevation snow frequency changed from {apr} in Apr-Sep to {may} in May-Sep ({chg}%).".format(
                thr=main_threshold,
                apr=fmt(float(apr["mean_snow_fraction_high_elevation"])),
                may=fmt(float(may["mean_snow_fraction_high_elevation"])),
                chg=fmt(snow_drop),
            )
        )
    if main:
        lines.append(
            "Under the selected window (start_month = {sm}) and NDWI > {thr:.1f}, mean retained land area was {land} km2 and mean masked water area was {water} km2.".format(
                sm=main_start_month,
                thr=main_threshold,
                land=fmt(float(main["mean_land_area_km2"])),
                water=fmt(float(main["mean_water_area_km2"])),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    csv_path = Path(args.csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(csv_path)
    summary = summarize(rows, args.main_threshold)
    write_csv(out_dir / "water_snow_sensitivity_summary.csv", summary)
    write_md(out_dir / "water_snow_sensitivity_summary.md", summary, args.main_start_month, args.main_threshold)


if __name__ == "__main__":
    main()
