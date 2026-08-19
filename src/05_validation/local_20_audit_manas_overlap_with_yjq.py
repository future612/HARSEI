#!/usr/bin/env python
"""Audit whether Manas field samples overlap valid YJQ HARSEI/RSEI pixels."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio


ROOT = Path(r"D:\Codex\260724 小论文")
NODATA = -9999.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=ROOT / "revise/manas_field_salinity_validation/manas_field_samples_with_srsi.csv")
    parser.add_argument("--year", type=int, default=2019)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "revise/manas_field_salinity_validation")
    return parser.parse_args()


def sample_raster(points: pd.DataFrame, raster_path: Path) -> tuple[int, int]:
    inside_bbox = 0
    finite = 0
    with rasterio.open(raster_path) as src:
        band = src.read(1)
        for lon, lat in zip(points["lon"], points["lat"]):
            if src.bounds.left <= lon <= src.bounds.right and src.bounds.bottom <= lat <= src.bounds.top:
                inside_bbox += 1
                row, col = src.index(lon, lat)
                if 0 <= row < src.height and 0 <= col < src.width:
                    val = band[row, col]
                    if src.nodata is not None and val == src.nodata:
                        continue
                    if val == NODATA:
                        continue
                    if np.isfinite(val):
                        finite += 1
    return inside_bbox, finite


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(args.samples)
    rows = []
    rasters = {
        "HARSEI": ROOT / f"revise/annual_harsei_outputs/rasters/HARSEI/YJQ_HARSEI_{args.year}.tif",
        "RSEI": ROOT / f"revise/annual_harsei_outputs/rasters/RSEI/YJQ_RSEI_{args.year}.tif",
        "SRSI_band": ROOT / f"revise/gee_downloads/YJQ_HARSEI_annual_inputs_2000_2024/YJQ_ecocomponents_{args.year}.tif",
    }
    for dataset in ["validation", "modeling", "all"]:
        pts = samples[samples["dataset"] == dataset]
        for name, path in rasters.items():
            inside, finite = sample_raster(pts, path)
            rows.append({
                "dataset": dataset,
                "year": args.year,
                "raster": name,
                "points_total": int(len(pts)),
                "points_inside_raster_bbox": inside,
                "points_on_valid_yjq_pixels": finite,
            })

    out_csv = args.out_dir / "manas_yjq_overlap_audit.csv"
    with out_csv.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    try:
        (args.out_dir / "manas_yjq_overlap_audit.md").write_text(pd.DataFrame(rows).to_markdown(index=False), encoding="utf-8")
    except Exception:
        pass
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
