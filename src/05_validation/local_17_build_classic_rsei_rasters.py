#!/usr/bin/env python
"""Build annual classic RSEI rasters for HARSEI-RSEI difference analysis.

The classic RSEI here uses NDVI, WET, NDBSI and LST only. Each component is
normalized with the same fixed 2000-2024 reference ranges used in the HARSEI
revision workflow, then projected onto a pooled PC1 model fitted from the
existing pooled sample table.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import rasterio


RSEI_BANDS = ("NDVI", "WET", "NDBSI", "LST")
NODATA = -9999.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(r"D:\Codex\260724 小论文")
    parser.add_argument("--eco-dir", type=Path, default=root / "revise/gee_downloads/YJQ_HARSEI_annual_inputs_2000_2024")
    parser.add_argument("--ranges", type=Path, default=root / "revise/annual_harsei_outputs/tables/fixed_normalization_ranges.csv")
    parser.add_argument("--pca-sample", type=Path, default=root / "revise/gee_downloads/YJQ_HARSEI_annual_inputs_2000_2024/YJQ_pooled_pca_sample_2000_2024.csv")
    parser.add_argument("--out-dir", type=Path, default=root / "revise/annual_harsei_outputs/rasters/RSEI")
    parser.add_argument("--years", default="2000-2024", help="Year range, e.g. 2000-2024, or comma list.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def parse_years(text: str) -> list[int]:
    if "-" in text:
        a, b = [int(v) for v in text.split("-", 1)]
        return list(range(a, b + 1))
    return [int(v.strip()) for v in text.split(",") if v.strip()]


def read_ranges(path: Path) -> dict[str, tuple[float, float, int]]:
    ranges: dict[str, tuple[float, float, int]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            ranges[row["variable"]] = (float(row["p_low"]), float(row["p_high"]), int(row["direction"]))
    return ranges


def normalize_band(arr: np.ndarray, lo: float, hi: float, direction: int) -> np.ndarray:
    out = (arr.astype(np.float32) - lo) / (hi - lo + 1e-6)
    out = np.clip(out, 0.0, 1.0)
    if direction < 0:
        out = 1.0 - out
    out[~np.isfinite(arr)] = np.nan
    return out.astype(np.float32)


def pca_from_matrix(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mean = np.nanmean(x, axis=0)
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


def fit_pca(sample_csv: Path, ranges: dict[str, tuple[float, float, int]]) -> dict[str, np.ndarray | float | int]:
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
            if ok and all(np.isfinite(vals)):
                rows.append(vals)
    x = np.asarray(rows, dtype=np.float64)
    mean, pc1, eigvals, explained = pca_from_matrix(x)
    scores = (x - mean) @ pc1
    return {
        "mean": mean,
        "pc1": pc1,
        "explained_pc1": float(explained[0]),
        "eigenvalue_pc1": float(eigvals[0]),
        "score_p1": float(np.nanpercentile(scores, 1)),
        "score_p99": float(np.nanpercentile(scores, 99)),
        "sample_size": int(x.shape[0]),
    }


def compute_rsei(stack: np.ndarray, ranges: dict[str, tuple[float, float, int]], pca: dict[str, np.ndarray | float | int]) -> np.ndarray:
    normalized = []
    for i, var in enumerate(RSEI_BANDS):
        lo, hi, direction = ranges[var]
        normalized.append(normalize_band(stack[i], lo, hi, direction))
    x = np.stack(normalized)
    valid = np.all(np.isfinite(x), axis=0)
    flat = x.reshape(x.shape[0], -1).T
    scores = (flat - pca["mean"]) @ pca["pc1"]  # type: ignore[operator]
    rsei = normalize_band(
        scores.reshape(x.shape[1], x.shape[2]),
        float(pca["score_p1"]),
        float(pca["score_p99"]),
        1,
    )
    rsei[~valid] = np.nan
    return rsei.astype(np.float32)


def main() -> None:
    args = parse_args()
    years = parse_years(args.years)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ranges = read_ranges(args.ranges)
    pca = fit_pca(args.pca_sample, ranges)

    summary_path = args.out_dir / "RSEI_pooled_PC1_model.csv"
    with summary_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["variable", "pc1_loading"])
        for var, loading in zip(RSEI_BANDS, pca["pc1"]):  # type: ignore[arg-type]
            writer.writerow([var, float(loading)])
        writer.writerow(["explained_pc1", float(pca["explained_pc1"])])
        writer.writerow(["sample_size", int(pca["sample_size"])])

    for year in years:
        eco_path = args.eco_dir / f"YJQ_ecocomponents_{year}.tif"
        out_path = args.out_dir / f"YJQ_RSEI_{year}.tif"
        if out_path.exists() and not args.overwrite:
            print(f"skip existing {out_path}")
            continue
        if not eco_path.exists():
            print(f"missing {eco_path}")
            continue
        with rasterio.open(eco_path) as src:
            stack = src.read([1, 2, 3, 4]).astype(np.float32)
            profile = src.profile.copy()
            nodata = src.nodata
        if nodata is not None:
            stack[stack == nodata] = np.nan
        stack[stack == NODATA] = np.nan
        rsei = compute_rsei(stack, ranges, pca)
        profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
        out = np.where(np.isfinite(rsei), rsei, NODATA).astype(np.float32)
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(out, 1)
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
