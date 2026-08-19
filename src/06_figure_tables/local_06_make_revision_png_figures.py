#!/usr/bin/env python
"""
Create PNG versions of the manuscript revision figures from HARSEI summary
tables. This script is intentionally separate from local_05 so the SVG/text
package remains usable with only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def save_trajectory(annual: list[dict[str, str]], fig_dir: Path, dpi: int) -> None:
    years = [int(row["year"]) for row in annual]
    awrsei = [f(row, "AWRSEI_mean") for row in annual]
    hai_reverse = [f(row, "HAI_reverse_mean") for row in annual]
    harsei = [f(row, "HARSEI_mean") for row in annual]

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    ax.plot(years, awrsei, color="#1b9e77", lw=2.0, marker="o", ms=3.2, label="AWRSEI")
    ax.plot(years, hai_reverse, color="#d95f02", lw=2.0, marker="s", ms=3.0, label="HAI reverse")
    ax.plot(years, harsei, color="#377eb8", lw=2.5, marker="^", ms=3.2, label="HARSEI")
    ax.set_xlim(min(years), max(years))
    ax.set_ylim(0, 1)
    ax.set_xticks([2000, 2005, 2010, 2015, 2020, 2024])
    ax.set_xlabel("Year")
    ax.set_ylabel("Index value")
    ax.set_title("Annual AWRSEI, HAI-reverse, and HARSEI Means (2000-2024)")
    ax.grid(True, color="#dddddd", linewidth=0.7, alpha=0.8)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    fig.savefig(fig_dir / "fig_annual_harsei_trajectory.png", dpi=dpi)
    plt.close(fig)


def save_sensitivity(
    snow_rows: list[dict[str, str]],
    threshold_rows: list[dict[str, str]],
    fig_dir: Path,
    dpi: int,
) -> None:
    snow = {}
    for row in snow_rows:
        if abs(f(row, "ndwi_threshold") - 0.2) < 1e-9:
            snow[int(f(row, "start_month"))] = f(row, "mean_snow_fraction_high_elevation")

    thresh = {}
    for row in threshold_rows:
        t = f(row, "ndwi_threshold")
        thresh.setdefault(t, []).append(f(row, "HARSEI_mean"))
    threshold_mean = {t: sum(vals) / len(vals) for t, vals in thresh.items()}

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2), constrained_layout=True)

    labels = ["Apr-Sep", "May-Sep"]
    vals = [snow.get(4, 0.0), snow.get(5, 0.0)]
    axes[0].bar(labels, vals, color=["#8da0cb", "#66c2a5"], edgecolor="#444444", linewidth=0.7)
    axes[0].set_ylim(0, 0.25)
    axes[0].set_ylabel("Snow/ice frequency")
    axes[0].set_title("High-elevation snow frequency")
    axes[0].grid(True, axis="y", color="#dddddd", linewidth=0.7, alpha=0.8)
    for i, val in enumerate(vals):
        axes[0].text(i, val + 0.006, f"{val:.4f}", ha="center", va="bottom", fontsize=9)

    thrs = sorted(threshold_mean)
    means = [threshold_mean[t] for t in thrs]
    axes[1].bar([f">{t:.1f}" for t in thrs], means, color="#fc8d62", edgecolor="#444444", linewidth=0.7)
    axes[1].set_ylim(0.5250, 0.5255)
    axes[1].set_ylabel("Mean HARSEI")
    axes[1].set_xlabel("NDWI water threshold")
    axes[1].set_title("Threshold sensitivity (zoomed y-axis)")
    axes[1].grid(True, axis="y", color="#dddddd", linewidth=0.7, alpha=0.8)
    for i, val in enumerate(means):
        axes[1].text(i, val + 0.00001, f"{val:.6f}", ha="center", va="bottom", fontsize=8)

    fig.suptitle("Water Mask and Snow Sensitivity", fontsize=13, fontweight="bold")
    fig.savefig(fig_dir / "fig_water_snow_sensitivity.png", dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    tables_dir = Path(args.tables_dir)
    fig_dir = Path(args.out_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    annual = read_csv(tables_dir / "annual_index_summary.csv")
    snow = read_csv(tables_dir / "water_snow_sensitivity_summary.csv")
    threshold = read_csv(tables_dir / "water_threshold_index_sensitivity_current_mask.csv")

    save_trajectory(annual, fig_dir, args.dpi)
    save_sensitivity(snow, threshold, fig_dir, args.dpi)


if __name__ == "__main__":
    main()
