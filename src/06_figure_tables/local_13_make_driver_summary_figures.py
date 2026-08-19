#!/usr/bin/env python
"""Create revised exogenous-driver summary figures for manuscript revision."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXCLUDE = {"LUCC", "HAI", "NDVI", "WET", "NDBSI", "SRSI", "LST", "AWRSEI", "POP", "LIGHT"}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--geodetector-dir", required=True)
    p.add_argument("--rf-dir", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def read_geodetector(path: Path):
    vals = defaultdict(list)
    for fpath in sorted(path.glob("GD_factor_*.csv")):
        with fpath.open("r", encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                if row["Factor"] not in EXCLUDE:
                    vals[row["Factor"]].append(float(row["q"]))
    return sorted(((k, statistics.mean(v)) for k, v in vals.items()), key=lambda x: x[1], reverse=True)


def read_rf(path: Path):
    vals = defaultdict(list)
    with (path / "importance_timeseries.csv").open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            vals[row["Variable"]].append(float(row["CompositeScore"]))
    return sorted(((k, statistics.mean(v)) for k, v in vals.items()), key=lambda x: x[1], reverse=True)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gd = read_geodetector(Path(args.geodetector_dir))[:8]
    rf = read_rf(Path(args.rf_dir))[:8]

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), constrained_layout=True)
    for ax, data, title, xlabel, color in [
        (axes[0], gd, "GeoDetector", "Mean q", "#4daf4a"),
        (axes[1], rf, "RF-SHAP", "Mean CompositeScore", "#377eb8"),
    ]:
        names = [x[0] for x in data][::-1]
        values = [x[1] for x in data][::-1]
        ax.barh(names, values, color=color, edgecolor="#444444", linewidth=0.5)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        for i, value in enumerate(values):
            ax.text(value + max(values) * 0.01, i, f"{value:.3f}", va="center", fontsize=8)
    fig.suptitle("Exogenous driver importance after excluding HARSEI components", fontsize=13, fontweight="bold")
    fig.savefig(out_dir / "fig_revised_exogenous_driver_summary.png", dpi=300)
    plt.close(fig)

    with (out_dir / "exogenous_driver_summary.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "variable", "mean_value"])
        for k, v in gd:
            w.writerow(["GeoDetector_q", k, v])
        for k, v in rf:
            w.writerow(["RF_SHAP_CompositeScore", k, v])
    print(out_dir / "fig_revised_exogenous_driver_summary.png")


if __name__ == "__main__":
    main()
