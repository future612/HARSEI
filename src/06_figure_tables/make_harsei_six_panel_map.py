from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Polygon, Rectangle
from osgeo import gdal
from PIL import Image


ROOT = Path(r"D:\Codex\260724 小论文")
RASTER_DIR = ROOT / "revise" / "annual_harsei_outputs" / "rasters" / "HARSEI"
OUT_DIR = ROOT / "意见修改" / "图表" / "HARSEI_六期空间分布图"

YEARS = [2000, 2005, 2010, 2015, 2020, 2024]
LETTERS = ["a", "b", "c", "d", "e", "f"]

# Fixed HARSEI grade thresholds used in the revised manuscript.
CLASS_BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.000001]
CLASS_LABELS = ["Worst", "Poor", "Moderate", "Good", "Excellent"]
CLASS_COLORS = ["#cf4b37", "#f4bd1f", "#7ee000", "#22a77b", "#0b2d6c"]


def read_harsei(year: int) -> np.ma.MaskedArray:
    path = RASTER_DIR / f"YJQ_HARSEI_{year}.tif"
    ds = gdal.Open(str(path))
    if ds is None:
        raise FileNotFoundError(path)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype("float32")
    nodata = band.GetNoDataValue()
    mask = ~np.isfinite(arr)
    if nodata is not None:
        mask |= arr == nodata
    mask |= arr < -100
    return np.ma.array(arr, mask=mask)


def classify(values: np.ma.MaskedArray) -> np.ma.MaskedArray:
    data = values.filled(np.nan)
    classes = np.digitize(data, CLASS_BINS[1:-1], right=False) + 1
    classes = classes.astype("float32")
    classes[values.mask] = np.nan
    return np.ma.array(classes, mask=values.mask)


def union_valid_bbox(arrays: dict[int, np.ma.MaskedArray], pad: int = 24) -> tuple[slice, slice]:
    valid = np.zeros(next(iter(arrays.values())).shape, dtype=bool)
    for arr in arrays.values():
        valid |= ~arr.mask
    rows, cols = np.where(valid)
    r0 = max(int(rows.min()) - pad, 0)
    r1 = min(int(rows.max()) + pad + 1, valid.shape[0])
    c0 = max(int(cols.min()) - pad, 0)
    c1 = min(int(cols.max()) + pad + 1, valid.shape[1])
    return slice(r0, r1), slice(c0, c1)


def add_compass(fig: plt.Figure) -> None:
    ax = fig.add_axes([0.890, 0.790, 0.075, 0.125])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)

    for i in range(16):
        angle = np.deg2rad(90 - i * 22.5)
        half_width = np.deg2rad(8)
        radius = 0.95 if i % 2 == 0 else 0.72
        base_radius = 0.15
        p_tip = (radius * np.cos(angle), radius * np.sin(angle))
        p_left = (base_radius * np.cos(angle + np.pi + half_width), base_radius * np.sin(angle + np.pi + half_width))
        p_right = (base_radius * np.cos(angle + np.pi - half_width), base_radius * np.sin(angle + np.pi - half_width))
        face = "black" if i % 2 == 0 else "white"
        ax.add_patch(Polygon([p_tip, p_left, (0, 0), p_right], closed=True, facecolor=face, edgecolor="black", linewidth=0.7))

    circle = plt.Circle((0, 0), 0.97, fill=False, color="black", linewidth=0.8)
    ax.add_patch(circle)
    ax.text(0, 1.17, "N", ha="center", va="center", fontsize=8.5, family="serif")
    ax.text(1.17, 0, "E", ha="center", va="center", fontsize=8.5, family="serif")
    ax.text(0, -1.17, "S", ha="center", va="center", fontsize=8.5, family="serif")
    ax.text(-1.17, 0, "W", ha="center", va="center", fontsize=8.5, family="serif")


def add_scale_bar(fig: plt.Figure) -> None:
    x0, y0 = 0.075, 0.079
    seg_w, h = 0.065, 0.015
    fig.patches.extend(
        [
            Rectangle((x0, y0), seg_w, h, transform=fig.transFigure, facecolor="black", edgecolor="black", linewidth=0.8),
            Rectangle((x0 + seg_w, y0), seg_w, h, transform=fig.transFigure, facecolor="white", edgecolor="black", linewidth=0.8),
        ]
    )
    for label, x in [("0", x0), ("250", x0 + seg_w), ("500", x0 + 2 * seg_w)]:
        fig.text(x, y0 + h + 0.012, label, ha="center", va="bottom", fontsize=10, family="serif")
    fig.text(x0 + seg_w, y0 - 0.017, "km", ha="center", va="top", fontsize=10, family="serif")


def add_legend(fig: plt.Figure) -> None:
    y = 0.072
    box_w, box_h = 0.052, 0.031
    x = 0.255
    gap = 0.145
    for color, label in zip(CLASS_COLORS, CLASS_LABELS):
        fig.patches.append(Rectangle((x, y), box_w, box_h, transform=fig.transFigure, facecolor=color, edgecolor="#666666", linewidth=0.6))
        fig.text(x + box_w + 0.010, y + box_h * 0.5, label, ha="left", va="center", fontsize=12, fontweight="bold", family="serif")
        x += gap


def write_class_area_table(arrays: dict[int, np.ma.MaskedArray], out_csv: Path) -> None:
    lines = ["year,class_code,class_label,pixels,percent"]
    for year, arr in arrays.items():
        cls = classify(arr)
        valid = ~cls.mask
        total = int(valid.sum())
        for code, label in enumerate(CLASS_LABELS, start=1):
            pixels = int(np.sum(cls.data[valid] == code))
            percent = pixels / total * 100 if total else 0
            lines.append(f"{year},{code},{label},{pixels},{percent:.6f}")
    out_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    mpl.rcParams["axes.unicode_minus"] = False

    arrays = {year: read_harsei(year) for year in YEARS}
    crop = union_valid_bbox(arrays)

    cmap = ListedColormap(CLASS_COLORS)
    cmap.set_bad("white")
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)

    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.0), dpi=600)
    fig.patch.set_facecolor("white")

    for ax, year, letter in zip(axes.ravel(), YEARS, LETTERS):
        cls = classify(arrays[year])[crop]
        ax.imshow(cls, cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_axis_off()
        ax.text(
            0.00,
            1.015,
            f"({letter}) {year}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=15,
            fontweight="bold",
            family="serif",
            clip_on=False,
        )

    plt.subplots_adjust(left=0.02, right=0.985, top=0.940, bottom=0.155, wspace=0.005, hspace=0.105)
    add_compass(fig)
    add_scale_bar(fig)
    add_legend(fig)

    png = OUT_DIR / "Fig_HARSEI_spatial_distribution_2000_2024.png"
    tif = OUT_DIR / "Fig_HARSEI_spatial_distribution_2000_2024.tif"
    pdf = OUT_DIR / "Fig_HARSEI_spatial_distribution_2000_2024.pdf"
    csv = OUT_DIR / "Fig_HARSEI_spatial_distribution_class_area.csv"
    note = OUT_DIR / "Fig_HARSEI_spatial_distribution_readme.txt"

    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(pdf, dpi=600, facecolor="white")
    plt.close(fig)

    with Image.open(png) as im:
        im.save(tif, compression="tiff_lzw", dpi=(600, 600))

    write_class_area_table(arrays, csv)
    note.write_text(
        "HARSEI six-period map generated from latest annual rasters in "
        f"{RASTER_DIR}\n"
        "Years: 2000, 2005, 2010, 2015, 2020, 2024\n"
        "Classification thresholds: Worst 0-0.2; Poor 0.2-0.4; "
        "Moderate 0.4-0.6; Good 0.6-0.8; Excellent 0.8-1.0.\n"
        "Output formats: PNG 600 dpi, TIFF LZW 600 dpi, PDF.\n",
        encoding="utf-8",
    )

    print(f"Saved: {png}")
    print(f"Saved: {tif}")
    print(f"Saved: {pdf}")
    print(f"Saved: {csv}")


if __name__ == "__main__":
    gdal.DontUseExceptions()
    main()
