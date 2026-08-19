from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Polygon, Rectangle
from osgeo import gdal
from PIL import Image


ROOT = Path(r"D:\Codex\260724 小论文")
RASTER_DIR = ROOT / "revise" / "annual_harsei_outputs" / "rasters" / "HARSEI"
OUT_DIR = ROOT / "意见修改" / "图表" / "HARSEI_转移矩阵_空间变化"
FIG_DIR = OUT_DIR / "figures"
RASTER_OUT = OUT_DIR / "rasters"
TABLE_DIR = OUT_DIR / "tables"

YEARS = [2000, 2005, 2010, 2015, 2020, 2024]
ADJACENT_PERIODS = [(2000, 2005), (2005, 2010), (2010, 2015), (2015, 2020), (2020, 2024)]
MAP_PERIODS = ADJACENT_PERIODS + [(2000, 2024)]
LETTERS = ["a", "b", "c", "d", "e", "f"]

GRADE_BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.000001]
GRADE_LABELS = ["Worst", "Poor", "Moderate", "Good", "Excellent"]
GRADE_LABELS_CN = ["Worst", "Poor", "Moderate", "Good", "Excellent"]
GRADE_COLORS = ["#cf4b37", "#f4bd1f", "#7ee000", "#22a77b", "#0b2d6c"]

CHANGE_LABELS = [
    "Significant deterioration",
    "Slight deterioration",
    "Unchanged",
    "Slight improvement",
    "Significant improvement",
]
CHANGE_RULES = ["delta_grade <= -2", "delta_grade = -1", "delta_grade = 0", "delta_grade = +1", "delta_grade >= +2"]
CHANGE_COLORS = ["#cf4b37", "#f4bd1f", "#7ee000", "#22a77b", "#0b2d6c"]


def setup_style() -> None:
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42


def read_raster(year: int) -> tuple[np.ma.MaskedArray, gdal.Dataset]:
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
    return np.ma.array(arr, mask=mask), ds


def classify_harsei(arr: np.ma.MaskedArray) -> np.ndarray:
    data = arr.filled(np.nan)
    classes = np.digitize(data, GRADE_BINS[1:-1], right=False) + 1
    classes = classes.astype("uint8")
    classes[arr.mask] = 0
    return classes


def change_class(start_grade: np.ndarray, end_grade: np.ndarray) -> np.ndarray:
    valid = (start_grade > 0) & (end_grade > 0)
    delta = end_grade.astype("int16") - start_grade.astype("int16")
    out = np.zeros(start_grade.shape, dtype="uint8")
    out[valid & (delta <= -2)] = 1
    out[valid & (delta == -1)] = 2
    out[valid & (delta == 0)] = 3
    out[valid & (delta == 1)] = 4
    out[valid & (delta >= 2)] = 5
    return out


def transition_matrix(start_grade: np.ndarray, end_grade: np.ndarray) -> tuple[np.ndarray, int]:
    valid = (start_grade > 0) & (end_grade > 0)
    mat = np.zeros((5, 5), dtype=np.int64)
    for i in range(1, 6):
        for j in range(1, 6):
            mat[i - 1, j - 1] = int(np.sum(valid & (start_grade == i) & (end_grade == j)))
    return mat, int(valid.sum())


def write_raster_like(ref_ds: gdal.Dataset, out_path: Path, arr: np.ndarray, dtype: int, nodata: float | int, options: list[str] | None = None) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    options = options or ["COMPRESS=LZW", "TILED=YES"]
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(str(out_path), ref_ds.RasterXSize, ref_ds.RasterYSize, 1, dtype, options)
    ds.SetGeoTransform(ref_ds.GetGeoTransform())
    ds.SetProjection(ref_ds.GetProjection())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(arr)
    band.FlushCache()
    ds.FlushCache()
    ds = None


def save_figure_bundle(fig: plt.Figure, stem: str) -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": FIG_DIR / f"{stem}.png",
        "tif": FIG_DIR / f"{stem}.tif",
        "pdf": FIG_DIR / f"{stem}.pdf",
    }
    fig.savefig(paths["png"], dpi=600, facecolor="white")
    fig.savefig(paths["pdf"], dpi=600, facecolor="white")
    with Image.open(paths["png"]) as im:
        im.save(paths["tif"], compression="tiff_lzw", dpi=(600, 600))
    return paths


def union_valid_bbox(arrays: list[np.ndarray], pad: int = 24) -> tuple[slice, slice]:
    valid = np.zeros(arrays[0].shape, dtype=bool)
    for arr in arrays:
        valid |= arr > 0
    rows, cols = np.where(valid)
    r0 = max(int(rows.min()) - pad, 0)
    r1 = min(int(rows.max()) + pad + 1, valid.shape[0])
    c0 = max(int(cols.min()) - pad, 0)
    c1 = min(int(cols.max()) + pad + 1, valid.shape[1])
    return slice(r0, r1), slice(c0, c1)


def add_compass(fig: plt.Figure, rect: list[float] = [0.890, 0.790, 0.075, 0.125]) -> None:
    ax = fig.add_axes(rect)
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
    ax.add_patch(plt.Circle((0, 0), 0.97, fill=False, color="black", linewidth=0.8))
    ax.text(0, 1.17, "N", ha="center", va="center", fontsize=8.5)
    ax.text(1.17, 0, "E", ha="center", va="center", fontsize=8.5)
    ax.text(0, -1.17, "S", ha="center", va="center", fontsize=8.5)
    ax.text(-1.17, 0, "W", ha="center", va="center", fontsize=8.5)


def add_scale_bar(fig: plt.Figure, x0: float = 0.785, y0: float = 0.145) -> None:
    seg_w, h = 0.046, 0.012
    fig.patches.extend(
        [
            Rectangle((x0, y0), seg_w, h, transform=fig.transFigure, facecolor="black", edgecolor="black", linewidth=0.8),
            Rectangle((x0 + seg_w, y0), seg_w, h, transform=fig.transFigure, facecolor="white", edgecolor="black", linewidth=0.8),
        ]
    )
    for label, x in [("0", x0), ("250", x0 + seg_w), ("500", x0 + 2 * seg_w)]:
        fig.text(x, y0 + h + 0.010, label, ha="center", va="bottom", fontsize=9)
    fig.text(x0 + 2 * seg_w + 0.008, y0 + h * 0.5, "km", ha="left", va="center", fontsize=9)


def add_change_legend(fig: plt.Figure) -> None:
    y = 0.055
    box_w, box_h = 0.045, 0.026
    starts = [0.075, 0.255, 0.425, 0.595, 0.775]
    for x, color, label in zip(starts, CHANGE_COLORS, CHANGE_LABELS):
        fig.patches.append(Rectangle((x, y), box_w, box_h, transform=fig.transFigure, facecolor=color, edgecolor="#666666", linewidth=0.6))
        text = label.replace(" deterioration", "\ndeterioration").replace(" improvement", "\nimprovement")
        fig.text(x + box_w + 0.008, y + box_h * 0.5, text, ha="left", va="center", fontsize=9.5, fontweight="bold")


def plot_spatial_change_maps(change_maps: dict[tuple[int, int], np.ndarray]) -> None:
    cmap = ListedColormap(CHANGE_COLORS)
    cmap.set_bad("white")
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    crop = union_valid_bbox(list(change_maps.values()))
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.0), dpi=600)
    fig.patch.set_facecolor("white")
    for ax, period, letter in zip(axes.ravel(), MAP_PERIODS, LETTERS):
        arr = change_maps[period][crop].astype("float32")
        arr[arr == 0] = np.nan
        ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_axis_off()
        ax.text(
            0.00,
            1.015,
            f"({letter}) {period[0]}-{period[1]}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
        )
    plt.subplots_adjust(left=0.02, right=0.985, top=0.940, bottom=0.170, wspace=0.005, hspace=0.105)
    add_compass(fig)
    add_scale_bar(fig)
    add_change_legend(fig)
    save_figure_bundle(fig, "Fig_HARSEI_spatial_transition_change_2000_2024")
    plt.close(fig)


def sample_values(values: np.ndarray, max_n: int = 140_000, seed: int = 20260810) -> np.ndarray:
    vals = values[np.isfinite(values)]
    if vals.size <= max_n:
        return vals
    rng = np.random.default_rng(seed)
    idx = rng.choice(vals.size, size=max_n, replace=False)
    return vals[idx]


def plot_violin_by_year(harsei_arrays: dict[int, np.ma.MaskedArray]) -> None:
    data = [sample_values(harsei_arrays[y].compressed(), seed=20260810 + y) for y in YEARS]
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=600)
    parts = ax.violinplot(data, positions=np.arange(len(YEARS)), widths=0.78, showmeans=True, showmedians=True, showextrema=False)
    for body, color in zip(parts["bodies"], ["#2f80ed"] * len(YEARS)):
        body.set_facecolor(color)
        body.set_alpha(0.45)
        body.set_edgecolor("#1b4f8a")
        body.set_linewidth(0.8)
    for key in ["cmeans", "cmedians"]:
        parts[key].set_color("#111111")
        parts[key].set_linewidth(1.0)
    ax.set_xticks(np.arange(len(YEARS)))
    ax.set_xticklabels([str(y) for y in YEARS], fontsize=11, fontweight="bold")
    ax.set_ylabel("HARSEI", fontsize=12, fontweight="bold")
    ax.set_ylim(0.2, 1.02)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.set_title("Distribution of HARSEI values in target years", fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure_bundle(fig, "Fig_HARSEI_violin_distribution_by_year")
    plt.close(fig)


def plot_violin_by_period(harsei_arrays: dict[int, np.ma.MaskedArray]) -> None:
    data = []
    labels = []
    for idx, (start, end) in enumerate(MAP_PERIODS):
        a = harsei_arrays[start]
        b = harsei_arrays[end]
        valid = (~a.mask) & (~b.mask)
        diff = np.full(a.shape, np.nan, dtype="float32")
        diff[valid] = b.data[valid] - a.data[valid]
        data.append(sample_values(diff[valid], seed=20260810 + idx))
        labels.append(f"{start}-{end}")
    fig, ax = plt.subplots(figsize=(9.4, 5.2), dpi=600)
    parts = ax.violinplot(data, positions=np.arange(len(labels)), widths=0.78, showmeans=True, showmedians=True, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("#22a77b")
        body.set_alpha(0.45)
        body.set_edgecolor("#0c6b51")
        body.set_linewidth(0.8)
    for key in ["cmeans", "cmedians"]:
        parts[key].set_color("#111111")
        parts[key].set_linewidth(1.0)
    ax.axhline(0, color="#555555", linestyle="--", linewidth=1.0)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=10.5, fontweight="bold")
    ax.set_ylabel("ΔHARSEI (end - start)", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
    ax.set_title("Distribution of HARSEI changes by period", fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure_bundle(fig, "Fig_HARSEI_change_violin_by_period")
    plt.close(fig)


def draw_ribbon(ax: plt.Axes, x0: float, x1: float, y0_low: float, y0_high: float, y1_low: float, y1_high: float, color: str, alpha: float = 0.34) -> None:
    dx = x1 - x0
    c0 = x0 + dx * 0.45
    c1 = x1 - dx * 0.45
    verts = [
        (x0, y0_high),
        (c0, y0_high),
        (c1, y1_high),
        (x1, y1_high),
        (x1, y1_low),
        (c1, y1_low),
        (c0, y0_low),
        (x0, y0_low),
        (x0, y0_high),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def plot_sankey(grade_arrays: dict[int, np.ndarray]) -> None:
    common = np.ones(next(iter(grade_arrays.values())).shape, dtype=bool)
    for year in YEARS:
        common &= grade_arrays[year] > 0
    total = int(common.sum())
    series = {year: grade_arrays[year][common] for year in YEARS}
    usable = 0.72
    bottom = 0.145
    gap = 0.012
    node_w = 0.022
    xs = np.linspace(0.07, 0.93, len(YEARS))
    node_pos: dict[tuple[int, int], tuple[float, float]] = {}

    fig, ax = plt.subplots(figsize=(11.0, 5.8), dpi=600)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    for stage_idx, year in enumerate(YEARS):
        y0 = bottom
        scale = usable - gap * 4
        for grade in range(1, 6):
            count = int(np.sum(series[year] == grade))
            h = count / total * scale
            node_pos[(year, grade)] = (y0, y0 + h)
            y0 += h + gap

    for stage_idx, (start, end) in enumerate(ADJACENT_PERIODS):
        src_offsets = {grade: node_pos[(start, grade)][0] for grade in range(1, 6)}
        dst_offsets = {grade: node_pos[(end, grade)][0] for grade in range(1, 6)}
        scale = usable - gap * 4
        for source_grade in range(1, 6):
            for target_grade in range(1, 6):
                count = int(np.sum((series[start] == source_grade) & (series[end] == target_grade)))
                if count == 0:
                    continue
                h = count / total * scale
                y0_low = src_offsets[source_grade]
                y0_high = y0_low + h
                y1_low = dst_offsets[target_grade]
                y1_high = y1_low + h
                src_offsets[source_grade] = y0_high
                dst_offsets[target_grade] = y1_high
                draw_ribbon(ax, xs[stage_idx] + node_w / 2, xs[stage_idx + 1] - node_w / 2, y0_low, y0_high, y1_low, y1_high, GRADE_COLORS[target_grade - 1])

    for year, x in zip(YEARS, xs):
        for grade in range(1, 6):
            y0, y1 = node_pos[(year, grade)]
            ax.add_patch(Rectangle((x - node_w / 2, y0), node_w, y1 - y0, facecolor=GRADE_COLORS[grade - 1], edgecolor="white", linewidth=0.8, zorder=3))
            if y1 - y0 > 0.045:
                pct = (y1 - y0) / (usable - gap * 4) * 100
                ax.text(x, (y0 + y1) / 2, f"{pct:.1f}%", ha="center", va="center", fontsize=7.5, color="white" if grade in [1, 5] else "black", zorder=4)
        ax.text(x, 0.905, str(year), ha="center", va="bottom", fontsize=12, fontweight="bold")

    legend_x = [0.135, 0.290, 0.445, 0.600, 0.755]
    legend_y = 0.055
    for x, label, color in zip(legend_x, GRADE_LABELS, GRADE_COLORS):
        ax.add_patch(Rectangle((x, legend_y - 0.014), 0.028, 0.028, facecolor=color, edgecolor="#666666", linewidth=0.5))
        ax.text(x + 0.036, legend_y, label, ha="left", va="center", fontsize=10.0, fontweight="bold")

    ax.set_title("HARSEI grade transitions across target years", fontsize=14, fontweight="bold", pad=12)
    fig.tight_layout()
    save_figure_bundle(fig, "Fig_HARSEI_grade_transition_sankey_2000_2024")
    plt.close(fig)


def plot_transition_heatmaps(matrices: dict[tuple[int, int], np.ndarray], totals: dict[tuple[int, int], int]) -> None:
    percent_mats = {period: matrices[period] / totals[period] * 100 for period in MAP_PERIODS}
    vmax = max(float(np.max(mat)) for mat in percent_mats.values())
    fig = plt.figure(figsize=(11.4, 7.0), dpi=600)
    gs = fig.add_gridspec(
        2,
        4,
        width_ratios=[1, 1, 1, 0.055],
        left=0.065,
        right=0.955,
        top=0.925,
        bottom=0.105,
        wspace=0.42,
        hspace=0.46,
    )
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)])
    cax = fig.add_subplot(gs[:, 3])
    fig.patch.set_facecolor("white")
    for ax, period, letter in zip(axes.ravel(), MAP_PERIODS, LETTERS):
        mat = percent_mats[period]
        im = ax.imshow(mat, cmap="YlGnBu", vmin=0, vmax=vmax)
        ax.set_title(f"({letter}) {period[0]}-{period[1]}", fontsize=12, fontweight="bold")
        ax.set_xticks(np.arange(5))
        ax.set_yticks(np.arange(5))
        ax.set_xticklabels(GRADE_LABELS, rotation=35, ha="right", fontsize=8)
        ax.set_yticklabels(GRADE_LABELS, fontsize=8)
        ax.set_xlabel("To grade", fontsize=9, fontweight="bold")
        ax.set_ylabel("From grade", fontsize=9, fontweight="bold")
        for i in range(5):
            for j in range(5):
                val = mat[i, j]
                if val == 0:
                    text = "0"
                elif val < 0.05:
                    text = "<0.1"
                else:
                    text = f"{val:.1f}"
                color = "white" if val > vmax * 0.45 else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=7.2, color=color)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Percent of valid pixels (%)", fontsize=10, fontweight="bold")
    save_figure_bundle(fig, "Fig_HARSEI_transition_matrix_heatmaps")
    plt.close(fig)


def summarize_distribution(vals: np.ndarray) -> list[float]:
    return [
        int(vals.size),
        float(np.nanmin(vals)),
        float(np.nanpercentile(vals, 25)),
        float(np.nanmedian(vals)),
        float(np.nanmean(vals)),
        float(np.nanpercentile(vals, 75)),
        float(np.nanmax(vals)),
    ]


def write_tables(
    harsei_arrays: dict[int, np.ma.MaskedArray],
    grade_arrays: dict[int, np.ndarray],
    change_maps: dict[tuple[int, int], np.ndarray],
    matrices: dict[tuple[int, int], np.ndarray],
    totals: dict[tuple[int, int], int],
) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    with (TABLE_DIR / "HARSEI_transition_matrix_long.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["period", "from_code", "from_grade", "to_code", "to_grade", "pixels", "percent_of_period_valid", "percent_of_from_grade"])
        for period in MAP_PERIODS:
            mat = matrices[period]
            row_sums = mat.sum(axis=1)
            for i in range(5):
                for j in range(5):
                    pixels = int(mat[i, j])
                    pct_total = pixels / totals[period] * 100 if totals[period] else 0
                    pct_from = pixels / row_sums[i] * 100 if row_sums[i] else 0
                    writer.writerow([f"{period[0]}-{period[1]}", i + 1, GRADE_LABELS[i], j + 1, GRADE_LABELS[j], pixels, f"{pct_total:.6f}", f"{pct_from:.6f}"])

    with (TABLE_DIR / "HARSEI_change_class_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["period", "change_code", "change_class", "rule", "pixels", "percent_of_valid_pixels"])
        for period in MAP_PERIODS:
            arr = change_maps[period]
            valid = arr > 0
            total = int(valid.sum())
            for code, (label, rule) in enumerate(zip(CHANGE_LABELS, CHANGE_RULES), start=1):
                pixels = int(np.sum(arr == code))
                pct = pixels / total * 100 if total else 0
                writer.writerow([f"{period[0]}-{period[1]}", code, label, rule, pixels, f"{pct:.6f}"])

    with (TABLE_DIR / "HARSEI_distribution_summary_by_year.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "n", "min", "q25", "median", "mean", "q75", "max"])
        for year in YEARS:
            vals = harsei_arrays[year].compressed()
            writer.writerow([year, *[f"{x:.6f}" if isinstance(x, float) else x for x in summarize_distribution(vals)]])

    with (TABLE_DIR / "HARSEI_change_distribution_summary_by_period.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["period", "n", "min", "q25", "median", "mean", "q75", "max"])
        for period in MAP_PERIODS:
            a = harsei_arrays[period[0]]
            b = harsei_arrays[period[1]]
            valid = (~a.mask) & (~b.mask)
            vals = b.data[valid] - a.data[valid]
            writer.writerow([f"{period[0]}-{period[1]}", *[f"{x:.6f}" if isinstance(x, float) else x for x in summarize_distribution(vals)]])

    common = np.ones(next(iter(grade_arrays.values())).shape, dtype=bool)
    for year in YEARS:
        common &= grade_arrays[year] > 0
    common_total = int(common.sum())
    with (TABLE_DIR / "HARSEI_sankey_links_common_valid.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["source_year", "source_code", "source_grade", "target_year", "target_code", "target_grade", "pixels", "percent_of_common_valid"])
        for start, end in ADJACENT_PERIODS:
            sg = grade_arrays[start][common]
            eg = grade_arrays[end][common]
            for i in range(1, 6):
                for j in range(1, 6):
                    pixels = int(np.sum((sg == i) & (eg == j)))
                    pct = pixels / common_total * 100 if common_total else 0
                    writer.writerow([start, i, GRADE_LABELS[i - 1], end, j, GRADE_LABELS[j - 1], pixels, f"{pct:.6f}"])

    readme = OUT_DIR / "README_HARSEI_transition_outputs.txt"
    readme.write_text(
        "HARSEI transition outputs generated from latest annual rasters in:\n"
        f"{RASTER_DIR}\n\n"
        "HARSEI grade thresholds:\n"
        "1 Worst: 0-0.2; 2 Poor: 0.2-0.4; 3 Moderate: 0.4-0.6; 4 Good: 0.6-0.8; 5 Excellent: 0.8-1.0.\n\n"
        "Spatial transition classes are based on grade difference (end grade - start grade):\n"
        "1 Significant deterioration: <= -2; 2 Slight deterioration: -1; 3 Unchanged: 0; "
        "4 Slight improvement: +1; 5 Significant improvement: >= +2.\n\n"
        "Adjacent transition periods: 2000-2005, 2005-2010, 2010-2015, 2015-2020 and 2020-2024. "
        "The 2000-2024 total change is also included to match the reference layout.\n",
        encoding="utf-8",
    )


def main() -> None:
    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RASTER_OUT.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    harsei_arrays: dict[int, np.ma.MaskedArray] = {}
    datasets: dict[int, gdal.Dataset] = {}
    grade_arrays: dict[int, np.ndarray] = {}
    for year in YEARS:
        arr, ds = read_raster(year)
        harsei_arrays[year] = arr
        datasets[year] = ds
        grade = classify_harsei(arr)
        grade_arrays[year] = grade
        write_raster_like(ds, RASTER_OUT / f"YJQ_HARSEI_grade_{year}.tif", grade, gdal.GDT_Byte, 0)

    change_maps: dict[tuple[int, int], np.ndarray] = {}
    matrices: dict[tuple[int, int], np.ndarray] = {}
    totals: dict[tuple[int, int], int] = {}

    for start, end in MAP_PERIODS:
        sg = grade_arrays[start]
        eg = grade_arrays[end]
        chg = change_class(sg, eg)
        mat, total = transition_matrix(sg, eg)
        change_maps[(start, end)] = chg
        matrices[(start, end)] = mat
        totals[(start, end)] = total
        write_raster_like(datasets[start], RASTER_OUT / f"YJQ_HARSEI_transition_change_class_{start}_{end}.tif", chg, gdal.GDT_Byte, 0)

        a = harsei_arrays[start]
        b = harsei_arrays[end]
        valid = (~a.mask) & (~b.mask)
        diff = np.full(a.shape, -9999.0, dtype="float32")
        diff[valid] = b.data[valid] - a.data[valid]
        write_raster_like(datasets[start], RASTER_OUT / f"YJQ_HARSEI_delta_{start}_{end}.tif", diff, gdal.GDT_Float32, -9999.0)

    plot_spatial_change_maps(change_maps)
    plot_violin_by_year(harsei_arrays)
    plot_violin_by_period(harsei_arrays)
    plot_sankey(grade_arrays)
    plot_transition_heatmaps(matrices, totals)
    write_tables(harsei_arrays, grade_arrays, change_maps, matrices, totals)

    print(f"Outputs saved to: {OUT_DIR}")
    print(f"Figures: {FIG_DIR}")
    print(f"Rasters: {RASTER_OUT}")
    print(f"Tables: {TABLE_DIR}")


if __name__ == "__main__":
    gdal.DontUseExceptions()
    main()
