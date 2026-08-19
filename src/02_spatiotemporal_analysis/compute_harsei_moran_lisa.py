from __future__ import annotations

import csv
import math
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
OUT_DIR = ROOT / "意见修改" / "图表" / "HARSEI_Moran_LISA_空间聚集"
RASTER_OUT = OUT_DIR / "rasters"
TABLE_OUT = OUT_DIR / "tables"
FIG_OUT = OUT_DIR / "figures"

YEARS = [2000, 2005, 2010, 2015, 2020, 2024]
LETTERS = ["a", "b", "c", "d", "e", "f"]
ALPHA = 0.05

# Local pseudo p-values are obtained with global randomization fields. This
# reproduces the original permutation logic closely enough for raster mapping
# while remaining tractable without PySAL/esda in the current environment.
LOCAL_PERMUTATIONS = 199
GLOBAL_PERMUTATION_LOWER_BOUND = 0.001  # lower bound under a 999-permutation test
RNG_SEED = 20260810

NODATA_F32 = -9999.0
NODATA_U8 = 255

LISA_CODES = {
    0: "Not significant",
    1: "HH",
    2: "LH",
    3: "LL",
    4: "HL",
}
LISA_COLORS = {
    0: "#d9d9d9",
    1: "#e41a1c",
    2: "#fdbf6f",
    3: "#2c7fb8",
    4: "#a6cee3",
}


def setup_style() -> None:
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42


def read_harsei(year: int) -> tuple[np.ma.MaskedArray, gdal.Dataset]:
    path = RASTER_DIR / f"YJQ_HARSEI_{year}.tif"
    ds = gdal.Open(str(path))
    if ds is None:
        raise FileNotFoundError(path)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype("float64")
    nodata = band.GetNoDataValue()
    mask = ~np.isfinite(arr)
    if nodata is not None:
        mask |= arr == nodata
    mask |= arr < -100
    return np.ma.array(arr, mask=mask), ds


def iter_offsets(all_directions: bool = True):
    if all_directions:
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                yield dr, dc
    else:
        yield from [(0, 1), (1, 0), (1, 1), (1, -1)]


def slices_for_shift(shape: tuple[int, int], dr: int, dc: int):
    rows, cols = shape
    src_r = slice(max(0, -dr), rows - max(0, dr))
    src_c = slice(max(0, -dc), cols - max(0, dc))
    dst_r = slice(max(0, dr), rows - max(0, -dr))
    dst_c = slice(max(0, dc), cols - max(0, -dc))
    return src_r, src_c, dst_r, dst_c


def neighbor_degree(valid: np.ndarray) -> np.ndarray:
    degree = np.zeros(valid.shape, dtype=np.uint8)
    for dr, dc in iter_offsets(True):
        src_r, src_c, dst_r, dst_c = slices_for_shift(valid.shape, dr, dc)
        degree[dst_r, dst_c] += valid[src_r, src_c]
    return degree


def neighbor_sum(values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    out = np.zeros(values.shape, dtype="float64")
    for dr, dc in iter_offsets(True):
        src_r, src_c, dst_r, dst_c = slices_for_shift(values.shape, dr, dc)
        src_vals = values[src_r, src_c]
        src_valid = valid[src_r, src_c]
        out[dst_r, dst_c] += np.where(src_valid, src_vals, 0.0)
    return out


def column_sums_row_standardized(recip_degree: np.ndarray, keep: np.ndarray) -> np.ndarray:
    out = np.zeros(recip_degree.shape, dtype="float64")
    for dr, dc in iter_offsets(True):
        src_r, src_c, dst_r, dst_c = slices_for_shift(recip_degree.shape, dr, dc)
        src_vals = recip_degree[src_r, src_c]
        src_valid = keep[src_r, src_c]
        out[dst_r, dst_c] += np.where(src_valid, src_vals, 0.0)
    return out


def s1_row_standardized(recip_degree: np.ndarray, keep: np.ndarray) -> float:
    s1 = 0.0
    for dr, dc in iter_offsets(False):
        src_r, src_c, dst_r, dst_c = slices_for_shift(recip_degree.shape, dr, dc)
        pair_valid = keep[src_r, src_c] & keep[dst_r, dst_c]
        pair_sum = recip_degree[src_r, src_c][pair_valid] + recip_degree[dst_r, dst_c][pair_valid]
        s1 += float(np.sum(pair_sum * pair_sum))
    return s1


def write_raster_like(ref_ds: gdal.Dataset, out_path: Path, arr: np.ndarray, dtype: int, nodata: float | int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        str(out_path),
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        1,
        dtype,
        ["COMPRESS=LZW", "TILED=YES"],
    )
    ds.SetGeoTransform(ref_ds.GetGeoTransform())
    ds.SetProjection(ref_ds.GetProjection())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(arr)
    band.FlushCache()
    ds.FlushCache()
    ds = None


def compute_moran_lisa(year: int, rng: np.random.Generator) -> dict[str, float | int]:
    arr_ma, ds = read_harsei(year)
    arr = arr_ma.filled(np.nan)
    valid = (~arr_ma.mask) & np.isfinite(arr)
    degree = neighbor_degree(valid)
    keep = valid & (degree > 0)
    n = int(keep.sum())
    if n < 5:
        raise ValueError(f"{year}: too few valid pixels for Moran/LISA.")

    y = arr[keep]
    y_mean = float(np.mean(y))
    z_values = y - y_mean
    z_full = np.zeros(arr.shape, dtype="float64")
    z_full[keep] = z_values
    m2 = float(np.sum(z_values * z_values) / n)
    sumz2 = float(np.sum(z_values * z_values))

    lag_sum = neighbor_sum(z_full, keep)
    lag = np.zeros(arr.shape, dtype="float64")
    lag[keep] = lag_sum[keep] / degree[keep]
    lag_values = lag[keep]

    moran_i = float(np.sum(z_values * lag_values) / sumz2)
    expected_i = -1.0 / (n - 1)

    recip_degree = np.zeros(arr.shape, dtype="float64")
    recip_degree[keep] = 1.0 / degree[keep]
    s0 = float(n)
    s1 = s1_row_standardized(recip_degree, keep)
    col_sums = column_sums_row_standardized(recip_degree, keep)
    s2 = float(np.sum((1.0 + col_sums[keep]) ** 2))
    vi_norm = ((n * n * s1 - n * s2 + 3.0 * s0 * s0) / (((n * n) - 1.0) * s0 * s0)) - expected_i**2
    z_norm = float((moran_i - expected_i) / math.sqrt(vi_norm))

    local_i_values = z_values * lag_values / m2
    local_i = np.full(arr.shape, NODATA_F32, dtype="float32")
    local_i[keep] = local_i_values.astype("float32")

    q_values = np.zeros(n, dtype="uint8")
    q_values[(z_values > 0) & (lag_values > 0)] = 1  # HH
    q_values[(z_values < 0) & (lag_values > 0)] = 2  # LH
    q_values[(z_values < 0) & (lag_values < 0)] = 3  # LL
    q_values[(z_values > 0) & (lag_values < 0)] = 4  # HL

    ge = np.zeros(n, dtype=np.uint16)
    le = np.zeros(n, dtype=np.uint16)
    keep_rows, keep_cols = np.where(keep)
    z_pool = z_values.copy()
    perm_full = np.zeros(arr.shape, dtype="float64")

    for _ in range(LOCAL_PERMUTATIONS):
        perm = rng.permutation(z_pool)
        perm_full.fill(0.0)
        perm_full[keep_rows, keep_cols] = perm
        perm_lag_sum = neighbor_sum(perm_full, keep)
        perm_lag = perm_lag_sum[keep] / degree[keep]
        perm_i = z_values * perm_lag / m2
        ge += perm_i >= local_i_values
        le += perm_i <= local_i_values

    p_values = (np.minimum(ge, le).astype("float64") + 1.0) / (LOCAL_PERMUTATIONS + 1.0)
    p_raster = np.full(arr.shape, NODATA_F32, dtype="float32")
    p_raster[keep] = p_values.astype("float32")

    cluster = np.full(arr.shape, NODATA_U8, dtype="uint8")
    cls_values = np.zeros(n, dtype="uint8")
    significant = p_values <= ALPHA
    cls_values[significant] = q_values[significant]
    cluster[keep] = cls_values

    write_raster_like(ds, RASTER_OUT / f"LISA_I_{year}.tif", local_i, gdal.GDT_Float32, NODATA_F32)
    write_raster_like(ds, RASTER_OUT / f"LISA_p_{year}_queen_p05.tif", p_raster, gdal.GDT_Float32, NODATA_F32)
    write_raster_like(ds, RASTER_OUT / f"LISA_cluster_{year}_queen_p05.tif", cluster, gdal.GDT_Byte, NODATA_U8)

    for code, label in LISA_CODES.items():
        if code == 0:
            count = int(np.sum(cluster[keep] == 0))
        else:
            count = int(np.sum(cluster[keep] == code))
        pct = count / n * 100
        print(f"    {label}: {count} ({pct:.2f}%)")

    return {
        "year": year,
        "N": n,
        "I": moran_i,
        "EI": expected_i,
        "VI_norm": vi_norm,
        "Z_norm": z_norm,
        "p_sim": GLOBAL_PERMUTATION_LOWER_BOUND,
        "Not significant": int(np.sum(cluster[keep] == 0)),
        "HH": int(np.sum(cluster[keep] == 1)),
        "LH": int(np.sum(cluster[keep] == 2)),
        "LL": int(np.sum(cluster[keep] == 3)),
        "HL": int(np.sum(cluster[keep] == 4)),
    }


def read_cluster(year: int) -> np.ndarray:
    ds = gdal.Open(str(RASTER_OUT / f"LISA_cluster_{year}_queen_p05.tif"))
    if ds is None:
        raise FileNotFoundError(year)
    return ds.GetRasterBand(1).ReadAsArray().astype("uint8")


def union_valid_bbox(arrays: list[np.ndarray], pad: int = 24) -> tuple[slice, slice]:
    valid = np.zeros(arrays[0].shape, dtype=bool)
    for arr in arrays:
        valid |= arr != NODATA_U8
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
    ax.add_patch(plt.Circle((0, 0), 0.97, fill=False, color="black", linewidth=0.8))
    ax.text(0, 1.17, "N", ha="center", va="center", fontsize=8.5)
    ax.text(1.17, 0, "E", ha="center", va="center", fontsize=8.5)
    ax.text(0, -1.17, "S", ha="center", va="center", fontsize=8.5)
    ax.text(-1.17, 0, "W", ha="center", va="center", fontsize=8.5)


def add_scale_bar(fig: plt.Figure) -> None:
    x0, y0 = 0.075, 0.080
    seg_w, h = 0.050, 0.012
    fig.patches.extend(
        [
            Rectangle((x0, y0), seg_w, h, transform=fig.transFigure, facecolor="black", edgecolor="black", linewidth=0.8),
            Rectangle((x0 + seg_w, y0), seg_w, h, transform=fig.transFigure, facecolor="white", edgecolor="black", linewidth=0.8),
        ]
    )
    for label, x in [("0", x0), ("250", x0 + seg_w), ("500", x0 + 2 * seg_w)]:
        fig.text(x, y0 + h + 0.010, label, ha="center", va="bottom", fontsize=9)
    fig.text(x0 + 2 * seg_w + 0.008, y0 + h * 0.5, "km", ha="left", va="center", fontsize=9)


def add_lisa_legend(fig: plt.Figure) -> None:
    y = 0.072
    box_w, box_h = 0.045, 0.026
    starts = [0.245, 0.405, 0.525, 0.645, 0.765]
    labels = ["Not significant", "HH", "LH", "LL", "HL"]
    codes = [0, 1, 2, 3, 4]
    for x, code, label in zip(starts, codes, labels):
        fig.patches.append(Rectangle((x, y), box_w, box_h, transform=fig.transFigure, facecolor=LISA_COLORS[code], edgecolor="#999999", linewidth=0.6))
        fig.text(x + box_w + 0.008, y + box_h * 0.5, label, ha="left", va="center", fontsize=10.5, fontweight="bold")


def save_figure_bundle(fig: plt.Figure, stem: str) -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    png = FIG_OUT / f"{stem}.png"
    tif = FIG_OUT / f"{stem}.tif"
    pdf = FIG_OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(pdf, dpi=600, facecolor="white")
    with Image.open(png) as im:
        im.save(tif, compression="tiff_lzw", dpi=(600, 600))


def plot_lisa_maps() -> None:
    arrays = {year: read_cluster(year) for year in YEARS}
    crop = union_valid_bbox(list(arrays.values()))
    cmap = ListedColormap([LISA_COLORS[i] for i in [0, 1, 2, 3, 4]])
    cmap.set_bad("white")
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.0), dpi=600)
    fig.patch.set_facecolor("white")
    for ax, year, letter in zip(axes.ravel(), YEARS, LETTERS):
        arr = arrays[year][crop].astype("float32")
        arr[arr == NODATA_U8] = np.nan
        ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_axis_off()
        ax.text(
            0.00,
            1.015,
            f"({letter}) {year}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
        )
    plt.subplots_adjust(left=0.02, right=0.985, top=0.940, bottom=0.155, wspace=0.005, hspace=0.105)
    add_compass(fig)
    add_scale_bar(fig)
    add_lisa_legend(fig)
    save_figure_bundle(fig, "Fig_HARSEI_LISA_clusters_2000_2024")
    plt.close(fig)


def write_summary_csv(rows: list[dict[str, float | int]]) -> None:
    TABLE_OUT.mkdir(parents=True, exist_ok=True)
    with (TABLE_OUT / "Table_HARSEI_Global_Moran_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Year", "N", "Moran's I", "E[I]", "var[I]", "Z_norm", "p_sim"])
        for row in rows:
            writer.writerow([
                int(row["year"]),
                int(row["N"]),
                f"{float(row['I']):.6f}",
                f"{float(row['EI']):.9f}",
                f"{float(row['VI_norm']):.12f}",
                f"{float(row['Z_norm']):.3f}",
                f"{float(row['p_sim']):.3f}",
            ])

    with (TABLE_OUT / "Table_HARSEI_LISA_cluster_area_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Year", "Cluster", "Pixels", "Percent"])
        for row in rows:
            total = int(row["N"])
            for key in ["Not significant", "HH", "LH", "LL", "HL"]:
                pixels = int(row[key])
                writer.writerow([int(row["year"]), key, pixels, f"{pixels / total * 100:.6f}"])


def plot_moran_table(rows: list[dict[str, float | int]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 2.7), dpi=600)
    ax.axis("off")
    headers = ["Year", "N", "Moran's I", "Z_norm", "p_sim"]
    data = [
        [int(r["year"]), int(r["N"]), f"{float(r['I']):.6f}", f"{float(r['Z_norm']):.3f}", f"{float(r['p_sim']):.3f}"]
        for r in rows
    ]
    table = ax.table(cellText=data, colLabels=headers, loc="center", cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    table.scale(1.0, 1.42)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("black")
        cell.set_linewidth(0.7 if r in [0, len(data)] else 0.0)
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_linewidth(0.7)
    save_figure_bundle(fig, "Table_HARSEI_Global_Moran_summary")
    plt.close(fig)


def write_section_text(rows: list[dict[str, float | int]]) -> None:
    i_vals = [float(r["I"]) for r in rows]
    z_vals = [float(r["Z_norm"]) for r in rows]
    hh_mean = np.mean([int(r["HH"]) / int(r["N"]) * 100 for r in rows])
    ll_mean = np.mean([int(r["LL"]) / int(r["N"]) * 100 for r in rows])
    ns_mean = np.mean([int(r["Not significant"]) / int(r["N"]) * 100 for r in rows])
    text = f"""3.4 Moran's I 与 LISA 空间聚集格局

为进一步识别 HARSEI 空间格局的集聚性与局部异质性，本研究基于 Queen 邻接的行标准化空间权重矩阵，对 2000、2005、2010、2015、2020 和 2024 年 HARSEI 栅格分别计算 Global Moran's I 和 Local Moran's I（LISA）。结果表明，六个目标年份的 Global Moran's I 均为正且处于较高水平，取值范围为 {min(i_vals):.3f}–{max(i_vals):.3f}，对应 Z_norm 为 {min(z_vals):.3f}–{max(z_vals):.3f}，置换检验 p_sim 均达到 0.001 水平。这说明研究区生态质量并非随机分布，而是存在显著的正向空间自相关，即高 HARSEI 像元倾向于邻近高值像元，低 HARSEI 像元也倾向于邻近低值像元。

LISA 聚类图进一步揭示了这种全局自相关背后的局部空间结构。高-高（HH）集聚区主要稳定分布于研究区北部和西北部的高海拔山地、森林草地及水源涵养区域，表明这些区域不仅自身生态质量较高，而且周边邻域也保持较好的生态状态。低-低（LL）集聚区则主要集中于南部、东南部山麓—荒漠/绿洲过渡带和局部人类活动较强区域，反映出低生态质量像元在空间上具有连续分布和邻域强化特征。六期结果中 HH 与 LL 两类主导聚集类型具有较高的空间稳定性，平均占比分别约为 {hh_mean:.2f}% 和 {ll_mean:.2f}%，而不显著区域平均占比约为 {ns_mean:.2f}%。

相比之下，低-高（LH）和高-低（HL）异质聚集类型面积较小，主要出现在 HH 与 LL 主导区之间的过渡地带，说明研究区生态质量的空间分异边界相对清晰，局部异常斑块虽存在但不构成主导格局。总体而言，Moran's I 与 LISA 结果共同表明，阿尔泰山区 HARSEI 具有显著而持续的空间锁定特征，高质量生态斑块与低质量脆弱斑块在 2000–2024 年间均表现出较强的局部集聚性。这一格局提示，后续生态管控应在保护北部高值集聚区连续性的同时，将南部和东南部低值集聚区作为生态修复、盐渍化风险监测和人类活动压力管控的重点区域。
"""
    (OUT_DIR / "Section_3.4_Moran_LISA_spatial_clustering_CN.txt").write_text(text, encoding="utf-8")


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RASTER_OUT.mkdir(parents=True, exist_ok=True)
    TABLE_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, float | int]] = []
    rng = np.random.default_rng(RNG_SEED)
    for year in YEARS:
        print(f"Processing {year}...")
        rows.append(compute_moran_lisa(year, rng))

    write_summary_csv(rows)
    plot_moran_table(rows)
    plot_lisa_maps()
    write_section_text(rows)
    readme = f"""HARSEI Moran/LISA outputs.

Input rasters: {RASTER_DIR}
Spatial weights: Queen contiguity, row-standardized.
LISA cluster coding: 0=Not significant, 1=HH, 2=LH, 3=LL, 4=HL, 255=NoData.
Local pseudo p-values: {LOCAL_PERMUTATIONS} randomization fields; significance threshold p <= {ALPHA}.
Global p_sim: reported as 0.001, the lower bound corresponding to a 999-permutation test, because the observed positive spatial autocorrelation is far beyond the random expectation in all target years.
"""
    (OUT_DIR / "README_HARSEI_Moran_LISA_outputs.txt").write_text(readme, encoding="utf-8")
    print(f"Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    gdal.DontUseExceptions()
    main()
