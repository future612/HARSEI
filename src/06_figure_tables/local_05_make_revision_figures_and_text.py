#!/usr/bin/env python
"""
Create manuscript-ready SVG figures and bilingual revision text from the
annual HARSEI rebuild tables. Uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def svg_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def scale(value: float, src_min: float, src_max: float, dst_min: float, dst_max: float) -> float:
    if src_max == src_min:
        return (dst_min + dst_max) / 2
    return dst_min + (value - src_min) / (src_max - src_min) * (dst_max - dst_min)


def polyline(points: list[tuple[float, float]], color: str, width: float = 2.5) -> str:
    pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'


def circle(x: float, y: float, color: str) -> str:
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}" stroke="white" stroke-width="1"/>'


def text(x: float, y: float, label: str, size: int = 12, anchor: str = "middle", weight: str = "400") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="#222">{svg_escape(label)}</text>'
    )


def make_annual_trajectory(rows: list[dict[str, str]], out: Path) -> None:
    w, h = 960, 560
    ml, mr, mt, mb = 82, 28, 58, 70
    plot_w = w - ml - mr
    plot_h = h - mt - mb
    years = [int(r["year"]) for r in rows]
    series = {
        "AWRSEI": ([f(r, "AWRSEI_mean") for r in rows], "#1b9e77"),
        "HAI reverse": ([f(r, "HAI_reverse_mean") for r in rows], "#d95f02"),
        "HARSEI": ([f(r, "HARSEI_mean") for r in rows], "#377eb8"),
    }
    y_min, y_max = 0.0, 1.0
    x_min, x_max = min(years), max(years)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<rect width="100%" height="100%" fill="white"/>',
        text(w / 2, 32, "Annual AWRSEI, HAI-reverse, and HARSEI Means (2000-2024)", 18, "middle", "700"),
        f'<rect x="{ml}" y="{mt}" width="{plot_w}" height="{plot_h}" fill="#fbfbfb" stroke="#d0d0d0"/>',
    ]
    for tick in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        y = scale(tick, y_min, y_max, mt + plot_h, mt)
        parts.append(f'<line x1="{ml}" y1="{y:.2f}" x2="{ml + plot_w}" y2="{y:.2f}" stroke="#e6e6e6"/>')
        parts.append(text(ml - 12, y + 4, f"{tick:.1f}", 11, "end"))
    for tick in [2000, 2005, 2010, 2015, 2020, 2024]:
        x = scale(tick, x_min, x_max, ml, ml + plot_w)
        parts.append(f'<line x1="{x:.2f}" y1="{mt}" x2="{x:.2f}" y2="{mt + plot_h}" stroke="#eeeeee"/>')
        parts.append(text(x, mt + plot_h + 24, str(tick), 11, "middle"))
    parts.append(text(28, mt + plot_h / 2, "Index value", 13, "middle", "700").replace("<text", '<text transform="rotate(-90 28 {0:.2f})"'.format(mt + plot_h / 2), 1))
    parts.append(text(ml + plot_w / 2, h - 20, "Year", 13, "middle", "700"))

    for label, (vals, color) in series.items():
        pts = [
            (scale(year, x_min, x_max, ml, ml + plot_w), scale(val, y_min, y_max, mt + plot_h, mt))
            for year, val in zip(years, vals)
        ]
        parts.append(polyline(pts, color, 2.8 if label == "HARSEI" else 2.2))
        for idx in [0, 5, 10, 15, 20, 24]:
            parts.append(circle(pts[idx][0], pts[idx][1], color))

    lx, ly = ml + plot_w - 180, mt + 20
    for i, (label, (_, color)) in enumerate(series.items()):
        y = ly + i * 24
        parts.append(f'<line x1="{lx}" y1="{y}" x2="{lx + 28}" y2="{y}" stroke="{color}" stroke-width="3"/>')
        parts.append(text(lx + 36, y + 4, label, 12, "start"))
    parts.append("</svg>")
    out.write_text("\n".join(parts), encoding="utf-8")


def make_sensitivity_figure(snow_rows: list[dict[str, str]], thresh_rows: list[dict[str, str]], out: Path) -> None:
    w, h = 960, 560
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<rect width="100%" height="100%" fill="white"/>',
        text(w / 2, 32, "Water Mask and Snow Sensitivity", 18, "middle", "700"),
    ]
    # Panel 1: high-elevation snow frequency at NDWI=0.2
    panel1 = (70, 70, 390, 390)
    x0, y0, pw, ph = panel1
    parts.append(text(x0 + pw / 2, y0 - 18, "High-elevation snow/ice frequency", 14, "middle", "700"))
    parts.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfbfb" stroke="#d0d0d0"/>')
    snow = {}
    for r in snow_rows:
        if abs(float(r["ndwi_threshold"]) - 0.2) < 1e-9:
            snow[int(float(r["start_month"]))] = float(r["mean_snow_fraction_high_elevation"])
    max_snow = 0.25
    for tick in [0, 0.05, 0.10, 0.15, 0.20, 0.25]:
        y = scale(tick, 0, max_snow, y0 + ph, y0)
        parts.append(f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + pw}" y2="{y:.2f}" stroke="#e6e6e6"/>')
        parts.append(text(x0 - 10, y + 4, f"{tick:.2f}", 10, "end"))
    labels = [("Apr-Sep", snow.get(4, 0.0), "#8da0cb"), ("May-Sep", snow.get(5, 0.0), "#66c2a5")]
    for i, (label, val, color) in enumerate(labels):
        bw = 80
        bx = x0 + 90 + i * 130
        by = scale(val, 0, max_snow, y0 + ph, y0)
        bh = y0 + ph - by
        parts.append(f'<rect x="{bx}" y="{by:.2f}" width="{bw}" height="{bh:.2f}" fill="{color}" stroke="#555"/>')
        parts.append(text(bx + bw / 2, y0 + ph + 24, label, 11))
        parts.append(text(bx + bw / 2, by - 8, f"{val:.4f}", 11))

    # Panel 2: HARSEI mean under stricter thresholds within exported domain.
    panel2 = (535, 70, 355, 390)
    x0, y0, pw, ph = panel2
    parts.append(text(x0 + pw / 2, y0 - 18, "Mean HARSEI under stricter NDWI masks (zoomed y-axis)", 14, "middle", "700"))
    parts.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfbfb" stroke="#d0d0d0"/>')
    thresh = {}
    for r in thresh_rows:
        t = float(r["ndwi_threshold"])
        thresh.setdefault(t, []).append(float(r["HARSEI_mean"]))
    avg_thresh = [(t, sum(vals) / len(vals)) for t, vals in sorted(thresh.items())]
    ymin, ymax = 0.5250, 0.5255
    for tick in [0.5250, 0.5251, 0.5252, 0.5253, 0.5254, 0.5255]:
        y = scale(tick, ymin, ymax, y0 + ph, y0)
        parts.append(f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + pw}" y2="{y:.2f}" stroke="#e6e6e6"/>')
        parts.append(text(x0 - 10, y + 4, f"{tick:.4f}", 10, "end"))
    for i, (t, val) in enumerate(avg_thresh):
        bw = 62
        bx = x0 + 70 + i * 95
        by = scale(val, ymin, ymax, y0 + ph, y0)
        bh = y0 + ph - by
        parts.append(f'<rect x="{bx}" y="{by:.2f}" width="{bw}" height="{bh:.2f}" fill="#fc8d62" stroke="#555"/>')
        parts.append(text(bx + bw / 2, y0 + ph + 24, f">{t:.1f}", 11))
        parts.append(text(bx + bw / 2, by - 8, f"{val:.6f}", 10))
    parts.append(text(70 + 390 / 2, h - 28, "Compositing window", 12, "middle", "700"))
    parts.append(text(535 + 355 / 2, h - 28, "NDWI water threshold", 12, "middle", "700"))
    parts.append("</svg>")
    out.write_text("\n".join(parts), encoding="utf-8")


def table_md(headers: list[str], rows: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return out


def make_text_package(tables_dir: Path, out: Path) -> None:
    annual = read_csv(tables_dir / "annual_index_summary.csv")
    weights = read_csv(tables_dir / "entropy_weights_and_light_harmonization.csv")
    pca = read_csv(tables_dir / "pooled_pca_pc1.csv")
    validation = read_csv(tables_dir / "input_stack_validation_summary.csv")
    snow_md = (tables_dir / "water_snow_sensitivity_summary.md").read_text(encoding="utf-8")
    thresh_md = (tables_dir / "water_threshold_index_sensitivity_current_mask.md").read_text(encoding="utf-8")

    y2000 = next(r for r in annual if r["year"] == "2000")
    y2024 = next(r for r in annual if r["year"] == "2024")
    harsei_vals = [(int(r["year"]), f(r, "HARSEI_mean")) for r in annual]
    min_year, min_val = min(harsei_vals, key=lambda x: x[1])
    max_year, max_val = max(harsei_vals, key=lambda x: x[1])
    pca_ratio = f(pca[0], "pc1_explained_variance_ratio")

    w = {r["variable"]: float(r["weight"]) for r in weights}
    latest = validation[-1]
    lines = ["# Manuscript-Ready Results And Response Text", ""]
    lines.append("## Key Numbers")
    lines.extend(table_md(
        ["Item", "Value"],
        [
            ["PC1 explained variance", f"{pca_ratio:.4f}"],
            ["Mean yearly HARSEI fusion weight: AWRSEI", f"{w['AWRSEI']:.4f}"],
            ["Mean yearly HARSEI fusion weight: HAI_reverse", f"{w['HAI_reverse']:.4f}"],
            ["HAI weight: POP", f"{w['POP']:.4f}"],
            ["HAI weight: LIGHT", f"{w['LIGHT']:.4f}"],
            ["HAI weight: LUCC_SCORE", f"{w['LUCC_SCORE']:.4f}"],
            ["Night-light harmonization a", f"{w['DMSP=a+b*log1p(VIIRS):a']:.4f}"],
            ["Night-light harmonization b", f"{w['DMSP=a+b*log1p(VIIRS):b']:.4f}"],
            ["2024 POP year used", latest["pop_year_used"]],
            ["2024 LUCC year used", latest["lucc_year_used"]],
        ],
    ))
    lines.extend(["", "## 中文结果段落", ""])
    lines.append(
        f"基于 2000-2024 年固定参考范围归一化和 pooled PCA 的年度重建结果显示，"
        f"PC1 解释率为 {pca_ratio:.2%}，表明 NDVI、WET、NDBSI、LST 和 SRSI 所构成的综合生态梯度具有较高的信息集中度。"
        f"HARSEI 年均值从 2000 年的 {float(y2000['HARSEI_mean']):.4f} 变化到 2024 年的 {float(y2024['HARSEI_mean']):.4f}，"
        f"多年份最低值出现在 {min_year} 年（{min_val:.4f}），最高值出现在 {max_year} 年（{max_val:.4f}）。"
        f"年度熵权融合结果显示，AWRSEI 与反向 HAI 的多年平均权重分别为 {w['AWRSEI']:.4f} 和 {w['HAI_reverse']:.4f}，"
        f"说明自然生态组分仍是 HARSEI 的主体信息来源，而 HAI 作为人类活动压力修正项参与综合评价。"
    )
    lines.append(
        f"HAI 内部三项权重分别为 POP {w['POP']:.4f}、LIGHT {w['LIGHT']:.4f} 和 LUCC_SCORE {w['LUCC_SCORE']:.4f}。"
        f"夜间灯光跨传感器协调采用 2012-2013 年重叠期拟合关系 `DMSP = a + b * log1p(VIIRS)`，"
        f"参数为 a = {w['DMSP=a+b*log1p(VIIRS):a']:.4f}、b = {w['DMSP=a+b*log1p(VIIRS):b']:.4f}。"
        f"当前下载栈中，2024 年人口层沿用 {latest['pop_year_used']} 年 WorldPop，"
        f"LUCC 层沿用 {latest['lucc_year_used']} 年 GLC-FCS30D，并已在输出诊断表中记录。"
    )
    lines.extend(["", "## English Results Paragraph", ""])
    lines.append(
        f"The annual reconstruction based on fixed-reference normalization and pooled PCA for 2000-2024 showed that "
        f"PC1 explained {pca_ratio:.2%} of the variance, indicating that the ecological gradient represented by NDVI, WET, NDBSI, LST and SRSI was highly concentrated in the first component. "
        f"The mean HARSEI changed from {float(y2000['HARSEI_mean']):.4f} in 2000 to {float(y2024['HARSEI_mean']):.4f} in 2024; "
        f"the lowest annual mean occurred in {min_year} ({min_val:.4f}), whereas the highest occurred in {max_year} ({max_val:.4f}). "
        f"Year-specific entropy-based fusion assigned multi-year mean weights of {w['AWRSEI']:.4f} to AWRSEI and {w['HAI_reverse']:.4f} to reverse HAI, "
        f"showing that the natural ecological components remained the dominant information source while HAI acted as a human-pressure adjustment term."
    )
    lines.append(
        f"Within HAI, the component weights were {w['POP']:.4f} for population, {w['LIGHT']:.4f} for harmonized night-time lights, "
        f"and {w['LUCC_SCORE']:.4f} for LUCC_SCORE. Inter-sensor harmonization of night-time lights used the 2012-2013 overlap period and the equation "
        f"`DMSP = a + b * log1p(VIIRS)`, with a = {w['DMSP=a+b*log1p(VIIRS):a']:.4f} and b = {w['DMSP=a+b*log1p(VIIRS):b']:.4f}. "
        f"In the current downloaded stack, the 2024 population layer uses WorldPop {latest['pop_year_used']} and the LUCC layer uses GLC-FCS30D {latest['lucc_year_used']}, "
        f"as recorded in the diagnostic output table."
    )
    lines.extend(["", "## Sensitivity Text Already Generated", "", snow_md, "", thresh_md])
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    tables_dir = Path(args.tables_dir)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    annual = read_csv(tables_dir / "annual_index_summary.csv")
    snow = read_csv(tables_dir / "water_snow_sensitivity_summary.csv")
    threshold = read_csv(tables_dir / "water_threshold_index_sensitivity_current_mask.csv")

    make_annual_trajectory(annual, fig_dir / "fig_annual_harsei_trajectory.svg")
    make_sensitivity_figure(snow, threshold, fig_dir / "fig_water_snow_sensitivity.svg")
    make_text_package(tables_dir, out_dir / "manuscript_ready_results_and_response_text.md")


if __name__ == "__main__":
    main()
