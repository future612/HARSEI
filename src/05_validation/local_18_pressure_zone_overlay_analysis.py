#!/usr/bin/env python
"""Quantitatively overlay HARSEI-RSEI differences with pressure-zone masks.

Supported zone inputs:
  - raster masks: .tif/.tiff, values > 0 are treated as zone pixels
  - vector masks: .shp/.gpkg/.geojson, all geometries are rasterized as zone pixels

Built-in zone inputs:
  - ISRIC salinity classes cropped in the previous validation workflow:
      class > 0 -> saline_any
      class >= 2 -> saline_moderate_high

Optional user/downloaded inputs:
  - irrigated-zone mask from GFSAD/GMIA/GEE
  - mining-zone mask from Global Mining Areas v2 or local mine polygons
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.vrt import WarpedVRT


NODATA = -9999.0


def parse_args() -> argparse.Namespace:
    root = Path(r"D:\Codex\260724 小论文")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--harsei-dir", type=Path, default=root / "revise/annual_harsei_outputs/rasters/HARSEI")
    parser.add_argument("--rsei-dir", type=Path, default=root / "revise/annual_harsei_outputs/rasters/RSEI")
    parser.add_argument("--salinity-dir", type=Path, default=root / "revise/external_validation_inputs/isric_global_soil_salinity_remote_crops")
    parser.add_argument("--extra-zone-dir", type=Path, default=root / "revise/external_validation_inputs/pressure_zone_masks")
    parser.add_argument("--irrigation-mask", type=Path, default=root / "revise/external_validation_inputs/pressure_zone_masks/YJQ_irrigated_area_mask.tif")
    parser.add_argument("--mining-mask", type=Path, default=root / "revise/external_validation_inputs/mining_maus_2022/global_miningarea_v2_30arcsecond.tif")
    parser.add_argument("--out-dir", type=Path, default=root / "revise/pressure_zone_overlay_validation")
    parser.add_argument("--years", default="2000-2024")
    return parser.parse_args()


def parse_years(text: str) -> list[int]:
    if "-" in text:
        a, b = [int(v) for v in text.split("-", 1)]
        return list(range(a, b + 1))
    return [int(v.strip()) for v in text.split(",") if v.strip()]


def read_raster(path: Path) -> tuple[np.ndarray, dict]:
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        nodata = src.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[arr == NODATA] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return arr, profile


def ref_profile(path: Path) -> dict:
    with rasterio.open(path) as src:
        return src.profile.copy()


def raster_mask_to_ref(path: Path, ref: dict) -> np.ndarray:
    with rasterio.open(path) as src:
        with WarpedVRT(
            src,
            crs=ref["crs"],
            transform=ref["transform"],
            width=ref["width"],
            height=ref["height"],
            resampling=Resampling.nearest,
        ) as vrt:
            arr = vrt.read(1)
            nodata = vrt.nodata
    if nodata is not None:
        arr = np.where(arr == nodata, 0, arr)
    return np.isfinite(arr) & (arr > 0)


def vector_mask_to_ref(path: Path, ref: dict) -> np.ndarray:
    import geopandas as gpd

    gdf = gpd.read_file(path)
    if gdf.empty:
        return np.zeros((ref["height"], ref["width"]), dtype=bool)
    if gdf.crs is None:
        raise ValueError(f"Vector has no CRS: {path}")
    gdf = gdf.to_crs(ref["crs"])
    shapes = [(geom, 1) for geom in gdf.geometry if geom is not None and not geom.is_empty]
    arr = rasterize(
        shapes,
        out_shape=(ref["height"], ref["width"]),
        transform=ref["transform"],
        fill=0,
        dtype="uint8",
    )
    return arr > 0


def load_external_mask(path: Path, ref: dict) -> np.ndarray | None:
    if not path.exists():
        return None
    suffix = path.suffix.lower()
    if suffix in {".tif", ".tiff"}:
        return raster_mask_to_ref(path, ref)
    if suffix in {".shp", ".gpkg", ".geojson", ".json"}:
        return vector_mask_to_ref(path, ref)
    raise ValueError(f"Unsupported mask type: {path}")


def salinity_union_masks(sal_dir: Path, ref: dict) -> dict[str, np.ndarray]:
    any_mask = np.zeros((ref["height"], ref["width"]), dtype=bool)
    moderate_mask = np.zeros_like(any_mask)
    if not sal_dir.exists():
        return {}
    for path in sorted(sal_dir.glob("*.tif")):
        with rasterio.open(path) as src:
            with WarpedVRT(
                src,
                crs=ref["crs"],
                transform=ref["transform"],
                width=ref["width"],
                height=ref["height"],
                resampling=Resampling.nearest,
            ) as vrt:
                arr = vrt.read(1).astype(np.float32)
        arr[~np.isfinite(arr)] = np.nan
        any_mask |= arr > 0
        moderate_mask |= arr >= 2
    masks = {}
    if any_mask.any():
        masks["ISRIC_saline_class_gt0_union"] = any_mask
    if moderate_mask.any():
        masks["ISRIC_saline_class_ge2_union"] = moderate_mask
    return masks


def finite_values(arr: np.ndarray) -> np.ndarray:
    return arr[np.isfinite(arr)]


def stats_for_mask(values: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    vals = finite_values(values[mask])
    if vals.size == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan, "std": np.nan, "p25": np.nan, "p75": np.nan}
    return {
        "n": int(vals.size),
        "mean": float(np.nanmean(vals)),
        "median": float(np.nanmedian(vals)),
        "std": float(np.nanstd(vals)),
        "p25": float(np.nanpercentile(vals, 25)),
        "p75": float(np.nanpercentile(vals, 75)),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_diff_raster(path: Path, diff: np.ndarray, profile: dict) -> None:
    out_profile = profile.copy()
    out_profile.update(count=1, dtype="float32", nodata=NODATA, compress="deflate")
    out = np.where(np.isfinite(diff), diff, NODATA).astype(np.float32)
    with rasterio.open(path, "w", **out_profile) as dst:
        dst.write(out, 1)


def try_make_plot(summary_csv: Path, out_dir: Path) -> None:
    try:
        import pandas as pd
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Plot skipped: {exc}")
        return
    df = pd.read_csv(summary_csv)
    if df.empty:
        return
    pivot = df.pivot_table(index="year", columns="zone", values="diff_median_zone_minus_background", aggfunc="first")
    fig, ax = plt.subplots(figsize=(9, 4.6), dpi=300)
    ax.axhline(0, color="0.25", lw=0.8)
    for col in pivot.columns:
        label = col.replace("ISRIC_saline_class_gt0_union", "ISRIC saline > 0")
        label = label.replace("ISRIC_saline_class_ge2_union", "ISRIC salinity >= 2")
        label = label.replace("_", " ")
        ax.plot(pivot.index, pivot[col], marker="o", lw=1.8, label=label)
    ax.set_xlabel("Year")
    ax.set_ylabel("Median(HARSEI - RSEI): zone - background")
    years = list(pivot.index)
    tick_step = 2 if len(years) <= 25 else 5
    ax.set_xticks(years[::tick_step])
    ax.tick_params(axis="x", labelrotation=45)
    ax.grid(True, color="0.88", lw=0.8)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pressure_zone_harsei_minus_rsei_overlay.png")
    plt.close(fig)

    agg = df.groupby("zone").agg(
        years=("year", "count"),
        mean_diff_median_zone_minus_background=("diff_median_zone_minus_background", "mean"),
        min_diff_median_zone_minus_background=("diff_median_zone_minus_background", "min"),
        max_diff_median_zone_minus_background=("diff_median_zone_minus_background", "max"),
        mean_HARSEI_zone_minus_background=("HARSEI_mean_zone_minus_background", "mean"),
        mean_RSEI_zone_minus_background=("RSEI_mean_zone_minus_background", "mean"),
        median_zone_pixels=("zone_pixels", "median"),
    ).reset_index()
    agg.to_csv(out_dir / "pressure_zone_overlay_cross_year_summary.csv", index=False, encoding="utf-8-sig")
    try:
        (out_dir / "pressure_zone_overlay_cross_year_summary.md").write_text(agg.to_markdown(index=False), encoding="utf-8")
    except Exception:
        pass


def main() -> None:
    args = parse_args()
    years = parse_years(args.years)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    diff_dir = args.out_dir / "difference_rasters"
    diff_dir.mkdir(parents=True, exist_ok=True)

    first_ref = args.harsei_dir / f"YJQ_HARSEI_{years[0]}.tif"
    if not first_ref.exists():
        raise FileNotFoundError(first_ref)
    ref = ref_profile(first_ref)

    zones = salinity_union_masks(args.salinity_dir, ref)
    if args.extra_zone_dir.exists():
        explicit_paths = {
            args.irrigation_mask.resolve() if args.irrigation_mask.exists() else None,
            args.mining_mask.resolve() if args.mining_mask.exists() else None,
        }
        for suffix in ("*.tif", "*.tiff", "*.shp", "*.gpkg", "*.geojson", "*.json"):
            for path in sorted(args.extra_zone_dir.glob(suffix)):
                if "_mask" not in path.stem.lower():
                    continue
                if path.resolve() in explicit_paths:
                    continue
                loaded = load_external_mask(path, ref)
                if loaded is not None and loaded.any():
                    zones[path.stem] = loaded
    irrigation = load_external_mask(args.irrigation_mask, ref)
    if irrigation is not None and irrigation.any():
        zones["irrigated_area"] = irrigation
    mining = load_external_mask(args.mining_mask, ref)
    if mining is not None and mining.any():
        zones["mining_area"] = mining

    manifest_rows = []
    for name, mask in zones.items():
        manifest_rows.append({"zone": name, "pixels_on_reference_grid": int(mask.sum())})
    write_csv(args.out_dir / "pressure_zone_mask_manifest.csv", manifest_rows, ["zone", "pixels_on_reference_grid"])

    if not zones:
        missing = [
            "No pressure-zone masks were available.",
            f"Expected irrigation mask: {args.irrigation_mask}",
            f"Expected mining mask: {args.mining_mask}",
            f"Expected ISRIC salinity crops under: {args.salinity_dir}",
        ]
        (args.out_dir / "missing_pressure_zone_inputs.md").write_text("\n".join(missing), encoding="utf-8")
        print("No masks available; wrote missing input note.")
        return

    rows = []
    for year in years:
        harsei_path = args.harsei_dir / f"YJQ_HARSEI_{year}.tif"
        rsei_path = args.rsei_dir / f"YJQ_RSEI_{year}.tif"
        if not harsei_path.exists() or not rsei_path.exists():
            print(f"skip {year}: missing HARSEI/RSEI")
            continue
        harsei, profile = read_raster(harsei_path)
        rsei, _ = read_raster(rsei_path)
        diff = harsei - rsei
        valid = np.isfinite(diff) & np.isfinite(harsei) & np.isfinite(rsei)
        save_diff_raster(diff_dir / f"YJQ_HARSEI_minus_RSEI_{year}.tif", diff, profile)

        for zone_name, zone_mask in zones.items():
            z = valid & zone_mask
            bg = valid & ~zone_mask
            diff_z = stats_for_mask(diff, z)
            diff_bg = stats_for_mask(diff, bg)
            harsei_z = stats_for_mask(harsei, z)
            harsei_bg = stats_for_mask(harsei, bg)
            rsei_z = stats_for_mask(rsei, z)
            rsei_bg = stats_for_mask(rsei, bg)
            rows.append({
                "year": year,
                "zone": zone_name,
                "zone_pixels": diff_z["n"],
                "background_pixels": diff_bg["n"],
                "diff_mean_zone": diff_z["mean"],
                "diff_mean_background": diff_bg["mean"],
                "diff_mean_zone_minus_background": diff_z["mean"] - diff_bg["mean"],
                "diff_median_zone": diff_z["median"],
                "diff_median_background": diff_bg["median"],
                "diff_median_zone_minus_background": diff_z["median"] - diff_bg["median"],
                "HARSEI_mean_zone": harsei_z["mean"],
                "HARSEI_mean_background": harsei_bg["mean"],
                "HARSEI_mean_zone_minus_background": harsei_z["mean"] - harsei_bg["mean"],
                "RSEI_mean_zone": rsei_z["mean"],
                "RSEI_mean_background": rsei_bg["mean"],
                "RSEI_mean_zone_minus_background": rsei_z["mean"] - rsei_bg["mean"],
            })

    summary_csv = args.out_dir / "pressure_zone_harsei_minus_rsei_overlay_summary.csv"
    fields = [
        "year", "zone", "zone_pixels", "background_pixels",
        "diff_mean_zone", "diff_mean_background", "diff_mean_zone_minus_background",
        "diff_median_zone", "diff_median_background", "diff_median_zone_minus_background",
        "HARSEI_mean_zone", "HARSEI_mean_background", "HARSEI_mean_zone_minus_background",
        "RSEI_mean_zone", "RSEI_mean_background", "RSEI_mean_zone_minus_background",
    ]
    write_csv(summary_csv, rows, fields)
    try_make_plot(summary_csv, args.out_dir)
    print(f"Wrote {summary_csv}")
    print(f"Wrote {args.out_dir / 'pressure_zone_mask_manifest.csv'}")


if __name__ == "__main__":
    main()
