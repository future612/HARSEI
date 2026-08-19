#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create Fig. 16 AWRSEI-HAI bivariate state-pressure zoning GeoTIFF."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from osgeo import gdal


ROOT = Path(r"D:\Codex\260724 小论文")
AWRSEI = ROOT / "revise" / "annual_harsei_outputs" / "rasters" / "AWRSEI" / "YJQ_AWRSEI_2024.tif"
HAI = ROOT / "revise" / "annual_harsei_outputs" / "rasters" / "HAI" / "YJQ_HAI_2024.tif"
OUT_DIR = ROOT / "意见修改" / "图表" / "地图_TIF_ArcGIS" / "Fig16_AWRSEI-HAI bivariate state-pressure zoning in 2024"
OUT_CLASS = OUT_DIR / "Fig16_AWRSEI_HAI_bivariate_state_pressure_2024_class.tif"
OUT_RGB = OUT_DIR / "Fig16_AWRSEI_HAI_bivariate_state_pressure_2024_RGB.tif"
OUT_TABLE = OUT_DIR / "Fig16_AWRSEI_HAI_bivariate_state_pressure_2024_legend_summary.csv"
NODATA_FLOAT = -9999.0
NODATA_CLASS = 0


CLASSES = {
    1: {
        "label": "Low state / low pressure",
        "label_cn": "低状态-低压力",
        "color": (253, 219, 199),
    },
    2: {
        "label": "Low state / high pressure",
        "label_cn": "低状态-高压力",
        "color": (178, 24, 43),
    },
    3: {
        "label": "High state / low pressure",
        "label_cn": "高状态-低压力",
        "color": (33, 102, 172),
    },
    4: {
        "label": "High state / high pressure",
        "label_cn": "高状态-高压力",
        "color": (146, 197, 222),
    },
}


def read_single(path: Path):
    ds = gdal.Open(str(path))
    if ds is None:
        raise RuntimeError(f"Cannot open raster: {path}")
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        arr[np.isclose(arr, nodata)] = np.nan
    arr[np.isclose(arr, NODATA_FLOAT)] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    arr[arr < -1.0e20] = np.nan
    arr[arr > 1.0e20] = np.nan
    return ds, arr


def check_grid(ref, other):
    if ref.RasterXSize != other.RasterXSize or ref.RasterYSize != other.RasterYSize:
        raise ValueError("AWRSEI and HAI raster sizes differ.")
    if tuple(round(v, 12) for v in ref.GetGeoTransform()) != tuple(round(v, 12) for v in other.GetGeoTransform()):
        raise ValueError("AWRSEI and HAI geotransforms differ.")
    if ref.GetProjection() != other.GetProjection():
        raise ValueError("AWRSEI and HAI projections differ.")


def write_class_tif(ref_ds, cls: np.ndarray):
    if OUT_CLASS.exists():
        OUT_CLASS.unlink()
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        str(OUT_CLASS),
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        1,
        gdal.GDT_Byte,
        options=["COMPRESS=LZW", "TILED=YES"],
    )
    out.SetGeoTransform(ref_ds.GetGeoTransform())
    out.SetProjection(ref_ds.GetProjection())
    band = out.GetRasterBand(1)
    color_table = gdal.ColorTable()
    color_table.SetColorEntry(0, (255, 255, 255, 0))
    for code, meta in CLASSES.items():
        r, g, b = meta["color"]
        color_table.SetColorEntry(code, (r, g, b, 255))
    band.SetRasterColorTable(color_table)
    band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
    band.SetNoDataValue(NODATA_CLASS)
    band.SetDescription("Fig16 AWRSEI-HAI bivariate state-pressure class")
    band.WriteArray(cls)
    band.FlushCache()
    out.FlushCache()
    out = None


def write_rgb_tif(ref_ds, cls: np.ndarray):
    if OUT_RGB.exists():
        OUT_RGB.unlink()
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        str(OUT_RGB),
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        3,
        gdal.GDT_Byte,
        options=["COMPRESS=LZW", "TILED=YES"],
    )
    out.SetGeoTransform(ref_ds.GetGeoTransform())
    out.SetProjection(ref_ds.GetProjection())
    rgb = np.full((3, cls.shape[0], cls.shape[1]), 255, dtype=np.uint8)
    for code, meta in CLASSES.items():
        mask = cls == code
        for band_idx, value in enumerate(meta["color"]):
            rgb[band_idx][mask] = value
    for i in range(3):
        band = out.GetRasterBand(i + 1)
        band.WriteArray(rgb[i])
        band.FlushCache()
    out.FlushCache()
    out = None


def write_summary(cls: np.ndarray, aw_thr: float, hai_thr: float, aw: np.ndarray, hai: np.ndarray):
    total = int(np.sum(cls > 0))
    pixel_area_km2 = 0.998877  # MODIS sinusoidal-like 1 km degree-grid approximation after export; pixel count is primary.
    with OUT_TABLE.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "class_code",
                "class_name_en",
                "class_name_cn",
                "AWRSEI_condition",
                "HAI_condition",
                "color_rgb",
                "pixels",
                "percent",
                "approx_area_km2",
                "AWRSEI_median_threshold",
                "HAI_median_threshold",
                "class_mean_AWRSEI",
                "class_mean_HAI",
            ]
        )
        for code, meta in CLASSES.items():
            mask = cls == code
            pixels = int(mask.sum())
            aw_cond = "AWRSEI < median" if code in (1, 2) else "AWRSEI >= median"
            hai_cond = "HAI < median" if code in (1, 3) else "HAI >= median"
            writer.writerow(
                [
                    code,
                    meta["label"],
                    meta["label_cn"],
                    aw_cond,
                    hai_cond,
                    "-".join(str(v) for v in meta["color"]),
                    pixels,
                    f"{100.0 * pixels / total:.6f}" if total else "",
                    f"{pixels * pixel_area_km2:.3f}",
                    f"{aw_thr:.10f}",
                    f"{hai_thr:.10f}",
                    f"{float(np.nanmean(aw[mask])):.10f}" if pixels else "",
                    f"{float(np.nanmean(hai[mask])):.10f}" if pixels else "",
                ]
            )


def main():
    gdal.UseExceptions()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    aw_ds, aw = read_single(AWRSEI)
    hai_ds, hai = read_single(HAI)
    check_grid(aw_ds, hai_ds)
    valid = np.isfinite(aw) & np.isfinite(hai)
    aw_thr = float(np.nanmedian(aw[valid]))
    hai_thr = float(np.nanmedian(hai[valid]))

    cls = np.zeros(aw.shape, dtype=np.uint8)
    cls[valid & (aw < aw_thr) & (hai < hai_thr)] = 1
    cls[valid & (aw < aw_thr) & (hai >= hai_thr)] = 2
    cls[valid & (aw >= aw_thr) & (hai < hai_thr)] = 3
    cls[valid & (aw >= aw_thr) & (hai >= hai_thr)] = 4

    write_class_tif(aw_ds, cls)
    write_rgb_tif(aw_ds, cls)
    write_summary(cls, aw_thr, hai_thr, aw, hai)

    print(f"AWRSEI median threshold: {aw_thr:.6f}")
    print(f"HAI median threshold: {hai_thr:.6f}")
    print(f"Saved class GeoTIFF: {OUT_CLASS}")
    print(f"Saved RGB GeoTIFF: {OUT_RGB}")
    print(f"Saved legend summary: {OUT_TABLE}")


if __name__ == "__main__":
    main()
