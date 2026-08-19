#!/usr/bin/env python
"""Crop the ISRIC global soil salinity tiles to the YJQ/HARSEI extent.

This helper avoids downloading full 32768 x 32768 global tiles. It reads the
single ISRIC tile covering the study area through GDAL's /vsicurl/ support and
writes local GeoTIFF crops that can be consumed by the validation script.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import rasterio
from rasterio.windows import from_bounds


DEFAULT_YEARS = (2000, 2002, 2005, 2009, 2016)
URL_TEMPLATE = (
    "https://files.isric.org/public/global_soil_salinity/"
    "salmap{year}/salMap{year}-0000000000-0000098304.tif"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--harsei-ref",
        type=Path,
        default=Path("revise/annual_harsei_outputs/rasters/HARSEI/YJQ_HARSEI_2000.tif"),
        help="Reference HARSEI raster used only for the study-area bounds.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("revise/external_validation_inputs/isric_global_soil_salinity_remote_crops"),
        help="Directory for cropped ISRIC salinity GeoTIFFs.",
    )
    parser.add_argument(
        "--years",
        default=",".join(str(y) for y in DEFAULT_YEARS),
        help="Comma-separated years to crop.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    years = [int(y.strip()) for y in args.years.split(",") if y.strip()]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Make remote reads robust on Windows/Schannel and avoid probing unrelated URLs.
    os.environ.setdefault("GDAL_HTTP_UNSAFESSL", "YES")
    os.environ.setdefault("GDAL_HTTP_SSL_VERIFYPEER", "NO")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")

    with rasterio.open(args.harsei_ref) as ref:
        bounds = ref.bounds

    for year in years:
        url = URL_TEMPLATE.format(year=year)
        out_path = args.out_dir / f"ISRIC_global_soil_salinity_{year}_YJQ_crop.tif"
        print(f"Cropping {year}: {url}")

        with rasterio.open(url) as src:
            window = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
            data = src.read(1, window=window)
            transform = src.window_transform(window)
            profile = src.profile.copy()
            profile.update(
                driver="GTiff",
                height=data.shape[0],
                width=data.shape[1],
                transform=transform,
                compress="deflate",
                tiled=True,
                blockxsize=256,
                blockysize=256,
                BIGTIFF="IF_SAFER",
            )
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(data, 1)

        print(f"  wrote {out_path} ({out_path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
