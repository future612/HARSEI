// Annual HARSEI input export for YJQ / Altai study area.
// Copy this file into the Google Earth Engine Code Editor and run it.
//
// Main purpose:
// 1) Export annual ecological component stacks for 2000-2024:
//    NDVI, WET, NDBSI, LST_C, SRSI, plus NDWI/snow diagnostics.
// 2) Export a pooled sample CSV for fixed-reference normalization and pooled PCA.
// 3) Optionally export annual HAI input layers: POP, night-time lights, LUCC score.
//
// Recommended scientific workflow:
// - Do NOT run yearly min-max or yearly PCA in GEE.
// - Use the downloaded annual inputs with local_02_build_harsei_from_downloads.py
//   to compute fixed normalization, pooled PCA, fixed entropy weights, AWRSEI,
//   HAI and HARSEI.

var ROI_ASSET = 'projects/jlu2024621038/assets/yjq';
var roiFc = ee.FeatureCollection(ROI_ASSET);
var roi = roiFc.geometry();

var START_YEAR = 2000;
var END_YEAR = 2024;
var START_MONTH = 5;  // May. Use 4 for Apr-Sep sensitivity.
var END_MONTH_EXCLUSIVE = 10;  // Oct 1 means May-Sep or Apr-Sep.
var SCALE = 1000;
var DRIVE_FOLDER = 'YJQ_HARSEI_annual_inputs_2000_2024';
var NDWI_WATER_THRESHOLD = 0.2;
var PCA_SAMPLE_PER_YEAR = 8000;

var EXPORT_ECO_COMPONENTS = true;
var EXPORT_POOLED_SAMPLE = true;
var EXPORT_DIAGNOSTIC_BANDS = true;

// Optional HAI inputs. The dedicated HAI export script gee_02 is preferred.
var EXPORT_HAI_INPUTS = false;
// Example only. Set these after uploading CLCD rasters as GEE image assets.
// If your assets are projects/jlu2024621038/assets/CLCD/CLCD_v01_2000,
// use CLCD_ASSET_PREFIX = 'projects/jlu2024621038/assets/CLCD/CLCD_v01_'
// and CLCD_ASSET_SUFFIX = ''.
var CLCD_ASSET_PREFIX = '';
var CLCD_ASSET_SUFFIX = '';
var GLC_FCS30D_LAST_YEAR = 2022;
var WORLDPOP_LAST_YEAR = 2020;

var ECO_BANDS = ['NDVI', 'WET', 'NDBSI', 'LST', 'SRSI'];

Map.centerObject(roiFc, 6);
Map.addLayer(roiFc, {color: 'red'}, 'ROI');

function dateStart(year) {
  return ee.Date.fromYMD(year, START_MONTH, 1);
}

function dateEnd(year) {
  return ee.Date.fromYMD(year, END_MONTH_EXCLUSIVE, 1);
}

function emptyBand(name) {
  return ee.Image.constant(0).rename(name).updateMask(ee.Image.constant(0));
}

function prepMod13(img) {
  var good = img.select('SummaryQA').lte(1);
  return img.select('NDVI')
    .multiply(0.0001)
    .rename('NDVI')
    .updateMask(good)
    .copyProperties(img, ['system:time_start']);
}

function prepMod09(img) {
  var qa = img.select('QA');
  var state = img.select('StateQA');
  // MODLAND QA bits 0-1: 0 ideal, 1 less than ideal but usable.
  var modlandOk = qa.bitwiseAnd(3).lte(1);
  // StateQA cloud state bits 0-1, cloud shadow bit 2, cirrus bits 8-9,
  // internal snow bit 12, and MOD35 snow/ice bit 15.
  var cloudClear = state.bitwiseAnd(3).eq(0);
  var noCloudShadow = state.rightShift(2).bitwiseAnd(1).eq(0);
  var cirrusOk = state.rightShift(8).bitwiseAnd(3).lte(1);
  var noInternalSnow = state.rightShift(12).bitwiseAnd(1).eq(0);
  var noMod35Snow = state.rightShift(15).bitwiseAnd(1).eq(0);
  var good = modlandOk
    .and(cloudClear)
    .and(noCloudShadow)
    .and(cirrusOk)
    .and(noInternalSnow)
    .and(noMod35Snow);
  var bands = [
    'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b04',
    'sur_refl_b05', 'sur_refl_b06', 'sur_refl_b07'
  ];
  var renamed = ['RED', 'NIR', 'BLUE', 'GREEN', 'SWIR1', 'SWIR2', 'SWIR3'];
  return img.select(bands, renamed)
    .multiply(0.0001)
    .updateMask(good)
    .copyProperties(img, ['system:time_start']);
}

function prepMod11(img) {
  var qc = img.select('QC_Day');
  var mandatoryGood = qc.bitwiseAnd(3).lte(1);
  var lstErrorOk = qc.rightShift(6).bitwiseAnd(3).lte(1);
  return img.select('LST_Day_1km')
    .multiply(0.02)
    .subtract(273.15)
    .rename('LST')
    .updateMask(mandatoryGood.and(lstErrorOk))
    .copyProperties(img, ['system:time_start']);
}

function annualNdvi(year) {
  var col = ee.ImageCollection('MODIS/061/MOD13A1')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year));
  return col.map(prepMod13).mean().rename('NDVI');
}

function annualSnowFraction(year) {
  var col = ee.ImageCollection('MODIS/061/MOD13A1')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year));
  return col.map(function(img) {
    return img.select('SummaryQA').eq(2).rename('snow');
  }).mean().rename('snow_frac');
}

function annualMod09SnowFraction(year) {
  var col = ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year));
  return col.map(function(img) {
    var state = img.select('StateQA');
    var snow = state.rightShift(12).bitwiseAnd(1).eq(1)
      .or(state.rightShift(15).bitwiseAnd(1).eq(1));
    return snow.rename('mod09_snow');
  }).mean().rename('mod09_snow_frac');
}

function annualReflectance(year) {
  return ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year))
    .map(prepMod09)
    .mean();
}

function annualLst(year) {
  return ee.ImageCollection('MODIS/061/MOD11A2')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year))
    .map(prepMod11)
    .mean()
    .rename('LST');
}

function computeNdbsi(sr) {
  var red = sr.select('RED');
  var nir = sr.select('NIR');
  var blue = sr.select('BLUE');
  var green = sr.select('GREEN');
  var swir1 = sr.select('SWIR1');

  var si = swir1.add(red).subtract(nir.add(blue))
    .divide(swir1.add(red).add(nir).add(blue))
    .rename('SI');

  var ibiNumerator = swir1.multiply(2).divide(swir1.add(nir))
    .subtract(nir.divide(nir.add(red)).add(green.divide(green.add(swir1))));
  var ibiDenominator = swir1.multiply(2).divide(swir1.add(nir))
    .add(nir.divide(nir.add(red)).add(green.divide(green.add(swir1))));
  var ibi = ibiNumerator.divide(ibiDenominator).rename('IBI');

  return si.add(ibi).divide(2).rename('NDBSI');
}

function computeWet(sr) {
  // MODIS tasseled-cap wetness coefficients from Lobser and Cohen (2007).
  // Band order: RED, NIR, BLUE, GREEN, SWIR1, SWIR2, SWIR3.
  return sr.select('RED').multiply(0.1147)
    .add(sr.select('NIR').multiply(0.2489))
    .add(sr.select('BLUE').multiply(0.2408))
    .add(sr.select('GREEN').multiply(0.3132))
    .add(sr.select('SWIR1').multiply(-0.3122))
    .add(sr.select('SWIR2').multiply(-0.6416))
    .add(sr.select('SWIR3').multiply(-0.5087))
    .rename('WET');
}

function computeSalinityDiagnostics(sr) {
  var red = sr.select('RED');
  var blue = sr.select('BLUE');
  var green = sr.select('GREEN');
  var nir = sr.select('NIR');
  var swir1 = sr.select('SWIR1');
  var eps = ee.Image.constant(1e-6);

  var ndviSr = nir.subtract(red).divide(nir.add(red)).rename('NDVI_SR');
  var si1 = green.multiply(red).sqrt().rename('SI1');
  var ndsiSal = red.subtract(nir).divide(red.add(nir)).rename('NDSI_SAL');
  var si3Grb = green.multiply(red).divide(blue.max(eps)).rename('SI3_GRB');
  var siS1 = swir1.subtract(blue).divide(swir1.add(blue)).rename('SI_SWIR_BLUE');

  // Candidate SRSI retained from the submitted manuscript:
  // SRSI = sqrt((NDVI_SR - 1)^2 + SI1^2), SI1 = sqrt(green * red).
  // The manuscript must cite its original source or replace this band with one
  // of the standard salinity diagnostics exported below.
  var srsi = ndviSr.subtract(1).pow(2).add(si1.pow(2)).sqrt().rename('SRSI');
  return ee.Image.cat([srsi, si1, ndsiSal, si3Grb, siS1, ndviSr]);
}

function annualEcoComponents(year) {
  year = ee.Number(year).toInt();
  var sr = annualReflectance(year);
  var ndvi = annualNdvi(year);
  var lst = annualLst(year);
  var wet = computeWet(sr);
  var ndbsi = computeNdbsi(sr);
  var sal = computeSalinityDiagnostics(sr);
  var srsi = sal.select('SRSI');
  var ndwi = sr.select('GREEN').subtract(sr.select('NIR'))
    .divide(sr.select('GREEN').add(sr.select('NIR')))
    .rename('NDWI');
  var waterMask = ndwi.lte(NDWI_WATER_THRESHOLD).rename('land_mask');
  var snowFrac = annualSnowFraction(year);
  var mod09SnowFrac = annualMod09SnowFraction(year);

  var base = ee.Image.cat([ndvi, wet, ndbsi, lst, srsi]);
  var diagnostics = ee.Image.cat([
    ndwi,
    waterMask,
    snowFrac,
    mod09SnowFrac,
    sal.select(['SI1', 'NDSI_SAL', 'SI3_GRB', 'SI_SWIR_BLUE', 'NDVI_SR'])
  ]);
  var eco = (EXPORT_DIAGNOSTIC_BANDS ? base.addBands(diagnostics) : base)
    .updateMask(waterMask)
    .clip(roi)
    .toFloat()
    .set('year', year)
    .set('start_month', START_MONTH)
    .set('end_month_exclusive', END_MONTH_EXCLUSIVE)
    .set('ndwi_water_threshold', NDWI_WATER_THRESHOLD)
    .set('wet_reference', 'Lobser and Cohen 2007 MODIS tasseled-cap wetness')
    .set('lst_source', 'MODIS/061/MOD11A2 LST_Day_1km scale 0.02 K, converted to Celsius')
    .set('srsi_formula', 'sqrt((NDVI_SR - 1)^2 + SI1^2); SI1=sqrt(green*red); verify manuscript citation');

  return eco;
}

function annualPop(year) {
  var popYear = Math.min(year, WORLDPOP_LAST_YEAR);
  return ee.ImageCollection('WorldPop/GP/100m/pop')
    .filter(ee.Filter.eq('year', popYear))
    .mosaic()
    .select('population')
    .rename('POP')
    .toFloat();
}

function annualLights(year) {
  var dmsp = (year <= 2013) ?
    ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')
      .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
      .select('stable_lights')
      .mean()
      .rename('LIGHT_DMSP') :
    emptyBand('LIGHT_DMSP');

  // Use VCMCFG instead of VCMSLCFG because VCMCFG starts in 2012 and overlaps
  // with DMSP-OLS in 2012-2013 for inter-sensor harmonization.
  var viirs = (year >= 2012) ?
    annualViirsWeighted(year) :
    emptyBand('LIGHT_VIIRS');

  var raw = (year <= 2013) ? dmsp.rename('LIGHT_RAW') : viirs.rename('LIGHT_RAW');
  return ee.Image.cat([raw, dmsp, viirs]).toFloat();
}

function annualViirsWeighted(year) {
  var col = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG')
    .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
    .map(function(img) {
      var rad = img.select('avg_rad').max(ee.Image.constant(0));
      var cov = img.select('cf_cvg');
      return rad.multiply(cov).rename('rad_x_cov')
        .addBands(cov.rename('cf_cvg'))
        .updateMask(cov.gt(0));
    });
  return col.select('rad_x_cov').sum()
    .divide(col.select('cf_cvg').sum())
    .rename('LIGHT_VIIRS');
}

function annualClcdScore(year) {
  if (CLCD_ASSET_PREFIX === '') {
    return emptyBand('LUCC_SCORE');
  }
  var clcd = ee.Image(CLCD_ASSET_PREFIX + year + CLCD_ASSET_SUFFIX).rename('CLCD');
  // Common CLCD codes: 1 cropland, 2 forest, 3 shrub, 4 grassland,
  // 5 water, 6 snow/ice, 7 barren, 8 impervious, 9 wetland.
  // HAI land-use score follows the manuscript's 0-10 rule:
  // built-up/impervious = 10, cropland = 7, grassland = 4, water = 1,
  // and all remaining land-cover types = 0.
  return clcd.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9],
      [7, 0, 0, 4, 1, 0, 0, 10, 0],
      0
    )
    .rename('LUCC_SCORE')
    .toFloat();
}

function annualGlcFcs30dScore(year) {
  var lcYear = Math.min(year, GLC_FCS30D_LAST_YEAR);
  var bandIndex = lcYear - 1999;
  var lc = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/annual')
    .mosaic()
    .select('b' + bandIndex)
    .rename('LUCC_CLASS');
  return lc.remap(
      [190, 10, 11, 12, 130, 210],
      [10, 7, 7, 7, 4, 1],
      0
    )
    .rename('LUCC_SCORE')
    .toFloat();
}

function annualHaiInputs(year) {
  var eco = annualEcoComponents(year);
  var mask = eco.select('land_mask');
  var popYear = Math.min(year, WORLDPOP_LAST_YEAR);
  return ee.Image.cat([
      annualPop(year),
      annualLights(year),
      (CLCD_ASSET_PREFIX === '' ? annualGlcFcs30dScore(year) : annualClcdScore(year)),
      ee.Image.constant(popYear).rename('POP_YEAR_USED'),
      ee.Image.constant(CLCD_ASSET_PREFIX === '' ? Math.min(year, GLC_FCS30D_LAST_YEAR) : year)
        .rename('LUCC_YEAR_USED')
    ])
    .updateMask(mask)
    .clip(roi)
    .toFloat()
    .set('year', year);
}

if (EXPORT_ECO_COMPONENTS) {
  for (var year = START_YEAR; year <= END_YEAR; year++) {
    var eco = annualEcoComponents(year);
    Export.image.toDrive({
      image: eco,
      description: 'YJQ_ecocomponents_' + year,
      folder: DRIVE_FOLDER,
      fileNamePrefix: 'YJQ_ecocomponents_' + year,
      region: roi,
      scale: SCALE,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF'
    });
  }
}

if (EXPORT_HAI_INPUTS) {
  for (var hy = START_YEAR; hy <= END_YEAR; hy++) {
    var hai = annualHaiInputs(hy);
    Export.image.toDrive({
      image: hai,
      description: 'YJQ_hai_inputs_' + hy,
      folder: DRIVE_FOLDER,
      fileNamePrefix: 'YJQ_hai_inputs_' + hy,
      region: roi,
      scale: SCALE,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF'
    });
  }
}

if (EXPORT_POOLED_SAMPLE) {
  var pooled = ee.FeatureCollection([]);
  for (var sy = START_YEAR; sy <= END_YEAR; sy++) {
    var sampleImg = annualEcoComponents(sy).select(ECO_BANDS);
    var sample = sampleImg.sample({
      region: roi,
      scale: SCALE,
      numPixels: PCA_SAMPLE_PER_YEAR,
      seed: sy,
      tileScale: 4,
      geometries: false
    }).map(function(f) {
      return f.set('year', sy);
    });
    pooled = pooled.merge(sample);
  }

  Export.table.toDrive({
    collection: pooled,
    description: 'YJQ_pooled_pca_sample_2000_2024',
    folder: DRIVE_FOLDER,
    fileNamePrefix: 'YJQ_pooled_pca_sample_2000_2024',
    fileFormat: 'CSV'
  });
}

// Quick map preview for one year.
var previewYear = 2024;
var preview = annualEcoComponents(previewYear);
Map.addLayer(preview.select('NDVI'), {min: 0, max: 0.8, palette: ['white', 'green']}, 'NDVI ' + previewYear);
Map.addLayer(preview.select('SRSI'), {min: 0, max: 1.5, palette: ['green', 'yellow', 'red']}, 'SRSI ' + previewYear);
