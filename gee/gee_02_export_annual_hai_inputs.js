// Annual HAI input export for YJQ / Altai study area.
// Copy this file into the Google Earth Engine Code Editor and run it after or
// alongside gee_01_export_annual_inputs.js.
//
// Exports:
//   YJQ_hai_inputs_2000.tif ... YJQ_hai_inputs_2024.tif
//
// First five bands are read by local_02_build_harsei_from_downloads.py:
//   1 POP
//   2 LIGHT_RAW
//   3 LIGHT_DMSP
//   4 LIGHT_VIIRS
//   5 LUCC_SCORE
//
// Default land-cover mode follows the submitted HAI land-use code:
// SAT-IO GLC_FCS30D annual land cover, scored on a 0-10 human-activity scale.
// If you upload CLCD as annual GEE assets, switch LANDCOVER_MODE to 'CLCD_ASSET'.

var ROI_ASSET = 'projects/jlu2024621038/assets/yjq';
var roiFc = ee.FeatureCollection(ROI_ASSET);
var roi = roiFc.geometry();

var START_YEAR = 2000;
var END_YEAR = 2024;
var SCALE = 1000;
var DRIVE_FOLDER = 'YJQ_HARSEI_annual_inputs_2000_2024';
var START_MONTH = 5;
var END_MONTH_EXCLUSIVE = 10;
var NDWI_WATER_THRESHOLD = 0.2;
var WORLDPOP_LAST_YEAR = 2020;

// 'GLC_FCS30D', 'CLCD_ASSET', or 'MCD12Q1'.
var LANDCOVER_MODE = 'GLC_FCS30D';

// GLC_FCS30D annual bands b1-b23 correspond to 2000-2022 in the SAT-IO asset.
// If no 2024 layer is provided, 2024 uses 2022 as the nearest available layer.
var GLC_FCS30D_LAST_YEAR = 2022;

// Use only if LANDCOVER_MODE = 'CLCD_ASSET'.
// Example: projects/jlu2024621038/assets/CLCD/CLCD_v01_2000
var CLCD_ASSET_PREFIX = 'projects/jlu2024621038/assets/CLCD/CLCD_v01_';
var CLCD_ASSET_SUFFIX = '';
var CLCD_LAST_YEAR = 2024;  // Current CLCD Zenodo latest record includes 2024.

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

function prepMod09(img) {
  var qa = img.select('QA');
  var state = img.select('StateQA');
  var modlandOk = qa.bitwiseAnd(3).lte(1);
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
  var bands = ['sur_refl_b02', 'sur_refl_b04'];  // NIR, GREEN
  var renamed = ['NIR', 'GREEN'];
  return img.select(bands, renamed)
    .multiply(0.0001)
    .updateMask(good)
    .copyProperties(img, ['system:time_start']);
}

function annualLandMask(year) {
  var sr = ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(dateStart(year), dateEnd(year))
    .map(prepMod09)
    .mean();
  var ndwi = sr.select('GREEN').subtract(sr.select('NIR'))
    .divide(sr.select('GREEN').add(sr.select('NIR')));
  return ndwi.lte(NDWI_WATER_THRESHOLD).rename('land_mask');
}

function annualPop(year) {
  // WorldPop/GP/100m/pop availability in the current GEE export is through 2020.
  // Carry 2020 forward for 2021-2024 and report this explicitly in the revision.
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

  // VCMCFG starts in 2012, giving 2012-2013 overlap with DMSP-OLS. The local
  // script harmonizes VIIRS to DMSP-like units using this overlap.
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

function mcd12q1Score(year) {
  // MCD12Q1 official availability is 2001-2024. Use 2001 for year 2000.
  var lcYear = Math.max(year, 2001);
  var lc = ee.ImageCollection('MODIS/061/MCD12Q1')
    .filterDate(ee.Date.fromYMD(lcYear, 1, 1), ee.Date.fromYMD(lcYear + 1, 1, 1))
    .first()
    .select('LC_Type1')
    .rename('LUCC_CLASS');

  // IGBP to HAI land-use score, following the manuscript's 0-10 rule:
  // built-up = 10, cropland = 7, grassland = 4, water = 1, other classes = 0.
  var score = lc.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 7, 10, 7, 0, 0, 1],
      0
    )
    .rename('LUCC_SCORE')
    .toFloat();
  return ee.Image.cat([
    score,
    lc.toFloat(),
    ee.Image.constant(lcYear).rename('LUCC_YEAR_USED')
  ]);
}

function glcFcs30dScore(year) {
  var lcYear = Math.min(year, GLC_FCS30D_LAST_YEAR);
  var bandIndex = lcYear - 1999;
  var lc = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/annual')
    .mosaic()
    .select('b' + bandIndex)
    .rename('LUCC_CLASS');

  // GLC_FCS30D to HAI land-use score, following the submitted code:
  // built-up (190) = 10; cropland (10, 11, 12) = 7;
  // grassland (130) = 4; water (210) = 1; all remaining classes = 0.
  var score = lc.remap(
      [190, 10, 11, 12, 130, 210],
      [10, 7, 7, 7, 4, 1],
      0
    )
    .rename('LUCC_SCORE')
    .toFloat();
  return ee.Image.cat([
    score,
    lc.toFloat(),
    ee.Image.constant(lcYear).rename('LUCC_YEAR_USED')
  ]);
}

function clcdScore(year) {
  var lcYear = Math.min(year, CLCD_LAST_YEAR);
  var clcd = ee.Image(CLCD_ASSET_PREFIX + lcYear + CLCD_ASSET_SUFFIX).rename('LUCC_CLASS');
  // Common CLCD codes: 1 cropland, 2 forest, 3 shrub, 4 grassland,
  // 5 water, 6 snow/ice, 7 barren, 8 impervious, 9 wetland.
  // HAI land-use score follows the manuscript's 0-10 rule:
  // built-up/impervious = 10, cropland = 7, grassland = 4, water = 1,
  // and all remaining land-cover types = 0.
  var score = clcd.remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9],
      [7, 0, 0, 4, 1, 0, 0, 10, 0],
      0
    )
    .rename('LUCC_SCORE')
    .toFloat();
  return ee.Image.cat([
    score,
    clcd.toFloat(),
    ee.Image.constant(lcYear).rename('LUCC_YEAR_USED')
  ]);
}

function annualLucc(year) {
  if (LANDCOVER_MODE === 'CLCD_ASSET') {
    return clcdScore(year);
  }
  if (LANDCOVER_MODE === 'MCD12Q1') {
    return mcd12q1Score(year);
  }
  return glcFcs30dScore(year);
}

function annualHaiInputs(year) {
  var mask = annualLandMask(year);
  var lucc = annualLucc(year);
  var popYear = Math.min(year, WORLDPOP_LAST_YEAR);
  return ee.Image.cat([
      annualPop(year),
      annualLights(year),
      lucc.select('LUCC_SCORE'),
      lucc.select('LUCC_CLASS'),
      mask,
      ee.Image.constant(popYear).rename('POP_YEAR_USED'),
      lucc.select('LUCC_YEAR_USED')
    ])
    .updateMask(mask)
    .clip(roi)
    .toFloat()
    .set('year', year)
    .set('landcover_mode', LANDCOVER_MODE);
}

for (var year = START_YEAR; year <= END_YEAR; year++) {
  var hai = annualHaiInputs(year);
  Export.image.toDrive({
    image: hai,
    description: 'YJQ_hai_inputs_' + year,
    folder: DRIVE_FOLDER,
    fileNamePrefix: 'YJQ_hai_inputs_' + year,
    region: roi,
    scale: SCALE,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

var previewYear = 2024;
var preview = annualHaiInputs(previewYear);
Map.addLayer(preview.select('POP'), {min: 0, max: 20, palette: ['white', 'orange', 'red']}, 'POP ' + previewYear);
Map.addLayer(preview.select('LIGHT_RAW'), {min: 0, max: 10, palette: ['black', 'yellow', 'white']}, 'LIGHT ' + previewYear);
Map.addLayer(preview.select('LUCC_SCORE'), {min: 0, max: 10, palette: ['green', 'yellow', 'red']}, 'LUCC_SCORE ' + previewYear);
