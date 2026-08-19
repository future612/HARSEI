// NDWI threshold and snow-contamination sensitivity for the annual HARSEI rebuild.
//
// Copy this file into the Google Earth Engine Code Editor and run it after or
// alongside gee_01_export_annual_inputs.js. It exports one CSV table:
//   YJQ_water_snow_sensitivity_2000_2024.csv
//
// The table compares:
//   - Apr-Sep versus May-Sep compositing windows
//   - NDWI water thresholds 0.0, 0.1, 0.2, 0.3
//   - all-area versus high-elevation MOD09A1 snow frequency

var ROI_ASSET = 'projects/jlu2024621038/assets/yjq';
var roiFc = ee.FeatureCollection(ROI_ASSET);
var roi = roiFc.geometry();

var START_YEAR = 2000;
var END_YEAR = 2024;
var START_MONTH_OPTIONS = [4, 5];
var END_MONTH_EXCLUSIVE = 10;
var NDWI_THRESHOLDS = [0.0, 0.1, 0.2, 0.3];
var HIGH_ELEVATION_M = 2500;
var SCALE = 1000;
var DRIVE_FOLDER = 'YJQ_HARSEI_annual_inputs_2000_2024';

var areaKm2 = ee.Image.pixelArea().divide(1e6).rename('area_km2');
var highElev = ee.Image('USGS/SRTMGL1_003')
  .select('elevation')
  .gte(HIGH_ELEVATION_M)
  .rename('high_elevation')
  .clip(roi);

Map.centerObject(roiFc, 6);
Map.addLayer(roiFc, {color: 'red'}, 'ROI');
Map.addLayer(highElev.updateMask(highElev), {palette: ['white']}, 'high elevation');

function dateStart(year, startMonth) {
  return ee.Date.fromYMD(year, startMonth, 1);
}

function dateEnd(year) {
  return ee.Date.fromYMD(year, END_MONTH_EXCLUSIVE, 1);
}

function cloudMaskNoSnow(img) {
  var qa = img.select('QA');
  var state = img.select('StateQA');
  var modlandOk = qa.bitwiseAnd(3).lte(1);
  var cloudClear = state.bitwiseAnd(3).eq(0);
  var noCloudShadow = state.rightShift(2).bitwiseAnd(1).eq(0);
  var cirrusOk = state.rightShift(8).bitwiseAnd(3).lte(1);
  return modlandOk.and(cloudClear).and(noCloudShadow).and(cirrusOk);
}

function prepNdwiReflectance(img) {
  return img.select(['sur_refl_b02', 'sur_refl_b04'], ['NIR', 'GREEN'])
    .multiply(0.0001)
    .updateMask(cloudMaskNoSnow(img))
    .copyProperties(img, ['system:time_start']);
}

function annualNdwi(year, startMonth) {
  var sr = ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(dateStart(year, startMonth), dateEnd(year))
    .map(prepNdwiReflectance)
    .mean();
  return sr.select('GREEN').subtract(sr.select('NIR'))
    .divide(sr.select('GREEN').add(sr.select('NIR')))
    .rename('NDWI')
    .clip(roi);
}

function annualSnowFraction(year, startMonth) {
  var col = ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(dateStart(year, startMonth), dateEnd(year));
  return col.map(function(img) {
    var state = img.select('StateQA');
    var snow = state.rightShift(12).bitwiseAnd(1).eq(1)
      .or(state.rightShift(15).bitwiseAnd(1).eq(1));
    return snow.rename('snow');
  }).mean().rename('mod09_snow_frac').clip(roi);
}

function getRegionValue(img, reducer, bandName) {
  return img.reduceRegion({
    reducer: reducer,
    geometry: roi,
    scale: SCALE,
    maxPixels: 1e13,
    tileScale: 4
  }).get(bandName);
}

function sensitivityFeature(year, startMonth, threshold) {
  var ndwi = annualNdwi(year, startMonth);
  var snowFrac = annualSnowFraction(year, startMonth);
  var waterMask = ndwi.gt(threshold);
  var landMask = ndwi.lte(threshold);

  var ndwiStats = ndwi.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: roi,
    scale: SCALE,
    maxPixels: 1e13,
    tileScale: 4
  });

  var props = ee.Dictionary({
    year: year,
    start_month: startMonth,
    end_month_exclusive: END_MONTH_EXCLUSIVE,
    ndwi_threshold: threshold,
    high_elevation_m: HIGH_ELEVATION_M,
    ndwi_mean: ndwiStats.get('NDWI_mean'),
    ndwi_stdDev: ndwiStats.get('NDWI_stdDev'),
    land_area_km2: getRegionValue(areaKm2.updateMask(landMask), ee.Reducer.sum(), 'area_km2'),
    water_area_km2: getRegionValue(areaKm2.updateMask(waterMask), ee.Reducer.sum(), 'area_km2'),
    high_elevation_area_km2: getRegionValue(areaKm2.updateMask(highElev), ee.Reducer.sum(), 'area_km2'),
    snow_fraction_all: getRegionValue(snowFrac, ee.Reducer.mean(), 'mod09_snow_frac'),
    snow_fraction_high_elevation: getRegionValue(
      snowFrac.updateMask(highElev), ee.Reducer.mean(), 'mod09_snow_frac'
    ),
    snow_area_fraction_gt_0_1: getRegionValue(
      snowFrac.gt(0.1), ee.Reducer.mean(), 'mod09_snow_frac'
    )
  });

  return ee.Feature(null, props);
}

var rows = ee.FeatureCollection([]);
for (var year = START_YEAR; year <= END_YEAR; year++) {
  for (var smi = 0; smi < START_MONTH_OPTIONS.length; smi++) {
    for (var ti = 0; ti < NDWI_THRESHOLDS.length; ti++) {
      rows = rows.merge(ee.FeatureCollection([
        sensitivityFeature(year, START_MONTH_OPTIONS[smi], NDWI_THRESHOLDS[ti])
      ]));
    }
  }
}

Export.table.toDrive({
  collection: rows,
  description: 'YJQ_water_snow_sensitivity_2000_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'YJQ_water_snow_sensitivity_2000_2024',
  fileFormat: 'CSV'
});

// Quick visual check for the most recent year and the submitted NDWI threshold.
var previewYear = 2024;
var previewNdwi = annualNdwi(previewYear, 5);
Map.addLayer(previewNdwi, {min: -0.5, max: 0.5, palette: ['brown', 'white', 'blue']}, 'NDWI May-Sep 2024');
Map.addLayer(previewNdwi.gt(0.2).selfMask(), {palette: ['blue']}, 'NDWI > 0.2 water 2024');
