// Export irrigation-zone masks for the YJQ pressure-zone validation.
//
// Run this in the Google Earth Engine Code Editor. After the Drive export
// finishes, download the file and place it at:
// D:\Codex\260724 小论文\revise\external_validation_inputs\pressure_zone_masks\YJQ_irrigated_area_mask.tif

// -------------------- 1. Study Area --------------------
var roi = ee.FeatureCollection('projects/jlu2024621038/assets/yjq').geometry();
Map.centerObject(roi, 6);

// -------------------- 2. GFSAD LGRIP Irrigated Cropland --------------------
// Class labels in projects/sat-io/open-datasets/GFSAD/LGRIP30:
// 0 Ocean / water bodies
// 1 Non-croplands
// 2 Irrigated croplands
// 3 Rainfed croplands
var lgripCol = ee.ImageCollection('projects/sat-io/open-datasets/GFSAD/LGRIP30');
var lgripFirst = ee.Image(lgripCol.first());
var lgripProj = lgripFirst.projection();

// ImageCollection.mosaic() can lose the default projection required by
// reduceResolution(), so we explicitly restore the native LGRIP projection.
var lgrip30 = lgripCol.mosaic()
  .rename('lgrip_class')
  .setDefaultProjection(lgripProj);

var irrigated30 = lgrip30.eq(2)
  .rename('irrigated_30m')
  .setDefaultProjection(lgripProj);

print('LGRIP first image', lgripFirst);
print('LGRIP native projection', lgripProj);

// Aggregate 30 m irrigated pixels to an approximately 1 km fraction.
// The binary mask keeps pixels where at least 25% of the 1 km cell is irrigated
// cropland. If this is too strict/loose, export the fraction raster as well and
// adjust the threshold locally.
var irrigatedFraction = irrigated30
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 4096
  })
  .reproject({
    crs: 'EPSG:4326',
    scale: 1000
  })
  .clip(roi)
  .rename('irrigated_fraction');

print('Irrigated fraction projection', irrigatedFraction.projection());

var irrigatedMask = irrigatedFraction.gte(0.25)
  .rename('irrigated_mask')
  .toUint8()
  .clip(roi);

// Optional cross-check: Deepak Nagaraj global irrigation maps, available for
// 2001-2015 in the GEE community catalog; value 2 marks highly irrigated areas.
var irrigationMaps = ee.ImageCollection('users/deepakna/global_irrigation_maps');
var deepak2010 = ee.Image(irrigationMaps
  .filter(ee.Filter.date('2010-01-01', '2010-12-31'))
  .first())
  .eq(2)
  .rename('deepak_highly_irrigated_2010')
  .toUint8()
  .clip(roi);

Map.addLayer(irrigatedMask.selfMask(), {palette: ['#1f78b4']}, 'GFSAD irrigated mask >=25%');
Map.addLayer(deepak2010.selfMask(), {palette: ['#33a02c']}, 'Deepak highly irrigated 2010');
Map.addLayer(roi, {}, 'ROI');

// -------------------- 3. Exports --------------------
Export.image.toDrive({
  image: irrigatedMask,
  description: 'YJQ_irrigated_area_mask',
  folder: 'YJQ_validation_zone_masks',
  fileNamePrefix: 'YJQ_irrigated_area_mask',
  region: roi,
  crs: 'EPSG:4326',
  scale: 1000,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: irrigatedFraction.toFloat(),
  description: 'YJQ_irrigated_fraction_1000m',
  folder: 'YJQ_validation_zone_masks',
  fileNamePrefix: 'YJQ_irrigated_fraction_1000m',
  region: roi,
  crs: 'EPSG:4326',
  scale: 1000,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: deepak2010,
  description: 'YJQ_deepak_highly_irrigated_2010_mask',
  folder: 'YJQ_validation_zone_masks',
  fileNamePrefix: 'YJQ_deepak_highly_irrigated_2010_mask',
  region: roi,
  crs: 'EPSG:4326',
  scale: 1000,
  maxPixels: 1e13
});
