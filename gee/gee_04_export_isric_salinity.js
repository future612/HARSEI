// Export ISRIC / Ivushkin et al. Global Soil Salinity Maps for external validation.
//
// Dataset:
//   ee.ImageCollection("projects/sat-io/open-datasets/global_soil_salinity")
// Years used for validation:
//   2000, 2002, 2005, 2009, 2016
//
// After export, download the GeoTIFFs from Google Drive and place them in:
//   D:\Codex\260724 小论文\revise\external_validation_inputs\isric_global_soil_salinity

// -------------------- 1. ROI --------------------
var roi = ee.FeatureCollection("projects/jlu2024621038/assets/yjq");
Map.centerObject(roi, 6);
Map.addLayer(roi, {}, "YJQ ROI");

// -------------------- 2. ISRIC Global Soil Salinity --------------------
var salCol = ee.ImageCollection("projects/sat-io/open-datasets/global_soil_salinity");

print("ISRIC Global Soil Salinity collection size", salCol.size());
print("First image metadata", salCol.first());
print("Available system:index values", salCol.aggregate_array("system:index"));
print("Available time_start values", salCol.aggregate_array("system:time_start"));

var years = [2000, 2002, 2005, 2009, 2016];
var driveFolder = "YJQ_ISRIC_salinity_validation";

function getSalinityImage(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, "year");
  var yearText = year.format();

  var byDate = salCol.filterDate(start, end);
  var byIndex = salCol.filter(ee.Filter.stringContains("system:index", yearText));
  var selected = ee.Image(ee.Algorithms.If(byDate.size().gt(0), byDate.first(), byIndex.first()));

  return selected
    .select([0])
    .rename("isric_salinity")
    .toFloat()
    .clip(roi)
    .set("year", year);
}

years.forEach(function(year) {
  var img = getSalinityImage(year);
  print("Selected salinity image " + year, img);
  Map.addLayer(img, {min: 0, max: 5, palette: ["2c7bb6", "abd9e9", "ffffbf", "fdae61", "d7191c"]}, "ISRIC salinity " + year, false);

  Export.image.toDrive({
    image: img,
    description: "ISRIC_global_soil_salinity_" + year + "_YJQ",
    folder: driveFolder,
    fileNamePrefix: "ISRIC_global_soil_salinity_" + year + "_YJQ",
    region: roi.geometry(),
    scale: 250,
    maxPixels: 1e13
  });
});

