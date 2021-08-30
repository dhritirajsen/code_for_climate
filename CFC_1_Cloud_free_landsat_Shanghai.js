// Create a cloud-free composite with default parameters.
var rawCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1") // L5- LANDSAT/LT05/C01/T1
  .filterBounds(geometry)
  .sort("CLOUD_COVER")
  .filterDate('2019-01-01', '2019-12-31');
var composite = ee.Algorithms.Landsat.simpleComposite(rawCollection);
//change bands to 5,3,2 for L8
Map.addLayer(composite.clip(geometry), {bands: ['B5', 'B3', 'B2'], max: 100}, 'Shanghai 2019', true);

// Create a cloud-free composite with default parameters.
var rawCollection = ee.ImageCollection("LANDSAT/LT05/C01/T1") // L5- LANDSAT/LT05/C01/T1
  .filterBounds(geometry)
  .sort("CLOUD_COVER")
  .filterDate('1990-01-01', '1990-12-31');
var composite5 = ee.Algorithms.Landsat.simpleComposite(rawCollection);
//change bands to 5,3,2 for L8
Map.addLayer(composite5.clip(geometry), {bands: ['B5', 'B3', 'B2'], max: 100}, 'Shanghai 1990', true);