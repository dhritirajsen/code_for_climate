/**
 * Examples of shoreline extraction
 *  
 * Quick example of some the approaches for detecting shorelines
 * in coastal areas
 *
 * @author Edward P. Morris (UCA, www.fast-space-project.eu) 
 * 
 */

Map.centerObject(table, 12);

////////////////////////////////////////////////////////////////////
// Digital elevation models
////////////////////////////////////////////////////////////////////
 
// Elevation viualisation
var imageVisParam = {"opacity":1,"min":-50,"max":30,"palette":['#0000ff','#6038ff','#875bff','#a67cff','#bf9cff','#d7bdff','#ecddff','#ffffe0','#dfedc0','#c1dba2','#a2c985','#83b767','#62a44a','#3f922c','#008000']};

// Elevation images
var gmted2010 = ee.Image("USGS/GMTED2010"),
    srtm30 = ee.Image("JAXA/ALOS/AW3D30_V1_1").select('AVE'),
    etopo1 = ee.Image("NOAA/NGDC/ETOPO1");

// Coastal elevation
var maxElev = ee.Number(5);
var minElev = ee.Number(-10);
var csrtm30 = srtm30.select(['AVE'],['e'])
.updateMask(srtm30.lte(maxElev).and(srtm30.gt(minElev)));
var elevMask = ee.Image.constant(1).updateMask(csrtm30.gt(0));

Map.addLayer(elevMask,imageVisParam, 'SRTM30');

////////////////////////////////////////////////////////////////////
// Optical OTSU segmentation
////////////////////////////////////////////////////////////////////

// Define process arguments
// temporal range
var sdate = ee.String('2015-01-01');
var edate = ee.String('2015-12-31');
// cloud cover (metadata)
var max_cloud_cover = ee.Number(30);
// per-pixel cloud exclusion; fixme this seems odd?
var cs_pct = ee.Number(90);

/**
 * Return system time of each image in collection
 *
 * @param {ee.ImageCollection} coll The collection to extract timestamps
 * 
 * @return {ee.List} A list of ee.Features, for export wrap this in an 
 *  ee.FeatureCollection
 * 
 * @export
 */
function getSysTime(coll) {
    return coll.toList(coll.size(), 0).map(function(im) {
      // Get the image system:time_start'
      return ee.Feature(null, 
        {'timestamp': ee.Image(im).get('system:time_start')});
    });
}

/**
  * Compute a cloud score for Sentinel 2.
  * 
  * This expects the input image to have the common
  * band names: ["R", "B", etc].
  * from Chris Herwig <herwig@google.com>
  * https://code.earthengine.google.com/3890321b866c5a24b7fde43256424885
  * 
  * @param {ee.Image} img A sentinel 2 image with bands B, G, R, NIR, SWIR1 and SWIR2.
  * 
  * @return An ee.Image where each pixel represents a potential cloud score, 
  *     higher numbers represent more likely a cloud. Note this is a simple indication
  *     prone to false positives.
  * @export
  */
  function s2cloudScore (img) {
  
  // A helper to apply an expression and linearly rescale the output.
  var rescale = function(img, exp, thresholds) {
    return img.expression(exp, {img: img})
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  };

  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1.0);
  // Clouds are reasonably bright in the blue band.
  score = score.min(rescale(img, 'img.B', [0.1, 0.3]));

  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.R + img.G + img.B', [0.2, 0.8]));

  // Clouds are reasonably bright in all infrared bands.
  score = score.min(
      rescale(img, 'img.NIR + img.SWIR1 + img.SWIR2', [0.3, 0.8]));

  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['G', 'SWIR1']);
  return score.min(rescale(ndsi, 'img', [0.8, 0.6]));
}

/**
 * Binary land mask based on OTSU method.
 * 
 * Finds optimum (dynamic per image) threshold to split image into 'Land' and 'Water'.
 * Image should have a substantial amount of water visible.
 * by JFreidman
 * 
 * @todo Add error mechanism when no clear threshold is present, i.e., OTSU fails.
 *
 * @param {ee.Image} im - A single band image to be thresholded, such as MNDWI.
 * 
 * @returns {ee.Image} An image where a value of 1 is land.
 */
function FindLandMask(im) {
  var thresh = ee.Number(0.75);
  var sigma = ee.Number(1);
  var scale = ee.Number(400);
  
  im = ee.Image(im).clip(table);
  var sf = im.get('system:footprint');
  var edge = ee.Algorithms.CannyEdgeDetector(im,thresh,sigma).focal_max(1);
  var hist_info = ee.Dictionary(im.mask(edge).reduceRegion(
      ee.Reducer.histogram(510), table, scale).get('wi'));
  var hist = hist_info.get('histogram');
  var buckets = hist_info.get('bucketMeans');
  var erf = ee.Array(hist).erf();
  var len = ee.List(hist).length();
  var counter = ee.List.sequence(1,len.subtract(1));
  
  // OTSU to build binary image
  function OTSU(i) {
      
      var len = ee.List(hist).length();
      var total = ee.Array(hist).accum(0).get([-1]);
      var q_L = ee.Array(ee.List(hist).slice(0,ee.Number(i))).divide(total).accum(0).get([-1]);
      var q_H = ee.Array(ee.List(hist).slice(ee.Number(i))).divide(total).accum(0).get([-1]);
      var miu_L = ee.Array(ee.List(hist).slice(0,ee.Number(i))).divide(total)
                    .multiply(ee.Array(ee.List.sequence(1,ee.Number(i))))
                    .accum(0).get([-1]).divide(q_L);
      var miu_H = ee.Array(ee.List(hist).slice(ee.Number(i))).divide(total)
                    .multiply(ee.Array(ee.List.sequence(ee.Number(i),ee.Number(len).subtract(1))))
                    .accum(0).get([-1]).divide(q_H);
      var sigma = q_L.multiply(q_H).multiply((miu_L.subtract(miu_H)).pow(2));
      return sigma;
    
  }
  var out = counter.map(OTSU);
  var maxxer = out.reduce(ee.Reducer.max());
  var threshold_index = ee.Number(out.indexOf(maxxer)).subtract(1);
  var threshold = ee.List(buckets).get(threshold_index.round()); 
  threshold = ee.Number(threshold);
  
  // Find thresholds for low and high ranges
  var hist_l = ee.List(hist).slice(threshold_index.round());
  var buckets_l = ee.List(buckets).slice(threshold_index.round());
  var threshold_l = ee.Number(ee.Array(buckets_l)
    .mask(ee.Array(hist_l).erf().gte(0.99)).reduce(ee.Reducer.median(), [0]).get([0]));
  var hist_w = ee.List(hist).slice(0, threshold_index.round());
  var buckets_w = ee.List(buckets).slice(0, threshold_index.round());
  var threshold_w = ee.Number(ee.Array(buckets_w)
    .mask(ee.Array(hist_w).erf().gte(0.99)).reduce(ee.Reducer.median(), [0]).get([0]));
  
  // Creat land-water image
  var land = im.where(im.gte(threshold),1).where(im.lt(threshold),0);
  
  return ee.Image(land)
    .set({
      'thresh_mid':threshold,
      'thresh_low':threshold_w,
      'thresh_high':threshold_l,
      'histogram':hist_info,
      //'erf':erf,
      'system:footprint': sf
      
    });
  
}

var s2 = ee.ImageCollection('COPERNICUS/S2')
        .filterBounds(table)
        .filterDate(sdate,edate)
        .filterMetadata(
          'CLOUDY_PIXEL_PERCENTAGE', 'less_than', max_cloud_cover
          )
        // Change to standard names
        .select(
          ee.List(['B2','B3','B4','B8','B11', 'B12', 'QA60']),
          ee.List(['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'QA60'])
          )
        // Cloud masked mndwi
        .map(function(im){
          var dater = im.get('system:time_start');
          var geo = im.get('system:footprint');
          im.updateMask(
            im.select('QA60').eq(1024).not()
            .focal_min(3.5));
          im.updateMask(
            im.select('QA60').eq(2048).not()
            .focal_min(3.5));
          //var buffer = im.geometry().buffer(-3500);
          return im.normalizedDifference(['SWIR1','G'])
            .select(['nd'],['wi'])
            .addBands(s2cloudScore(im.divide(10000)).select(['constant'],['cs']))
            //.clip(buffer)
            .clip(table)
            .set({
              'system:time_start':dater,
              'system:footprint':geo
            });
        });

var s2t = s2.select('cs').reduce(ee.Reducer.percentile([cs_pct]));
var s2_filt= s2.map(
  function(im){
    return im.updateMask(im.select('cs').lt(s2t));
    
  });

// Merge the collections, mean and scale
var ls = ee.ImageCollection(
  s2.select('wi')
  );

print('Images available:', ls.size());  

// Mean water index  
var mls = ls.mean().clip(table);
print('Histogram of water index:', ui.Chart.image.histogram(mls, table, 400));

// Create land-water-mask using fixed threshold
// Note this is different for time-ensemble-averaged images!!
var mndwiThresh = -0.2;
var lwfix = ee.Image.constant(1).updateMask(mls.gt(mndwiThresh));
Map.addLayer(lwfix.updateMask(lwfix), {"opacity":0.5,"palette":["black"]}, 'lw mask-TEA-S2-FIX');

// Create land-water-mask and get thresholds
var lw = FindLandMask(mls);
Map.addLayer(lw.updateMask(lw), {"opacity":0.5,"bands":["wi"],"palette":["blue"]}, 'lw mask-TEA-S2-OTSU');
//print(lw);
var low = ee.Number(lw.get('thresh_low'));
print('Lower threshold', low);
print('Mid threshold', ee.Number(lw.get('thresh_mid')));
var high = ee.Number(lw.get('thresh_high'));
print('High threshold', high);
