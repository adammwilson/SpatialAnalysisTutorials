# Introduction to Raster Package
Adam M. Wilson  
February 23, 2015  




----

This script is available:

  * [SpatialAnalysisTutorials repository](http://github.com/adammwilson/SpatialAnalysisTutorials/blob/master/R_RasterIntroduction)
  * Plain text (.R) with commented text 
  [here](https://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/master/R_RasterIntroduction/R_CallingOutsideFunctions.R)
 

## Starting R on Omega

Remember to `source` the .bashrc file at the `$` prompt:
```{}
source .bashrc
```

And load the raster package (either from your own privaite library or from mine).

```r
library(raster)  # OR, 
library(raster,lib.loc="/lustre/scratch/client/fas/geodata/aw524/R/")
```

## Work with climate data

First set the path to the data directory.  You'll need to uncomment the line setting the directory to `lustre/...`.


```r
datadir="~/Downloads/bio1-9_30s_bil"
#datadir="/lustre/scratch/client/fas/geodata/aw524/data/"

outputdir="~/scratch/data"
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)
```

This time, let's look at BIO9, which is the [Mean Temperature of Driest Quarter](http://www.worldclim.org/bioclim):

```r
## first define the path to the data:
file=paste0(datadir,"/bio_9.bil") 
## now use that to load the raster dataset:
data=raster(file)
data
```

```
## class       : RasterLayer 
## dimensions  : 18000, 43200, 777600000  (nrow, ncol, ncell)
## resolution  : 0.008333333, 0.008333333  (x, y)
## extent      : -180, 180, -60, 90  (xmin, xmax, ymin, ymax)
## coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs 
## data source : /Users/adamw/Downloads/bio1-9_30s_bil/bio_9.bil 
## names       : bio_9 
## values      : -521, 366  (min, max)
```


## Calling outside functions
R can do almost anything, but it isn't always the fastest at it...  For example, let's compare cropping a tif file using the `raster` package and `gdal_translate`:

First let's do it using the `crop` function:

```r
r1 = crop(data, extent(-77,-74.5,6,7))
r1
```

```
## class       : RasterLayer 
## dimensions  : 120, 300, 36000  (nrow, ncol, ncell)
## resolution  : 0.008333333, 0.008333333  (x, y)
## extent      : -77, -74.5, 6, 7  (xmin, xmax, ymin, ymax)
## coord. ref. : +proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs 
## data source : in memory
## names       : bio_9 
## values      : 58, 281  (min, max)
```


To call functions using the `system` command, you first need to _assemble_ the text string equivalent to what you would run on the command line outside of R.  For example:

```r
paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif")
```

```
## [1] "gdal_translate -projwin -77 7 -74.5 6 ~/Downloads/bio1-9_30s_bil/bio_9.bil ~/scratch/data/cropped.tif"
```

We can save that as an object and _wrap_ it in with a `system()` call to actually run it:

```r
cmd=paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif")
system(cmd)
```

Or even do it all at once.  The depth of nested commands is a matter of personal taste (being consise vs. being clear).  

```r
system(paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif"))
```


### Processing speed

To compare processing speed, let's use the `system.time` function to record how long the operation takes.  To make it fair, we'll write the output tif to disk (rather than keeping it in RAM).  


```r
t1=system.time(
  r1 <<- crop(data,
              extent(-77,-74.5,6,7),
              file=paste0(outputdir,"/cropped.tif"),
              overwrite=T))
```

And let's run it one more time and time it to compare:

```r
t2=system.time(
  system(
    paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif")
    )
  )
```



The `crop` command took 0.15 seconds, while `gdal_translate` took 0.019.  That's a 87.3X speedup!  Whether this matters to you depends on the scale of your project.


## Another example: gdalwarp
But there are programs that go beyond the existing functionality of R, so `system()` can be a useful way to keep your entire workflow in R and call other programs as needed.  For example, it's common to aquire data in different spatial _projections_ that need to be reprojected to a common format prior to analysis.  A great tool for this is `gdalwarp`, which can also mosaic and perform some simple mosaicing procedures (such as taking the mean of multiple layers).


```r
command=paste0("gdalwarp -r cubic -t_srs '+proj=utm +zone=11 +datum=WGS84' ",outputdir,"/cropped.tif", " ",outputdir,"/reprojected.tif ")
system(command)
```

