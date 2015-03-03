---
title: "Introduction to Raster Package"
author: "Adam M. Wilson"
date: "February 23, 2015"
output: 
  html_document:
    toc: true
    keep_md: true
---




----

This script is available:

  * [SpatialAnalysisTutorials repository](http://github.com/adammwilson/SpatialAnalysisTutorials/blob/master/R_RasterIntroduction)
  * Plain text (.R) with commented text 
  [here](https://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/master/R_RasterIntroduction/R_CallingOutsideFunctions.R)
 

## Starting R on Omega

Remember to `source` the .bashrc file:
```{}
source .bashrc
```



```r
library(raster)
```

```
## Error in library(raster): there is no package called 'raster'
```

```r
library(raster,lib.loc="/lustre/scratch/client/fas/geodata/aw524/R/")
```

```
## Loading required package: sp
```

## Work with climate data

First set the path to the data directory.  You'll need to uncomment the line setting the directory to `lustre/...`.


```r
datadir="~/Downloads/bio1-9_30s_bil"
#datadir="/lustre/scratch/client/fas/geodata/aw524/data/"

outputdir="~/scratch/data"
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)
```



```r
file=paste0(datadir,"/bio_9.bil") 
## now use that to load the raster dataset:
tmean=raster(file)
```

```
## Loading required namespace: rgdal
```

```
## Error in .rasterObjectFromFile(x, band = band, objecttype = "RasterLayer", : Cannot create a RasterLayer object from this file. (file does not exist)
```


## Calling outside functions
R can do almost anything, but it isn't always the fastest at it...  For example, let's compare cropping a tif file using the `raster` package and `gdal_translate`:


```r
system.time(r1 <<- crop(tmean, extent(-77,-74.5,6,7)))
```

```
## Error in crop(tmean, extent(-77, -74.5, 6, 7)): error in evaluating the argument 'x' in selecting a method for function 'crop': Error: object 'tmean' not found
```

```
## Timing stopped at: 0.001 0 0.001
```

```r
system.time(system(paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif")))
```

```
##    user  system elapsed 
##   0.001   0.009   0.011
```


But there are programs that go beyond the existing functionality of R, so `system()` can be a useful way to keep your entire workflow in R and call other programs as needed.  For example, pkfilter 

```r
system.time(system(paste0("pkfilter -i ",file," -o test.tif -f stdev")))
```

```
##    user  system elapsed 
##   0.002   0.007   0.010
```
