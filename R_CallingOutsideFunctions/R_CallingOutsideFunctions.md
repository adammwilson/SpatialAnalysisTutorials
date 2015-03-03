# Introduction to Raster Package
Adam M. Wilson  
February 23, 2015  




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
## Loading required package: sp
```

```r
library(raster,lib.loc="/lustre/scratch/client/fas/geodata/aw524/R/")
```

## Work with climate data

First set the path to the data directory.  You'll need to uncomment the line setting the directory to `lustre/...`.


```r
datadir="~/Downloads/bio1-9_30s_bil"
#datadir="/lustre/scratch/client/fas/geodata/aw524/data/"
```



```r
file=paste0(datadir,"/bio_9.bil")
## now use that to load the raster dataset:
tmean=raster(file)
```

```
## rgdal: version: 0.8-16, (SVN revision 498)
## Geospatial Data Abstraction Library extensions to R successfully loaded
## Loaded GDAL runtime: GDAL 1.11.0, released 2014/04/16
## Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/3.1/Resources/library/rgdal/gdal
## Loaded PROJ.4 runtime: Rel. 4.8.0, 6 March 2012, [PJ_VERSION: 480]
## Path to PROJ.4 shared files: /Library/Frameworks/R.framework/Versions/3.1/Resources/library/rgdal/proj
```


## Calling outside functions
R can do almost anything, but it isn't always the fastest at it...  For example, let's compare cropping a tif file using the `raster` package and `gdal_translate`:

`

