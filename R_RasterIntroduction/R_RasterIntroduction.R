#' ---
#' title: "Introduction to Raster Package"
#' author: "Adam M. Wilson"
#' date: "February 23, 2015"
#' output: 
#'   html_document:
#'     toc: true
#'     keep_md: true
#' ---
#' 
#' 

#' 
#' ----
#' 
#' This script is available:
#' 
#'   * [SpatialAnalysisTutorials repository](`r paste0("http://github.com/",repo)`)
#'   * Plain text (.R) with commented text 
#'   [here](`r paste0("https://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/master/R_RasterIntroduction/",output)`)
#'  
#' 
#' ## Starting R
#' 
#' On the cluster, it's a little more complicated because the system is set up to host many different versions of software simutaneously and thus you need to load the correct `modules` before starting the software.  If the following lines are not in your .bashrc file, you'll need to run them before starting R.   
#' 
## ----,eval=F-------------------------------------------------------------
## module load Apps/R/3.0.2
## module load Rpkgs/RGDAL
## module load Rpkgs/GGPLOT2/1.0.0

#' 
#' We've added those modules to your `.bashrc` so you shouldn't have to run them separately.   To start R, simply type `R`:
## ----, eval=FALSE--------------------------------------------------------
## [aw524@login-0-0 ~]$ R
## 
## R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"
## Copyright (C) 2013 The R Foundation for Statistical Computing
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## R is free software and comes with ABSOLUTELY NO WARRANTY.
## You are welcome to redistribute it under certain conditions.
## Type 'license()' or 'licence()' for distribution details.
## 
##   Natural language support but running in an English locale
## 
## R is a collaborative project with many contributors.
## Type 'contributors()' for more information and
## 'citation()' on how to cite R or R packages in publications.
## 
## Type 'demo()' for some demos, 'help()' for on-line help, or
## 'help.start()' for an HTML browser interface to help.
## Type 'q()' to quit R.
## >

#' 
#' 
#' 
#' ## Loading packages
#' Loading packages works pretty much the same way on the cluster as it does on your computer (e.g. `library(...)`), but you will need to specify the location of the library on disk if there is no R package module available.    However, this only works for packages that you already have installed on your system.  To install new packages, you can use `install.packages()` (or use the package manager if you are working in RStudio). 
#' 
#' R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options or just wait until it asks you.
#' 
#' If you run `install.packages(...)` on Omega, you will probably see something like this:
#' 
## ----, eval=FALSE--------------------------------------------------------
## install.packages("raster")
## Installing package into ‘/lustre/home/client/apps/fas/Rpkgs/GGPLOT2/1.0.0/3.0’
## (as ‘lib’ is unspecified)
## Warning in install.packages("raster") :
##   'lib = "/lustre/home/client/apps/fas/Rpkgs/GGPLOT2/1.0.0/3.0"' is not writable
## Would you like to use a personal library instead?  (y/n) y
## Would you like to create a personal library
## ~/R/x86_64-unknown-linux-gnu-library/3.0
## to install packages into?  (y/n) y
## --- Please select a CRAN mirror for use in this session ---

#' 
#' Which is fine, you can install a personal copy of a library in any directory you choose (the default location is usually fine).  Once the package is installed, you can simply type `library(package)` where `package` is the name of the package you want to load.  
#' 
## ----message=F,warning=FALSE---------------------------------------------
library(raster)

#' 
#' Alternatively, you can point to Adam's personal library to avoid this extra step:
#' 
## ------------------------------------------------------------------------
library(raster,lib.loc="/lustre/scratch/client/fas/geodata/aw524/R/")

#' 
#' Below we'll run though several of the examples from the excellant [Raster vignette](http://cran.r-project.org/web/packages/raster/vignettes/Raster.pdf) by Robert J. Hijmans.
#' 
## ------------------------------------------------------------------------
x <- raster()
x

#' 
## ------------------------------------------------------------------------
x <- raster(ncol=36, nrow=18, xmn=-1000, xmx=1000, ymn=-100, ymx=900)
res(x)
res(x) <- 100
res(x)
ncol(x)

#' 
#' 
#' 
## ------------------------------------------------------------------------
# change the numer of columns (affects resolution)
ncol(x) <- 18
ncol(x)
res(x)

#' 
#' ## Spatial Projections
#' Raster package uses the standard [coordinate reference system (CRS)](http://www.spatialreference.org).  For example, see the projection format for the [_standard_ WGS84](http://www.spatialreference.org/ref/epsg/4326/).
#' ``{r}
#' projection(x) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
#' x
#' ```
#' 
#' 
## ------------------------------------------------------------------------
r <- raster(ncol=10, nrow=10)
ncell(r)
hasValues(r)

#' 
#' 
#' Use the 'values' function > # e.g.,
## ------------------------------------------------------------------------
values(r) <- 1:ncell(r)
hasValues(r)
inMemory(r)
values(r)[1:10]

plot(r, main='Raster with 100 cells')

#' 
#' > You can change the memory options using the `maxmemory` option in `rasterOptions()` 
#' 
#' ## Work with climate data
#' 
#' First set the path to the data directory.  You'll need to uncomment the line setting the directory to `lustre/...`.
#' 
## ------------------------------------------------------------------------
#datadir="data/"
datadir="/lustre/scratch/client/fas/geodata/aw524/data/"

#' 
#' 
## ------------------------------------------------------------------------
## we can use the paste0() command to construct a filename
paste0(datadir,"worldclim/bio_1.bil")

## now use that to load the raster dataset:
tmean=raster(paste0(datadir,"worldclim/bio_1.bil"))
## and inspect the object
tmean

#' 
#' Note the min/max of the raster.  What are the units?  It's always important to read over the documentation for any dataset to be sure you are using it correctly.  In this case, the WorldClim temperature dataset has a `gain` of 0.01, meaning that to work in degrees Celsius, you need to multiply by 0.01.  You can do this in the raster package with `gain()` 
## ------------------------------------------------------------------------
gain(tmean)=0.01
plot(tmean)

#' 
#' Let's dig a little deeper into the data object:
#' 
## ------------------------------------------------------------------------
## is it held in RAM?
inMemory(tmean)
## How big is it?
object.size(tmean)

## can we work with it directly in RAM?
canProcessInMemory(tmean)

## And the full data structure:
str(tmean)

#' 
#' 
#' # Modifying raster* objects
#' 
#' ## Spatial cropping
## ------------------------------------------------------------------------
## crop to a latitude/longitude box
r1 <- crop(tmean, extent(-77,-74.5,6,7))
r1
plot(r1)

#' 
#' ## Spatial aggregation
## ------------------------------------------------------------------------
## aggregate using a function
ra <- aggregate(r1, 20, fun=mean)
plot(ra)

#' > Now try to aggregate to the minimum (`min`) value with each 10 pixel window
#' 
#' ## Focal ("moving window") operation
## ------------------------------------------------------------------------
## apply a function over a moving window
rf <- focal(r1, w=matrix(1,11,11), fun=mean)
plot(rf)

#' 
## ------------------------------------------------------------------------
## apply a function over a moving window
rf_min <- focal(r1, w=matrix(1,11,11), fun=min)
rf_max <- focal(r1, w=matrix(1,11,11), fun=max)
rf_range=rf_max-rf_min

## or just use the range function
rf_range2 <- focal(r1, w=matrix(1,11,11), fun=range)
plot(rf_range2)

#' 
#' 
#' ## Raster calculations
#' 
#' the `raster` package has many options for _raster algebra_, including `+`, `-`, `*`, `/`, logical operators such as `>`, `>=`, `<`, `==`, `!` and functions such as `abs`, `round`, `ceiling`, `floor`, `trunc`, `sqrt`, `log`, `log10`, `exp`, `cos`, `sin`, `max`, `min`, `range`, `prod`, `sum`, `any`, `all`.
#' 
#' So, for example, you can 
## ------------------------------------------------------------------------
cellStats(r1,range)

## add 10
s = r1 + 10
cellStats(s,range)

## take the square root
s = sqrt(r1)
cellStats(s,range)

# round values
r = round(r1)
cellStats(r,range)

# find cells with values less than 2 degrees C
r = r1 < 2
cellStats(r,range)
plot(r)

# multiply s time r and add 5
s = s * r1 + 5
cellStats(s,range)


#' 
#' 
