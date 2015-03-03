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
#' ## Starting R on Omega
#' 
#' Remember to `source` the .bashrc file:
#' ```{}
#' source .bashrc
#' ```
#' 
#' 
## ------------------------------------------------------------------------
library(raster)
library(raster,lib.loc="/lustre/scratch/client/fas/geodata/aw524/R/")

#' 
#' ## Work with climate data
#' 
#' First set the path to the data directory.  You'll need to uncomment the line setting the directory to `lustre/...`.
#' 
## ------------------------------------------------------------------------
datadir="~/Downloads/bio1-9_30s_bil"
#datadir="/lustre/scratch/client/fas/geodata/aw524/data/"

outputdir="~/scratch/data"
if(!file.exists(outputdir)) dir.create(outputdir,recursive=T)

#' 
#' 
## ------------------------------------------------------------------------
file=paste0(datadir,"/bio_9.bil") 
## now use that to load the raster dataset:
tmean=raster(file)

#' 
#' 
#' ## Calling outside functions
#' R can do almost anything, but it isn't always the fastest at it...  For example, let's compare cropping a tif file using the `raster` package and `gdal_translate`:
#' 
## ------------------------------------------------------------------------
system.time(r1 <<- crop(tmean, extent(-77,-74.5,6,7)))
system.time(system(paste0("gdal_translate -projwin -77 7 -74.5 6 ",file, " ",outputdir,"/cropped.tif")))

#' 
#' 
#' But there are programs that go beyond the existing functionality of R, so `system()` can be a useful way to keep your entire workflow in R and call other programs as needed.  For example, pkfilter 
## ------------------------------------------------------------------------
system.time(system(paste0("pkfilter -i ",file," -o test.tif -f stdev")))

