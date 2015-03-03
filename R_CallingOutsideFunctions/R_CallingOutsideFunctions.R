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
#' `
#' 
