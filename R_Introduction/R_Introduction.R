#' ---
#' title: "Introduction to R"
#' author: "Adam M. Wilson"
#' date: "February 11, 2015"
#' output: 
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     keep_md: true
#' ---
#' 
#' 

#' 
#' 
#' #  Introduction
#' This script is available:
#' 
#'   * [SpatialAnalysisTutorials repository](`r paste0("http://github.com/",repo)`)
#'   * Plain text (.R) with commented text [here](`r paste0("http://raw.githubusercontent.com/",repo,"/",output)`)
#'  
#' 
#' ## Vectors
#' 
## ----echo=FALSE----------------------------------------------------------
1:10

#' 
#' 
#' ## Matrices
#' 
#' 
#' 
#' ## Loading Packages
#' 
#' To load a package, you can simply type `library(package)` where `package` is the name of the package you want to load.  However, this only works for packages that you already have installed on your system.  To install new packages, you can use `install.packages()` or use the package manager. 
#' 
#' > R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options.
#' 
#' 
## ----message=F,warning=FALSE---------------------------------------------
library(ggplot2)

presentation_theme <- theme_grey()+
  theme(text = element_text(size = 25, colour = "black"))


#' 
#' 
#' If you don't have the packages above, install them in the package manager or by running `install.packages("ggplot2")`. 
#' 
#' 
#' ## Working with Raster Data
## ------------------------------------------------------------------------
library(raster)

#' 
#' 
#' 
#' # Coda
#' 
#' 
#' 

#' 
#' 
#' 
