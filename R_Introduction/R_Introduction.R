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
#' We can also assign values to variables:
## ------------------------------------------------------------------------
x=1:10
x

#' 
#' Subset the vector using `x[ ]` notation
## ------------------------------------------------------------------------
x[5]
x[1:5]

#' 
#' 
#' And do simple arithmetic:
## ------------------------------------------------------------------------
x+2

#' 
#' Or use a function to summarize:
## ------------------------------------------------------------------------
mean(x)

#' 
#' > Calculate the standard deviation of `x`
#' 
#' ## Matrices
#' You can also store matrices:
## ------------------------------------------------------------------------
y=matrix(1:30,ncol=5)
y

#' 
#' Which behave much like vectors:
## ------------------------------------------------------------------------
y+2

#' 
#' 
#' ## Plotting
#' Plotting in R can be as simple as typing `plot(x)`:
## ------------------------------------------------------------------------
plot(x)

#' 
#' But of course that's not such a meaningful plot (it's just a sequence of numbers and is rather dull).  Fortunately, though, R has some really flexible graphics packages such as `ggplot2`.  
#' 
#' ## Loading Packages
#' 
#' To load a package, you can simply type `library(package)` where `package` is the name of the package you want to load.  However, this only works for packages that you already have installed on your system.  To install new packages, you can use `install.packages()` or use the package manager. 
#' 
#' > R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options.
#' 
#' 
## ----message=F,warning=FALSE---------------------------------------------


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
