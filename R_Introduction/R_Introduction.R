#' ---
#' title: "Introduction to R"
#' author: "Adam M. Wilson"
#' date: "February 11, 2015"
#' output: 
#'   html_document:
#'     toc: true
#'     keep_md: true
#' ---
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
#' 
#' ## Variables
#' 
## ----echo=FALSE----------------------------------------------------------
x=1
x

#' 
#' We can also assign a vector to a variable:
## ------------------------------------------------------------------------
x=c(5,8,14,91,3,36,14,30)
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
#' To calculate the mean, you could do it _manually_ like this
#' 
## ------------------------------------------------------------------------
(5+8+14+91+3+36+14+30)/8

#' 
#' Or use a function:
## ------------------------------------------------------------------------
mean(x)

#' 
#' Type `?functionname` to learn more about a function, e.g. `?mean`.  In RStudio, you can also search in the help panel.  There are other arguments too: `mean(x, trim = 0, na.rm = FALSE, ...)`
#' 
#' If you press `TAB` after a function name (such as `mean(`), it will show function arguments.
#' 
#' > Try to calculate the mean of `c(3,6,12,89)`.    
#' 
#' 
#' There are many ways to generate data in R such as sequences:
## ------------------------------------------------------------------------
seq(from=0, to=1, by=0.25)

#' and random numbers that follow a statistical distribution (such as the normal):
#' 
## ------------------------------------------------------------------------
a=rnorm(100,mean=0,sd=10)
hist(a)

#' 
#' 
#' ## Matrices
#' You can also use matrices (2-dimensional arrays of numbers):
## ------------------------------------------------------------------------
y=matrix(1:30,ncol=5)
y

#' 
#' Which behave much like vectors:
## ------------------------------------------------------------------------
y+2

#' 
#' and have 2-dimensional indexing:
## ------------------------------------------------------------------------
y[2,3]

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
## ----message=F,warning=FALSE---------------------------------------------
library(raster)

#' 
#' > R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options or just wait until it asks you.
#' 
#' 
#' If you don't have the packages above, install them in the package manager or by running `install.packages("raster")` (or use the package manager). 
#' 
#' ## Working with Raster Data
## ------------------------------------------------------------------------

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
