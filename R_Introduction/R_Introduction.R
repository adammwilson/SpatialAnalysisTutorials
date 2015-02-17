#' ---
#' title: "Introduction to R"
#' author: "Adam M. Wilson"
#' date: "February 17, 2015"
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
#'   [here](`r paste0("https://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/master/R_Introduction/",output)`)
#'  
#' 
#' 
#' ## Variables
#' 
## ------------------------------------------------------------------------
x=1
x

#' 
#' We can also assign a vector to a variable:
## ------------------------------------------------------------------------
x=c(5,8,14,91,3,36,14,30)
x

#' 
#' And do simple arithmetic:
## ------------------------------------------------------------------------
x+2

#' 
#' Note that `R` is case sensitive, if you ask for `X` instead of `x`, you will get an error
## ----,eval=FALSE---------------------------------------------------------
## X

#' 
#' `Error: object 'X' not found`
#' 
#' ### Variable naming conventions
#' Naming your variables is your business, but there are [5 conventions](http://www.r-bloggers.com/consistent-naming-conventions-in-r/) to be aware of:
#' 
#' * alllowercase: _e.g._ `adjustcolor`
#' * period.separated: _e.g._ `plot.new`
#' * underscore_separated: _e.g._ `numeric_version`
#' * lowerCamelCase: _e.g._ `addTaskCallback`
#' * UpperCamelCase: _e.g._ `SignatureMethod`
#' 
#' ### Subsetting
#' Subset the vector using `x[ ]` notation
## ------------------------------------------------------------------------
x[5]

#' 
#' You can use a `:` to quickly generate a sequence:
## ------------------------------------------------------------------------
1:5

#' 
#' and use that to subset as well:
## ------------------------------------------------------------------------
x[1:5]

#' 
#' 
#' ### Using Functions
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
#' In RStudio, if you press `TAB` after a function name (such as `mean( `), it will show function arguments.
#' 
#' ![Autocomplete screenshot][pic1]
#' 
#' [pic1]: img/autocomplete.png "Autocomplete Screenshot"
#' 
#' 
#' > Try to calculate the mean of `c(3,6,12,89)`.    
#' 
#' Writing functions in R is pretty easy.  First think about how to break it down into parts:
## ------------------------------------------------------------------------
x1= sum(x)
x2=length(x)
x1/x2

#' 
#' Then put it all back together and create a new function called `mymean`:
## ------------------------------------------------------------------------
mymean=function(f){
  sum(f)/length(f)
}

mymean(f=x)

#' 
#' ### Missing data:  dealing with `NA` values
#' But, that simple function would be vulnerable to potential problems (such as missing data).  
#' 
## ------------------------------------------------------------------------
x3=c(5,8,NA,91,3,NA,14,30,100)

## Calculate the mean using the new function
mymean(x3)

## Use the built-in function (with and without na.rm=T)
mean(x3)
mean(x3,na.rm=T)

#' Writing simple functions is easy, writing robust, reliable functions can be hard...
#' 
#' ### Logical values
#' 
#' R also has standard conditional tests to generate `TRUE` or `FALSE` values (which also behave as `0`s and `1`s.  These are often useful for filtering data (e.g. identify all values greater than 5).  The logical operators are `<`, `<=`, `>`, `>=`, `==` for exact equality and `!=` for inequality.
#' 
## ------------------------------------------------------------------------
  x3 > 75
 
  x3 == 40
 
  x3 >   15

sum(x3>15,na.rm=T)

#' 
#' And of course you can save the results as variables:
## ------------------------------------------------------------------------
result =  x3 >  3
result

#' 
#' > Try this:  define a function that counts how many values in a vector are less than or equal to 12
#' 
#' ### Generating Data
#' 
#' There are many ways to generate data in R such as sequences:
## ------------------------------------------------------------------------
seq(from=0, to=1, by=0.25)

#' and random numbers that follow a statistical distribution (such as the normal):
#' 
## ------------------------------------------------------------------------
a=rnorm(100,mean=0,sd=10)
## illustrate with a histogram
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
#' ## Data Frames
#' Data frames are similar to matrices, but more flexible.  Matrices must be all the same type (e.g. all numbers), while a data frame can include multiple data types (e.g. text, factors, numbers). Dataframes are commonly used when doing statistical modeling in R.  
#' 
## ------------------------------------------------------------------------
data = data.frame( x = c(11,12,14),
                   y = c(19,20,21),
                   z = c(10,9,7))
data

mean(data$z)

mean(data[["z"]])

mean(data[,3])

#' 
#' ## Loading Packages
#' 
#' One of the best things about 
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
