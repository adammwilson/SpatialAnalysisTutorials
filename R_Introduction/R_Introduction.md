# Introduction to R
Adam M. Wilson  
February 11, 2015  




#  Introduction
This script is available:

  * [SpatialAnalysisTutorials repository](http://github.com/adammwilson/SpatialAnalysisTutorials/blob/master/R_Introduction)
  * Plain text (.R) with commented text 
  [here](https://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/master/R_Introduction/R_Introduction.R)
 


## Variables


```
## [1] 1
```

We can also assign a vector to a variable:

```r
x=c(5,8,14,91,3,36,14,30)
x
```

```
## [1]  5  8 14 91  3 36 14 30
```

Subset the vector using `x[ ]` notation

```r
x[5]
```

```
## [1] 5
```

```r
x[1:5]
```

```
## [1] 1 2 3 4 5
```


And do simple arithmetic:

```r
x+2
```

```
##  [1]  3  4  5  6  7  8  9 10 11 12
```

To calculate the mean, you could do it _manually_ like this


```r
(5+8+14+91+3+36+14+30)/8
```

```
## [1] 25.125
```

Or use a function:

```r
mean(x)
```

```
## [1] 5.5
```

Type `?functionname` to learn more about a function, e.g. `?mean`.  In RStudio, you can also search in the help panel.  There are other arguments too: `mean(x, trim = 0, na.rm = FALSE, ...)`

If you press `TAB` after a function name (such as `mean(`), it will show function arguments.

> Try to calculate the mean of `c(3,6,12,89)`.    


There are many ways to generate data in R such as sequences:

```r
seq(from=0, to=1, by=0.25)
```

```
## [1] 0.00 0.25 0.50 0.75 1.00
```
and random numbers that follow a statistical distribution (such as the normal):


```r
a=rnorm(100,mean=0,sd=10)
## illustrate with a histogram
hist(a)
```

![](R_Introduction_files/figure-html/unnamed-chunk-9-1.png) 


## Matrices
You can also use matrices (2-dimensional arrays of numbers):

```r
y=matrix(1:30,ncol=5)
y
```

```
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    7   13   19   25
## [2,]    2    8   14   20   26
## [3,]    3    9   15   21   27
## [4,]    4   10   16   22   28
## [5,]    5   11   17   23   29
## [6,]    6   12   18   24   30
```

Which behave much like vectors:

```r
y+2
```

```
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    3    9   15   21   27
## [2,]    4   10   16   22   28
## [3,]    5   11   17   23   29
## [4,]    6   12   18   24   30
## [5,]    7   13   19   25   31
## [6,]    8   14   20   26   32
```

and have 2-dimensional indexing:

```r
y[2,3]
```

```
## [1] 14
```

## Data Frames
Data frames are similar to matrices, but more flexible.  Matrices must be all the same type (e.g. all numbers), while a data frame can include multiple data types (e.g. text, factors, numbers).  


```r
data = data.frame( x = c(11,12,14),
                   y = c(19,20,21),
                   z = c(10,9,7))
data
```

```
##    x  y  z
## 1 11 19 10
## 2 12 20  9
## 3 14 21  7
```

```r
mean(data$z)
```

```
## [1] 8.666667
```

```r
mean(data[["z"]])
```

```
## [1] 8.666667
```

```r
mean(data[,3])
```

```
## [1] 8.666667
```

## Loading Packages

To load a package, you can simply type `library(package)` where `package` is the name of the package you want to load.  However, this only works for packages that you already have installed on your system.  To install new packages, you can use `install.packages()` or use the package manager. 


```r
library(raster)
```

> R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options or just wait until it asks you.


If you don't have the packages above, install them in the package manager or by running `install.packages("raster")` (or use the package manager). 

## Working with Raster Data




# Coda







