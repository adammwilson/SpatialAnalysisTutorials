# Introduction to R
Adam M. Wilson  
February 11, 2015  





#  Introduction
This script is available:

  * [SpatialAnalysisTutorials repository](http://github.com/adammwilson/SpatialAnalysisTutorials/blob/master/R_Introduction)
  * Plain text (.R) with commented text [here](http://raw.githubusercontent.com/adammwilson/SpatialAnalysisTutorials/blob/master/R_Introduction/R_Introduction.R)
 

## Vectors


```
##  [1]  1  2  3  4  5  6  7  8  9 10
```

We can also assign values to variables:

```r
x=1:10
x
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10
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

Or use a function to summarize:

```r
mean(x)
```

```
## [1] 5.5
```

<span style="color:red">
Calculate the standard deviation of `x`
</span>

## Matrices
You can also store matrices:

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


## Plotting
Plotting in R can be as simple as typing `plot(x)`:

```r
plot(x)
```

![](R_Introduction_files/figure-html/unnamed-chunk-9-1.png) 

But of course that's not such a meaningful plot (it's just a sequence of numbers and is rather dull).  Fortunately, though, R has some really flexible graphics packages such as `ggplot2`.  

## Loading Packages

To load a package, you can simply type `library(package)` where `package` is the name of the package you want to load.  However, this only works for packages that you already have installed on your system.  To install new packages, you can use `install.packages()` or use the package manager. 

> R may ask you to choose a CRAN mirror. CRAN is the distributed network of servers that provides access to R's software.  It doesn't really matter which you chose, but closer ones are likely to be faster.  From RStudio, you can select the mirror under Tools→Options.





If you don't have the packages above, install them in the package manager or by running `install.packages("ggplot2")`. 


## Working with Raster Data

```r
library(raster)
```



# Coda







