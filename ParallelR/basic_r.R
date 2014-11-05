#######################################
#######################################
######################################
# Stefano Casalegno
# Introduction to R

# The object of this document is to help you starting to use the R environment for statistical analysis and graphics.
# You can read and follow the the text, meanwhile you can copy the commands included into the frames part of this document and paste them in an interactive R session.
# Once you have familiarity with the general functioning of R and of R's objects you can further advance in learning R with online manuals and guides. There is a large variety of documentation available at

#     *      www.r-project.org.
#     *      http://www.statmethods.net/index.html
#     *      http://wiki.r-project.org/rwiki/doku.php

# as well efficacious learning tools. I would recommend that the user experiment with commands by, for example, trying different options to those stated. This experimentation is an important part of learning R using this synthetic document.

#############################################################
#### Starting R, install packages, getting help, stopping R
## Start R

#from a shell window type
R

#In the bash terminal will appear the following text:

# R version 2.10.0 (2009-10-26)
# Copyright (C) 2009 The R Foundation for Statistical Computing
# ISBN 3-900051-07-0
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# Natural language support but running in an English locale\\
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.

# the > sign and the following blinking cursor is advising you are in the R environment. If you like to enter in administrative mode type sudo R and you will be able to install packages

stefano@stefano-linux:~\$ sudo R



## Install a new package

# install.packages(pkgs=“PACKAGE NAME”, lib= “library directories where to install the packages”, dependencies=TRUE)

  install.packages(pkgs="randomForest", lib= "/usr/lib/R/site-library", dependencies=TRUE)

# to verify wich packages you have installed and where the library directories are type :

library()

# You will see the list of your available library as follow:

#  class                   Functions for Classification
#  colorspace              Color Space Manipulation
# gdata                   Various R programming tools for data
#                          manipulation
# gtools                  Various R programming tools
# maptools                Tools for reading and handling spatial objects
# MASS                    Main Package of Venables and Ripleys MASS
# mda                     Mixture and flexible discriminant analysis
# nnet                    Feed-forward Neural Networks and Multinomial
#                         Log-Linear Models
# randomForest            Breiman and Cutlers random forests for
#                         classification and regression
# rgdal                   Bindings for the Geospatial Data Abstraction
#                         Library
# sp                      classes and methods for spatial data
# spatial                 Functions for Kriging and Point Pattern
#                         Analysis
# vcd                     Visualizing Categorical Data
# xtable                  Export tables to LaTeX or HTML
# Packages in library '/usr/lib/R/site-library':
# base                    The R Base Package
# boot                    Bootstrap R (S-Plus) Functions (Canty)
# class                   Functions for Classification
# grid                    The Grid Graphics Package
# KernSmooth              Functions for kernel smoothing for Wand & Jones
#                         (1995)
# lattice                 Lattice Graphics
# MASS                    Main Package of Venables and Ripleys MASS
# Matrix                  Sparse and Dense Matrix Classes and Methods
# methods                 Formal Methods and Classes
# mgcv                    GAMs with GCV/AIC/REML smoothness estimation
#                         and GAMMs by PQL
# nlme                    Linear and Nonlinear Mixed Effects Models
# nnet                    Feed-forward Neural Networks and Multinomial
#                         Log-Linear Models
# rpart                   Recursive Partitioning
# spatial                 Functions for Kriging and Point Pattern
#                         Analysis
# splines                 Regression Spline Functions and Classes
# stats                   The R Stats Package
# stats4                  Statistical Functions using S4 Classes
#  survival                Survival analysis, including penalised
#                         likelihood.
# tcltk                   Tcl/Tk Interface
# tools                   Tools for Package Development
# utils                   The R Utils Package
# (END)

# to exit the list of library mode type q

q



## Getting help

# R gives help on function and commands. On-line help gives useful information as well. Getting used to R help is a key to successful statistical modelling. The on-line help can be accessed in HTML format by typing:

  help.start()

# A keyword search is possible using the Search Engine and Keywords link. You can use as well the help() or ? functions. For example, if we want to know how to use the matrix() function, the following two commands are equivalents:

  help(matrix)
  ? matrix

# Once you are in the help section, you can use “enter” to scroll down the menu and q to go back at the R prompt.

# The str(object.name) command is used to display the internal structure of an R object
# The summary(object.name) command gives a summary of a a object, usually a statistical summary but it is generic meaning it has different operations for different classes of a
# dir() show files in the current directory
# ls.str() str() for each variable in the search path
# getwd() is asking for the current working directory


## Stop R

  q()

# When you quit, R will ask you if you want to save the workspace (that is, all of the variables you have defined in this session); say “no” in order to avoid clutter.

# Should an R command seem to be stuck or take longer than you’re willing to wait, type Control-C.

# Calling linux shell scripting commands

# system(”…”) is used to call any linux scripting commands within the R environment> system(“pwd”) is equal to write getwd()

#####################
## Inputs and outputs

# Once you opened an R session and eventually loaded the library you need you can start exploring your data

# Loading data

# load(file.name) function, loads R type datasets written with the save function load(file.name)

##  Saving data

# save(object.name.1, object.name.2, … ) function save the specified object in XDR platform independent binary format

## Reading tables

# read.table(“filename”) Reads a file in table format and creates a data frame from it, with cases corresponding to lines and variables to fields in the file. The default separator sep=”” is any whitespace. You might need sep=”,” or ”;” and so on. Use header=TRUE to read the first line as a header of column names. The as.is=TRUE specification is used to prevent character vectors from being converted to factors. The comment.char=”“ specification is used to prevent ”#” from being interpreted as a comment and use skip=n to skip n lines before reading data;
# For more datails

  ? read.table
 
landuse04=read.csv("~/ost4sem/exercise/basic_adv_r/inputs/2004_landuse.csv", header=TRUE, sep=",", dec=".", na.string=":")

# read.csv(“filename”) is set to read comma separated files. read.csv(file.name, header = TRUE, sep = ”,”, quote=“\””, dec=”.”, fill =TRUE, comment.char=””, …)

# read.delim(“filename”) is set for reading tab-delimited files

# read.fwf() reads a table of fixed width formatted data into a ’data.frame’; widths is an integer vector, giving the widths of the fixed-width fields



# Show variables and data in your workspace

# the list function ls() outputs a list of existing R objects

  ls()
#   > [1] "landuse04"

# the structure function str(object.name) inform you on the structure of a specific object the summary function summary(object.name) inform on basic statistics of a specific object

# Save and remove data or R objects

# save(file, …) saves the specified objects (…) in the XDR platform-independent binary format

  save(landuse04, file="~/ost4sem/exercise/basic_adv_r/outputs/landuse2004")

# save.image(file) saves all objects

  save(file="~/ost4sem/exercise/basic_adv_r/outputs/landuse2004_and_more")

# rm(file, …) removes object you created or data you uploaded

  rm(landuse04)

# No object are present in memory now, use ls function to check it

  ls()
#   character(0)

# But since you saved the landuse2004 data yo can reload it using the load() function and check its structrure using the str function

  load("~/ost4sem/exercise/basic_adv_r/outputs/landuse2004")
  str(landuse04)

# as a result you will see the following informations:

# 'data.frame':   72 obs. of  13 variables:                    
# $ geographic.Unit                     : Factor w/ 72 levels "be10 Région de Bruxelles-Capitale/Brussels Hoofdstedelijk Gewest",..: 20 17 16 22 19 18 15 2 1 8 ...                                                                                                                                            
#  $ landuse                             : logi  NA NA NA NA NA NA ...                                                                                   
#  $ total.Total.area                    : num  NA NA NA NA NA ...                                                                                       
#  $ agriarea.Utilized.agricultural.area : num  NA NA NA NA NA ...                                                                                       
#  $ arabland.Arable.land                : num  NA NA NA NA NA ...                                                                                       
#  $ forest.Wooded.area                  : num  NA NA NA NA NA ...
#  $ garden.Private.gardens              : num  NA NA NA NA NA NA 0.2 0 0 0.1 ...
#  $ grasland.Permanent.grassland        : num  NA NA NA NA NA ...
#  $ greenfod.Green.fodder.on.arable.land: num  NA NA NA NA NA ...
#  $ fallow.Fallow                       : num  NA NA NA NA NA NA 23.6 0 0 7.4 ...
#  $ olivepl.Olive.plantations           : int  NA NA NA NA NA NA 0 0 0 0 ...
#  $ permcrop.Permanent.crops            : num  NA NA NA NA NA NA 21.4 0 0 19.2 ...
#  $ vineyard.Vineyards                  : num  NA NA NA NA NA NA 0 0 0 0 ...
# >


# Variables and calculations

# R has an interactive calculations functioning. The command is executed and results are dysplayed.
# R uses: +, -, /, and ^ for addition, subtraction, multiplication, division and exponentiation, respectively

2+2 

# R will execute and dysplay
# 
# > [1] 4
# 
# the [1] at the beginning of the line is just R printing an index of element numbers. If you print a result that displays on multiple lines, R will put an index at the beginning of each line.

2*5
# [1] 10

10/2
# [1] 5

2^3
# [1] 8



# variable settings

# you can simply create a variable typing: variable name = function, constant or calculation.
# Example:

x =3*2

# The results of 3*2 is not displayed. In fact, the x variable value is stored in memory without printing it. To display the x value you can use:

  print(x) 
#   [1] 6

# or

  x
#   [1] 6

# Most users apply a similar syntax using the ← character string instead of the = character.

  x <- 3
  x
#   [1] 3

# Also remember that R is case sensitive, print(X) or X is different from x. For instance

  a <- 3
  a
#   [1] 3
  A
#   Error: object 'A' not found

# Variable names in R must begin with a letter, followed by alphanumeric characters.

  3e = 2
#   Error: unexpected input in "3e "

# In long names you can use ”.” or “_” as in

# very.long.variable.name.X or very_long_variable_name_Y but you can’t use blank spaces in variable names. Avoid single letter names such us: c, l, q, t, C, D, F, I, and T, which are either built-in R functions or hard to tell apart.

  very.long.variable.name.X3 = 3
  very.long.variable.name.X3
#   [1] 3



# Interactive calculations
# 
# Once defined variables you can use them in interactive calculations :

  b = 2*2
  a = 2*3
  a*b
#   [1] 24

# And you can use variables in formulas :

  c = 60 /(a+b)
  c
#  [1] 6

typing a;b you can display a and b variables at the same time:

 a;b
#   [1] 6
#   [1] 4

# If you omit to close a parenthesis, R will display a *+* sign.

 c = 60 /(a+b
  +

# In this case you can either close the parenthesis in the next line or type ctrl + c to go back at a new starting prompt.

## Operators order

# When using more complex formulas be aware of the importance of the order of operators. Parenthesis have priority on exponentiation, or powers, then comes multiplication and division, finally addition and subtraction. The following command:

  C = ((a + 2 * sqrt(b))/(a + 8 * sqrt(b)))/2
  C
#   [1] 0.2272727

# is different from:

  C = a + 2 * sqrt(b) / a + 8 * sqrt(b) / 2
  C
#   [1] 14.66667

# as well as

  100-40/2^4
#   [1] 97.5

# is different from:

  (100-40)/2^4 
#   [1] 3.75

# and

    -2^4
#     [1] -16

is different from:

   (-2)^4
#    [1] 16



## Logical values

# R can perform conditional tests and generate True or False values as results.
# The logical operators are <, ⇐, >, >=, == for exact equality and != for inequality.

  x = 60
  x > 100
#   [1] FALSE
 
  x == 70
#   [1] FALSE
 
  x >   3
#   [1] TRUE
 
  x = 100
#   [1] TRUE

# Logical values can be stored as variables:

  x = 60
  logical.value =  x >  3
  logical.value
#   [1] TRUE

# In addition if c1 and c2 are logical expressions, then c1 & c2 is their intersection (“and”), c1 | c2 is their union (“or”), and ! !c1 is the negation of c1.
#################
#################
##  R objects

# The entities R operates on are technically known as objects. Examples are vectors of numeric (real) or complex values, vectors of logical values and vectors of character strings. These are known as “atomic” structures since their components are all of the same type, or mode, namely numeric, complex, logical, character and raw. R also operates on objects called lists, which are of mode list. These are ordered sequences of objects which individually can be of any mode. lists are known as “recursive” rather than atomic structures since their components can themselves be lists in their own right.
# 
# The other recursive structures are those of mode function and expression. Functions are the objects that form part of the R system along with similar user written functions, which we discuss in some detail later. Expressions as objects form an advanced part of R which will not be discussed in this guide, except indirectly when we discuss formulae used with modeling in R.
# 
# By the mode of an object we mean the basic type of its fundamental constituents. This is a special case of a “property” of an object. Another property of every object is its length. The functions mode(object) and length(object) can be used to find out the mode and length of any defined structure 10.
# 
# Further properties of an object are usually provided by attributes(object), see Getting and setting attributes. Because of this, mode and length are also called “intrinsic attributes” of an object. For example, if z is a complex vector of length 100, then in an expression mode(z) is the character string “complex” and length(z) is 100.

### Vectors

# Vectors are combinations of scalars in a string structure. Vectors must have their values all of the same mode. Thus any given vector must be unambiguously either logical, numeric, complex, character or raw. (The only apparent exception to this rule is the special “value” listed as NA for quantities not available, but in fact there are several types of NA). Note that a vector can be empty and still have a mode. For example the empty character string vector is listed as character(0) and the empty numeric vector as numeric(0).
# 
# c(…) is the generic function to combine arguments with the default forming a vector; with RECURSIVE=TRUE descends through lists combining all elements into one vector.
# To see details for the generic function c(…) and combine arguments forming a vector.

? c

# As an example we can create a simple vector of seven values typing:

  c(2,3,4,5,10,5,8)
#   [1]  2  3  4  5 10  5  8

# We can generate a sequence using the syntax :

  1:10
#   [1]  1  2  3  4  5  6  7  8  9 10

seq()
# We can generate the same sequence of 1:10 command using the seq() function. The syntax will be :

 seq(1,10)
#   [1]  1  2  3  4  5  6  7  8  9 10

# The seq() function “seq(from = number, to = number, by = number)” allow to create a vector starting from a value to another by a defined increment:

  seq(1,10, 0.25)
#   [1]  1.00  1.25  1.50  1.75  2.00  2.25  2.50  2.75  3.00  3.25  3.50  3.75
#   [13]  4.00  4.25  4.50  4.75  5.00  5.25  5.50  5.75  6.00  6.25  6.50  6.75
#   [25]  7.00 7.25  7.50  7.75  8.00  8.25  8.50  8.75  9.00  9.25  9.50  9.75
#   [37]  10.0
 
  seq(from = 1, to = 10, by =  0.25)
#   [1]  1.00  1.25  1.50  1.75  2.00  2.25  2.50  2.75  3.00  3.25  3.50  3.75
#   [13]  4.00  4.25  4.50  4.75  5.00  5.25  5.50  5.75  6.00  6.25  6.50  6.75
#   [25]  7.00  7.25  7.50  7.75  8.00  8.25  8.50  8.75  9.00  9.25  9.50  9.75
#   [37] 10.00

# rep() The replicate function “rep(x,times)” allow to replicate a vector several times in a more complex vector. Calculations can be included to form vectors as well and functions can be combined in the same command:

  one2tree = 1:3 
  rep(one2tree,10) 
#  [1] 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3
 
  c(10*0:10)
#   [1]   0  10  20  30  40  50  60  70  80  90 100
 
  rep(c (5*40:1, 5*1:40, 5, 6,7,8, 3, 2001:2014), 2)
#   [1]  200  195  190  185  180  175  170  165  160  155  150  145  140  135  130
#   [16]  125  120  115  110  105  100   95   90   85   80   75   70   65   60   55
#   [31]   50   45   40   35   30   25   20   15   10    5    5   10   15   20   25
#   [46]   30   35   40   45   50   55   60   65   70   75   80   85   90   95  100
#   [61]  105  110  115  120  125  130  135  140  145  150  155  160  165  170  175
#   [76]  180  185  190  195  200    5    6    7    8    3 2001 2002 2003 2004 2005
#   [91] 2006 2007 2008 2009 2010 2011 2012 2013 2014  200  195  190  185  180  175
#   [106]  170  165  160  155  150  145  140  135  130  125  120  115  110  105  100
#   [121]   95   90   85   80   75   70   65   60   55   50   45   40   35   30   25
#   [136]   20   15   10    5    5   10   15   20   25   30   35   40   45   50   55
#   [151]   60   65   70   75   80   85   90   95  100  105  110  115  120  125  130
#   [166]  135  140  145  150  155  160  165  170  175  180  185  190  195  200    5
#   [181]    6    7    8    3 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011
#   [196] 2012 2013 2014
#  
  rep(seq(1,3,0.5),3)
#   [1] 1.0 1.5 2.0 2.5 3.0 1.0 1.5 2.0 2.5 3.0 1.0 1.5 2.0 2.5 3.0



## Missing Values

# In some cases the components of a vector or of an R object more in general, may not be completely known. When an element or value is “not available” or a “missing value” in the statistical sense, a place within a vector may be reserved for it by assigning it the special value NA.
# Any operation on an NA becomes an NA.
# 
# The function is.na(x) gives a logical vector of the same size as x with value TRUE if and only if the corresponding element in x is NA.

   > z <- c(1:3,NA);  ind <- is.na(z)
   > ind
#    [1] FALSE FALSE FALSE  TRUE

# There is a second kind of “missing” values which are produced by numerical computation, the so-called Not a Number, NaN, values. Examples are 0/0 or Inf - Inf which both give NaN since the result cannot be defined sensibly.

  > Inf-Inf
#   [1] NaN
  > 0/0
#   [1] NaN

# In summary, is.na(xx) is TRUE both for NA and NaN values. To differentiate these, is.nan(xx) is only TRUE for NaNs. Missing values are sometimes printed as <NA> when character vectors are printed without quotes.

   z <- c(1:3,NA);  is.not.available <- is.na(z) ; is.not.a.number <-is.nan(z)
   is.not.a.number
#   [1] FALSE FALSE FALSE FALSE
  is.not.available 
#   [1] FALSE FALSE FALSE  TRUE

###########
##  Matrix

# Matrices or more generally arrays are multi-dimensional generalizations of vectors. In fact, they are vectors that can be indexed by two or more indices and will be printed in special ways. See Arrays and matrices.

#     *      factors provide compact ways to handle categorical data. See Factors.
#     *      lists are a general form of vector in which the various elements need not be of the same type, and are often themselves vectors or lists. Lists provide a convenient way to return the results of a statistical computation. See Lists.

# matrix() function creates a matrix from the given set of values. We use the matrix(x, nrow=, ncol=) function to set the matrix cells values, the number of raws and the number of columns. We can use the colnames() and rawnames() functions to set the column and raw names of the matrix-like object.

  matrix(data = NA, nrow = 2, ncol = 3) 
  example.matrix = matrix(0,2,3)
  example.matrix
#        [,1] [,2] [,3]
#   [1,]    0    0    0
#   [2,]    0    0    0
 
  example.matrix[1,]
#   [1] 0 0 0
 
  example.matrix[,2]
#   [1] 0 0
 
  example.matrix[1,] = 1:3
  example.matrix[2,] = c(5,10,4)
  example.matrix
#        [,1] [,2] [,3]
#   [1,]    1    2    3
#   [2,]    5   10    4
#  
 
  matrix.head = c("col a","col b","column c")
  matrix.side = c("first raw","second raw")
  str(matrix.side)
#   chr [1:2] "first raw" "second raw"

# When using ” ” whe create and refer to a character type “chr” input

  numeric.vector = c(rep(c (5*10:1, 5, 6), 2))
  character.vector  = as.character(numeric.vector)
  str(character.vector)
#   chr [1:24] "50" "45" "40" "35" "30" "25" "20" "15" "10" ...
 
 
  colnames(example.matrix) = matrix.head
  rownames(example.matrix) = matrix.side
  example.matrix
#              col a col b column c
#   first raw      1     2        3
#   second raw     5    10        4
#  
  str(example.matrix)
#   num [1:2, 1:3] 1 5 2 10 3 4
#    - attr(*, "dimnames")=List of 2
#     ..$ : chr [1:2] "first raw" "second raw"
#     ..$ : chr [1:3] "col a" "col b" "column c"



##########
##  Array

# An array can be considered as a multiply subscripted collection of data entries, for example numeric. R allows simple facilities for creating and handling arrays, and in particular the special case of matrices.
# 
# As well as giving a vector structure a dim attribute, arrays can be constructed from vectors by the array function, which has the form
# array(data_vector, dim_vector)

   Z <- array(1:24, c(3,4,2))
   Z
# #    , , 1
#         [,1] [,2] [,3] [,4]
#    [1,]    1    4    7   10
#    [2,]    2    5    8   11
#    [3,]    3    6    9   12
#    , , 2
#         [,1] [,2] [,3] [,4]
#    [1,]   13   16   19   22
#    [2,]   14   17   20   23
#    [3,]   15   18   21   24



#############
## Data Frames

# Data frames are matrix-like structures, in which the columns can be of different types. Think of data frames as data matrices with one row per observational unit but with (possibly) both numerical and categorical variables. Many experiments are best described by data frames: the treatments are categorical but the response is numeric.
# As a result R dataframes are tightly coupled collections of variables which share many of the properties of matrices and of lists. Data frames are used as the fundamental data structure by most of R's modeling software.
# A data frame is a list with class “data.frame”. There are restrictions on lists that may be made into data frames, namely :

#     *      The components must be vectors (numeric, character, or logical), factors, numeric matrices, lists, or other data frames.
#     *      Matrices, lists, and data frames provide as many variables to the new data frame as they have columns, elements, or variables, respectively.
#     *      Numeric vectors, logicals and factors are included as is, and character vectors are coerced to be factors, whose levels are the unique values appearing in the vector.
#     *      Vector structures appearing as variables of the data frame must all have the same length, and matrix structures must all have the same row size.

# see

   ? data.frame

# to learn how to construct a dataframe

   my.data.frame = data.frame(v=1:4,ch=c("a","b","c","d"),n=10)
   my.data.frame
#     v ch  n
#    1 1  a 10
#    2 2  b 10
#    3 3  c 10
#    4 4  d 10

# or

  my.data.frame = data.frame(vector=1:4,character=c("a","b","c","d"),const.vector=10,row.names =c("data1","data2","data3","data4"))
  my.data.frame
#         vector character const.vector
#   data1      1         a           10
#   data2      2         b           10
#   data3      3         c           10
#   data4      4         d           10

# data selection and manipulation
# 
# You can extract data from dataframes using the start and $ sign:

  my.data.frame[["character"]]
#   [1] a b c d
#   Levels: a b c d
 
  my.data.frame[[2]]
#   [1] a b c d
#   Levels: a b c d

# call the 3rd value of the character vector :

  my.data.frame[[2]][3]
#   [1] c
#   Levels: a b c d

# or using the $ syntax:

my.data.frame$vector
#   [1] 1 2 3 4
 
  my.data.frame$character[2:3]
#   [1] b c
#   Levels: a b c d

# You can add single arguments to a data frame, query informations, select and manipulate arguments or single values from a dataframe

my.data.frame$new
#   NULL
 
  my.data.frame$new = c(10,11,20,40)
  my.data.frame
#         vector character const.vector new
#   data1      1         a           10  10
#   data2      2         b           10  11
#   data3      3         c           10  20
#   data4      4         d           10  40
  str(my.data.frame)
#   'data.frame':	4 obs. of  4 variables:
#   $ vector      : int  1 2 3 4
#   $ character   : Factor w/ 4 levels "a","b","c","d": 1 2 3 4
#   $ const.vector: num  10 10 10 10
#   $ new         : num  10 11 20 40


# length(object.name) number of elements in an object such as matrix vector or dataframes:

  length(my.data.frame$new) 
#   [1] 4


# which(object.name) and which.max(object.name) return the index of a specific or of the greatest element of an object

  which.max(my.data.frame$new) 
#   [1] 4
 
  which(my.data.frame$new==20) 
#   [1] 3


# max(object.name) return the value of the greatest element

  max(my.data.frame$new) 
#   [1] 40


# sort(object.name) sort from small to big

  sort(my.data.frame$new) 
#   [1] 10 11 20 40

# rev(object.name) from big to small

  rev(sort(my.data.frame$new)) # 
#   [1] 40 20 11 10


# subset(object.name, …) subset returns a selection of an R-object with respect to criteria ( typically comparisons: x$V1 < 10); if the R-object is a data frame, the option select gives the variables to be kept or dropped using a minus sign

  subset(my.data.frame, my.data.frame$new == 20)
#      v ch  n new
#    3 3  c 10  20
 
sample(my.data.frame$new, 3)
# [1] 10 11 20
 
sample(my.data.frame$new, 3)
# [1] 40 10 11
 
sample(my.data.frame$new, 3)
# [1] 11 10 40


# More examples
# 
# The following R command give an example of simple procedure of importing data, cleaning a table by extracting relevant information, checking the presence of missing data.

  landuse04=read.csv("~/ost4sem/studycase/Lab_scripts/inputs/2004_landuse.csv", header=TRUE, sep=",", dec=".", na.string=":")
 
  forests04 = subset(landuse04, landuse04$forest.Wooded.area >= 0 )
  forests04$landuse = NULL
  forests04.check=na.fail(forests04)
  forests04$total.Total.area[1] = NA
  forests04.check=na.fail(forests04)
#   Error in na.fail.default(forests04) : missing values in object

# now recover the situation at the begin with no NA

  forests04 = subset(landuse04, landuse04$forest.Wooded.area >=0 )
  forests04$landuse = NULL
  forests04.check=na.fail(forests04)
  str(forests04)
#   'data.frame':   15 obs. of  12 variables:
#   data.frame':   15 obs. of  12 variables:
#   $ geographic.Unit                     : Factor w/ 72 levels "be10¬†R√©gion de Brxelles-Capitale/Brussels Hoofdstedelijk Gewest",..: 15 2 1 8 3 4 5 6 7 14 ...
#   $ total.Total.area                    : num  3052.8 16.1 16.1 1352.2 286.7 ...
#   $ agriarea.Utilized.agricultural.area : num  1393.8 0.2 0.2 633.8 91.3 ...
#   $ arabland.Arable.land                : num  840 0.2 0.2 430.8 60.2 ...
#   $ forest.Wooded.area                  : num  606.5 1.8 1.8 108.1 34 ...
#   $ garden.Private.gardens              : num  0.2 0 0 0.1 0 0 0 0 0 0.1 ...
#   $ grasland.Permanent.grassland        : num  530 0.1 0.1 181.4 28.7 ...
#   $ greenfod.Green.fodder.on.arable.land: num  248.6 0 0 161.4 41.7 ...
#   $ fallow.Fallow                       : num  23.6 0 0 7.4 1.1 1.2 1.2 2.3 1.7 16.2 ...
#   $ olivepl.Olive.plantations           : int  0 0 0 0 0 0 0 0 0 0 ...
#   $ permcrop.Permanent.crops            : num  21.4 0 0 19.2 1.5 9.7 3.1 4.1 0.9 2.1 ...
#   $ vineyard.Vineyards                  : num  0 0 0 0 0 0 0 0 0 0 ...

# Do you see something strange? look at the forests04$geographic.Unit level of factors and the dataframe number of variables !

# Let's fix it now!

  library(gdata)
  forests04 = drop.levels(forests04)
  str(forests04)
 
#   'data.frame':	15 obs. of  12 variables:
#   $ geographic.Unit                     : Factor w/ 15 levels "be1 Région de Bruxelles-Capitale/Brussels Hoofdstedelijk Gewest",..: 15 1 4 2 5 6 7 8 9 3 ...
# ...



###  Functions

# Functions are themselves objects in R which can be stored in the project's workspace. This provides a simple and convenient way to extend R.
# Usage: in writing your own function you provide one or more argument or names for the function, an expression (or body of the function) and a value is produced equal to the output function result

# function(arglist) expr function definition
return(value)

# example

 
myfunction <- function(x) x^5
myfunction(3)
# 243
 
body(myfunction) <- quote(5^x)
## or equivalently  body(myfunction) <- expression(5^x)
myfunction(3) 
# [1] 125
body(myfunction)
5^x
myfunction
function (x) 
5^x

