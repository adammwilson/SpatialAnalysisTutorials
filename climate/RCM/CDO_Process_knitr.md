Example integration of R and CDO to process global climate model output
======

Explore the data for this exercise
-----------------

First, set the working directory and load libraries.

```r
setwd("/home/user/ost4sem/exercise/SpatialAnalysisTutorials/RCM")
```

```
## Error: cannot change working directory
```

```r
library(ncdf)
```

```
## Error: there is no package called 'ncdf'
```

```r
library(rasterVis)
```

```
## Error: there is no package called 'rasterVis'
```

```r
library(sp)
library(rgdal)
```

```
## Error: there is no package called 'rgdal'
```

```r
library(reshape)
```

```
## Error: there is no package called 'reshape'
```

```r
library(lattice)
library(xtable)
```

```
## Error: there is no package called 'xtable'
```

If you get the error: Error in library(xx): there is no package called 'xx'
 then run the following for the missing package (replace ncdf with the package that is missing):
 install.packages("ncdf")
 and then re-run the lines above to load the newly installed libraries

Look at the files in the data directory:

```r
list.files("../data")
```

```
## [1] "GFDL_Current.nc"      "GFDL_Future.nc"       "gfdl_RCM3_Current.nc"
## [4] "gfdl_RCM3_Future.nc"  "NewEngland.dbf"       "NewEngland.prj"      
## [7] "NewEngland.shp"       "NewEngland.shx"       "README.md"
```


To run a CDO command from R, we have to 'wrap' it with a "system" call.  
Try this:

```r
system("cdo -V")
```


That command tells R to run the command "cdo -V" at the system level (at the linux prompt rather than R).
 "-V says you want to see what version of CDO you are have installed.
 If it worked, you should see a list of the various components of CDO
 and their version numbers and dates in the R console.

Let's take a closer look at one of those files using the "variable description" (VARDES) command:

```r
system("cdo vardes ../data/gfdl_RCM3_Future.nc")
```

note that there is a "../data/" command before the filename. This tells cdo to go up one level
in the directory tree, then go into the data directory to find the file.  Using relative paths like this allow the top-level folder to be moved without breaking links (even on other computers).

 What is the output?  What variables are in that file?
 You can also explore the file with panoply with this command:

```r
system("/usr/local/PanoplyJ/panoply.sh ../data/gfdl_RCM3_Future.nc &")
```

 
Adding the "&" on the end tells linux to start panoply and not wait for it to close/finish, but it remains a 'sub' process of R here, so if you close R, the panoply window will also die. See nohup (http://en.wikipedia.org/wiki/Nohup) for ways to spawn an independant process.

 CDO commands generally have the following syntax:  cdo command inputfile outputfile 
 and create a new file (outputfile) that is the result of running the command on the 
 inputfile.  Pretty easy, right?  

 What if you wanted to select only one variable (selvar) out of the netcdf file?  
 Try this:

```r
system("cdo selvar,tmean ../data/gfdl_RCM3_Future.nc gfdl_RCM3_Future_tmean.nc")
```


 Now inspect the resulting file (gfdl_RCM3_Future_tmean.nc) with Panoply or using the vardes command above.  You can open additional files using the panoply menus).  What's different about it?  Can you guess what the selvar command does?  

 One powerful feature of the CDO tools is the ability to string together commands to perform multiple operations at once.  You just have to add a dash ("-") before any subsequent commands.  The inner-most commands are run first. If you have a multi-core processor, these separate commands will be run in parallel if possible (up to 128 threads).  This comes in 
 handy when you are working with hundreds of GB of data.  Now convert the daily data  to mean annual temperature with the following command:

```r
system("cdo timmean -selvar,tmean ../data/gfdl_RCM3_Future.nc gfdl_RCM3_Future_mat.nc")
```


 That command:
 * first uses selvar to subset only tmean (mean temperatures), 
 * uses timmean to take the average over the entire 'future' timeseries.

Now explore the new file (gfdl_RCM3_Future_mat.nc) with Panoply.
* How many time steps does it have?
* What do the values represent?

 Now calculate the same thing (mean annual temperature) from the model output for the "Current" 
 time period:

```r
system("cdo timmean -selvar,tmean ../data/gfdl_RCM3_Current.nc gfdl_RCM3_Current_mat.nc")
```


 And calculate the difference between them:

```r
system("cdo sub gfdl_RCM3_Future_mat.nc gfdl_RCM3_Current_mat.nc gfdl_RCM3_mat_dif.nc")
```

Look at this file in Panoply.  
 * How much (mean annual) temperature change is projected for this model?

Or you can string the commands above together and do all of the previous  operations in one step with the following command  (and if you have a dual core - or more - processor it will use all of them, isn't that cool!):

```r
system("cdo sub -timmean -selvar,tmean ../data/gfdl_RCM3_Future.nc -timmean -selvar,tmean ../data/gfdl_RCM3_Current.nc gfdl_RCM3_mat_dif.nc")
```


Now load these summarized data into R as a 'raster' class (defined by the raster library) 
 using the raster package which allows you to read directly from a netcdf file.
 You could also use the rgdal package if it has netcdf support (which it doesn't by default in windows)

```r
mat_dif = raster("gfdl_RCM3_mat_dif.nc", varname = "tmean", band = 1)
```

```
## Error: could not find function "raster"
```


Also load a polygon (shapefile) of New England to overlay on the grid so you know what you are looking at
 the shapefile has polygons, but since we just want this for an overlay, let's convert to lines as well

```r
ne = as(readOGR("../data/NewEngland.shp", "NewEngland"), "SpatialLines")
```

```
## Error: could not find function "readOGR"
```


Take a look at the data:

```r
levelplot(mat_dif, margin = F) + layer(sp.lines(ne))
```

```
## Error: object 'mat_dif' not found
```


Questions:
* So how much does that model say it will warm by mid-century?  
* What patterns do you see at this grain?
* Any difference between water and land?
* Why would this be?

 The previous plot shows the information over space, but how about a distribution of values:

```r
hist(mat_dif, col = "grey", main = "Change (future-current) \n Mean Annual Temperature", 
    xlab = "Change (Degrees C)", xlim = c(1.3, 2.5))
```

```
## Error: object 'mat_dif' not found
```


 We now need to calculate the change in mean annual temp and precipitation in the GCM output.

```r
system("cdo sub -timmean ../data/GFDL_Future.nc -timmean ../data/GFDL_Current.nc GFDL_dif.nc")
```

You might see some warnings about missing "scanVarAttributes", this is due to how IRI writes it's netCDF files
 in this example we don't need to worry about them


## Comparison of RCM & GCM output
_____________________________________________

Calculate the mean change (future-current) for all the variables in the GFDL_RCM3 dataset (tmin, tmax, tmean, and precipitation) like we did for just temperature above (Hint: do them one by one or simply don't subset a variable). After you run the command, check the output file with Panoply to check that it has the correct number of variables and only one time step.


```r
system("cdo sub -timmean ../data/gfdl_RCM3_Future.nc -timmean ../data/gfdl_RCM3_Current.nc GFDL_RCM3_dif.nc")
```


Let's compare the GCM and RCM data and resolutions by plotting them.  First we need to read them in.  Let's use the raster package again:

```r
gcm = brick(raster("GFDL_dif.nc", varname = "tmean"), raster("GFDL_dif.nc", 
    varname = "ptot"))
```

```
## Error: could not find function "brick"
```

```r
rcm = brick(raster("GFDL_RCM3_dif.nc", varname = "tmean"), raster("GFDL_RCM3_dif.nc", 
    varname = "ptot"))
```

```
## Error: could not find function "brick"
```


It would be nice if we could combine the two model types in a single two-panel plot with common color-bar to facilitate comparison.  For understandable reasons, it's not possible to merge rasters of different resolutions as a single raster stack or brick object.  But you can combine the plots easily with functions in the rasterVis/latticeExtra packages like this:

```r
c(levelplot(gcm[[1]], margin = F, main = "Projected change in \n mean annual temperature (C)"), 
    levelplot(rcm[[1]], margin = F)) + layer(sp.lines(ne))
```

```
## Error: object 'gcm' not found
```

```r

c(levelplot(gcm[[2]], margin = F, main = "Projected change in \n mean annual precipitation (mm)"), 
    levelplot(rcm[[2]], margin = F)) + layer(sp.lines(ne))
```

```
## Error: object 'gcm' not found
```


* Do you see differences in patterns or absolute values?
* How do you think the differences in resolution would affect analysis of biological/anthropological data?
* Is the increased resolution useful?

## Calculate the mean change BY SEASON for all the variables in the GFDL_RCM3 dataset
The desired output will be a netcdf file with 4 time steps, the mean change for each of the 4 seasons.
  

```r
system("cdo sub -yseasmean ../data/GFDL_Future.nc -yseasmean ../data/GFDL_Current.nc  gfdl_seasdif.nc")
```


Read the file into R using the 'brick' command to read in all 4 time steps from tmean as separate raster images.  Then update update the column names (to spring, summer, etc.) and plot the change in each of the four seasons.


```r
difseas = brick("gfdl_seasdif.nc", varname = "tmean", band = 1:4)
```

```
## Error: could not find function "brick"
```

```r
names(difseas) = c("Winter", "Spring", "Summer", "Fall")
```

```
## Error: object 'difseas' not found
```

```r
levelplot(difseas, col.regions = rev(heat.colors(20))) + layer(sp.lines(ne))
```

```
## Error: object 'difseas' not found
```


Which season is going to warm the most?
Let's marginalize across space and look at the data in a few different ways:
* density plot (like a histogram) with the four seasons.
* "violin" plot

```r
densityplot(difseas, auto.key = T)
```

```
## Error: object 'difseas' not found
```

```r
bwplot(difseas, ylab = "Temperature Change", xlab = "Season", scales = list(y = list(lim = c(1, 
    5))))
```

```
## Error: object 'difseas' not found
```


We can also look at the relationship between variables with a scatterplot matrix:

```r
splom(difseas, colramp = colorRampPalette("blue"))
```

```
## Error: object 'difseas' not found
```


Or with a tabular correlation matrix:

```r
print(xtable(layerStats(difseas, stat = "pearson")[[1]], caption = "Pearson Correlation between seasonal change", 
    digits = 3), type = "html")
```

```
## Error: could not find function "xtable"
```


# Additional excercises
If you have time, try to:
* make these same plots for precipitation (ptot) by copying the code above and changing the variable name (varname="ptot")
* develop some quantiative summaries of future change using summary() on the difseas object


