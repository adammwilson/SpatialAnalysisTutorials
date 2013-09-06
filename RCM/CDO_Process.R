
### set the working directory
setwd("/home/user/ClimateLabs/RCM")

## load some libraries that we will need
library(ncdf)
library(rasterVis)
library(sp)
library(rgdal)
library(reshape)
library(lattice)
##  If you get the error: Error in library(xx): there is no package called 'xx'
## then run the following for the missing package (replace ncdf with the package that is missing):
## install.packages("ncdf")
## and then re-run the lines above to load the newly installed libraries

## look at what files are in the data directory
list.files("../data")

## to run a CDO command from R, we have to 'wrap' it with a "system" call.  
## Try this:
system("cdo -V")

## that command tells R to run the command "cdo -V" at the system level.
## -V says you want to see what version of CDO you are have installed.
## If it worked, you should see a list of the various components of CDO
## and their version numbers and dates in the R console.  If you don't see this output, let us know.

## Let's take a closer look at one of those files using the "variable description" (VARDES) command:
system("cdo vardes ../data/gfdl_RCM3_Future.nc")
## note that there is a "../data/" command before the filename. This tells cdo to go up one level
## in the directory tree, then go into the data directory to find the file.  
## using relative paths like this allow the top-level folder to be moved
## without breaking links (even on other computers)

## What is the output?  What variables are in that file?
## You can also explore the file with panoply
## adding the "&" on the end tells linux to start panoply and not wait for it to close/finish
## But it remains a 'sub' process of R here, so if you close R, the panoply window will also die.
## See nohup (http://en.wikipedia.org/wiki/Nohup) for ways to spawn an independant process.
system("/usr/local/PanoplyJ/panoply.sh ../data/gfdl_RCM3_Future.nc &")

## CDO commands generally have the following syntax:  cdo command inputfile outputfile 
## and create a new file (outputfile) that is the result of running the command on the 
## inputfile.  Pretty easy, right?  

## What if you wanted to select only one variable (selvar) out of the netcdf file?  
## Try this:
system("cdo selvar,tmean ../data/gfdl_RCM3_Future.nc gfdl_RCM3_Future_tmean.nc")

## Now inspect the resulting file (gfdl_RCM3_Future_tmean.nc) with Panoply
## (you can open additional files using the panoply menus)
## or using the vardes command above.  
## What's different about it?  Do you understand what the selvar command does?  

## One powerful feature of the CDO tools is the ability to string together commands 
## to perform multiple operations at once.  You just have to add a dash ("-") before 
## any subsequent commands.  The inner-most commands are run first. If you have a multi-core processor, these separate 
## commands will be run in parallel if possible (up to 128 threads).  This comes in 
## handy when you are working with hundreds of GB of data.  Now convert the daily data 
## to mean annual temperature with the following command:
system("cdo timmean -selvar,tmean ../data/gfdl_RCM3_Future.nc gfdl_RCM3_Future_mat.nc")

## That command first uses selvar to subset only tmean (mean temperatures), 
## then uses timmean to take the average over the entire 'future' timeseries.
## Now explore the new file (gfdl_RCM3_Future_mat.nc) with Panoply.  
## How many time steps does it have?  What do the values represent?

## Now calculate the same thing (mean annual temperature) from the model output for the "Current" 
## time period:
system("cdo timmean -selvar,tmean ../data/gfdl_RCM3_Current.nc gfdl_RCM3_Current_mat.nc")

## And calculate the difference between them:
system("cdo sub gfdl_RCM3_Future_mat.nc gfdl_RCM3_Current_mat.nc gfdl_RCM3_mat_dif.nc")
## Look at this file in Panoply.  
## How much (mean annual) temperature change is projected for this model?

## Or you can string the commands above together and do all of the previous 
## operations in one step with the following command 
## (and if you have a dual core - or more - processor it will use all of them, isn't that cool!):
system("cdo sub -timmean -selvar,tmean ../data/gfdl_RCM3_Future.nc -timmean -selvar,tmean ../data/gfdl_RCM3_Current.nc gfdl_RCM3_mat_dif.nc")

## How does that command differ (in speed and ease) to calculating mean annual temperature like we did it last week? 

##  Now load these summarized data into R as a 'raster' class (defined by the raster library) 
## using the raster package which allows you to read directly from a netcdf file.
## You could also use the rgdal package if it has netcdf support (which it doesn't by default in windows)
mat_dif=raster("gfdl_RCM3_mat_dif.nc",varname="tmean",band=1)

##load a polygon (shapefile) of New England to overlay on the grid so you know what you are looking at
## the shapefile has polygons, but since we just want this for an overlay, let's convert to lines as well
ne=as(readOGR("NewEngland.shp","NewEngland"),"SpatialLines")

## take a look at the data
levelplot(mat_dif,margin=F)+layer(sp.lines(ne))

## So how much does that model say it will warm by mid-century?  
## What patterns do you see at this scale?  An difference between water and land?  Why would this be?

## The previous plot shows the information over space, but how about a distribution of values:
hist(mat_dif,col="grey",main="Change (future-current) \n Mean Annual Temperature",
     xlab="Change (Degrees C)",xlim=c(1.3,2.5))

######################################################################
######################################################################
## Now let's compare the RCM output with the GCM output

## We now need to calculate the change in mean annual temp and precipitation in the GCM output.
## You will have to edit the file names unless you saved them with exactly these names:
system("cdo sub -timmean ../data/GFDL_Future.nc -timmean ../data/GFDL_Current.nc GFDL_dif.nc")
## you might see some warnings about missing variables, this is due to how IRI writes it's netCDF files
## in this example we don't need to worry about them

## You can also calculate the seasonal mean change from this file like we did above:
system("cdo sub -yseasmean ../data/GFDL_Future.nc -yseasmean ../data/GFDL_Current.nc  GFDL_difseas.nc")


## To compare the GCM and RCM data and resolutions, 
## open the mean annual change (future-current) in Panoply: GFDL_dif.nc and GFDL_RCM3_mat_dif.nc.  What do you see?
## Set the scale bars to be the same (perhaps min: 1 and max:3) so you can compare the colors.
## Do you see differences in patterns or absolute values?
## How do you think the differences in resolution would affect analysis of biological/anthropological data?
## Is the increased resolution useful?

## OR we can plot a comparison via R
gcm=brick(raster("GFDL_dif.nc",varname="tmean"),raster("GFDL_dif.nc",varname="ptot"))
rcm=brick(raster("GFDL_RCM3_dif.nc",varname="tmean"),raster("GFDL_RCM3_dif.nc",varname="ptot"))

## it's not possible to merge rasters of different resolutions as a single raster or brick object
## but you can combine the plots easily with some functions in the rasterVis/latticeExtra packages like this:
c(levelplot(gcm[[1]],margin=F,main="Projected change in \n mean annual temperature (C)"),
  levelplot(rcm[[1]],margin=F))+layer(sp.lines(ne))

c(levelplot(gcm[[2]],margin=F,main="Projected change in \n mean annual precipitation (mm)"),
  levelplot(rcm[[2]],margin=F))+layer(sp.lines(ne))

######################################################################
######################################################################
######################################################################
## Below are a list of two data products that we want you to develop using the skills you just learned. 
## This will require you to run commands similar to the code above to create new summaries of the daily data.

## 1) Calculate the mean change (future-current) for all the variables in the GFDL_RCM3 dataset (tmin, tmax, tmean, and precipitation) like we did for just temperature above (Hint: do them one by one or simply don't subset a variable). After you run the command, check the output file with Panoply to check that it has the correct number of variables and only one time step.

system("cdo sub -timmean ../data/gfdl_RCM3_Future.nc -timmean ../data/gfdl_RCM3_Current.nc GFDL_RCM3_dif.nc")

  ## 2) Calculate the mean change BY SEASON for all the variables in the GFDL_RCM3 dataset.
  ## The desired output will be a netcdf file with 4 time steps, the mean change for the 4 seasons.
  
system("cdo sub -yseasmean ../data/GFDL_Future.nc -yseasmean ../data/GFDL_Current.nc  gfdl_seasdif.nc")

  ## Read the file into R with the following command
      difseas=brick("gfdl_seasdif.nc",varname="tmean",band=1:4)
      ## that uses the 'brick' command to read in all 4 time steps from tmean as separate raster images
      ## now update the column names
      names(difseas)=c("Winter","Spring","Summer","Fall") 

      ## then plot the change in each of the four seasons
      levelplot(difseas,col.regions=rev(heat.colors(20)))+layer(sp.lines(ne))

      ## Which season is going to warm the most?
      ## let's look at that another way, as a density plot (like a histogram) with the four seasons.
      densityplot(difseas,auto.key=T)
      ## or the same data again as a "violin" plot
      bwplot(difseas,col="grey",ylab="Temperature Change",xlab="Season")
      ##  Oops, default y limits are a bit too tight, let's fix that...
      bwplot(difseas,col="grey",ylab="Temperature Change",xlab="Season",scales=list(y=list(lim=c(1,5))))

##############################################################################
### if you have time, make these same plots for precipitation (ptot) by copying the code above and changing the variable name (varname="ptot")...

