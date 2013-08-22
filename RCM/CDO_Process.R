
### set the working directory
setwd("/media/Data/Work/Courses/CyberInfrastructure/ClimateModule/Courses/Day3")
setwd("C:/Documents and Settings/EEB/Desktop/junk")

## load some libraries that we will need
library(ncdf)
library(raster)
library(sp)
library(rgdal)
library(reshape)
library(lattice)
##  If you get the error: Error in library(xx): there is no package called 'xx'
## then run the following for the missing package (replace ncdf with the package that is missing):
## install.packages("ncdf")
## and then re-run the lines above to load the newly installed libraries

## look at what files are there
list.files()
## do you see the two .nc files you downloaded?  
## If not, you either did not set the working directory correctly 
## or did not save the files in the right place.

## to run a CDO command from R, we have to 'wrap' it with a "system" call.  
## Try this:
system("cdo -V")

## that command tells R to run the command "cdo -V" at the system level.
## -V says you want to see what version of CDO you are have installed.
## If it worked, you should see a list of the various components of CDO
## and their version numbers and dates in the R console.  If you don't see this output, let Adam know.

## Let's take a closer look at one of those files using the "variable description" (VARDES) command:
system("cdo vardes gfdl_RCM3_Future.nc")

## What is the output?  What variables are in that file?


## CDO commands generally have the following syntax:  cdo command inputfile outputfile 
## and create a new file (outputfile) that is the result of running the command on the 
## inputfile.  Pretty easy, right?  

## What if you wanted to select only one variable (selvar) out of the netcdf file?  
## Try this:
system("cdo selvar,tmean gfdl_RCM3_Future.nc gfdl_RCM3_Future_tmean.nc")

## Now inspect the resulting file (gfdl_RCM3_Future_tmean.nc) with Panoply 
## or using the vardes command above.  
## What's different about it?  Do you understand what the selvar command does?  

## One powerful feature of the CDO tools is the ability to string together commands 
## to perform multiple operations at once.  You just have to add a dash ("-") before 
## any subsequent commands.  The inner-most commands are run first. If you have a multi-core processor, these separate 
## commands will be run in parallel if possible (up to 128 threads).  This comes in 
## handy when you are working with hundreds of GB of data.  Now convert the daily data 
## to mean annual temperature with the following command:
system("cdo timmean -selvar,tmean gfdl_RCM3_Future.nc gfdl_RCM3_Future_mat.nc")

## That command first uses selvar to subset only tmean (mean temperatures), 
## then uses timmean to take the average over the entire 'future' timeseries.
## Now explore the new file (gfdl_RCM3_Future_mat.nc) with Panoply.  
## How many time steps does it have?  What do the values represent?

## Now calculate the same thing (mean annual temperature) from the model output for the "Current" 
## time period:
system("cdo timmean -selvar,tmean gfdl_RCM3_Current.nc gfdl_RCM3_Current_mat.nc")

## And calculate the difference between them:
system("cdo sub gfdl_RCM3_Future_mat.nc gfdl_RCM3_Current_mat.nc gfdl_RCM3_mat_dif.nc")
## Look at this file in Panoply.  
## How much (mean annual) temperature change is projected for this model?

## Or you can string the commands above together and do all of the previous 
## operations in one step with the following command 
## (and if you have a dual core - or more - processor it will use all of them, isn't that cool!):
system("cdo sub -timmean -selvar,tmean gfdl_RCM3_Future.nc -timmean -selvar,tmean gfdl_RCM3_Current.nc gfdl_RCM3_mat_dif.nc")

## How does that command differ (in speed and ease) to calculating mean annual temperature like we did it last week? 

##  Now load these summarized data into R as a "SpatialGridDataFrame" 
## using the raster package which allows you to read directly from a netcdf file.
## You could also use the rgdal package if it has netcdf support (which it doesn't by default in windows)
mat_dif=as(raster("gfdl_RCM3_mat_dif.nc",varname="tmean",band=1),"SpatialGridDataFrame")
## That command first loads the data as a 'raster' object then converts it to a "SpatialGridDataFrame" 
## for easy plotting with spplot.

##load a polygon (shapefile) of New England to overlay on the grid so you know what you are looking at
ne=as(readOGR("NewEngland.shp","NewEngland"),"SpatialLines")

## take a look at the data
spplot(mat_dif,sp.layout=list("sp.lines",ne))

## if you don't like those colors, try some others (and add axis):
spplot(mat_dif,sp.layout=list("sp.lines",ne),col.regions=rev(heat.colors(20)),
scales=list(draw=T))

## So how much does that model say it will warm by mid-century?  
## What patterns do you see at this scale?  An difference between water and land?  Why would this be?

## The previous plot shows the information over space, but how about a distribution of values:
hist(mat_dif$values,col="grey",main="Change (future-current) Mean Annual Temperature",xlab="Change (Degrees C)")

######################################################################
######################################################################
## Now let's compare the RCM output with the GCM output

## We now need to calculate the change in mean annual temp and precipitation in the GCM output.
## You will have to edit the file names unless you saved them with exactly these names:
system("cdo sub -timmean GFDL_Future.nc -timmean GFDL_Current.nc GFDL_dif.nc")
## you might see some warnings about missing variables, this is due to how IRI writes it's netCDF files
## in this example we don't need to worry about them

## You can also calculate the seasonal mean change from this file like we did above:
system("cdo sub -yseasmean GFDL_Future.nc -yseasmean GFDL_Current.nc  GFDL_difseas.nc")


## To compare the GCM and RCM data and resolutions, 
## open the mean annual change (future-current) in Panoply: GFDL_dif.nc and GFDL_RCM3_mat_dif.nc.  What do you see?
## Set the scale bars to be the same (perhaps min: 1 and max:3) so you can compare the colors.
## Do you see differences in patterns or absolute values?
## How do you think the differences in resolution would affect analysis of biological/anthropological data?
## Is the increased resolution useful?

######################################################################
######################################################################
## Ok, now you have to start thinking (a little).  
######################################################################
## Below are a list of two data products that we want you to develop using the skills you just learned. 
## This will require you to run commands similar to the code above to create new summaries of the daily data.
## But you will have to make small changes to the code above to calculate the products below.
## I suggest that you copy the commands you need from above (don't edit the above code) and paste them below each question below and edit them there.  

## When you finish, ask Adam to come over and check (or you can email him the two lines of code if you finish after class).

## 1) Calculate the mean change (future-current) for all the variables in the GFDL_RCM3 dataset (tmin, tmax, tmean, and precipitation) like we did for just temperature above (Hint: do them one by one or simply don't subset a variable). After you run the command, check the output file with Panoply to check that it has the correct number of variables and only one time step.


INSERT YOUR CODE HERE******************

  ## 2) Calculate the mean change BY SEASON for all the variables in the GFDL_RCM3 dataset.
  ## The desired output will be a netcdf file with 4 time steps, the mean change for the 4 seasons.
  ## Hint: you will substitute "timmean" used above for something else.
  ## You may have to look in the CDO documentation to find the correct command.  Name the output file: gfdl_seasdif.nc
  

INSERT YOUR CODE HERE******************


      ## After you have finished number 2, read the file into R with the following command
      ## If you named your file something other than gfdl_seasdif.nc, you'll have to edit the file name below
      ## to match what you called your output above.
      difseas=as(brick("gfdl_seasdif.nc",varname="tmean",band=1:4),"SpatialGridDataFrame")
      ## that uses the 'brick' command to read in all 4 time steps from tmean as separate raster images
      ## now update the column names
      colnames(difseas@data)=c("Winter","Spring","Summer","Fall") 

      ## then plot the change in each of the four seasons
      spplot(difseas,sp.layout=list("sp.lines",ne),col.regions=rev(heat.colors(20)),scales=list(draw=T))

      ## Which season is going to warm the most?
      ## let's look at that another way, as a density plot (like a histogram) with the four seasons.
      densityplot(~value,groups=variable,data=melt(difseas@data),auto.key=T)
      ## or the same data again as a boxplot
      boxplot(value~variable,data=melt(difseas@data),col="grey",ylab="Temperature Change",xlab="Season")

##############################################################################
### if you have time, make these same plots for precipitation (ptot) by copying the code above and changing the variable name (varname="ptot")...

