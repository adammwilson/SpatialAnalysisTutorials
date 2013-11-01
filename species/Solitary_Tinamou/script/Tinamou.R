library(raster)
library (dismo)
library (XML)
library (rgeos)
library (maptools)
library(sp)
library(rasterVis)
library(rgdal)
library(car)

## download point data of occurrences from the Global Biodiversity Information Facility (GBIF) dataset 
gbif_points = gbif('Tinamus' , 'solitarius' , download=T , geo=T)
gbif_points=gbif_points[!is.na(gbif_points$lat),]

## import the ebird points
ebird = read.table("/home/user/Dropbox/EEB713/week7/Solitary_Tinamou/pointdata/lat_long_ebd.txt" ,header = TRUE  )

## import a presence-absence dataset from parks
parks = read.table ("/home/user/ost4sem/exercise/Solitary_Tinamou/pointdata/Tinamus_solitarius_rowids.asc" , header = TRUE  )

## import expert range
tin_range=readOGR("/home/user/ost4sem/exercise/Solitary_Tinamou/shp/","iucn_birds_proj")
tin_range=spTransform(tin_range,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

## Make some pseudo-absences
pa=as.data.frame(coordinates(spsample(tin_range,70,type="random")))

## Build a combined dataset (with source, presence/absence, and coordinates)
points=rbind.data.frame(
  data.frame(src="gbif",obs=1,lat=gbif_points$lat,lon=gbif_points$lon),
  data.frame(src="ebird",obs=1,lat=ebird$LATITUDE,lon=ebird$LONGITUDE),
  data.frame(src="parks",obs=parks$presence,lat=parks$st_y,lon=parks$st_x),  
  data.frame(src="psuedo",obs=0,lat=pa$y,lon=pa$x)  
  )
## turn it into a spatial dataframe and define projection
coordinates(points)=c('lon','lat')
projection(points)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## Create a combined src_presence field for easy plotting
points$type=paste(points$src,points$obs,sep="_")

##import a world country boundary to ground the map
World  = readShapePoly ("/home/user/ost4sem/exercise/Solitary_Tinamou/shp/world_country_admin_boundary_shapefile_with_fips_codes.shp")
projection(World)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## check out the data
spplot(points,zcol="type",pch=1:5,col.regions=c("red","red","black","red","black"))+
  layer(sp.polygons( World))# ,xlim=c(-80 , -15) , ylim=c(-40,-5) ) ;

### uh oh, there are some records out in the ocean (typical GBIF!)
### There are several ways to handle this, one is to clip by the 'expert' range of the species.

## plot it again with the range
spplot(points,zcol="type",pch=1:5,col.regions=c("red","red","black","red"))+
  layer(sp.polygons( World))+
  layer(sp.polygons(tin_range,fill="grey"),under=T)

## So there are a few points just outside the range, but those in the ocean are most likely wrong.  
## Let's add the distance to the range polygon as a way to filter the observations:
## First we need a equidistant projection to do the calculation
dproj=CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs") 
points$dist=gDistance(spTransform(points,dproj),spTransform(tin_range,dproj),byid=T)[1,]
## that adds 'distance' (in meters) from each point to the polygon, check out the values:
hist(points$dist/1000,xlab="km to expert range")
## so some points are > 2000km from the range, let's drop any more than 10km
points=points[points$dist<10000,]

## plot it again with the range
spplot(points,zcol="type",pch=1:5,col.regions=c("red","red","black","red","black"))+
  layer(sp.polygons( World))+
  layer(sp.polygons(tin_range,fill="grey"),under=T)
## that looks better.



## plot it again adding psuedo-absences
spplot(points,zcol="type",pch=1:5,col.regions=c("red","red","black","red","black"))+
  layer(sp.polygons( World))+
  layer(sp.polygons(tin_range,fill="grey"),under=T)+
  layer(sp.points(pa,col="purple"))

## that looks better.


################################################################################
## Import some evironmental data (Climate & NPP) and align it to a common grid
env=stack(list.files(path = "/home/user/ost4sem/exercise/Solitary_Tinamou/worldclim/", pattern="*.tif$" , full.names = TRUE ))
## do some renaming for convenience
names(env)=sub("_34","",names(env))
names(env)[names(env)=="MOD17A3_Science_NPP_mean_00_12_clip"]="npp"
## set missing value in npp
NAvalue(env[["npp"]])=65535

## get total % forest
forest=sum(env[[grep("consensus",names(env))]])
names(forest)="forest"
## add forest into the env stack
env=stack(env,forest)

##  List all available environmental data
names(env)

## See how the points compare to altitude
var="prec12"
levelplot(env[[var]],col.regions=terrain.colors(20),margin=F)+
  layer(sp.polygons(tin_range))+
  layer(sp.points(points[points$obs==0,],col="black"))+ #add absences
  layer(sp.points(points[points$obs==1,],col="red"))    #add presences


## variable selection is tricky business and we're not going to dwell on it here...  
## we'll use the following variables
vars=c("tmean12","prec12","npp","alt","forest")
## look at the environmental variability
plot(env[[vars]])

## add the environmental data to each point
pointsd=extract(env[[vars]],points,sp=T)

## now check out some bivariate plots
## feel free to explore other variables

xyplot(npp~alt,groups=type,data=pointsd@data,auto.key=T)+
  layer(panel.ellipse(x,y,groups=pointsd$type,subscripts=T,level=.68))

xyplot(tmean12~prec12,groups=type,data=pointsd@data,auto.key=T)+
  layer(panel.ellipse(x,y,groups=pointsd$type,subscripts=T,level=.68))

xyplot(npp~forest,groups=type,data=pointsd@data,auto.key=T)+
  layer(panel.ellipse(x,y,groups=pointsd$type,subscripts=T,level=.68))

## Fit a very simple GLM to the data
m1=glm(obs~tmean12+prec12+npp+alt+forest,data=pointsd@data)
## Show the model summary
summary(m1)

## Use the model to predict p(presence) for all sites
pred=predict(env[[vars]],m1)

## Plot the predictions
levelplot(pred,col.regions=rainbow(20,start=.2,end=.9),margin=F)+
  layer(sp.polygons(tin_range,lwd=2))+
  layer(sp.points(points[points$obs==0,],col="black",cex=2,lwd=4))+ #add absences
  layer(sp.points(points[points$obs==1,],col="red",cex=2,lwd=4))    #add presences


## perhaps we should mask the predictions to the expert range
pred2=mask(pred,tin_range)

levelplot(pred2,col.regions=rainbow(20,start=.2,end=.9),margin=F)+
  layer(sp.polygons( World))+
  layer(sp.polygons(tin_range,lwd=2))+
  layer(sp.points(points[points$obs==0,],col="black",cex=2,lwd=4))+ #add absences
  layer(sp.points(points[points$obs==1,],col="red",cex=2,lwd=4))    #add presences

## quick accuracy assessment
predpoints=extract(pred,points,sp=T)
boxplot(layer~obs,data=predpoints@data,ylab="Estimated Propability of Occurrence",xlab="Observed Presence/Absence")

## to exclude the 'psuedo-absence' points, add this: [pointsd$type!="psuedo",]
## after pointsd@data

## what if we only used percent forest?
boxplot(forest~obs,data=pointsd@data,ylab="Percent Forest",xlab="Observed Presence/Absence")

##########################
## the polygons of the parks (and the associate presence/absence status) is available in 
## /home/user/ost4sem/exercise/Solitary_Tinamou/shp/protected_areas
## How would you use these polygons rather than the point data?
## check out the help for the extract function.
