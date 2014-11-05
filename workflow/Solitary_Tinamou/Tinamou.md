
# Simple species distribution model workflow

In this session we will perform a simple species distribution model workflow for the Solitary Tinamou (Tinamus solitarius).  
![Tinamus](figure/TinamusSolitariusSmit.jpg)
Illustration by Joseph Smit, 1895


## Objectives

In this session we will:

 1. Download and process some raster environmental data
 2. Process occurrence data from various sources
 3. Fit a Bayesian species distribution model using the observations and environmental data
 4. Predict across the landscape and write the results to disk as a geotif (for use in GIS, etc.)
 
 

```
## Loading required package: sp Checking rgeos availability: TRUE Loading
## required package: lattice Loading required package: latticeExtra Loading
## required package: RColorBrewer Loading required package: hexbin Loading
## required package: grid rgeos version: 0.3-2, (SVN revision 413M) GEOS
## runtime version: 3.3.3-CAPI-1.7.4 Polygon checking: TRUE
## 
## rgdal: version: 0.8-11, (SVN revision 479M) Geospatial Data Abstraction
## Library extensions to R successfully loaded Loaded GDAL runtime: GDAL
## 1.9.2, released 2012/10/08 Path to GDAL shared files:
## /Library/Frameworks/R.framework/Versions/3.0/Resources/library/rgdal/gdal
## Loaded PROJ.4 runtime: Rel. 4.8.0, 6 March 2012, [PJ_VERSION: 480] Path to
## PROJ.4 shared files:
## /Library/Frameworks/R.framework/Versions/3.0/Resources/library/rgdal/proj
## Loading required package: coda Linked to JAGS 3.3.0 Loaded modules:
## basemod,bugs
```

```
## Error: cannot change working directory
```



## Data processing

### Import Gridded Environmental Data
Import some evironmental data (Climate, NPP, & Forest) and align it to a common grid

```r
# Check if data already exists
if (!file.exists("data/bio14_34_clip.tif")) system("bash DataPrep.sh")
```


Read them in as a raster stack

```r
env = stack(list.files(path = "data/", pattern = "*_clip.tif$", full.names = TRUE))
```

```
## Error: subscript out of bounds
```

```r
## do some renaming for convenience
names(env) = sub("_34", "", names(env))
```

```
## Error: object 'env' not found
```

```r
names(env) = sub("_clip", "", names(env))
```

```
## Error: object 'env' not found
```

```r
names(env)[names(env) == "MOD17A3_Science_NPP_mean_00_12"] = "npp"
```

```
## Error: object 'env' not found
```

```r
## set missing value in npp
NAvalue(env[["npp"]]) = 65535
```

```
## Error: object 'env' not found
```

```r
## get total % forest
forest = sum(env[[grep("consensus", names(env))]])
```

```
## Error: object 'env' not found
```

```r
names(forest) = "forest"
```

```
## Error: object 'forest' not found
```

```r
## add forest into the env stack
env = stack(env, forest)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'stack': Error: object 'env' not found
```

```r
## List all available environmental data
names(env)
```

```
## Error: object 'env' not found
```


### Import point observations
Download point data of occurrences from the Global Biodiversity Information Facility (GBIF) dataset 

```r
gbif_points = gbif("Tinamus", "solitarius", download = T, geo = T)
```

```
## Loading required package: XML
```

```
## http://data.gbif.org/ws/rest/occurrence/count?scientificname=Tinamus+solitarius&coordinatestatus=true
```

```
## Error: invalid request
```

```r
gbif_points = gbif_points[!is.na(gbif_points$lat), ]
```

```
## Error: object 'gbif_points' not found
```

Import the ebird points

```r
ebird = read.table("data/lat_long_ebd.txt", header = TRUE)
```

```
## Warning: cannot open file 'data/lat_long_ebd.txt': No such file or
## directory
```

```
## Error: cannot open the connection
```


Import a presence-absence shapefile from park checklists.

```r
parks = readOGR("data/", "protected_areas")
```

```
## Error: Cannot open file
```

```r
## Many of the parks with no observered presences were recorded as NA in the
## 'Presence' Column. Replace them with 0s.
parks$Presence[is.na(parks$Presence)] = 0
```

```
## Error: object 'parks' not found
```

```r
## generate an 'absence' dataset by sampling from the parks with no observed
## presences
nulls = coordinates(spsample(parks[parks$Presence == 0, ], 25, type = "stratified"))
```

```
## Error: error in evaluating the argument 'obj' in selecting a method for
## function 'coordinates': Error in spsample(parks[parks$Presence == 0, ],
## 25, type = "stratified") : error in evaluating the argument 'x' in
## selecting a method for function 'spsample': Error: object 'parks' not
## found
```



Import IUCN expert range map

```r
tin_range = readOGR("data/", "iucn_birds_proj")
```

```
## Error: Cannot open file
```

```r
tin_range = spTransform(tin_range, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'spTransform': Error: object 'tin_range' not found
```


Build a combined dataset (with source, presence/absence, and coordinates)

```r
points = rbind.data.frame(data.frame(src = "gbif", obs = 1, lat = gbif_points$lat, 
    lon = gbif_points$lon), data.frame(src = "ebird", obs = 1, lat = ebird$LATITUDE, 
    lon = ebird$LONGITUDE), data.frame(src = "parks", obs = 0, lat = nulls[, 
    "x2"], lon = nulls[, "x1"]))
```

```
## Error: object 'gbif_points' not found
```

```r
## turn it into a spatial dataframe and define projection
coordinates(points) = c("lon", "lat")
```

```
## Error: unable to find an inherited method for function 'coordinates<-' for
## signature '"standardGeneric"'
```

```r
projection(points) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
```

```
## Error: 'crs' is not a slot in class "standardGeneric"
```

```r

## Create a combined src_presence field for easy plotting
points$type = paste(points$src, points$obs, sep = "_")
```

```
## Error: object of type 'closure' is not subsettable
```



Import a world country boundary to ground the map

```r
World = readShapePoly("data/world_country_admin_boundary_shapefile_with_fips_codes.shp")
```

```
## Error: Error opening SHP file
```

```r
projection(World) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
```

```
## Error: object 'World' not found
```


As we saw before, there are a few points just outside the range, but those in the ocean are most likely wrong.  Let's add the distance to the range polygon as a way to filter the observations.  First we need a equidistant projection to do the calculation

```r
dproj = CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs")
points$dist = gDistance(spTransform(points, dproj), spTransform(tin_range, dproj), 
    byid = T)[1, ]
```

```
## Error: error in evaluating the argument 'obj' in selecting a method for
## function 'is.projected': Error in spTransform(points, dproj) : load
## package rgdal for spTransform methods Calls: spTransform -> spTransform
```

```r
## that adds 'distance' (in meters) from each point to the polygon so some
## points are > 2000km from the range, let's drop any more than 10km
points = points[points$dist < 10000, ]
```

```
## Error: error in evaluating the argument 'i' in selecting a method for
## function '[': Error in points$dist : object of type 'closure' is not
## subsettable
```


Check out the data in a plot

```r
spplot(points, zcol = "type", pch = 1:3, col.regions = c("red", "red", "black")) + 
    layer(sp.polygons(parks, col = NA, fill = ifelse(parks$Presence == 0, "black", 
        "red")), under = T) + layer(sp.polygons(World)) + layer(sp.polygons(tin_range, 
    fill = "grey"), under = T)
```

```
## Error: unable to find an inherited method for function 'spplot' for
## signature '"standardGeneric"'
```




Variable selection is tricky business and we're not going to dwell on it here... We'll use the following variables

```r
vars = c("bio5", "bio6", "bio13", "bio14", "npp", "forest")
```

[Worldclim "BIO" variables](http://www.worldclim.org/bioclim)

 * BIO5 = Max Temperature of Warmest Month
 * BIO6 = Min Temperature of Coldest Month
 * BIO13 = Precipitation of Wettest Month
 * BIO14 = Precipitation of Driest Month

To faciliate model fitting and interpretation, let's scale the environmental data

```r
senv = scale(env[[vars]])
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'scale': Error: object 'env' not found
```

```r
## Make a plot to explore the data
levelplot(senv, col.regions = rainbow(100, start = 0.2, end = 0.9), cuts = 99)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'levelplot': Error: object 'senv' not found
```


Add the (scaled) environmental data to each point

```r
pointsd = extract(senv, points, sp = F)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'extract': Error: object 'senv' not found
```

```r
## create single data.frame to hold all data for modelling
pointsd2 = data.frame(obs = points$obs, pointsd)
```

```
## Error: object of type 'closure' is not subsettable
```


## Fit a simple GLM to the data

```r
m1 = glm(obs ~ bio5 + bio6 + bio13 + bio14 + npp + forest, data = pointsd2, 
    family = "binomial")
```

```
## Error: object 'pointsd2' not found
```

```r
summary(m1)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for
## function 'summary': Error: object 'm1' not found
```


## Simple Bayesian Distribution Model

```r
## create the data object 
jags.data <- list(N.cells = nrow(pointsd2),
                  obs=points$obs,
                  X=data.frame(1,pointsd2[,vars]),
                  nBeta=length(vars)+1)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'nrow': Error: object 'pointsd2' not found
```

```r
 
# define the model
cat("
    model
{
    # priors
    for (l in 1:nBeta) {
    beta[l] ~ dnorm(0,0.01)
    }
    
    # likelihood
    for(i in 1:N.cells)
{
    # The observation as the result of a bernoulli outcome
    obs[i] ~ dbern(p[i])
    # logit transformation
    p[i]<-1/(1+exp(-lp.lim[i]))
    # Alternatively, could use the built-in function
    # logit(p[i])<-lp.lim[i]
    # 'stabilize' the logit to prevent hitting size limits
    lp.lim[i]<-min(999,max(-999,lp[i])) 
    }
    # The regression 
    # (using matrix notation rather than lp<-beta1+beta2*X1, etc)
    lp <- X%*%beta
    }
    ", file="model.txt")

params <- c("beta","p")

jm <- jags.model("model.txt",
                 data = jags.data,
                 n.chains = 3,
                 n.adapt = 2000)
```

```
## Error: object 'jags.data' not found
```


The model has been defined and an initial adaptive run of 2000 iterations complete.  Let's take some samples.

```r
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 5000, thin = 5)
```

```
## Error: object 'jm' not found
```


Extract the posterior samples and convert to `mcmc.list` objects for plotting/summaries

```r
ps.beta = as.mcmc.list(jm.sample$beta)
```

```
## Error: object 'jm.sample' not found
```

```r
ps.p = as.mcmc.list(jm.sample$p)
```

```
## Error: object 'jm.sample' not found
```


### Check Convergence

```r
xyplot(ps.beta, main = "Beta", strip = strip.custom(factor.levels = c("intercept", 
    vars)))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'xyplot': Error: object 'ps.beta' not found
```

```r
gelman.diag(ps.beta, confidence = 0.95, autoburnin = F, multivariate = T)
```

```
## Error: object 'ps.beta' not found
```



## Summarize the posterior betas

```r
densityplot(ps.beta, main = "Posterior Distributions", strip = strip.custom(factor.levels = c("intercept", 
    vars)), scales = list(relation = "same"), layout = c(1, 7)) + layer(panel.abline(v = 0))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'densityplot': Error: object 'ps.beta' not found
```

```r
HPDinterval(ps.beta[[1]], prob = 0.95)
```

```
## Error: object 'ps.beta' not found
```


## Predict model to the grid


```r
## First subset area to speed up predictions
pext = extent(c(-50, -48, -26.5, -24))
penv = crop(senv, pext)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'crop': Error: object 'senv' not found
```

```r

## if you want to make predictions for the full grid, run this line:
## penv=senv

## Calculate posterior estimates of p(occurrence) for each cell This extracts
## the posterior coefficients, performs the regression, calculates the
## quantiles, and takes the inverse logit to get p(occurrence)

## niter will use a reduced number of posterior samples to generate the
## summaries
pred = calc(penv, function(x, niter = 30) {
    mu1 = apply(apply(ps.beta[[1]][1:niter, ], 1, function(y) y * c(1, x)), 
        2, sum, na.rm = T)
    mu2 = quantile(mu1, c(0.025, 0.5, 0.975), na.rm = T)
    p = 1/(1 + exp(-mu2))
    return(p)
})
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'calc': Error: object 'penv' not found
```

```r
names(pred) = c("Lower_CI_2.5", "Median", "Upper_CI_97.5")
```

```
## Error: object 'pred' not found
```

```r
## Write out the predictions
writeRaster(pred, file = "Prediction.tif", overwrite = T)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'writeRaster': Error: object 'pred' not found
```


Plot the predictions

```r
levelplot(pred,col.regions=rainbow(100,start=.2,end=.9),cuts=99,margin=F)+
  layer(sp.polygons(tin_range,lwd=2))+
  layer(sp.points(points[points$obs==0,],pch="-",col="black",cex=8,lwd=4))+ #add absences
  layer(sp.points(points[points$obs==1,],pch="+",col="black",cex=4,lwd=4))    #add presences
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'levelplot': Error: object 'pred' not found
```


# Summary

In this script we have illustrated a complete workflow, including:

 1. Calling a BASH script (including GDAL Functions) from R to perform data pre-processing
 2. Running a (simple) Bayesian Species Distribution Model using rjags
 3. Making spatial predictions from model posteriors
 4. Writing results to disk as a geotif (for use in GIS, etc.)
 
