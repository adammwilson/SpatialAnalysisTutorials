Spatially-explicit logistic regression
========================================================


```r
library(rgdal)
library(raster)
library(sp)
library(gam)
library(colorRamps)
library(ncf)
```


Load a vector dataset (shapefile) representing the San Diego bird atlas data for Purple finch:

```r
finch <- readOGR("finch", layer = "finch")
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: "finch", layer: "finch"
## with 414 features and 26 fields
## Feature type: wkbPolygon with 2 dimensions
```

```r
# Check how R handles vector data using the structure ('str') command to
# look at one of the polygons
str(finch[1, ])
```

```
## Formal class 'SpatialPolygonsDataFrame' [package "sp"] with 5 slots
##   ..@ data       :'data.frame':	1 obs. of  26 variables:
##   .. ..$ RC        : Factor w/ 414 levels "C-10","C-11",..: 27
##   .. ..$ FID_ATLAS1: num 0
##   .. ..$ AREA      : num 2.14e+08
##   .. ..$ ID_COARSE : Factor w/ 46 levels "c1","c10","c11",..: 1
##   .. ..$ FID_ATLA_1: num 11
##   .. ..$ AREA_1    : num 24584734
##   .. ..$ PERIMETER : num 19800
##   .. ..$ BIRD_ATLAS: num 13
##   .. ..$ BIRD_ATL_1: num 13
##   .. ..$ TRQ       : Factor w/ 396 levels "T09R1ENE","T09R1ESE",..: 19
##   .. ..$ BLOCKNAME : Factor w/ 395 levels "Agra","Agua Caliente",..: 286
##   .. ..$ area_mdm  : num 24584493
##   .. ..$ X_CEN     : num -117
##   .. ..$ Y_CEN     : num 33.4
##   .. ..$ ndvi      : num 129
##   .. ..$ meanelev  : num 337
##   .. ..$ minelev   : num 156
##   .. ..$ maxelev   : num 560
##   .. ..$ vegtypes  : int 19
##   .. ..$ maxtmp    : num 32.4
##   .. ..$ mintmp    : num 4.21
##   .. ..$ meanppt   : num 319
##   .. ..$ summert   : num 23.3
##   .. ..$ wintert   : num 11.7
##   .. ..$ urban     : num 10.5
##   .. ..$ present   : num 1
##   ..@ polygons   :List of 1
##   .. ..$ :Formal class 'Polygons' [package "sp"] with 5 slots
##   .. .. .. ..@ Polygons :List of 1
##   .. .. .. .. ..$ :Formal class 'Polygon' [package "sp"] with 5 slots
##   .. .. .. .. .. .. ..@ labpt  : num [1:2] 484697 3696562
##   .. .. .. .. .. .. ..@ area   : num 24584493
##   .. .. .. .. .. .. ..@ hole   : logi FALSE
##   .. .. .. .. .. .. ..@ ringDir: int 1
##   .. .. .. .. .. .. ..@ coords : num [1:37, 1:2] 487191 487167 487149 487153 486763 ...
##   .. .. .. ..@ plotOrder: int 1
##   .. .. .. ..@ labpt    : num [1:2] 484697 3696562
##   .. .. .. ..@ ID       : chr "0"
##   .. .. .. ..@ area     : num 24584493
##   ..@ plotOrder  : int 1
##   ..@ bbox       : num [1:2, 1:2] 482186 3694052 487191 3699070
##   .. ..- attr(*, "dimnames")=List of 2
##   .. .. ..$ : chr [1:2] "x" "y"
##   .. .. ..$ : chr [1:2] "min" "max"
##   ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slots
##   .. .. ..@ projargs: chr NA
```

```r
# Now look at the associated data frame (analogous to the *.dbf file that
# accompanied the shapefile)
head(finch@data)
```

```
##     RC FID_ATLAS1      AREA ID_COARSE FID_ATLA_1   AREA_1 PERIMETER
## 0  C-9          0 214130000        c1         11 24584734     19800
## 1 C-10          0 214130000        c1         12 23355580     19452
## 2 C-11          0 214130000        c1         22 23534602     19301
## 3 D-10          0 214130000        c1         43 24535492     19735
## 4 D-11          0 214130000        c1         45 23853714     19508
## 5  D-9          0 214130000        c1         46 23849098     19561
##   BIRD_ATLAS BIRD_ATL_1      TRQ      BLOCKNAME area_mdm  X_CEN Y_CEN
## 0         13         13 T09R3WNE        Rainbow 24584493 -117.2 33.41
## 1         14         14 T09R2WNW    Mt. Olympus 23355572 -117.1 33.40
## 2         24         24 T09R2WNE Trujillo Creek 23534551 -117.0 33.40
## 3         45         45 T09R2WSW    Gomez Creek 24535492 -117.1 33.36
## 4         47         47 T09R2WSE           Pala 23853715 -117.0 33.36
## 5         48         48 T09R3WSE  Monserate Mt. 23849098 -117.2 33.36
##    ndvi meanelev minelev maxelev vegtypes maxtmp mintmp meanppt summert
## 0 129.4    336.6  155.84   560.3       19  32.35   4.21   318.7   23.32
## 1 123.2    470.0  173.03   675.8       16  32.30   3.64   346.1   23.35
## 2 123.2    498.7  159.22  1008.6       22  32.88   3.55   342.9   23.58
## 3 128.0    244.2   87.45   552.7       25  32.88   4.57   301.7   23.58
## 4 112.5    228.6  110.70   611.5       27  33.55   4.75   292.6   24.08
## 5 142.1    208.3   82.69   471.4       23  32.47   4.93   288.1   23.17
##   wintert urban present
## 0   11.68 10.48       1
## 1   10.95  1.98       0
## 2   10.94  0.83       1
## 3   12.19  1.17       0
## 4   12.36  8.55       1
## 5   12.37 12.26       1
```


Scaling and centering the environmental variables:

```r
envi <- finch@data[, 15:25]
envi.scaled <- as.numeric(scale(envi))
finch@data[, 15:25] <- envi.scaled
```


Plotting the response (presence/absence data) and the predictor (NDVI):

```r
spplot(finch, zcol = c("present"), col.regions = c("white", "black"), colorkey = FALSE)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
spplot(finch, zcol = c("ndvi"), col.regions = topo.colors(20))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 


Now we will do the actual modelling. The first simple model links the presences and absences to NDVI.

This is the model that we will fit:

$\log ( \frac{p_i}{1-p_i} ) = \beta_0 + \beta_1 NDVI_i$

$o_i \sim Bernoulli(p_i)$

It can be fitted by simple glm() in R:

```r
ndvi.only <- glm(present ~ ndvi, data = finch@data, family = "binomial")
summary(ndvi.only)
```

```
## 
## Call:
## glm(formula = present ~ ndvi, family = "binomial", data = finch@data)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.8376  -0.4842  -0.2227  -0.0439   2.6604  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -2.939      0.296   -9.93   <2e-16 ***
## ndvi           2.652      0.322    8.23   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 385.67  on 413  degrees of freedom
## Residual deviance: 228.61  on 412  degrees of freedom
## AIC: 232.6
## 
## Number of Fisher Scoring iterations: 6
```

```r
## and let's extract predictions and residuals:
preds.ndvi.only <- predict(ndvi.only, type = "response")
resid.ndvi.only <- residuals(ndvi.only)
```


Now let's plot the logistic curve:

```r
newx <- data.frame(ndvi = seq(-2, 3, by = 0.1))
newy <- predict(ndvi.only, newdata = newx, type = "response")
plot(newx[, 1], newy, type = "l", xlab = "(Scaled) NDVI", ylab = "P of presence", 
    col = "red")
points(finch@data$ndvi, finch@data$present)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


The second model fits only the spatial trend in the data (using GAM and splines):

```r
space.only <- gam(present ~ s(X_CEN, 5) + s(Y_CEN, 5), data = finch@data, family = "binomial")
summary(space.only)
```

```
## 
## Call: gam(formula = present ~ s(X_CEN, 5) + s(Y_CEN, 5), family = "binomial", 
##     data = finch@data)
## Deviance Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.7776 -0.5682 -0.3386 -0.0396  2.6354 
## 
## (Dispersion Parameter for binomial family taken to be 1)
## 
##     Null Deviance: 385.7 on 413 degrees of freedom
## Residual Deviance: 270.6 on 403 degrees of freedom
## AIC: 292.6 
## 
## Number of Local Scoring Iterations: 12 
## 
## Anova for Parametric Effects
##              Df Sum Sq Mean Sq F value  Pr(>F)    
## s(X_CEN, 5)   1    0.2     0.2    0.22    0.64    
## s(Y_CEN, 5)   1   36.7    36.7   50.63 5.1e-12 ***
## Residuals   403  292.2     0.7                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Anova for Nonparametric Effects
##             Npar Df Npar Chisq  P(Chi)    
## (Intercept)                               
## s(X_CEN, 5)       4       37.2 1.6e-07 ***
## s(Y_CEN, 5)       4       11.1   0.025 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
preds.space.only <- predict(space.only, type = "response")
```


The third model uses both the NDVI and spatial trends to explain the finch's occurrences:

```r
space.and.ndvi <- gam(present ~ ndvi + s(X_CEN, 5) + s(Y_CEN, 5), data = finch@data, 
    family = "binomial")
summary(space.and.ndvi)
```

```
## 
## Call: gam(formula = present ~ ndvi + s(X_CEN, 5) + s(Y_CEN, 5), family = "binomial", 
##     data = finch@data)
## Deviance Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.2184 -0.3275 -0.1376 -0.0373  3.3394 
## 
## (Dispersion Parameter for binomial family taken to be 1)
## 
##     Null Deviance: 385.7 on 413 degrees of freedom
## Residual Deviance: 175.3 on 402 degrees of freedom
## AIC: 199.3 
## 
## Number of Local Scoring Iterations: 12 
## 
## Anova for Parametric Effects
##              Df Sum Sq Mean Sq F value  Pr(>F)    
## ndvi          1     54    53.8   48.44 1.4e-11 ***
## s(X_CEN, 5)   1      1     1.2    1.04    0.31    
## s(Y_CEN, 5)   1     30    29.8   26.87 3.4e-07 ***
## Residuals   402    446     1.1                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Anova for Nonparametric Effects
##             Npar Df Npar Chisq P(Chi)  
## (Intercept)                            
## ndvi                                   
## s(X_CEN, 5)       4      10.98  0.027 *
## s(Y_CEN, 5)       4       6.81  0.147  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
preds.space.and.ndvi <- predict(space.and.ndvi, type = "response")
resid.space.and.ndvi <- residuals(space.and.ndvi)
```


Now let's put all of the predictions and residuals together. We will use this for plotting:

```r
predictions <- data.frame(RC = finch@data$RC, preds.ndvi.only, resid.ndvi.only, 
    preds.space.only, preds.space.and.ndvi, resid.space.and.ndvi)
finch.preds <- merge(finch, predictions, by = "RC")
```


Here we will plot the predictions of the three models, together with the actual observed presences and absences:

```r
spplot(finch.preds, zcol = c("present", "preds.ndvi.only", "preds.space.only", 
    "preds.space.and.ndvi"), col.regions = matlab.like2(50))
```

```
## Error: unable to find an inherited method for function 'spplot' for
## signature '"data.frame"'
```


It is always useful to check the magnitude of spatial correlation in residuals:

```r
ndvi.only.cor <- correlog(finch.preds@data$X_CEN, finch.preds@data$Y_CEN, finch.preds@data$resid.ndvi.only, 
    increment = 0.2, resamp = 1)
```

```
## Error: trying to get slot "data" from an object (class "data.frame") that
## is not an S4 object
```

```r
space.and.envi.cor <- correlog(finch.preds@data$X_CEN, finch.preds@data$Y_CEN, 
    finch.preds@data$resid.space.and.ndvi, increment = 0.2, resamp = 1)
```

```
## Error: trying to get slot "data" from an object (class "data.frame") that
## is not an S4 object
```


And we can plot the correlograms:

```r
plot(ndvi.only.cor$mean.of.class, ndvi.only.cor$correlation, type = "b", xlab = "Distance class", 
    ylab = "Moran's I", main = "Residual correlograms")
```

```
## Error: error in evaluating the argument 'x' in selecting a method for
## function 'plot': Error: object 'ndvi.only.cor' not found
```

```r
points(space.and.envi.cor$mean.of.class, space.and.envi.cor$correlation, col = "red", 
    type = "b")
```

```
## Error: object 'space.and.envi.cor' not found
```

```r
abline(h = 0, lty = 2)
```

```
## Error: plot.new has not been called yet
```

```r
legend("topright", legend = c("ndvi.only", "space.and.ndvi"), col = c("black", 
    "red"), lwd = c(2, 2))
```

```
## Error: plot.new has not been called yet
```




