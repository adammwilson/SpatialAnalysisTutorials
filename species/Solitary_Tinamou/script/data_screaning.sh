
# see the data
cd ~/Solitary_Tinamou
more ebd_soltin1_relAug-2013.txt

# we are interested in lat long 
awk -F  "\t" '{  print $21 ,  $22 }'   ebd_soltin1_relAug-2013.txt > lat_long_ebd.txt



wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/tmean_34_tif.zip 
wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/tmax_34_tif.zip 
wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/tmin_34_tif.zip 
wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/alt_34_tif.zip 
wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/prec_34_tif.zip 



sudo R
install.packages("dismo")
install.packages("XML")
install.packages("maptools")
install.packages("rgeos")

library(raster)
library (dismo)
library (XML)
library (rgeos)
library (maptools)
library(sp)

gbif_data = gbif('Tinamus' , 'solitarius' , download=T , geo=T )
gbif_points = SpatialPoints( cbind ( gbif_data$lon[!is.na( gbif_data$lon)]  , gbif_data$lat[!is.na( gbif_data$lat)] ) )

ebird = read.table ("lat_long_ebd.txt" ,header = TRUE  )
ebird_points =  SpatialPoints(cbind (ebird$LONGITUDE ,  ebird$LATITUDE ))

field = read.table ("Tinamus_solitarius_rowids.asc" , header = TRUE  )
field_points =  SpatialPoints(cbind (field$st_x , field$st_y))

field_points_p =  SpatialPoints(cbind(subset(field, field$presence == 1 )$st_x, subset(field, field$presence == 1 )$st_y ))
field_points_a =  SpatialPoints(cbind(subset(field, field$presence == 0 )$st_x, subset(field, field$presence == 0 )$st_y)) 


World  = readShapePoly ("world_country_admin_boundary_shapefile_with_fips_codes.shp")

plot( World ,xlim=c(-80 , -15) , ylim=c(-40,-5) ) ;
points(field_points) 
points(ebird_points)
points(gbif_points)

# presence absence 

plot( World ,xlim=c(-65 , -20) , ylim=c(-20,-19) ) ;
points(field_points_p,col='red') ; 
points (ebird_points,col='red') ;
points (gbif_points,col='red') ; 

points(field_points_a,col='blue') ; 

presence=overlay(field_points_p,ebird_points,gbif_points)


w = getData('worldclim', var='tmin', res=0.5, lon=-50, lat=-20)

