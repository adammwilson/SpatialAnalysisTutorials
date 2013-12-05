#! /bin/bash

## Script to download data, unzip it, and clip to domain

## go up one directory from master script
cd ..

## Set some variables

datasrc="dropbox"  #this term can be either 'dropbox' or 'worldclim' to download the data directly

datadir=data       # this is where the data will go
mkdir $datadir

## download worldclim data directly from the server

if [ $datasrc = "worldclim" ] ; then 
    wget -P $datadir http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/bio_34_tif.zip 
fi 


if [ $datasrc = "dropbox" ]  ; then 
    cp /home/user/Dropbox/eeb713/week11/data/* $datadir
fi 

## loop through zipped files and uncompress them. Remove file that are not needed 
## keep
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month

for file in $datadir/*.zip ; do 
    unzip -o -d $datadir $file
    rm -f $datadir/bio[1-4]_34.tif    $datadir/bio[7-9]_34.tif  $datadir/bio1[0-2]_34.tif    $datadir/bio1[5-9]_34.tif    
done

# clip the data 
for file  in $datadir/*.tif ;  do 
  filename=`basename $file .tif`
  gdal_translate -projwin -56.500000 -12.0000000 -38.0000000 -30.0000000 $file $datadir/$filename"_clip.tif"
  rm -f $file 
done

# remove zipped and full rasters
