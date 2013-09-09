#Working with Climate Model Output
__________________

## Objective
To familiarize you with accessing and working with climate model output, in this exercise we'll compare the output from global climate models (GCMs) and regional climate models (RCMs). First, we'll download monthly data from the Geophysical Fluid Dynamics Laboratory (GFDL) GCM output using the IRI data library. Then we download and process daily RCM output from NARCCAP and compare them.  We’ve selected the GFDL model as an example for today. 

[![GFDL](http://www.gfdl.noaa.gov/pix/user_images/wga/cm2.6_sst.png)](http://www.gfdl.noaa.gov/brief-history-of-global-atmospheric-modeling-at-gfdl)

Here are a few relevant links:

* [Overview of the GFDL model](http://www.gfdl.noaa.gov/brief-history-of-global-atmospheric-modeling-at-gfdl)
* [NARCCAP website](http://www.gfdl.noaa.gov/brief-history-of-global-atmospheric-modeling-at-gfdl)

## Instructions
------------
### Global Climate Model Output

The IRI Data library at Columbia University provides convenient access to a variety of data sets especially climate and remotely sensed.  It has unusual flexibility in asking for *server side* processing of datasets.  This means you can ask for spatial or temporal subsets or even do some calculations with the data before you download it.  For more information, [check out their tutorials](http://iridl.ldeo.columbia.edu/dochelp/Tutorial/).

In this section we will download mean annual temperature and precipitation from the [IRI/LDEO site](http://iridl.ldeo.columbia.edu/).

#### Future Climate data

* Open the [IRI/LDEO website](http://iridl.ldeo.columbia.edu/) in a browser
* Click on the “expert” button on the left-hand column of the screen
* Paste the following code, which specifies the GFDL CM2.0 sresA2 run (2038-2070), into the expert box and click “OK”:

``` 
expert
SOURCES .WCRP .CMIP3
   (ipcc4/sresa2) @@
   (ipcc4/sresa2/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.sresa2.run1.atm.mo.xml .tas
  time (2038) (2070) RANGEEDGES
  lat (40N) (49N) RANGEEDGES
  lon (77W) (65W) RANGEEDGES
SOURCES .WCRP .CMIP3
   (ipcc4/sresa2) @@
   (ipcc4/sresa2/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.sresa2.run1.atm.mo.xml .pr
  c: 0.001 (m3 kg-1) :c
  mul
  c: 1000 (mm m-1) :c
  mulc: 86400 (s day-1) :c
  mul
  time (2038) (2070) RANGEEDGES
  lat (40N) (49N) RANGEEDGES
  lon (77W) (65W) RANGEEDGES
{tmean ptot}ds
```

>Do you understand what’s happing there? We’re asking for two sources (temp -tas- and precip -pr-) and subsetting both to a region and time. In addition, the lines that look like “ c:0.001(m3 kg-1):c“ convert the units from kg m-2 day-1 to mm/day (a more familiar unit). The new units (mm/ day) will be documented in the data you download (nice!).

Next:

1. Click on the red “Data Downloads & Files” button next to the views
* Scroll down and download the link that says “netCDF file.” 
* Save the file with the name: GFDL_Future.nc

#### Historical Data
For the past monthly data for precip and temperature copy and paste the following script for the GFDL CM2.0 20c run (1968-2000) into the expert box and click “OK”:

```
expert
SOURCES .WCRP .CMIP3
   (ipcc4/20c3m) @@
   (ipcc4/20c3m/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.20c3m.run2.atm.mo.xml .tas
  time (1968) (2000) RANGEEDGES
  lat (40N) (49N) RANGEEDGES
  lon (77W) (65W) RANGEEDGES
SOURCES .WCRP .CMIP3
   (ipcc4/20c3m) @@
   (ipcc4/20c3m/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.20c3m.run2.atm.mo.xml .pr
  c: 0.001 (m3 kg-1) :c
  mul
  c: 1000 (mm m-1) :c
  mul
  c: 86400 (s day-1) :c
  mul
  time (1968) (2000) RANGEEDGES
  lat (40N) (49N) RANGEEDGES
  lon (77W) (65W) RANGEEDGES
{tmean ptot}ds
```

Then:

1. Click on the red “Data Downloads & Files” button next to the views
2. Scroll down and download the link that says “netCDF file”
3. Save the file with the name: GFDL_Current.nc



>**OPTIONAL:** You could also calculate the difference (future-current) directly on the IRI webpage if you wanted, but we’re going to do it later using a new tool explained below. But, for the curious, here is how you would do it. For the future minus present data (in a single step) repeat step 1, and in the “expert” box, copy and paste the following script for the GFDL CM2.0 sresA2 run (yrs 2038-2070) minus the 20c run (yrs 1968-2000) monthly mean precipitation. Again, this is optional, you don’t need to run this for the rest of the exercise.

>```
expert
SOURCES .WCRP .CMIP3
   (ipcc4/sresa2) @@
   (ipcc4/sresa2/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.sresa2.run1.atm.mo.xml .pr
  time 0.0 monthlyAverage
  c: 0.001 (m3 kg-1) :c
  mul
  c: 1000 (mm m-1) :c
  mul
  c: 86400 (s day-1) :c
  mul
  time (Jan 2038) (Dec 2070) RANGE
  lon (85W) (65W) RANGEEDGES
  lat (30N) (50N) RANGEEDGES
SOURCES .WCRP .CMIP3
   (ipcc4/20c3m) @@
   (ipcc4/20c3m/gfdl_cm2_0) @@
   .pcmdi.ipcc4.gfdl_cm2_0.20c3m.run1.atm.mo.xml .pr
  time 0.0 monthlyAverage
  time (Jan 1968) (Dec 2000) RANGE
  lon (85W) (65W) RANGEEDGES
  lat (30N) (50N) RANGEEDGES
  time 12 splitstreamgrid
[time2]average
sub
```

>Then hit “ok”. Next click on the red “Data Files” button next to the views. Then, at the bottom of the page click on red “netCDF” .
Run the script again and download a netCDF for temperature change but replace “.pr” with “.tas”

### Regional Climate Model Output
Now we’ll work with output from a regional climate model (RCM) that was forced by the same GCM. We’ll get these data from the North American Regional Climate Change Assessment Program [http://www.narccap.ucar.edu/](http://www.narccap.ucar.edu/). Explore the site a little to see what they offer. Here’s a figure of their "region of interest."

![NARCCAP Domain](http://narccap.ucar.edu/img/narccap-domain.png)

I’ve already downloaded the data from NARCCAP you need for this exercise and put it the 'climate/data' folder, but it will be good to explore the site and see how the data are organized. 

1. Click on the “output data catalog” in the left column.
2. Choose the folder of data from the RCM3 regional climate model embedded in the GFDL global climate model (Row 5, Column 4)
3. Select [NARCCAP RCM3 gfdl-future Table 2](http://www.earthsystemgrid.org/dataset/narccap.rcm3.gfdl-future.table2.html)
4. Explore the “Variables” tab to see what variables are available.  For example:

Select [NARCCAP RCM3 gfdl-future Table 2](http://www.earthsystemgrid.org/dataset/narccap.rcm3.gfdl-future.table2.html) and explore the “Variables” tab to see what variables are available

```
huss
Description:	Surface_Specific_Humidity
Units:	 (kg kg-1)
Standard Name:	specific_humidity
Description:	"specific" means per unit mass. Specific humidity is the mass fraction of water vapor in (moist) air.
Units:	 (1)
Type:	CF
vas
Description:	Meridional_Surface_Wind_Speed
Units:	 (m s-1)
Standard Name:	northward_wind
Description:	"Northward" indicates a vector component which is positive when directed northward (negative southward). Wind is defined as a two-dimensional (horizontal) air velocity vector, with no vertical component. (Vertical motion in the atmosphere has the standard name upward_air_velocity.)
Units:	 (m s-1)
Type:	CF
ps
pr
Description:	Precipitation
Units:	 (kg m-2 s-1)
Standard Name:	precipitation_flux
Description:	In accordance with common usage in geophysical disciplines, "flux" implies per unit area, called "flux density" in physics.
Units:	 (kg m-2 s-1)
Type:	CF
tas
Description:	Surface_Air_Temperature
Units:	 (K)
Standard Name:	air_temperature
Description:	Air temperature is the bulk temperature of the air, not the surface (skin) temperature.
Units:	 (K)
Type:	CF
uas
Description:	Zonal_Surface_Wind_Speed
Units:	 (m s-1)
Standard Name:	eastward_wind
Description:	"Eastward" indicates a vector component which is positive when directed eastward (negative westward). Wind is defined as a two-dimensional (horizontal) air velocity vector, with no vertical component. (Vertical motion in the atmosphere has the standard name upward_air_velocity.)
Units:	 (m s-1)
Type:	CF
```
**You do NOT need to actually download any data from here.** These are the 3-hourly data, note the size of the files and time periods available. Do you see why we didn’t want all of you to download a dozen or more files from here at the same time? NARCCAP has not set up a mechanism to subset the data (like what is available via IRI Data Library) so the entire file has to be downloaded even if only a small portion is required. 

> **Optional:** To download data from NARCCAP you would do the following: 

>1. Register with an openID.
2. On the “Summary” tab, click: “Download files for this collection.”
3. To select a specific variable, check the box on the left and click “Sub-Select”

To avoid everyone downloading the full dataset, I’ve already downloaded the complete 3-hourly precipitation and temperature data from the NARCCAP site for several GCM-RCM combinations:

1. Historic (1968-2000)
2. Future (2038-2070) time periods
3. Historic *Observed* (1979-2004) from the NCEP reanalysis dataset

I then subsetted all of them them to the New England region and reprojected them to WGS84 latitude-longitude grid. I then took used the 3-hourly data to generate the daily total precipitation and min/mean/max temperature for each day. This reduced the size of the data from ~221GB to 1.2GB. I’ve posted the code I used to do this in climate/RCM/ as DownloadData.sh and ClipDaily.R if you are interested...

### Processing Climate data with CDO and R
___________________
To illustrate processing these data with R and other tools, we will use the Climate Data Operators ([CDO Tools](https://code.zmaw.de/projects/cdo)) developed by the Max-Plank Institute for Meteorology.  To continue, open the script climate/RCM/CDO_Process.Rmd in RStudio.

