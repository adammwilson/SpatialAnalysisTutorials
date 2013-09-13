#### This script downloads the complete matrix of GCM-RCM temperature and precipitation from NARCCAP.

### choose a destination for the raw, 3-hour files
setwd("/media/Data2/Data/narccap")

### Download the data
### Be careful, running this script will download ~220GB to your drive
system("sh DownloadData.sh")

### load libraries
library(ncdf4)

### make metadata table
fs=data.frame(file=list.files(pattern="*.nc"),path=list.files(pattern="*.nc",full=T),stringsAsFactors=F)
fs[,c("var","rcm","gcm","startdate")]=do.call(rbind,strsplit(fs$file,"_"))
fs$startdate=as.Date(do.call(rbind,strsplit(fs$startdate,"[.]"))[,1],"%Y%m%d")
fs$timeperiod=ifelse(fs$startdate<as.Date("2010-01-01"),"Current","Future")
fs$chunk=paste(fs$gcm,fs$rcm,fs$timeperiod,sep="_")

### first check to see if temp and precip exist for all model-times
table(fs$chunk,fs$var)

### loop through chunks to 1) subset to region 2) take daily sum/mean/min/max 3) merge to continuous timeseries
targetdir="/media/Data/Work/Regional/NewEngland/ClimateProjections/NARCCAP/2_DailyClip/"  #destination of daily files
ROI="-77,-65,40,49"  # Region to clip

for(i in unique(fs$chunk)){
  ## Precip
  pfs=fs$path[fs$chunk==i&grepl("pr",fs$path)]    #get precip files
  pclip=paste("-daysum -expr,'pr=pr*3600*3;' -shifttime,-5hour -sellonlatbox,",ROI," ",pfs,sep="",collapse=" ")
  system(paste("cdo mergetime ",pclip," ",targetdir,"totpr_",i,".nc",sep=""))
  ## Temp
  tfs=fs$path[fs$chunk==i&grepl("tas",fs$path)]  #get temp files
  tclipmax=paste("-daymax -shifttime,-5hour -sellonlatbox,",ROI," ",tfs,sep="",collapse=" ")
  tclipmin=paste("-daymin -shifttime,-5hour -sellonlatbox,",ROI," ",tfs,sep="",collapse=" ")
  tclipmean=paste("-daymean -shifttime,-5hour -sellonlatbox,",ROI," ",tfs,sep="",collapse=" ")
  system(paste("cdo mergetime ",tclipmax," ",targetdir,"maxtas_",i,".nc",sep=""))
  system(paste("cdo mergetime ",tclipmin," ",targetdir,"mintas_",i,".nc",sep=""))
  system(paste("cdo mergetime ",tclipmean," ",targetdir,"meantas_",i,".nc",sep=""))
  ## Merge them all together
  system(paste("cdo merge -chname,tas,tmean ", targetdir,"meantas_",i,
               ".nc -chname,tas,tmax ",targetdir,"maxtas_",i,
               ".nc -chname,tas,tmin ",targetdir,"mintas_",i,
               ".nc -chname,pr,ptot ",targetdir,"totpr_",i,
               ".nc ",targetdir,i,".nc",sep=""))
  ## Update attributes
  system(paste("ncatted -O -a units,ptot,c,c,\"mm/day\" ",targetdir,i,".nc",sep="")) 
  system(paste("ncatted -O -a long_name,ptot,o,c,\"Total Daily Precipitation\" ",targetdir,i,".nc",sep="")) 
  system(paste("ncatted -O -a long_name,tmin,o,c,\"Minimum Daily Temperature\" ",targetdir,i,".nc",sep="")) 
  system(paste("ncatted -O -a long_name,tmax,o,c,\"Maximum Daily Temperature\" ",targetdir,i,".nc",sep="")) 
  system(paste("ncatted -O -a long_name,tmean,o,c,\"Mean Daily Temperature\" ",targetdir,i,".nc",sep="")) 
  ## Project to WGS84
  gridfile="/media/Data/Work/Regional/NewEngland/ClimateProjections/NARCCAP/CDO_grid_lonlat"
  system(paste("cdo remapbic,",gridfile," ",targetdir,i,".nc ",targetdir,"wgs84_",i,".nc",sep=""))
  ## Remove temporary files
  sepfiles=c(paste(targetdir,"meantas_",i,".nc",sep=""),paste(targetdir,"maxtas_",i,".nc",sep=""),
    paste(targetdir,"mintas_",i,".nc",sep=""),paste(targetdir,"totpr_",i,".nc",sep=""))
    file.remove(sepfiles)
  print(paste("################################        Finshed ",i,"    ################################################"))
}

