#CALCULATE DD

count= function(x) length(na.omit(x))

fdir= "C:\\Users\\Buckley\\Google Drive\\Seasonality\\"

library(stringr)
library(dismo)
library(rnoaa) #http://recology.info/2015/07/weather-data-with-rnoaa/
library(zoo)
library(tidyr)
library(ggplot2)
library(tidyr)
library(maptools)
#for mapping
library(ggmap)
library(maps)
library(mapdata)
library(colorRamps)     # for matlab.like(...)
library(grid)
library(plyr); library(dplyr)
#devtools::install_github("hrbrmstr/nominatim")
library(nominatim) #for osm geocode

#Make wrapper for Julian function
julian.wrap=function(date1){ 
  #date1= as.Date(date1, "%Y-%m-$d")
  origin1 = as.Date(paste("1-1-",as.numeric(format(date1, "%Y")), sep=""), format="%m-%d-%Y")
  julian(date1, origin=origin1)
}

##source functions
#setwd(paste(fdir,"analysis\\",sep="") )
#source("DDFunctions.R")

#calculate slope
calcm= function(x, ys=ys){
  yd= which( !is.na(x) )
  if(length(yd)>0) mod1= lm( as.numeric(x[yd])~ys[yd] )   
  return(tryCatch(coef(summary(mod1))[2, ], error=function(e) c(NA,NA,NA,NA)))
}

calcm.poly= function(x, ys=ys){
  yd= which( !is.na(x) )
  if(length(yd)>4) mod1= lm( as.numeric(x[yd])~ poly(ys[yd],2) )   
  return(tryCatch( c(coef(summary(mod1))[2:3, ]), error=function(e) c(NA,NA,NA,NA,NA,NA,NA,NA)))
}
#===========================================
#LOAD AND CLEAN DATA

setwd(paste(fdir,"Database\\",sep="") )

dat= read.csv("ThermalDatabase_Nov2016.csv")

#restrict to sites with LDT (BDT.C) and dd threshold (EADDC) 
dat$BDT.C= as.numeric(dat$BDT.C)
dat$EADDC= as.numeric(dat$EADDC)
dat= dat[!is.na(dat$BDT.C) & !is.na(dat$EADDC),]

### NO LONGER USE
##Find sites without Location_new
##Use initial location information
#inds= which(is.na(dat$Location_new))

#dat$Location_new[inds]= dat$Location[inds]
#dat$lon_new[inds]= dat$lon[inds]
#dat$lat_new[inds]= dat$lat[inds]
#dat$quality[inds]= 4

#=============================================
## GEOREFERENCE
dat$Location_new= as.character(dat$Location_new)

#Find locations needing georeference
inds= which(is.na(dat$lon_new))
locs= as.character(dat$Location_new[inds])

#source from current database
inds2= which(!is.na(dat$lon_new))
locs2= dat[inds2, c("Location_new","lon_new","lat_new")]
nodups= which(duplicated(locs2$Location_new)==FALSE)
locs2= locs2[nodups,]

match1= match(locs, locs2$Location_new)
matched= inds[which(!is.na(match1))]
#check match
check= cbind(dat$Location_new[matched], locs2$Location_new[na.omit(match1)])

dat$lat_new[matched]= locs2$lat_new[na.omit(match1)]
dat$lon_new[matched]= locs2$lon_new[na.omit(match1)]

#--------------------------------------------
#Find locations needing georeference
inds= which(is.na(dat$lon_new))
locs= as.character(dat$Location_new[inds])

#source from previous georeference
setwd(paste(fdir,"out\\",sep="") )

geos= read.csv("georef_14Nov2016.csv")
geos$originalPlace= as.character(geos$originalPlace)
nodups= which(duplicated(geos$originalPlace)==FALSE)
geos= geos[nodups,]

match1= match(locs, geos$originalPlace)
matched= inds[which(!is.na(match1))]
#check match
check= cbind(dat$Location_new[matched], geos$originalPlace[na.omit(match1)])

dat$lat_new[matched]= geos$latitude[na.omit(match1)]
dat$lon_new[matched]= geos$longitude[na.omit(match1)]

#------------------------------  
#GEOREFERENCE ##NOTE DAILY MAX CALLS

#Find sites without lat lon
inds= which(is.na(dat$lon_new))

locs= unique(as.character(dat$Location_new[inds]))
## GEOCODING
#GEOCODE USING GOOGLE API
#locs.geo2<- dismo::geocode(locs, oneRecord=TRUE) ##USE
#locs.geo2<- ggmap::geocode(locs, output = "latlon")
#GEOCODE USING OPEN STREET MAPS #FAILS IF DOESN'T FIND MATCH
#locs.geo <- osm_geocode(locs, key="nG1nq56AVcJD3lFSQbfHxJDJzJP1HwVa") #key for huckley@uw.edu, login:huckley, password:colias

## FOR GEOCODE
match1= match(as.character(dat$Location_new[inds]), locs.geo2$originalPlace )
locs.geo2$originalPlace= as.character(locs.geo2$originalPlace)
#check match
check= cbind(dat$Location_new[inds], locs.geo2$originalPlace[match1])

dat$lon_new[inds]= locs.geo2$longitude[match1]
dat$lat_new[inds]= locs.geo2$latitude[match1]

#rerun to deal with nas
#Find sites without lat lon
inds= which(is.na(dat$lon_new))

locs2= unique(as.character(dat$Location_new[inds]))
#locs.geo3<- dismo::geocode(locs2, oneRecord=TRUE)
match1= match(as.character(dat$Location_new[inds]), locs.geo3$originalPlace )
dat$lon_new[inds]= locs.geo3$longitude[match1]
dat$lat_new[inds]= locs.geo3$latitude[match1]

#---------------
##FOR OSM
#match1= match(as.character(dat$Location_new[inds]), locs)
#dat$lon_new[inds]= locs.geo$lon[match1]
#dat$lat_new[inds]= locs.geo$lat[match1]

##update old location info with new values
#dat$Location= dat$Location_new
##make numeric, issue ?? with Chinese data as factor
#dat$lon= as.numeric(dat$lon)
#dat$lat=as.numeric(dat$lat)
#dat$lon= dat$lon_new
#dat$lat= dat$lat_new

#---------------------------------
dddat=dat
#issue with Chinese data ???
#dddat[which(dddat$lon==round(dddat$lon) ),]

##WRITE OUT
#setwd(paste(fdir,"out\\",sep="") )
#write.csv(dddat, "dddat_14Nov2016.csv" )

#write out locations
#write.csv(locs.geo2, "georef_15Nov2016.csv" )

#=========================================================
##UPDATE DATABASE

setwd(paste(fdir,"out\\",sep="") )

dat= read.csv("SeasonalityDatabase_Nov2016.csv")

#AGGREGATE ACROSS REPLICATES, also use bootstrapping?
##omit
#dat= dat[-which(dat$omit=="y"),]

#add index
ind_new= paste(dat$"Species",dat$"Author",dat$"Year",dat$"Location_new",sep="_")
ind= unique(ind_new)
dat$index= match(ind_new, ind)

#add parasitoid data
parasite= read.csv("Parasitoids_JGK.csv")

match1= match(parasite$DatabaseOrder, dat$DatabaseOrder)
check= cbind(dat$DatabaseOrder[match1], parasite$DatabaseOrder)
dat$parasitoid[match1]= parasite$parasitoid

## WRITE OUT
#write.csv(dat, "SeasonalityDatabase_MASTER_Nov2016.csv")

#=========================================================
