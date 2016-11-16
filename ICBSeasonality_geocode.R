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
#--------------------------------
#LOAD AND CLEAN DATA

setwd(paste(fdir,"Database\\",sep="") )

dat= read.csv("ThermalDatabase_updatelocs.csv")
#dat= read.csv("Thermal_database_for_PRATIQUE_20June2016.csv")
#dat= dat[,1:48]
#dat= dat[1:100,] ##subset for testing

#restrict to sites with LDT (BDT.C) and dd threshold (EADDC) 
dat$BDT.C= as.numeric(dat$BDT.C)
dat$EADDC= as.numeric(dat$EADDC)
dat= dat[!is.na(dat$BDT.C) & !is.na(dat$EADDC),]

#locs<- str_split_fixed(dat$Location, ",", 3)
#locs<- apply(locs, MARGIN=2, FUN="str_trim", side = c("both"))
#-------------------------------
#AGGREGATE ACROSS REPLICATES, also use bootstrapping?

#omit
dat= dat[-which(dat$omit=="y"),]

#group by index
dat= dat %>% group_by(index) %>% summarise(Species=head(Species)[1],Order=head(Order)[1],Family=head(Family)[1],Genus=head(Genus)[1],Species.1=head(Species.1)[1],BDT.C=mean(BDT.C), UDT.C=mean(UDT.C),EADDC=mean(EADDC),EEDDC=mean(EEDDC), pest=mean(pest), aquatic=mean(aquatic),pupal=mean(pupal),Location=head(Location)[1],lon=mean(lon),lat=mean(lat), Location_new=head(Location_new)[1],lon_new=mean(lon_new),lat_new=mean(lat_new), colony=head(colony)[1], quality=head(quality)[1])
                                            
dat= as.data.frame(dat)

#------------------------------  
#GEOREFERENCE ##NOTE DAILY MAX CALLS

#Find sites without Location_new
#Use initial location information
inds= which(is.na(dat$Location_new))

dat$Location_new[inds]= dat$Location[inds]
dat$lon_new[inds]= dat$lon[inds]
dat$lat_new[inds]= dat$lat[inds]
dat$quality[inds]= 4

#Find sites without lat lon
inds= which(is.na(dat$lon_new))

locs= unique(as.character(dat$Location_new[inds]))
## GEOCODING
#GEOCODE USING GOOGLE API
#locs.geo2<- dismo::geocode(locs, oneRecord=TRUE)
#locs.geo2<- ggmap::geocode(locs, output = "latlon")
#GEOCODE USING OPEN STREET MAPS #FAILS IF DOESN'T FIND MATCH
#locs.geo <- osm_geocode(locs, key="nG1nq56AVcJD3lFSQbfHxJDJzJP1HwVa") #key for huckley@uw.edu, login:huckley, password:colias

## FOR GEOCODE
match1= match(as.character(dat$Location_new[inds]), locs.geo2$originalPlace )
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

#update old location info with new values
dat$Location= dat$Location_new
#make numeric, issue ?? with Chinese data as factor
dat$lon= as.numeric(dat$lon)
dat$lat=as.numeric(dat$lat)
dat$lon= dat$lon_new
dat$lat= dat$lat_new

#---------------------------------
dddat=dat
#issue with Chinese data ???
#dddat[which(dddat$lon==round(dddat$lon) ),]

##WRITE OUT
#setwd(paste(fdir,"out\\",sep="") )
#write.csv(dddat, "dddat_14Nov2017.csv" )

#write out locations
#write.csv(locs.geo2, "georef_14Nov2017.csv" )

