fdir= "C:\\Users\\Buckley\\Google Drive\\Seasonality\\"

library(ncdf)
library(abind) #bind arrays

##source functions
setwd(paste(fdir,"analysis\\",sep="") )
source("DDFunctions.R")

#find developmental traits to use
setwd(paste(fdir,"out\\",sep="") )
dddat= read.csv("SeasonalityDatabase_MASTER_Nov2016.csv")

#set to rounded median
T0= 10
DDD= 300

#-------------------------
#LOAD NETCDF DATA
#DATA FROM HERE: http://www.metoffice.gov.uk/hadobs/hadghcnd/download.html

setwd("C:\\Users\\Buckley\\Google Drive\\Seasonality\\data\\TminTmaxGridded\\HadGHCND_TXTN_acts_1950-2014_15102015.nc\\")

txn<-  open.ncdf("HadGHCND_TXTN_acts_1971-1980_15102015.nc") # opens netcdf file example.nc as R object
print (txn) # show information about the structure of the data

lat = get.var.ncdf(txn, "latitude")
lon = get.var.ncdf(txn,"longitude")
time = get.var.ncdf(txn,"time") 

#make vector of years and calendar days corresponding to time dimersions
years= c(rep(1971,365),rep(1972,366),rep(1973,365),rep(1974,365),rep(1975,365),rep(1976,366),rep(1977,365),rep(1978,365),rep(1979,365),rep(1980,365) )
js= c(1:365,1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:365 )

tmax = get.var.ncdf(txn,varid="tmax") #dim:  96   73 3652
tmin = get.var.ncdf(txn,varid="tmin")

#remove ncdf
rm(txn)

#PLOT #flipped
image(tmax[,,1])

#can use vector notation to access netcdf: e.g., tmax[,,1]

#----------------------------
#CALCULATE DEGREE DAYS

#combine tmin and tmax
txn= abind(tmin, tmax, along=4)
dim(txn) #  96   73 3652    2

#change NaN to NA
txn[is.nan(txn)] = NA  

#test on subset without NAs
dd= apply(txn[10:15,11:15,1,], MARGIN=c(1,2), FUN=degree.days.mat, LDT=-5 )
#test with NAs 
dd= apply(txn[5:12,10:15,1,], MARGIN=c(1,2), FUN=degree.days.mat, LDT=-5 )

#run on all
#!NEED TO MODIFY FUNCTION TO DEAL WITH NAs
dd= apply(txn, MARGIN=c(1,2,3), FUN=degree.days.mat, LDT=T0 )

#Calculate cumulative sum cumsum()
#code template: dat = dat %>% group_by(year) %>% arrange(date) %>% mutate(cs = cumsum(dd)) 

#Find phenology of 1st generation

#Find number generations



