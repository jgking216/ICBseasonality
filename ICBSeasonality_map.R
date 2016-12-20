fdir= "C:\\Users\\Buckley\\Google Drive\\Seasonality\\"

library(ncdf)
library(abind) #bind arrays
library(dplyr)

##source functions
setwd(paste(fdir,"analysis\\",sep="") )
source("DDFunctions.R")

#find developmental traits to use
setwd(paste(fdir,"out\\",sep="") )
dddat= read.csv("SeasonalityDatabase_MASTER_Nov2016.csv")

#set to rounded median
T0= 10
DDD= 300

#----------------------------
degree.days.mat.na= function(Tdat,LDT){
  dd=NA
  
  Tmin=Tdat[1]
  Tmax=Tdat[2]
  
  if(!is.na(Tmin) &!is.na(Tmax) ){
  # entirely above LDT
  if(Tmin>=LDT) {dd=(Tmax+Tmin)/2-LDT}
  
  # intercepted by LDT
  ## for single sine wave approximation
  if(Tmin<LDT && Tmax>LDT){
    alpha=(Tmax-Tmin)/2
    theta1=asin(((LDT-(Tmax+Tmin))/alpha)*pi/180)
    dd=1/pi*(((Tmax+Tmin)/2-LDT)*(pi/2-theta1)+alpha*cos(theta1))
    if(!is.na(dd))if(dd<0){dd=0}
  } #matches online calculation
  
  # entirely below LDT
  if(Tmax<=LDT){ dd=0}
  } #check NA
    
  return(dd)
}

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
#PLOT SEASONALITY

tmax1= tmax[19,15,1:365]
tmin1= tmin[19,15,1:365]
tmean1= colMeans( rbind(tmin1,tmax1))

tmean1= cbind(tmean1, 1:365)

plot(tmean1,type="l", ylab="Temperature (°C)", cex.lab=1.2, cex=1.2)

#----------------------------
#CALCULATE DEGREE DAYS

#combine tmin and tmax
txn= abind(tmin, tmax, along=4)
dim(txn) #  96   73 3652    2

#change NaN to NA
txn[is.nan(txn)] = NA  

#test on subset without NAs
dd= apply(txn[10:15,11:15,1,], MARGIN=c(1,2), FUN=degree.days.mat.na, LDT=-5 )
#test with NAs 
dd= apply(txn[5:12,10:15,1,], MARGIN=c(1,2), FUN=degree.days.mat.na, LDT=-5 )

#run on all
dd= apply(txn, MARGIN=c(1,2,3), FUN=degree.days.mat.na, LDT=T0 )

#slow loop style 96, 73, 365
phen.dd= array(NA, dim=c(dim(dd)[1:2],20))
ngen.dd= array(NA, dim=dim(dd)[1:2])

for(rowk in 1:nrow(dd)){
  for(colk in 1:ncol(dd)){
dd1= dd[rowk,colk,]

if(!all(is.na(dd1))){

dd1= as.data.frame(dd1)
dd1$year= years
dd1$j= js

#revise to handle S hemisphere
if(lat[colk]<0) 
{inds.jd= which(dat$j>181)
inds.jj= which(dat$j<182)

dat[inds.jd, "j"] = dat$j[inds.jd] -181
dat[inds.jj, "j"] = dat$j[inds.jj] +184
} #end fix s hemi

dat = dd1 %>% group_by(year) %>% arrange(j)  %>%  mutate(cs = cumsum(dd1))
#fix for S hemisphere

#Egg to adult DD, First date beyond threshold
for(genk in 1:20){
  
  #Find phenology of generations
  dat1 = dat %>%  group_by(year) %>% slice(which.max(cs> (genk* DDD) ))
  phen.dd[rowk,colk,genk]= mean(dat1$j, na.rm=TRUE)
  
} #end gen loop

#number generations
phen.dd[rowk,colk,which(phen.dd[rowk,colk,]==1)]= NA
ngen.dd[rowk,colk]= length(na.omit(phen.dd[rowk,colk,] ))

} #end check na
} #end loop coloumns
} #end loop rows
  

#Calculate cumulative sum cumsum()
#code template: dat = dat %>% group_by(year) %>% arrange(date) %>% mutate(cs = cumsum(dd)) 
##dd.cumsum = dd1 %>% group_by(year) %>% arrange(date) %>% mutate_each( cumsum(dd)) 

#-----------------------------------------------------
#PLOT
library(lattice)
grid <- expand.grid(lon = lon, lat = lat)

#plot phenology
image(lon, rev(lat), phen.dd[,ncol(phen.dd):1,1] )
levelplot(phen.dd[,,1] ~ lon * lat, data=grid, col.regions = terrain.colors(100) )

#plot number generations
levelplot(ngen.dd ~ lon * lat, data=grid, col.regions =  terrain.colors(100) )

