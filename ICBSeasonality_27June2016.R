#CALCULATE DD

count= function(x) length(na.omit(x))

fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ICBseasonality\\"

library(stringr)
library(dismo)
library(rnoaa) #http://recology.info/2015/07/weather-data-with-rnoaa/
library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyr)
library(maptools)

#Make wrapper for Julian function
julian.wrap=function(date1){ 
  #date1= as.Date(date1, "%Y-%m-$d")
  origin1 = as.Date(paste("1-1-",as.numeric(format(date1, "%Y")), sep=""), format="%m-%d-%Y")
  julian(date1, origin=origin1)
}

#source functions
setwd(paste(fdir,"analysis\\",sep="") )
source("DDFunctions.R")
#--------------------------------
#LOAD AND CLEAN DATA

setwd(paste(fdir,"data\\",sep="") )

dat= read.csv("Thermal_database_for_PRATIQUE_20June2016.csv")
dat= dat[,1:48]
#dat= dat[1:100,] ##subset for testing

#restrict to sites with LDT (BDT.C) and dd threshold (EADDC) 
dat$BDT.c= as.numeric(dat$BDT.C)
dat$EADDC= as.numeric(dat$EADDC)
dat= dat[!is.na(dat$BDT.c) & !is.na(dat$EADDC),]

#locs<- str_split_fixed(dat$Location, ",", 3)
#locs<- apply(locs, MARGIN=2, FUN="str_trim", side = c("both"))
  
#georeference ##NOTE DAILY MAX CALLS
locs= unique(as.character(dat$Location))
locs.geo<- geocode(locs, oneRecord=FALSE)

match1= match(as.character(dat$Location), locs.geo$originalPlace )
dat$lon= locs.geo$longitude[match1]
dat$lat= locs.geo$latitude[match1]

dddat=dat

##WRITE OUT
#setwd(paste(fdir,"out\\",sep="") )
#write.csv(dddat, "dddat.csv" )

##READ BACK IN
#dddat= read.csv("dddat.csv")

#Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]

#-----------------------------
#Plot locations
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-180,180), ylim=c(-60,60))
points(dddat$lon, dddat$lat, col="red", cex=1)

#--------------------------------
#FIND CLOSEST GHCN STATIONS

#Read GGCN database inventory http://www.ncdc.noaa.gov/oa/climate/ghcn-daily/
setwd(paste(fdir,"data\\",sep=""))
stations=read.table("ghcnd-inventory.txt")
names(stations)=c("Id","Lat","Lon","Var","FirstYear","LastYear")
stations.un= unique(stations$Id)
stats.elev= read.csv("ghcnd-stations.csv")

#Restrict stations to those with recent data and at least 20 years of data
stations= stations[which(stations$LastYear>2010 & (stations$LastYear-stations$FirstYear)>20 ),] 

#Check have TMIN and TMAX dat
stations= stations[which(stations$Var=="TMAX" | stations$Var=="TMIN"),]
stations= spread(stations,Var, Var)
stations= stations[which(stations$TMAX=="TMAX" & stations$TMIN=="TMIN"),]

stat.coords= cbind(stations$Lon, stations$Lat)

#--------------------------------
#SET UP DATA STORAGE
phen.dat= array(NA, dim=c(nrow(dddat), 101, 20, 7), dimnames=list(NULL,as.character(1915:2015),as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q")) )
ngens= array(NA, dim=c(nrow(dddat), 101), dimnames=list(NULL,as.character(1915:2015)))
#----------------------------------
#ANALYSIS

for(stat.k in 1:nrow(dddat) ){  #1:nrow(sites)
  if( all(!is.na(dddat[stat.k,c("lon","lat")] )) ){
  min.dist<- order(spDistsN1(stat.coords, as.numeric(dddat[stat.k,c("lon","lat")]), longlat = TRUE))[1:100]
  min.site= stations$Id[min.dist]
  
  ind=0
  years=NA
  
  while(length(years)<10 & ind<100){
    ind= ind + 1
    tmax=try( ghcnd_search(min.site[ind], var = "TMAX") )
    tmin=try( ghcnd_search(min.site[ind], var = "TMIN") )
    while( (is.null(nrow(tmax$tmax)) | is.null(nrow(tmin$tmin))) & ind<100){ ind=ind+1
      tmax=try( ghcnd_search(min.site[ind], var = "TMAX") )
      tmin=try( ghcnd_search(min.site[ind], var = "TMIN") )
    }
    
  if(!is.null(nrow(tmax$tmax)>0) ){ #CHECK DATA
    #combine tmin and tmax
    match1= match(tmax$tmax$date, tmin$tmin$date)
    is.matched= !is.na(match1)
    
    dat= tmax$tmax
    dat$tmin[is.matched]= tmin$tmin$tmin[match1[is.matched]]
    
    #split date 
    date= as.Date(dat$date, "%Y-%m-$d")
    dat$year=as.numeric(format(date, "%Y"))
    dat$month=as.numeric(format(date, "%m"))
    dat$day=as.numeric(format(date, "%d"))
    dat$j= unlist(lapply(date, FUN="julian.wrap"))+1
    
    #Format data
    dat$tmax[which(dat$tmax==-9999)]= NA
    dat$tmax= dat$tmax/10 #correct for tenths of degrees or mm
    dat$tmin[which(dat$tmin==-9999)]= NA
    dat$tmin= dat$tmin/10 #correct for tenths of degrees or mm
    
    ## FIND YEARS WITH NEARLY COMPLETE DATA
    dat.agg= aggregate(dat, list(dat$year),FUN=count)
    years= dat.agg$Group.1[which(dat.agg$tmax>300)]
    dat= dat[which(dat$year %in% years),]
    
  } #END CHECK NULL
  } #END WHILE YEARS
    
  #RECORD SITE DATA
  dddat$Id[stat.k]= as.character(min.site[ind])
  dddat$st.lat[stat.k]= stations$Lat[min.dist[ind]]
  dddat$st.lon[stat.k]= stations$Lon[min.dist[ind]]
  dddat$st.elev[stat.k]= stations$elev[min.dist[ind]]
    
    #--------------------------------------
    #INTERPOLATE MISSING DATA
    dat$tmin= na.approx(dat$tmin, maxgap=7, na.rm = FALSE)
    dat$tmax= na.approx(dat$tmax, maxgap=7, na.rm = FALSE)
    #Tmean
    dat$tmean= (dat$tmin+dat$tmax)/2
    
    #cut missing data
    dat= dat[which(!is.na(dat$tmin) & !is.na(dat$tmax)),]
    
    #CALCULATE DEGREE DAYS
    dat$dd= apply( dat[,c("tmin","tmax")], MARGIN=1, FUN=degree.days.mat, LDT=dddat$BDT.C[stat.k] )
    #! ADD UDT
    
    #Estimate phenology across years based on DD accumulations
    #cumsum within groups
    dat = dat %>% group_by(year) %>% arrange(date) %>% mutate(cs = cumsum(dd))
    
    #Egg to adult DD, First date beyond threshold
    for(genk in 1:20){
    
    phen= dat %>%  group_by(year) %>% slice(which.max(cs> genk* dddat$EADDC[stat.k]))
    
    years.ind=1915:2015
    year.loop= unique(phen$year)
    year.loop= year.loop[which(year.loop>1914) ]
    
    for(yeark in year.loop){
    
      j= as.numeric(phen[phen$year==yeark,"j"])
     
       dat.yr= dat[dat$year==yeark ,]
    #TEMPS: 2 WEEKS AT EMERGENCE
      temps= as.numeric(unlist(dat.yr[dat.yr$j %in% (j-7):(j+7),"tmean"]))
      phen.dat[stat.k,which(years.ind==yeark),genk,"phen"]= j
      #MEAN
      phen.dat[stat.k,which(years.ind==yeark),genk,"Tmean.e"]= mean(temps)
      #STDEV
      phen.dat[stat.k,which(years.ind==yeark),genk,"Tsd.e"]= sd(temps)
      #10% QUANTILE
      phen.dat[stat.k,which(years.ind==yeark),genk,"T10q.e"]= quantile(temps, 0.1)
        
    #TEMPS: ACROSS GENERATIONS
      #FIND GEN TIME
      j.gs= as.numeric(ifelse(genk==1, dat.yr[which.min(dat.yr$dd>0),"j"], phen.dat[stat.k,which(years.ind==yeark),(genk-1),"phen"] ))
      temps= as.numeric(unlist(dat.yr[dat.yr$j %in% j.gs:j,"tmean"]))
      #MEAN
      phen.dat[stat.k,which(years.ind==yeark),genk,"Tmean"]= mean(temps)
      #STDEV
      phen.dat[stat.k,which(years.ind==yeark),genk,"Tsd"]= sd(temps)
      #10% QUANTILE
      phen.dat[stat.k,which(years.ind==yeark),genk,"T10q"]= quantile(temps, 0.1)
    
    } #end check data again
    } #end GEN LOOP
      
    
      #Number generations by year
    phen.dat2= as.data.frame(phen.dat[stat.k,,,"phen"])
    ngens[stat.k,]= apply(phen.dat2, MARGIN=1, FUN=function(x)length(which(x>1)) )  
    
    print(stat.k)
    
  } #check coordinates
} #end site (stat.k) loop

##SAVE OUTPUT
#setwd(paste(fdir,"out\\",sep="") )
#save(phen.dat, file = "phendat.rda")

##READ BACK IN
#load(file = "phendat.rda")

#--------------------------------
#PLOTS
#1. How do phenology shifts vary across phylogeny and latitude?
#2. Differences in temperature across generations?
#3. Do phenological advancements result in cooler late in development (pupal?) temperatures?
#4. Do phenological shifts result in insects experiencing more variable environments? Cool extremes?
#5. How does the potential for additional generations vary across latitudes?

#VAR: "phen"    "Tmean.e" "Tsd.e"   "T10q.e"  "Tmean"   "Tsd"     "T10q"   

#---------------------
#Fig 1. D0, DD reg, Number generations
#Slopes across generations vs latitude: temps, generation lengths




#---------------------
#Fig 2. Slopes through time vs latitude: adult phenology, number generations; dev temperature, dev temp sd; adult temperature, adult temp sd (only phenological advancements). [dev 10th quantile]

calcm= function(x, ys=ys){
  yd= which( !is.na(x) )
  if(length(yd)>0) mod1= lm( as.numeric(x[yd])~ys[yd] )   
  return(tryCatch(coef(summary(mod1))[2, ], error=function(e) c(NA,NA,NA,NA)))
}

ylabs= c("Adult phenology (J)","Number generations", "Developmental temperature mean (°C)","Developmental temperature sd (°C)", "Adult temperature mean (°C)", "Adult temperature sd (°C)")

setwd(paste(fdir,"figures\\",sep="") )
pdf("Phen_latPlots.pdf", height = 14, width = 10)
par(mfrow=c(3,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

for(i in 1:6){

#Calculate phenology shifts across years
if(i==1) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"phen"]))

#Number generations
if(i==2) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],ngens))
  
#Generation temps across years for ? generation
if(i==3) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tmean"]))

if(i==4) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tsd"]))

#Adult temps across years for ? generation
if(i==5) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tmean.e"]))

if(i==6) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tsd.e"]))

colnames(phen.dat2)[1]="siteID"

# CALCULATE SLOPES
ys= as.numeric( colnames(phen.dat2)[6:ncol(phen.dat2)])

slopes= apply(phen.dat2[6:ncol(phen.dat2)], MARGIN=1, FUN=calcm, ys=ys)
phen.dat3= cbind(phen.dat2[,1:5],t(slopes))
names(phen.dat3)[6:9]=c("Estimate","Std.Error","t value","P")

#restrict to significant shifts
phen.sig= phen.dat3[which(phen.dat3$P<0.05),]

plot(abs(phen.sig$lat), phen.sig$Estimate, ylab=ylabs[i],xlab="Absolute latitude (°)")
abline(h=0)

#------------------------------
#Store plots over time

} #end loop metrics

dev.off()
#plot(density(phen.dat3$Estimate))

#======================================
#PLOT TRENDS OVER TIME
#drop rows with all NAs
drop.row=apply(phen.dat2, MARGIN=1, FUN=function(x)all(is.na(x[6:length(x)])) )
if(!all(drop.row==FALSE)) phen.dat2= phen.dat2[-which(drop.row==TRUE),]

phen.dat3= gather(phen.dat2, "year", "phen",6:106)
phen.dat4= phen.dat3[which(!is.na(phen.dat3$phen)) ,]
#--------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("Phen_yearPlots.pdf", height = 14, width = 14)
par(mfrow=c(2,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

plot1 = ggplot(phen.dat3, aes(x=year, y=phen, group=siteID, color=lat)) +geom_line() #+ylim(0, 5)
plot1 = ggplot(phen.dat3, aes(x=year, y=phen, group=siteID, color=abs(lat) )) +geom_smooth(method=lm, se=FALSE)+ labs(y = ylab) #+ylim(0, 5)
plot(plot1)

dev.off()


#==================================================================






