#CALCULATE DD

count= function(x) length(na.omit(x))

fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ICBseasonality\\"

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

#Make wrapper for Julian function
julian.wrap=function(date1){ 
  #date1= as.Date(date1, "%Y-%m-$d")
  origin1 = as.Date(paste("1-1-",as.numeric(format(date1, "%Y")), sep=""), format="%m-%d-%Y")
  julian(date1, origin=origin1)
}

#source functions
setwd(paste(fdir,"analysis\\",sep="") )
source("DDFunctions.R")

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

setwd(paste(fdir,"data\\",sep="") )

dat= read.csv("Thermal_database_for_PRATIQUE_20June2016.csv")
dat= dat[,1:48]
#dat= dat[1:100,] ##subset for testing

#restrict to sites with LDT (BDT.C) and dd threshold (EADDC) 
dat$BDT.C= as.numeric(dat$BDT.C)
dat$EADDC= as.numeric(dat$EADDC)
dat= dat[!is.na(dat$BDT.C) & !is.na(dat$EADDC),]

#locs<- str_split_fixed(dat$Location, ",", 3)
#locs<- apply(locs, MARGIN=2, FUN="str_trim", side = c("both"))
  
#georeference ##NOTE DAILY MAX CALLS
locs= unique(as.character(dat$Location))
#locs.geo<- geocode(locs, oneRecord=FALSE)

match1= match(as.character(dat$Location), locs.geo$originalPlace )
dat$lon= locs.geo$longitude[match1]
dat$lat= locs.geo$latitude[match1]

dddat=dat

##WRITE OUT
#setwd(paste(fdir,"out\\",sep="") )
#write.csv(dddat, "dddat.csv" )

##READ BACK IN
#dddat= read.csv("dddat.csv")

##Restrict to dat with lat / lon
#dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]

#Add pest data
setwd(paste(fdir,"data\\PlantWisePests\\",sep="") )
#sdir= paste(fdir,"data\\PlantWisePests\\",sep="")
#files<-list.files(sdir,pattern="\\.csv$")

#pest= lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE))
#pest= unique(unlist(pest, recursive=TRUE, use.names = FALSE))
#write.csv(pest, "pest_global.csv" )
pest= read.csv("pest_global.csv" )

pests= paste(pest$gen, pest$spec, sep=" ")
match1= match(as.character(dddat$Species), pests )
matched= which(!is.na(match1))
dddat$pest=0
dddat$pest[matched]=1

#-----------------------
#REMOVE MITES
dddat= dddat[-which(dddat$Order=="Acari"),]

#-----------------------
#Code Aquatic
dddat$aquatic=0

dddat[which(dddat$Order %in% c("Ephemeroptera","Odonata", "Plecoptera") ),"aquatic"]=1

dddat[which(dddat$Order=="Hemiptera" & dddat$Family=="Gerridae"),"aquatic"]=1

dddat[which(dddat$Order=="Diptera" & dddat$Family %in% c("Chironomidae","Culicidae", "Simulidae","Tabanidae")),"aquatic"]=1

dddat[which(dddat$Order=="Diptera" & dddat$Family=="Sciomyzidae" & dddat$Genus=="Tabanus"),"aquatic"]=1
              
#-----------------------            
#Code Holometabolous
dddat$pupal=0

dddat[which(dddat$Order %in% c("Megaloptera","Raphidioptera", "Neuroptera","Coleoptera","Strepsiptera","Diptera","Mecoptera","Siphonaptera","Trichoptera","Lepidoptera","Hymenoptera","Miomoptera")),"pupal"]=1

##WRITE OUT
#setwd(paste(fdir,"out\\",sep="") )
#write.csv(dddat, "dddat.csv" )

#READ BACK IN
#setwd(paste(fdir,"out\\",sep="") )
#dddat= read.csv("dddat.csv" )

#-----------------------
#TABLE DATA
#GROUP BY ORDER, FAMILY
with(dddat, table(Order))

fam= with(dddat, table(Family))
fam[fam>30]
# Aphelinidae: tiny parasitic wasps
# Aphididae: aphids
# Braconidae: parasitoid wasps
# Coccinellidae: lady bugs
# Drosophilidae
# Noctuidae: owlet moths
# Pyralidae: snout moths
# Tetranychidae: spider moths

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
phen.dat= array(NA, dim=c(nrow(dddat), 101, 20, 10), dimnames=list(NULL,as.character(1915:2015),as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen", "Tmean.e.fixed", "Tmean.fixed")) )
ngens= array(NA, dim=c(nrow(dddat), 101), dimnames=list(NULL,as.character(1915:2015)))
phen.fixed= array(NA, dim=c(nrow(dddat), 20, 8), dimnames=list(NULL,as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen")) )
#----------------------------------
#ANALYSIS

#for(stat.k in 1:nrow(dddat) ){ 
for(stat.k in 1:nrow(dddat) ){
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
    dat.agg= aggregate(dat, list(dat$year),FUN=count)  ### PROBLEM IF RASTER LOADED
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
    #replace eroneous values
    phen[which(phen$j==1),"j"]=NA
    #drop years without generation
    phen= phen[!is.na(phen$j),]
    
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
    
      #duration of generation
      phen.dat[stat.k,which(years.ind==yeark),genk,"DaysGen"]= j-j.gs
      
    } #end year loop
    } #end GEN LOOP
      
      #Number generations by year
    phen.dat2= as.data.frame(phen.dat[stat.k,,,"phen"])
    ngens[stat.k,]= apply(phen.dat2, MARGIN=1, FUN=function(x)length(which(x>1)) )
    
    #Generation j and temp average across years
    for(genk in 1:20){
      j= round( mean(phen.dat[stat.k,,genk,"phen"], na.rm=TRUE) )
      
      if(!is.nan(j)){
      phen.fixed[stat.k,genk,"phen"]= j
      dat.gen= dat[dat$year==yeark ,]
      
      temps= as.numeric(unlist(dat[dat$j %in% (j-7):(j+7),"tmean"]))
      
      #MEAN
      phen.fixed[stat.k,genk,"Tmean.e"]= mean(temps)
      #STDEV
      phen.fixed[stat.k,genk,"Tsd.e"]= sd(temps)
      #10% QUANTILE
      phen.fixed[stat.k,genk,"T10q.e"]= quantile(temps, 0.1)
      
      #---
      #Yearly temp at fixed date
      temps= dat[dat$j %in% (j-7):(j+7),]
      temps1= temps %>% group_by(year) %>% summarise(mean(tmean, na.rm=TRUE))
      temps1= temps1[which(temps1$year %in% 1915:2015),]
      
      match1= match(as.character(unlist(temps1[,1])), names(phen.dat[stat.k,,genk,"Tmean.e.fixed"]))
      phen.dat[stat.k,match1,genk,"Tmean.e.fixed"]= unlist(temps1[,2])
      
      #---
      if(genk>1){
        j_1= round( mean(phen.dat[stat.k,,genk-1,"phen"], na.rm=TRUE) )
        
      temps= as.numeric(unlist(dat[dat$j %in% j_1:j,"tmean"]))
      
      #MEAN
      phen.fixed[stat.k,genk-1,"Tmean"]= mean(temps)
      #STDEV
      phen.fixed[stat.k,genk-1,"Tsd"]= sd(temps)
      #10% QUANTILE
      phen.fixed[stat.k,genk-1,"T10q"]= quantile(temps, 0.1)
      
      #gen length
      phen.fixed[stat.k,genk-1,"DaysGen"]= j-j_1+1
      
      #---
      #Yearly temp at fixed date
      temps= dat[dat$j %in% j_1:j,]
      temps1= temps %>% group_by(year) %>% summarise(mean(tmean, na.rm=TRUE))
      temps1= temps1[which(temps1$year %in% 1915:2015),]
      
      match1= match(as.character(unlist(temps1[,1])), names(phen.dat[stat.k,,genk,"Tmean.fixed"]))
      phen.dat[stat.k,match1,genk-1,"Tmean.fixed"]= unlist(temps1[,2])
      
      } #end gen loop
      #---
      
      } #end check gen exists
      } #end GEN LOOP
    
    print(stat.k)
    
  } #check coordinates
} #end site (stat.k) loop
 
##SAVE OUTPUT
#setwd(paste(fdir,"out\\",sep="") )
#saveRDS(phen.dat, "phendat.rds")
#saveRDS(phen.fixed, "phenfix.rds")
#saveRDS(ngens, "ngens.rds")
#saveRDS(dddat, "dddat_media.rds")

##READ BACK IN
setwd(paste(fdir,"out\\",sep="") )
phen.dat= readRDS("phendat.rds")
phen.fixed= readRDS("phenfix.rds")
ngens= readRDS("ngens.rds")
dddat= readRDS("dddat_media.rds")

#============================================================
#STATISTICS 
#Duration of generations
#Developmental temperatures across generations
#Adult temperatures across generations

#Account for negative generation durations
days.gen= phen.dat[,,,"DaysGen"]
days.gen2= days.gen
days.gen[days.gen<0] <- NA
phen.dat[,,,"DaysGen"]= days.gen

days.gen= phen.fixed[,,"DaysGen"]
days.gen3= days.gen
days.gen[days.gen<0] <- NA
phen.fixed[,,"DaysGen"]= days.gen

#--------------------------------
#PLOTS
#1. How do phenology shifts vary across phylogeny and latitude?
#2. Differences in temperature across generations?
#3. Do phenological advancements result in cooler late in development (pupal?) temperatures?
#4. Do phenological shifts result in insects experiencing more variable environments? Cool extremes?
#5. How does the potential for additional generations vary across latitudes?

#VAR: "phen"    "Tmean.e" "Tsd.e"   "T10q.e"  "Tmean"   "Tsd"     "T10q"   

#-----------------------------
#Fig 0. map, D0, DD reg, Number generations

#Plot locations
#http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

world <- map_data("world")

w1= ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill=NA, color="black") + 
  coord_fixed(1.3) +ylim(c(-55,80))+xlim(c(-170,195))
w2= w1 +xlab("Longitude (°)")+ylab("Latitude (°)")
w3= w2+ geom_point(data = dddat, aes(x = lon, y = lat, color= BDT.C, size=log(EADDC))) + scale_color_gradientn(colours=matlab.like(10))

#-------------------
#D0, DD reg, Number generations

ngens.ave= rowMeans(ngens, na.rm=TRUE)
phen.dat2= as.data.frame( cbind(dddat[,c(1:9,51:52)],ngens.ave))

s1= ggplot()+geom_point(data = phen.dat2, aes(x = lat, y = ngens.ave))+theme_bw()+geom_smooth(data = phen.dat2, formula=y~x, aes(x = lat, y = ngens.ave), method=loess, se=TRUE) + coord_flip()

s2= ggplot()+geom_point(data = phen.dat2, aes(x = lat, y = BDT.C))+theme_bw()+geom_smooth(data = phen.dat2, aes(x = lat, y = BDT.C), method=loess, se=TRUE)+ylim(c(-5,21.2))+ coord_flip()

s3= ggplot()+geom_point(data = phen.dat2, aes(x = lat, y = EADDC))+theme_bw()+geom_smooth(data = phen.dat2, aes(x = lat, y = EADDC), method=loess, se=TRUE)+ylim(c(0,2000))+ coord_flip()

#-----------
setwd(paste(fdir,"figures\\",sep="") )
pdf("DevMap.pdf", height = 5, width =15)
par(mfrow=c(1,1), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,4, widths=c(8,2,2,2) )))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(w3,vp=vplayout(1,1))
print(s2,vp=vplayout(1,2))
print(s3,vp=vplayout(1,3))
print(s1,vp=vplayout(1,4))

dev.off()

#---------------------
#Fig 1. Slopes across generations vs latitude: temps, generation lengths
# POLYNOMIAL FITS

ylabs= c("Developmental T mean (°C)","Developmental T sd (°C)", "Adult T mean (°C)", "Adult T sd (°C)","Generation length (days)")

setwd(paste(fdir,"figures\\",sep="") )
pdf("FixedPhen_latPlots.pdf", height = 14, width = 10)
par(mfrow=c(5,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

for(i in 1:5){
  
  #Generation temps across years for ? generation
  if(i==1) phen.dat2= as.data.frame( cbind(1:nrow(dddat), dddat[,c("Species","Order","lon","lat")],phen.fixed[,,"Tmean"]))
  
  if(i==2) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.fixed[,,"Tsd"]))
  
  #Adult temps across years for ? generation
  if(i==3) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.fixed[,,"Tmean.e"]))
  
  if(i==4) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.fixed[,,"Tsd.e"]))
  
  #Generation length
  if(i==5) phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.fixed[,,"DaysGen"]))
  
  colnames(phen.dat2)[1]="siteID"
  
  #---------------------------------
#!  ## CLUST BY LAT
  #elevation aggregate
  phen.dat2$abslat= abs(phen.dat2$lat)
  phen.dat2$lcut= cut(phen.dat2$abslat, breaks=c(0,25,35,45,90) )
  
  phen.dat4= aggregate(phen.dat2, list(phen.dat2$lcut), FUN=mean, na.rm=TRUE)
  names(phen.dat4)[1]= "latgroup"
  
  #convert to long format
  phen.dat4= gather(phen.dat4, "gen", "T",7:28)
  phen.dat4$gen= as.numeric(phen.dat4$gen)  
  
  p1 = ggplot(phen.dat4, aes(x=gen, y=T, group= latgroup, color=latgroup)) +geom_line() +theme_bw()
  
  #----------------------------------
  
  # CALCULATE SLOPES
  ys= as.numeric( colnames(phen.dat2)[6:ncol(phen.dat2)])
  
  slopes= apply(phen.dat2[6:ncol(phen.dat2)], MARGIN=1, FUN=calcm.poly, ys=ys)
  phen.dat3= cbind(phen.dat2[,1:5],t(slopes))
  names(phen.dat3)[6:13]=c("Estimate1","Estimate2","Std.Error1","Std.Error2","t value1","t value2","P1","P2")
  
  #restrict to significant shifts
  phen.sig= phen.dat3[which(phen.dat3$P2<0.05),]
  
  plot(abs(phen.sig$lat), phen.sig$Estimate2, ylab=ylabs[i],xlab="Absolute latitude (°)")
  abline(h=0)
  
  #curvature
  phen.sig= phen.dat3[which(phen.dat3$P1<0.05),]
  
  plot(abs(phen.sig$lat), phen.sig$Estimate1, ylab=ylabs[i],xlab="Absolute latitude (°)")
  abline(h=0)
  
  #------------------------------
  #Store plots over time
  #PLOT TRENDS OVER TIME
  
  #drop non-significant trends
  phen.dat2= phen.dat2[which(phen.dat3$"P1"<0.05),]
  
  #drop rows with all NAs
  drop.row=apply(phen.dat2, MARGIN=1, FUN=function(x)all(is.na(x[6:length(x)])) )
  if(!all(drop.row==FALSE)) phen.dat2= phen.dat2[-which(drop.row==TRUE),]
  
  phen.dat3= gather(phen.dat2, "gen", "phen",6:25)
  phen.dat3$gen= as.numeric(phen.dat3$gen)
  phen.dat4= phen.dat3[which(!is.na(phen.dat3$phen)) ,]
  
  plot1 = ggplot(phen.dat3, aes(x=gen, y=phen, group=siteID, color=abs(lat))) +geom_line() + labs(y = ylabs[i]) #+ylim(0, 5)
  #plot1 = ggplot(phen.dat3, aes(x=gen, y=phen, group=siteID, color=abs(lat) )) +geom_smooth(method=lm, se=FALSE)+ labs(y = ylab) #+ylim(0, 5)
  
  if(i==1) p1=plot1
  if(i==2) p2=plot1
  if(i==3) p3=plot1
  if(i==4) p4=plot1
  if(i==5) p5=plot1
  
} #end loop metrics

dev.off()

#-----------------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("FixedPhen_genPlots.pdf", height = 14, width = 14)
par(mfrow=c(3,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,2)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(2,1))
print(p4,vp=vplayout(2,2))
print(p5,vp=vplayout(3,1))

dev.off()

#for(i in 1:5){
#  px= paste("p",i,sep="")
#  plot(get(px))
#} #end i loop
#---------------------
#Fig 2. Slopes through time vs latitude: adult phenology, number generations; dev temperature, dev temp sd; adult temperature, adult temp sd (only phenological advancements). [dev 10th quantile]

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
#PLOT TRENDS OVER TIME

#drop non-significant trends
phen.dat2= phen.dat2[which(phen.dat3$"P"<0.05),]

#drop rows with all NAs
drop.row=apply(phen.dat2, MARGIN=1, FUN=function(x)all(is.na(x[6:length(x)])) )
if(!all(drop.row==FALSE)) phen.dat2= phen.dat2[-which(drop.row==TRUE),]

phen.dat3= gather(phen.dat2, "year", "phen",6:106)
phen.dat4= phen.dat3[which(!is.na(phen.dat3$phen)) ,]

#plot1 = ggplot(phen.dat3, aes(x=year, y=phen, group=siteID, color=lat)) +geom_line() #+ylim(0, 5)
plot1 = ggplot(phen.dat3, aes(x=year, y=phen, group=siteID, color=abs(lat) )) +geom_smooth(method=lm, se=FALSE)+ labs(y = ylabs[i]) #+ylim(0, 5)

if(i==1) p1=plot1
if(i==2) p2=plot1
if(i==3) p3=plot1
if(i==4) p4=plot1
if(i==5) p5=plot1
if(i==6) p6=plot1

} #end loop metrics

dev.off()
#plot(density(phen.dat3$Estimate))

#-----------------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("Phen_yearPlots.pdf", height = 14, width = 14)
par(mfrow=c(3,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,2)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(2,1))
print(p4,vp=vplayout(2,2))
print(p5,vp=vplayout(3,1))
print(p6,vp=vplayout(3,2))

dev.off()

#==================================================================
#Fig X. Plot shift
#x: Tpupal anomaly at fixed
#y: T pupal shift due to phenology

ylabs= c("Adult temperature shift (°C)", "Developmental temperature shift (°C)")

setwd(paste(fdir,"figures\\",sep="") )
pdf("Phen_anom.pdf", height = 5, width = 10)
par(mfrow=c(1,2), cex=1.2, mar=c(3, 3, 0.5, 0.5), oma=c(0,0,0,0), lwd=1, bty="o", tck=0.02, mgp=c(1, 0, 0))

for(i in 1:2){
    
    #Generation temps across years for ? generation
    if(i==1) {phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tmean"]))
    phen.normal= rowMeans(phen.dat[,,1,"Tmean.fixed"], na.rm=TRUE)
    phen.fix= phen.dat[,,1,"Tmean.fixed"] } 
   
    #Adult temps across years for ? generation
    if(i==2) {phen.dat2= as.data.frame( cbind(1:nrow(phen.dat), dddat[,c("Species","Order","lon","lat")],phen.dat[,,1,"Tmean.e"]))
    phen.normal= rowMeans(phen.dat[,,1,"Tmean.e.fixed"], na.rm=TRUE)
    phen.normal[is.nan(phen.normal)] = NA
    phen.fix= phen.dat[,,1,"Tmean.e.fixed"] } 
    
  colnames(phen.dat2)[1]="siteID"
  
  slopes1= matrix(NA, nrow=nrow(dddat), ncol=4)
  colnames(slopes1)=c("Estimate","Std.Error","t value","P")
  
    for(stat.k in 1:nrow(dddat) ){  #1:nrow(sites)
      
      #metrics
      Tanom= unlist(phen.dat2[stat.k,6:106] -phen.normal[stat.k])
      Tshift= unlist(phen.dat2[stat.k,6:106] -phen.fix[stat.k,])
      
      #Tdat=na.omit(cbind(Tanom, Tshift))
      #plot(Tdat[,1], Tdat[,2] )
      
      #Calculate slope
      slopes1[stat.k,]= calcm(Tanom, Tshift)
    } #end loop stats
    
  #plot
  phen.dat3= cbind(phen.dat2[,1:5],slopes1)
  names(phen.dat3)[6:9]=c("Estimate","Std.Error","t value","P")
  
  #restrict to significant shifts
  phen.sig= phen.dat3[which(phen.dat3$P<0.05),]
  
  plot(abs(phen.sig$lat), phen.sig$Estimate, ylab=ylabs[i],xlab="Absolute latitude (°)")
  abline(h=0)
  
} #end i loop 

dev.off()

#==================================================================
#ANALYZE DEVELOPMENTAL CONSTRAINTS

#source functions
setwd(paste(fdir,"out\\",sep="") )
dat= read.csv("dddat.csv")

p<- ggplot(data=dat, aes(x=abs(lat), y = BDT.C, shape=as.factor(pest), color=as.factor(pupal), size= as.factor(aquatic))) + scale_shape_manual(values = c(1,19)) + scale_size_manual(values = c(4,8)) 
# +xlab("Physiological temperature limit (°C)")+ylab("Cold range boundary temperature (°C)") 
p + geom_point()

p<- ggplot(data=dat, aes(x=abs(lat), y = EADDC, shape=as.factor(pest), color=as.factor(pupal), size= as.factor(aquatic))) + scale_shape_manual(values = c(1,19)) + scale_size_manual(values = c(4,8)) +ylim(0,1600)
p + geom_point()

p<- ggplot(data=dat, aes(x=BDT.C, y = EADDC, shape=as.factor(pupal), color=abs(lat), size= as.factor(pest))) + scale_shape_manual(values = c(1,19)) + scale_size_manual(values = c(4,8)) +ylim(0,1600) +xlim(-5,24)
p + geom_point()

p<- ggplot(data=dat, aes(x=BDT.C, y = EADDC, shape=as.factor(pupal), color=abs(lat) )) + scale_shape_manual(values = c(1,19)) +ylim(0,1600) +xlim(-5,24)
p + geom_point()

#---------------------
#plots by order

#restrict to orders with data
dat2=dat[dat$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]
#Neuroptera   Thysanoptera

setwd(paste(fdir,"figures\\",sep="") )
pdf("LDT_DDbyOrder.pdf", height = 4, width = 12)

p<- ggplot(data=dat2, aes(x=BDT.C, y = EADDC, shape=as.factor(pupal), color=abs(lat) ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(0,1600) +xlim(-5,24)
print(p + geom_point() )

dev.off()

#by latitude

pdf("DDD_ByLat.pdf", height = 4, width = 12)
p<- ggplot(data=dat2, aes(x=abs(lat), y = EADDC, shape=as.factor(pupal) ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(0,1600)
p + geom_point()
dev.off()

pdf("LDT_ByLat.pdf", height = 4, width = 12)
p<- ggplot(data=dat2, aes(x=abs(lat), y = BDT.C, shape=as.factor(pupal) ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(-5,24)
p + geom_point()
dev.off()

#number generations
ngens.ave= rowMeans(ngens, na.rm=TRUE)
#phen.dat2= as.data.frame( cbind(dddat[,c(1:9,51:52)],ngens.ave))
phen.dat2= as.data.frame( cbind(dddat[,c(1:6,7,9,51:52)],ngens.ave))
phen.dat2= phen.dat2[phen.dat2$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]

pdf("Ngens_ByLat.pdf", height = 4, width = 12)
s1= ggplot()+geom_point(data = phen.dat2, aes(x = abs(lat), y = ngens.ave))+facet_grid(.~Order)+theme_bw()+geom_smooth(data = phen.dat2, formula=y~x, aes(x = lat, y = ngens.ave), method=loess, se=TRUE) + coord_flip() +xlim(0,60)+ylim(0,20)
s1
dev.off()
#----------------------------
#model

mod1= lm(dat$EADDC~ poly(dat$BDT.C) *dat$pupal * abs(dat$lat) *dat$Order)
mod1= lm(dat$EADDC~ dat$Order * abs(dat$lat))

dat1= dat[!is.na(dat$lat),]
mod1= lm(dat$BDT.C~ dat$pupal + abs(dat$lat) + abs(dat$lat)^2)
mod1= lm(dat$BDT.C~ dat$Order * abs(dat$lat) )

#===============================
#Plot fitness surface

ngens[is.nan(ngens)] = NA 

phen.dat2= as.data.frame( cbind(1:nrow(dddat), dddat[,c("Species","Order","lon","lat","BDT.C")],ngens))
phen.dat3= as.data.frame( cbind(1:nrow(dddat), dddat[,c("Species","Order","lon","lat","BDT.C")],rowMeans(ngens, na.rm=T) ))
colnames(phen.dat3)[7]= "Ngen"

#restrict to orders with data
phen.dat3=phen.dat3[phen.dat3$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]

#ave gens across years
setwd(paste(fdir,"figures\\",sep="") )
pdf("Ngen_byLDT.pdf", height = 4, width = 12)

p<- ggplot(data=phen.dat3, aes(x=BDT.C, y = Ngen, color=abs(lat) ))+facet_grid(.~Order)  +ylim(0,20) +xlim(-5,25) +geom_smooth(data = phen.dat3, formula=y~x, aes(x=BDT.C, y = Ngen), method=loess, se=TRUE)
p + geom_point()
dev.off()

#plot Ngen by latitude
setwd(paste(fdir,"figures\\",sep="") )
pdf("Ngen_byLat.pdf", height = 4, width = 12)
p<- ggplot(data=phen.dat3, aes(x=abs(lat), y = Ngen, color=BDT.C ))+facet_grid(.~Order)  +ylim(0,20) +xlim(0,60) +geom_smooth(data = phen.dat3, formula=y~x, aes(x=abs(lat), y = Ngen), method=loess, se=TRUE)
p + geom_point()
dev.off()

#subset to Ngen data
phen.dat3= phen.dat3[which(!is.nan(phen.dat3$Ngen) ),]

#aggregate into latitudinal bins
phen.dat3$latbin= cut(phen.dat3$lat,breaks= c(-45,-20, 0,10, 20,30,40,50,61) ) #c(-45,0,20,40,61)

phen.w= phen.dat3 %>%
  group_by( latbin ) %>% 
  do(mod = lm(Ngen ~ BDT.C,na.action=na.omit,  data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2]) %>%
  select(-mod)


