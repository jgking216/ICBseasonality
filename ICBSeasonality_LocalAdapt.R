fdir= "/Volumes/GoogleDrive/My Drive/Seasonality/"

#LOAD libraries
library(ggplot2)
library(dplyr)
library(colorRamps)     # for matlab.like(...)
library(akima) #for interpolations
library(tidyr)
library(sp)
library(rnoaa) #http://recology.info/2015/07/weather-data-with-rnoaa/
library(zoo)

day_of_year<- function(day, format="%Y-%m-%d"){
  day=  as.POSIXlt(day, format=format)
  return(as.numeric(strftime(day, format = "%j")))
}

#LOAD DATA
setwd(paste(fdir,"out/",sep="") )

##READ BACK IN
dddat= read.csv("SeasonalityDatabase_MASTER.csv")

##Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]

#remove phys outliers
dddat$omit=NA
dddat[which(dddat$BDT.C< (-7)),"omit"]="y" #drops 3
dddat[which(dddat$EADDC>2000),"omit"]="y"  #drops 9

#drop outliers
dddat= dddat[which(is.na(dddat$omit)),]

#group by index
dddat$index= as.factor(dddat$index)

dddat= dddat %>% group_by(index) %>% summarise(Species=head(Species)[1],Order=head(Order)[1],Family=head(Family)[1],Genus=head(Genus)[1],Species.1=head(Species.1)[1],BDT.C=mean(BDT.C), UDT.C=mean(UDT.C),EADDC=mean(EADDC),EEDDC=mean(EEDDC), pest=mean(pest), aquatic=mean(aquatic),pupal=mean(pupal),Location=head(Location)[1],lon=mean(lon),lat=mean(lat), quality=head(quality)[1], parasitoid=head(parasitoid)[1])

dddat= as.data.frame(dddat)

#change trait names
names(dddat)[which(names(dddat)=="BDT.C")]<-"T0"
names(dddat)[which(names(dddat)=="EADDC")]<-"G"

#------------------------
#Count of locations by species

la.spec= dddat %>% group_by(Species) %>% summarise(Npop= length(lon)  )
la.spec= subset(la.spec, la.spec$Npop>4)

la.dat= subset(dddat, dddat$Species %in% la.spec$Species)

p<- ggplot(data=la.dat, aes(x=T0, y = G, color=Species ))+facet_grid(.~Order) +theme_bw()+ geom_point()  +geom_smooth(method=lm, se=FALSE)

#by latitude
p2<- ggplot(data=la.dat, aes(x=abs(lat), y = T0, color=Species ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +xlab("Absolute latitude (°)")+theme_bw()+ geom_point() +geom_smooth(method=lm, se=FALSE)
#+ scale_color_gradientn(colours=matlab.like(20)) #+scale_shape_discrete(guide = FALSE)

p3<- ggplot(data=la.dat, aes(x=abs(lat), y = G, color=Species ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +xlab("Absolute latitude (°)") +theme_bw() + geom_point() +geom_smooth(method=lm, se=FALSE)

#------------------------
# BY GENUS

la.gen= dddat %>% group_by(Genus) %>% summarise(Npop= length(lon)  )
la.gen= subset(la.gen, la.gen$Npop>9)

la.dat= subset(dddat, dddat$Genus %in% la.gen$Genus)

p<- ggplot(data=la.dat, aes(x=T0, y = G, color=Genus, size=abs(lat) ))+facet_grid(.~Order) +theme_bw()+ geom_point()  +geom_smooth(method=lm, se=FALSE)+scale_size(range = c(1, 4))

#by latitude
p2<- ggplot(data=la.dat, aes(x=abs(lat), y = T0, color=Genus ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +xlab("Absolute latitude (°)")+theme_bw()+ geom_point() +geom_smooth(method=lm, se=FALSE)
#+ scale_color_gradientn(colours=matlab.like(20)) #+scale_shape_discrete(guide = FALSE)

p3<- ggplot(data=la.dat, aes(x=abs(lat), y = G, color=Genus ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +xlab("Absolute latitude (°)") +theme_bw() + geom_point() +geom_smooth(method=lm, se=FALSE)

#plot out
setwd("/Volumes/GoogleDrive/My Drive/Buckley/work/ICBSeasonality/figures/LocalAdaptation/") 
pdf("GeneraToG.pdf",height = 6, width = 10)
p2
p3
dev.off()

#-----------------------
#FITNESS IMPLICATIONS

#LOAD DATA
setwd(paste(fdir,"out/",sep="") )
phen.dat= readRDS("phendat.rds")
phen.fixed= readRDS("phenfix.rds")
ngens= readRDS("ngens.rds")
dddat= readRDS("dddat_media.rds")

#drop omit data
dropi= which( dddat$omit=="y" )
phen.dat= phen.dat[-dropi,,,]
phen.fixed= phen.fixed[-dropi,,]
ngens= ngens[-dropi,]
dddat= dddat[-dropi,]

#Restrict to focal genera
inds= which(dddat$Genus %in% la.gen$Genus)
phen.dat= phen.dat[inds,,,]
phen.fixed= phen.fixed[inds,,]
ngens= ngens[inds,]
dddat= dddat[inds,]

#==================================================================
#FITNESS
#number generations
ngens.ave= rowMeans(ngens, na.rm=TRUE)
#phen.dat2= as.data.frame( cbind(dddat[,c(1:9)],ngens.ave))
phen.dat2= as.data.frame( cbind(la.dat[,c(1:6,7,9)],ngens.ave))
phen.dat2= phen.dat2[phen.dat2$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]

#Plot fitness surface
ngens[is.nan(ngens)] = NA 

phen.dat2= as.data.frame( cbind(1:nrow(la.dat), la.dat[,c("Species","Genus","Order","lon","lat","T0","G")],ngens))
phen.dat3= as.data.frame( cbind(1:nrow(la.dat), la.dat[,c("Species","Genus","Order","lon","lat","T0","G")],rowMeans(ngens, na.rm=T) ))
colnames(phen.dat3)[9]= "Ngen"

#subset to Ngen data
phen.dat3= phen.dat3[which(!is.nan(phen.dat3$Ngen) ),]

#--------------
#Fitness across latitude
f1=ggplot(data=phen.dat3, aes(x=abs(lat), y = Ngen, color=Genus ))+facet_grid(.~Order) +theme_bw() + geom_point() +geom_smooth(method=lm, se=FALSE)

#==========================================================
#PHENOLOGY

phen.dat2= as.data.frame( cbind(1:nrow(la.dat), la.dat[,c("Species","Genus","Order","lon","lat", "T0","G")],phen.dat[,,1,"phen"]))
colnames(phen.dat2)[1]="siteID"

#mean phenology 1970 to 2015
phen.dat2$phen= rowMeans(phen.dat2[,9:54], na.rm=TRUE)

#--------------
#Phenology across latitude
f2=ggplot(data=phen.dat2, aes(x=abs(lat), y = phen, color=Genus ))+facet_grid(.~Order) +theme_bw() + geom_point() +geom_smooth(method=lm, se=FALSE)

#plot out
setwd("/Volumes/GoogleDrive/My Drive/Buckley/work/ICBSeasonality/figures/LocalAdaptation/") 
pdf("Genera_ToGNgenPhenology.pdf",height = 6, width = 10)
p2
p3
f1
f2
dev.off()

#==========================================================
#VIRTUAL RECIPROCAL TRANSPLANTS

#--------------------------------
#FIND CLOSEST GHCN STATIONS

#Read GGCN database inventory http://www.ncdc.noaa.gov/oa/climate/ghcn-daily/
setwd(paste(fdir,"data/",sep=""))
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
# phen.dat= array(NA, dim=c(nrow(la.dat), length(1970:2015), 20, 10), dimnames=list(NULL,as.character(1970:2015),as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen", "Tmean.e.fixed", "Tmean.fixed")) )
# ngens= array(NA, dim=c(nrow(la.dat), length(1970:2015)), dimnames=list(NULL,as.character(1970:2015)))
# phen.fixed= array(NA, dim=c(nrow(la.dat), 20, 8), dimnames=list(NULL,as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen")) )

#----------------------------------
#ANALYSIS

genera= unique(la.dat$Genus)
#or switch to order?

j1.all= array(NA, dim=c(length(genera),length(1970:2015),100,100) )
t1.all= array(NA, dim=c(length(genera),length(1970:2015),100,100) )
ngen.all= array(NA, dim=c(length(genera),length(1970:2015),100,100) )

for(genus.k in 1:length(genera)){
  
  gen.dat= subset(la.dat, la.dat$Genus==genera[genus.k])
  stat.inds= order(gen.dat$lat)[!duplicated(sort(gen.dat$lat))]

for(stat.k in 1:length(stat.inds) ){  
    min.dist<- order(spDistsN1(stat.coords, as.numeric(gen.dat[stat.k,c("lon","lat")]), longlat = TRUE))[1:100]
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
        dat$tmin= NA
        dat$tmin[is.matched]= tmin$tmin$tmin[match1[is.matched]]
        
        #split date 
        date= as.Date(dat$date, "%Y-%m-$d")
        dat$year=as.numeric(format(date, "%Y"))
        dat$month=as.numeric(format(date, "%m"))
        dat$day=as.numeric(format(date, "%d"))
        dat$j= unlist(lapply(date, FUN="day_of_year"))
        
        #Format data
        dat$tmax= as.numeric(dat$tmax)
        dat$tmin= as.numeric(dat$tmin)
        
        dat$tmax[which(dat$tmax==-9999)]= NA
        dat$tmax= dat$tmax/10 #correct for tenths of degrees or mm
        dat$tmin[which(dat$tmin==-9999)]= NA
        dat$tmin= dat$tmin/10 #correct for tenths of degrees or mm
        
        #Catch other NA values
        dat$tmax[which(dat$tmax>200)]= NA
        dat$tmin[which(dat$tmin>200)]= NA
        dat$tmax[which(dat$tmax< -200)]= NA
        dat$tmin[which(dat$tmin< -200)]= NA
        
        ## FIND YEARS WITH NEARLY COMPLETE DATA
        dat.agg= aggregate(dat, list(dat$year),FUN=function(x)length(na.omit(x))  )  
        years= dat.agg$Group.1[which(dat.agg$tmax>300)]
        dat= dat[which(dat$year %in% years),]
        
      } #END CHECK NULL
    } #END WHILE YEARS
    
    #RECORD SITE DATA
    gen.dat$Id[stat.k]= as.character(min.site[ind])
    gen.dat$st.lat[stat.k]= stations$Lat[min.dist[ind]]
    gen.dat$st.lon[stat.k]= stations$Lon[min.dist[ind]]
    gen.dat$st.elev[stat.k]= stations$elev[min.dist[ind]]
    
    #--------------------------------------
    #INTERPOLATE MISSING DATA
    dat$tmin= na.approx(dat$tmin, maxgap=7, na.rm = FALSE)
    dat$tmax= na.approx(dat$tmax, maxgap=7, na.rm = FALSE)
    #Tmean
    dat$tmean= (dat$tmin+dat$tmax)/2
    
    #cut missing data
    dat= dat[which(!is.na(dat$tmin) & !is.na(dat$tmax)),]
    
    #---------------------------
    #TRANSPLANT
    
    dds= array(NA, dim=c( nrow(gen.dat), length(stat.inds), 20, 2) )
    
    #normalize S hemisphere to N hemisphere 
    if(gen.dat[stat.k,"lat"]<0) 
    {inds.jd= which(dat$j>181)
    inds.jj= which(dat$j<182)
    
    dat[inds.jd, "j"] = dat$j[inds.jd] -181
    dat[inds.jj, "j"] = dat$j[inds.jj] +184
    } #end fix s hemi
    
    #CALCULATE DEGREE DAYS
    for(pop.k in 1:nrow(gen.dat) ){
    
      #dds[,pop.k,stat.k]
      dat$dd= apply( dat[,c("tmin","tmax")], MARGIN=1, FUN=degree.days.mat, LDT=gen.dat$T0[pop.k] )
    
    dat.dd = dat %>% group_by(year) %>% arrange(j) %>% mutate(cs = cumsum(dd))
    
    #Estimate phenology across years based on DD accumulations
    #cumsum within groups
    #Egg to adult DD, First date beyond threshold
    
    js= matrix(NA, 20, length(1970:2015) )
    ts=  matrix(NA, 20, length(1970:2015) )
    
    for(genk in 1:20){
      
      phen= dat.dd %>%  group_by(year) %>% slice(which.max(cs> (genk* gen.dat$G[pop.k]) ))
     
      #replace eroneous values
      phen[which(phen$j==1),"j"]=NA #& phen$cs==0
      #drop years without generation
      phen= phen[!is.na(phen$j),]
      
      #use only years starting 1970
      years.ind=1970:2015
      year.loop= unique(phen$year)
      year.loop= year.loop[which(year.loop>1969 & year.loop<2016) ]
      
     if(length(year.loop)>0) for(yeark in 1: length(year.loop)){
        
        year1= year.loop[yeark]
        j= as.numeric(phen[phen$year==year1,"j"])
        
        dat.yr= dat.dd[dat.dd$year==year1,]
        
        #FIND GEN TIME
        j.gs= as.numeric(ifelse(genk==1, dat.yr[which.min(dat.yr$dd>0),"j"], js[genk-1,yeark] ))
        
        #TEMPS: ACROSS GENERATIONS
        ts[genk,yeark]= mean(as.numeric(unlist(dat.yr[dat.yr$j %in% j.gs:j,"tmean"])))
      
        js[genk,yeark]= j
        
      } #end year loop
    } #end GEN LOOP
    
    #Number generations by year
    
    ##DROP?
    #years.match= match(year.loop, years.ind)
    
    all.na= apply(js, MARGIN=1, function(x)all(is.na(x)) )
    
    ngen.all[genus.k,,stat.k,pop.k]= apply(js, MARGIN=2, FUN=function(x)length(which(x>1)) )
    #correct for years without data
    ngen.all[genus.k,which(all.na==TRUE),stat.k,pop.k]=NA
    
    j1.all[genus.k,,stat.k,pop.k]= js[1,]
    t1.all[genus.k,,stat.k,pop.k]= ts[1,]
    
    } #end population loop
    
  } #loop station
} #loop genus

##SAVE OUTPUT
setwd(paste(fdir,"out/" ,sep=""))
saveRDS(phen.dat, "phendat_gen.rds")
saveRDS(phen.fixed, "phenfix_gen.rds")
saveRDS(ngens, "ngens_gen.rds")
saveRDS(dddat, "dddat_media_gen.rds")

##READ BACK IN
# setwd(paste(fdir,"out/",sep="") )
# phen.dat= readRDS("phendat_gen.rds")
# phen.fixed= readRDS("phenfix_gen.rds")
# ngens= readRDS("ngens_gen.rds")
# dddat= readRDS("dddat_media_gen.rds")

#drop omit data
dropi= which( dddat$omit=="y" )
phen.dat= phen.dat[-dropi,,,]
phen.fixed= phen.fixed[-dropi,,]
ngens= ngens[-dropi,]
dddat= dddat[-dropi,]

#============================================================
#PLOT RECIPROCAL TRANPLANT RESULTS

require(reshape2); require(ggplot2)
library(cowplot)

#Does local adaptation enable more generations, earlier phenology, diff temperature?

#surface plots
#for each genus
#across years, for each station, plot data for all populations

setwd("/Volumes/GoogleDrive/My Drive/Buckley/work/ICBSeasonality/figures/LocalAdaptation/") 
pdf("RecipTran.pdf",height = 5, width = 10)

for(genus.k in 1:length(genera)){
  
  gen.dat= subset(la.dat, la.dat$Genus==genera[genus.k])
  lat.ord= order(abs(gen.dat$lat))

  ngen.g= ngen.all[genus.k,,lat.ord,lat.ord]
  j1.g= j1.all[genus.k,,lat.ord,lat.ord]
  t1.g= t1.all[genus.k,,lat.ord,lat.ord]

  #melt
  ngen.m= melt(ngen.g, varnames=c("year","station","pop") )
  j1.m= melt(j1.g, varnames=c("year","station","pop") )
  t1.m= melt(t1.g, varnames=c("year","station","pop") )
  
  #average across years
  ngen.y= aggregate(ngen.m, by= list(ngen.m$station, ngen.m$pop), FUN=mean, na.rm=TRUE )
  j1.y= aggregate(j1.m, by= list(j1.m$station, j1.m$pop), FUN=mean, na.rm=TRUE )
  t1.y= aggregate(t1.m, by= list(t1.m$station, t1.m$pop), FUN=mean, na.rm=TRUE )
  
  #surface plots
  
  #ngen       
  n.plot= ggplot(ngen.y) + 
    aes(x = station, y = pop, z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="number\ngenerations") #, breaks=c(1,2, 5,10,20)) 
#  +theme_bw(base_size=18)+xlab("T0 (°C)")+ylab("G")+ggtitle(lats[i])+ theme(legend.position="right")+ coord_fixed(ratio = 0.01) 
  
  #j1
  j.plot= ggplot(j1.y) + 
    aes(x = station, y = pop, z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="j") 
  
  #t1
  t.plot= ggplot(t1.y) + 
    aes(x = station, y = pop, z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="temperature") 

  plot_grid(n.plot, j.plot, t.plot, labels = c("A", "B","C"), ncol=3)
  
  } #end loop genera

dev.off()
