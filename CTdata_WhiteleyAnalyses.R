library(ggplot2)
library(reshape2)
library(nlme)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Whiteley2019/data/")

dat1= read.csv("Hoffmannetal2012_1.csv")
dat2= read.csv("Hoffmannetal2012_2.csv")

#dat1, CTs, only 12 larvae
#dat2, LTs, 316 adults, 8 eggs, 68 larvae, 3 nymphs, 21 pupae

dat2.a= dat2[dat2$Life.stage=="Adults",]
dat2.e= dat2[dat2$Life.stage=="Eggs",]
dat2.l= dat2[dat2$Life.stage=="Larvae",]
dat2.n= dat2[dat2$Life.stage=="Nymphs",]
dat2.p= dat2[dat2$Life.stage=="Pupae",]

match(dat2.e$Species.name, dat2.a$Species.name)
match(dat2.l$Species.name, dat2.a$Species.name)
match(dat2.n$Species.name, dat2.a$Species.name)
match(dat2.p$Species.name, dat2.a$Species.name)

dat2.a[c(117,137,149,263,264),]

#----------------------------------
#DEVEL DATA
dddat= read.csv("SeasonalityDatabase_MASTER.csv")
##Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]
#remove phys outliers
dddat$omit=NA
dddat[which(dddat$BDT.C< (-7)),"omit"]="y" #drops 3
dddat[which(dddat$EADDC>2000),"omit"]="y"  #drops 9

#match
length(unique(na.omit(match(dddat$Species, dat1$Species.name)))) #13
length(unique(na.omit(match(dddat$Species, dat2$Species.name)))) #29

#ULT
dat.sub=dat2[which(dat2$Metric=="ULT"),c("Species.name","LT50","LT100")]
names(dat.sub)[2:3]<-c("ULT50","ULT100")
#dddat=merge(dddat, dat.sub[,c()], by.x = "Species", by.y ="Species.name", all.x = TRUE, all.y = FALSE)
match1= match(dddat$Species, dat.sub$Species.name)
dddat$ULT50=NA
dddat$ULT100=NA
dddat$ULT50[!is.na(match1)]=dat.sub$ULT50[na.omit(match1)]
dddat$ULT100[!is.na(match1)]=dat.sub$ULT100[na.omit(match1)]
#LLT
dat.sub=dat2[which(dat2$Metric=="LLT"),c("Species.name","LT50","LT100")]
names(dat.sub)[2:3]<-c("LLT50","LLT100")
match1= match(dddat$Species, dat.sub$Species.name)
dddat$LLT50=NA
dddat$LLT100=NA
dddat$LLT50[!is.na(match1)]=dat.sub$LLT50[na.omit(match1)]
dddat$LLT100[!is.na(match1)]=dat.sub$LLT100[na.omit(match1)]

dddat[which(!is.na(dddat$ULT50) | !is.na(dddat$ULT100) |  !is.na(dddat$LLT50) | !is.na(dddat$LLT100)),]
#Plot T0 by LT
plot(dddat$LLT50, dddat$BDT.C, type="p", ylim=c(-10,20))
points(dddat$LLT100, dddat$BDT.C, type="p", pch="*")
plot(dddat$LLT50, dddat$EADDC, type="p", ylim=c(0,1800))
points(dddat$LLT100, dddat$EADDC, type="p", pch="*")

#----------------------
##ACROSS STAGES
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

#T0
dddat$ET[which(dddat$ET=="-")]<-NA
dddat$ET= as.numeric(as.character(dddat$ET))
plot(dddat$ET, dddat$LT, xlim=c(-5,25), ylim=c(-5,25))
abline(a=0,b=1)
mod1=lm(dddat$LT~dddat$ET)
abline(mod1, lty="dashed")

plot(dddat$ET, dddat$TP, xlim=c(-5,25), ylim=c(-5,25))
abline(a=0,b=1)
mod1=lm(dddat$TP~dddat$ET)
abline(mod1, lty="dashed")

ggplot(data=dddat, aes(x=ET, y=LT, color=abs(lat)))+geom_point()+geom_smooth(method="lm",se=FALSE)+ylim(-5,25)+xlim(-5,25)+geom_abline(intercept=0, slope=1)
ggplot(data=dddat, aes(x=ET, y=TP, color=abs(lat)))+geom_point()+geom_smooth(method="lm",se=FALSE)+ylim(-5,25)+xlim(-5,25)+geom_abline(intercept=0, slope=1)

#G
dddat$DE[which(dddat$DE=="-")]<-NA
dddat$DE= as.numeric(as.character(dddat$DE))
plot(dddat$DE, dddat$DLT, xlim=c(0,300),ylim=c(0,1000))
plot(dddat$DLT, dddat$DP,xlim=c(0,1000))

#PLOT AS PROPORTION
plot(dddat$DE/dddat$EADDC, dddat$DLT/dddat$EADDC, xlim=c(0,1),ylim=c(0,1))
plot(dddat$DLT/dddat$EADDC, dddat$DP/dddat$EADDC,xlim=c(0,1))

#STATS
dddat.o= subset(dddat, dddat$Order %in% c("Coleoptera", "Lepidoptera","Diptera","Hymenoptera") )
dat.t0= na.omit(dddat.o[,c("Species","Order","ET","LT","TP")])
#convert to long format
dat.t0= melt(dat.t0, id = c("Species","Order"))
colnames(dat.t0)[c(3,4)]= c("stage","T0")

mod1= lme(T0 ~ Order*stage, random =~stage|Species, data=dat.t0)

#ACROSS LATITUDE
#restrict to orders with sufficient data
dddat.o= subset(dddat, dddat$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") )

#T0
ggplot(data=dddat.o, aes(x=abs(lat), y=LT, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-5,25)
ggplot(data=dddat.o, aes(x=abs(lat), y=ET, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-5,25)
ggplot(data=dddat.o, aes(x=abs(lat), y=TP, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-5,25)

#subset of orders
dddat.o= subset(dddat, dddat$Order %in% c("Coleoptera", "Lepidoptera","Diptera","Hymenoptera") )
#melt
dat2= melt(dddat.o, id.vars=c("index","Species","lat","Order") , measure.vars=c("ET","LT","TP","BDT.C"))

ggplot(data=dat2, aes(x=abs(lat), y=value))+geom_point()+facet_grid(Order~variable)+geom_smooth(method="lm",se=TRUE, alpha=0.4)+ylim(0,25)

#stage specific variation vs latitude
dddat.o$T0mean= rowMeans(dddat.o[,c("ET","LT","TP")])
dddat.o$T0sd= abs(dddat.o$ET-dddat.o$T0mean)+abs(dddat.o$LT-dddat.o$T0mean)+abs(dddat.o$TP-dddat.o$T0mean)
ggplot(data=dddat.o, aes(x=abs(lat), y=T0sd))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-5,25)

#G
ggplot(data=dddat.o, aes(x=abs(lat), y=DE, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(0,300)
ggplot(data=dddat.o, aes(x=abs(lat), y=DLT, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(0,500)
ggplot(data=dddat.o, aes(x=abs(lat), y=DP, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(0,400)

#---
#PLOT RESIDUALS ACROSS LATITUDE
dddat1= dddat.o[which(!is.na(dddat.o$ET)&!is.na(dddat.o$LT)&!is.na(dddat.o$TP)),]
#just coleoptera and lepidoptera
#dddat1= dddat1[which(dddat1$Order %in% c("Coleoptera","Lepidoptera")),]

#larval vs egg
mod1=lm(LT~ET, data=dddat1)
dddat1$resid.EL= mod1$residuals
#pupal vs larval
mod1=lm(TP~LT, data=dddat1)
dddat1$resid.LP= mod1$residuals
#pupal vs egg
mod1=lm(TP~ET, data=dddat1)
dddat1$resid.EP= mod1$residuals

#melt
dat2= melt(dddat1, id.vars=c("index","Species","lat","Order") , measure.vars=c("resid.EL","resid.LP","resid.EP"))

ggplot(data=dat2, aes(x=abs(lat), y=value))+geom_point()+facet_grid(Order~variable)+geom_smooth(method="lm",se=TRUE, alpha=0.4)+ylim(-10,10)

#----
#pupal vs larval
dat.sub2= dddat.o[which(!is.na(dddat.o$LT) & !is.na(dddat.o$TP)),]
mod1=lm(dat.sub2$TP~dat.sub2$LT)
dat.sub2$resid= mod1$residuals
ggplot(data=dat.sub2, aes(x=abs(lat), y=resid, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-10,10)
#1:1
ggplot(data=dat.sub2, aes(x=abs(lat), y=TP-LT, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-10,10)

#pupal vs egg
dat.sub2= dddat.o[which(!is.na(dddat.o$TP) & !is.na(dddat.o$ET)),]
mod1=lm(dat.sub2$TP~dat.sub2$ET)
dat.sub2$resid= mod1$residuals
ggplot(data=dat.sub2, aes(x=abs(lat), y=resid, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-10,10)
#1:1
ggplot(data=dat.sub2, aes(x=abs(lat), y=TP-ET, color=BDT.C))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-10,10)

#--------------------
#PLOTTING BY STAGE
library(reshape2)
library(ggplot2)
library(nlme)

#find species with all data
dddat1= dddat.o[which(!is.na(dddat.o$ET)&!is.na(dddat.o$LT)&!is.na(dddat.o$TP)),]
#just coleoptera and lepidoptera
#dddat1= dddat1[which(dddat1$Order %in% c("Coleoptera","Lepidoptera")),]
dddat1$ind=1:nrow(dddat1)
#mean across index
### Passing further arguments through ...
##dcast(ffm, treatment ~ ., sum)
dddat1= dcast(dddat1, index ~., fun.aggregate = mean, value.var=c("ET","LT","TP") )
dddat1= aggregate(dddat1, list(dddat1$index,dddat1$Species,dddat1$lat,dddat1$Order), FUN="mean")
names(dddat1)[1:4]= c("index","Species","lat","Order")

#melt
dat2= melt(dddat1, id.vars=c("index","Species","lat","Order") , measure.vars=c("ET","LT","TP"))

ggplot(data=dat2, aes(x=variable, y=value, group=index,color=abs(lat)))+geom_line(lwd=0.5)+facet_wrap(~Order)+ylim(0,22) #+geom_smooth(method="lm",se=FALSE, alpha=0.4)

#-----
#STATS
#leps
dat3= dat2[which(dat2$Order=="Lepidoptera"),]
mod1= lme(value~variable+abs(lat), random=~1|index, data= dat3)
#T0: larval lower, pupal (NS) higher than egg
#T0 decreases with latitude

#coleoptera
dat3= dat2[which(dat2$Order=="Coleoptera"),]
mod1= lme(value~variable*abs(lat), random=~1|index, data= dat3)
mod1= lme(value~variable, random=~1|index, data= dat3)

#T0
#T0 decreases with latitude

#by species
mod1= lme(value~variable+abs(lat), random=~1|Species, data= dat3)
#standard deviations
#COLEOPTERA: 2.02 intercept, 

#across orders
mod1= lme(value~variable*abs(lat)*Order, random=~1|Species, data= dat2)
#model selection without random than random effect

#------------------------
#plot out develop time in different stages across temperature
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

temps=15:40

dt= function(temp, T0, G){
dev=NA
if(temp>T0)  dev= G/(temp-T0)
return(dev)
}

#PLOT DT BY STAGE
dat.sub=dddat.o[dddat.o$Genus=="Liriomyza",]
#dat.sub=dddat.o[dddat.o$Species=="Pieris brassicae",]

dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]
dat.sub1= dat.sub1[c(5,8,10),]

par(mfrow=c(1,3))

for(k in 1:nrow(dat.sub1)){
dat.sub= dat.sub1[k,]
dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
#if(k==1) 
plot(temps, dt.e, type="l", ylim=c(0,30), col=heat.colors(nrow(dat.sub1))[k])
#if(k>1)points(temps, dt.e, type="l", col=heat.colors(nrow(dat.sub1))[k])
dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
points(temps, dt.e+dt.l, type="l", col=heat.colors(nrow(dat.sub1))[k])
dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
points(temps, dt.e+dt.l+dt.p, type="l", col=heat.colors(nrow(dat.sub1))[k])
}

#PLOT G BY T
dat.sub=dddat.o[dddat.o$Order=="Lepidoptera",]
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]

par(mfrow=c(1,2))

#Gen time across species
for(k in 1:nrow(dat.sub1)){
  dat.sub= dat.sub1[k,]
  dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
  dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
  dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
  if(k==1) plot(temps, dt.e+dt.l+dt.p, type="l", col=heat.colors(nrow(dat.sub1))[k])
  if(k>1) points(temps, dt.e+dt.l+dt.p, type="l", col=heat.colors(nrow(dat.sub1))[k])
}

temps=seq(10,25,0.3)
#plot larval as proportion of total
for(k in 1:nrow(dat.sub1)){
  dat.sub= dat.sub1[k,]
  dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
  dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
  dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
  if(k==1) plot(temps, dt.l/(dt.e+dt.l+dt.p), type="l", col=heat.colors(nrow(dat.sub1))[k], ylim=c(0,1))
  if(k>1) points(temps, dt.l/(dt.e+dt.l+dt.p), type="l", col=heat.colors(nrow(dat.sub1))[k])
}

#LOOK FOR FECUNDITY DATA
#dddat.o[dddat.o$Genus=="Pieris",]
dat.sub= dddat.o
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]
fec=grep("fecundity",dat.sub1$Title)
dat.sub1[fec,c("Author","Year","index","Species")]
dat.sub1[fec[21],]

#load fecundity data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ICBseasonality/data/")
fec= read.csv('FecundTemp.csv')

fec2= melt(fec, id.vars=c("Species","Temp","index") , measure.vars=c("S_adult","Long","F"))
#Multiply S and L by 100 for plotting
fec2[which(fec2$variable=="S_adult"),"value"]=fec2[which(fec2$variable=="S_adult"),"value"]*100

ggplot(data=fec2, aes(x=Temp, y=value, color=variable))+facet_wrap(~index)+geom_line(lwd=0.5)

#species with good fecundity data
dat.sub1[dat.sub1$Species=="Platyptilia carduidactyla",]
dat.sub1[dat.sub1$Species=="Plutella xylostella",]

#-----
#add development data for Iran population
dddat.add= read.csv("SeasonalityDatabase_MASTER_add2019.csv")

dat.sub= dddat.add
dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
plot(temps, dt.e, type="l", ylim=c(0,30))
dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
points(temps, dt.e+dt.l, type="l")
dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
points(temps, dt.e+dt.l+dt.p, type="l")

#=================
#CHECK GLOBTHERM FOR ONTOGENETIC DATA

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")

tol.gt= read.csv('GlobalTherm_upload_10_11_17.csv')
#ALL ADULTS







