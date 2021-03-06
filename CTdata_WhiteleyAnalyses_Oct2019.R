library(ggplot2)
library(reshape2)
library(nlme)
library(viridis)

#-----------------------
wd=getwd()

#DEVELOPMENT DATA
setwd("./COIS_data/")
dddat= read.csv("SeasonalityDatabase_MASTER.csv")
##Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]
#remove phys outliers
dddat$omit=NA
dddat[which(dddat$BDT.C< (-7)),"omit"]="y" #drops 3
dddat[which(dddat$EADDC>2000),"omit"]="y"  #drops 9

dddat$DE= as.numeric(as.character(dddat$DE))

#summarize
dddat$index= factor(dddat$index)

hemi= dddat[dddat$pupal==0,]
#find species with all data
hemi= hemi[which(!is.na(hemi$ET)&!is.na(hemi$LT)&!is.na(hemi$DE)&!is.na(hemi$DLT) ),]
length(unique(hemi$index))
length(unique(hemi$Species))
#hemi: 101 populations, 80 species

holo= dddat[dddat$pupal==1,]
holo$countT0= rowSums(cbind(!is.na(holo$ET),!is.na(holo$LT),!is.na(holo$TP)))
holo$countG= rowSums(cbind(!is.na(holo$DE),!is.na(holo$DLT),!is.na(holo$DP)))
holo= holo[which(holo$countT0>=2 & holo$countG>=2 ),]
length(unique(holo$index))
length(unique(holo$Species))
#holo: 406 pop, 317 species

#=======================
#SUMMARIZE OTHER DATA
#CTmax and min
#CHECK OUT HOFFMANDATA FOR MATCHES; FEW
dat1= read.csv("Hoffmannetal2012_1.csv")
dat2= read.csv("Hoffmannetal2012_2.csv")

#dat1, CTs, only 12 larvae
dat1.l= dat1[dat1$Life.stage=="Larvae",]
dat1.a= dat1[dat1$Life.stage=="Adults",]

matched= match(dat1.l$Species.name, dat1.a$Species.name)
unique(dat1$Species.name)

#dat2, LTs, 316 adults, 8 eggs, 68 larvae, 3 nymphs, 21 pupae

dat2.a= dat2[dat2$Life.stage=="Adults",]
dat2.e= dat2[dat2$Life.stage=="Eggs",]
dat2.l= dat2[dat2$Life.stage=="Larvae",]
dat2.n= dat2[dat2$Life.stage=="Nymphs",]
dat2.p= dat2[dat2$Life.stage=="Pupae",]

#cases where LT is available for adult plut one juvenile stage
matched= c(match(dat2.e$Species.name, dat2.a$Species.name),match(dat2.l$Species.name, dat2.a$Species.name),match(dat2.n$Species.name, dat2.a$Species.name),match(dat2.p$Species.name, dat2.a$Species.name))
matched= na.omit(unique(matched))

#----
#Temperature size rule
#tsr= read.csv("tsr_KlokHarrison.csv")
#matched= na.omit(match(unique(tsr$Species), dddat$Species))
#dddat[matched,]

#-------
#Dell
#setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
#dat.full<-read.csv('Delletal2013.csv')

#dat.agg= aggregate(dat.full, by=list(dat.full$DataSeriesID, dat.full$ResStage), FUN=mean)
#names(dat.agg)[1:2]=c("ID","stage")
#dat.agg=dat.agg[,1:2]
#dat.agg[dat.agg$ID==37409 ,]

library(plyr)
#dat.agg1=count(dat.agg, 'ID')
#dat.agg1[dat.agg1$freq>1,]

#=======================

##ACROSS STAGES
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

#subset of orders
dddat.o= subset(dddat, dddat$Order %in% c("Coleoptera", "Lepidoptera","Diptera","Hymenoptera") )

#PLOTTING BY STAGE

dddat.o$ET= as.numeric(as.character(dddat.o$ET))
#find species with all data
dddat1= dddat.o[which(!is.na(dddat.o$ET)&!is.na(dddat.o$LT)&!is.na(dddat.o$TP)),]
#just coleoptera and lepidoptera
#dddat1= dddat1[which(dddat1$Order %in% c("Coleoptera","Lepidoptera")),]
dddat1$ind=1:nrow(dddat1)
#mean across index
dddat1= aggregate(dddat1, list(dddat1$index,dddat1$Species,dddat1$lat,dddat1$Order), FUN="mean", na.rm=TRUE)
names(dddat1)[1:4]= c("index","Species","lat","Order")

#melt
dat2= melt(dddat1, id.vars=c("index","Species","lat","Order") , measure.vars=c("ET","LT","TP"))

#FIGURE 1
dat2$stage<- "egg"
dat2$stage[dat2$variable=="LT"]<- "larvae"
dat2$stage[dat2$variable=="TP"]<- "pupae"
colnames(dat2)[3]="latitude"
colnames(dat2)[6]="T0"

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Whiteley2019/figures/")
pdf("Fig1_To_by_stage.pdf",height = 6, width = 6)
ggplot(data=dat2, aes(x=stage, y=T0, group=index,color=abs(latitude)))+geom_line(lwd=0.5)+facet_wrap(~Order)+ylim(0,22) +
  ylab("lower developmental temperature, T0 (°C)")+ theme_bw()+ theme(legend.position="bottom")+
  scale_colour_gradientn(colours = rev(viridis(20)), name="absolute latitude (°)" )
dev.off()

#-----
#STATS
#leps
dat3= dat2[which(dat2$Order=="Lepidoptera"),]
mod1= lme(T0~stage+abs(latitude), random=~1|index, data= dat3)
#T0: larval lower, pupal (NS) higher than egg
#T0 decreases with latitude

#coleoptera
dat3= dat2[which(dat2$Order=="Coleoptera"),]
mod1= lme(T0~stage+abs(latitude), random=~1|index, data= dat3)

#T0
#T0 decreases with latitude

#by species
mod1= lme(T0~variable+abs(latitude), random=~1|Species, data= dat3)
#standard deviations
#COLEOPTERA: 2.02 intercept, 

#across orders
mod1= lme(T0~variable+abs(latitude)*Order, random=~1|Species, data= dat2)
#model selection without random than random effect

#-----------------
#PARTITION VARIATION
#index is population 
mod1= lme(T0 ~ Order*stage, random =~stage|index, data=dat2)

#stage specific variation vs latitude
dddat1$T0mean= rowMeans(dddat1[,c("ET","LT","TP")])
dddat1$T0sd= abs(dddat1$ET-dddat1$T0mean)+abs(dddat1$LT-dddat1$T0mean)+abs(dddat1$TP-dddat1$T0mean)
dddat2= dddat1[,c("Order","lat","T0sd","T0mean")]
ggplot(data=dddat2, aes(x=abs(lat), y=T0sd))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Order)+ylim(-5,25)


#------------------------
#DEVELOPMENT TIME
#plot out develop time in different stages across temperature
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

temps=seq(10,40,0.1)

dt= function(temp, T0, G){
dev=NA
if(temp>T0)  dev= G/(temp-T0)
return(dev)
}

#====================
#FIGURE 2

dat.sub=dddat.o[dddat.o$Order=="Lepidoptera",]
dat.sub$z= dat.sub$LT - (dat.sub$ET+dat.sub$TP)/2

range(dat.sub$z, na.rm=T)
cols= heat.colors(15)

#compare T0s
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]
dat.sub1$label= paste(dat.sub1$Species, " \n eggs: T0=", round(dat.sub1$ET,1), ", G=", round(dat.sub1$DE,1), ", \n larvae: T0=", round(dat.sub1$LT,1), ", G=", round(dat.sub1$DLT,1), ", \n pupae: T0=", round(dat.sub1$TP,1), ", G=", round(dat.sub1$DP,1), sep="") #, ", z=", round(dat.sub1$z,1), ", ind=", dat.sub1$index

#find good examples
dat.sub1= dat.sub1[order(dat.sub1$z),]

dat.sub1= dat.sub1[which(dat.sub1$index %in% c(805,896,946)),]
dat.sub1= dat.sub1[c(1,2,4),]

#dat.sub2= dat.sub1[,c("Species","ET","DE","LT","DLT","TP","DP")]

for(k in 1:nrow(dat.sub1) ){  #nrow(dat.sub1)
  dat.sub.sel= dat.sub1[k,]
  
  dt.e1=unlist(lapply(temps,FUN=dt, T0=dat.sub.sel$ET, G=dat.sub.sel$DE))
  dt.e= as.data.frame(cbind(temps, dt.e1))
  dt.e$species= dat.sub.sel$Species
  dt.e$index= dat.sub.sel$index
  dt.e$stage= "larvae"
  colnames(dt.e)[2]<-"devtime"
  
  dt.l1=unlist(lapply(temps,FUN=dt, T0=dat.sub.sel$LT, G=dat.sub.sel$DLT))
  dt.l= as.data.frame(cbind(temps, dt.l1+dt.e1))
  dt.l$species= dat.sub.sel$Species
  dt.l$index= dat.sub.sel$index
  dt.l$stage= "pupae"
  colnames(dt.l)[2]<-"devtime"
  dt.l$devtime[dt.l$devtime>100]=100
  
  dt.p1=unlist(lapply(temps,FUN=dt, T0=dat.sub.sel$TP, G=dat.sub.sel$DP))
  dt.p= as.data.frame(cbind(temps, dt.p1+dt.l1+dt.e1))
  dt.p$species= dat.sub.sel$Species
  dt.p$index= dat.sub.sel$index
  dt.p$stage= "adult"
  colnames(dt.p)[2]<-"devtime"
  dt.p$devtime[dt.p$devtime>100]=100
  
  if(k==1) {dt.dat= cbind(rbind(dt.e, dt.l, dt.p),dat.sub.sel[,"z"]) }
  if(k>1) dt.dat= rbind(dt.dat, cbind(rbind(dt.e, dt.l, dt.p),dat.sub.sel[,"z"]) )

  }
dt.dat$stage= factor(dt.dat$stage, levels=c("larvae","pupae","adult") )

#make labels
dt.dat$lab=NA
dt.dat$lab= dat.sub1$label[match(dt.dat$index, dat.sub1$index)]
dt.dat$Species= dat.sub1$Species[match(dt.dat$index, dat.sub1$index)]
colnames(dt.dat)[6]<-"z"

#gather total dt
library(tidyr)
dt2= spread(dt.dat, stage, devtime)

#fill in for Chilo
#dt2[which(dt2$Species=="Chilo auricilius" & dt2$temps<22 & is.na(dt2$pupae) ),"pupae"]<-100
#dt2[which(dt2$Species=="Chilo auricilius" & dt2$temps<22 & is.na(dt2$adult) ),"adult"]<-100

#plot
plot1= ggplot(data=dt2, aes(x=temps))+facet_wrap(~Species) +
  geom_ribbon(aes(ymin = 0, ymax = larvae), fill = "cadetblue1") +
  geom_ribbon(aes(ymin = larvae, ymax = pupae ), fill = "cadetblue3") +
  geom_ribbon(aes(ymin = larvae, ymax = 100 ), fill = "cadetblue3") +
  geom_ribbon(aes(ymin = pupae, ymax = adult ), fill = "cadetblue4") +
  geom_ribbon(aes(ymin = adult, ymax = 100 ), fill = "white") +
  xlim(15,25)+ ylim(0,100)+
  ylab("Development time (days)")+xlab("temperature (°C)")+ theme_bw()+ theme(legend.position="bottom", strip.text = element_text(face = "italic")) 

#------
#PLOT G BY T
library(directlabels)
library(cowplot)

dat.sub=dddat.o[dddat.o$Order=="Lepidoptera",]
dat.sub$z= dat.sub$LT - (dat.sub$ET+dat.sub$TP)/2

range(dat.sub$z, na.rm=T)
cols= heat.colors(15)

#compare T0s
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]

temps=seq(10,25,0.3)

for(k in 1:nrow(dat.sub1)){
  dat.sub= dat.sub1[k,]
  dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
  dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
  dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
  z= dat.sub$LT - (dat.sub$ET+dat.sub$TP)/2
  y= dt.l/(dt.e+dt.l+dt.p)
  
  dat.k=as.data.frame(cbind(k,temps,y,z))
  dat.k$Species= dat.sub$Species
  dat.k$index= dat.sub$index
  if(k==1) dat.kall=dat.k
  if(k>1) dat.kall= rbind(dat.kall, dat.k)
}
 
  #make species labels
  dat.kall$sp.lab= dat.kall$Species
  dat.kall$sp.lab[!dat.kall$index %in% c(805,896,946)]<- NA
  
#ggplot version
  plot2l= ggplot(data=dat.kall, aes(x=temps, y=y, col=z, group=k))+geom_line(lwd=0.5)+
  ylab("Proportion development time")+xlab("temperature (°C)")+    
    geom_dl(data=dat.kall, aes(label = sp.lab), method = list(dl.combine("first.points"), hjust=0.2)) +
    scale_colour_gradientn(colours = (viridis(20)) )

  setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Whiteley2019/figures/")
pdf("Fig2_PropLarvalDevTime.pdf",height = 6, width = 6)  
  plot_grid(plot1, plot2l, ncol=1,
          labels = 'AUTO')
dev.off()
  
#==============================
#LOOK FOR FECUNDITY DATA

dat.sub=dddat[dddat$Species=="Plutella xylostella" | dddat$Species=="Platyptilia carduidactyla",]
#write out
#write.csv(dat.sub, "SeasonalityDatabase_sub.csv" )
#add development data for Iran population
#dddat.add= read.csv("SeasonalityDatabase_MASTER_add2019.csv")
#read back in
setwd(wd)
setwd("./COIS_data/")
dat.pp= read.csv("SeasonalityDatabase_sub.csv")

#Look for papers with fecundity in title
dat.sub= dddat.o
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]
fec=grep("fecundity",dat.sub1$Title)
dat.sub1[fec,c("Author","Year","index","Species")]
dat.sub1[fec[21],]

#load fecundity data
fec= read.csv('FecundTemp.csv')

fec2= melt(fec, id.vars=c("Species","Temp","index") , measure.vars=c("S_adult","Long","F"))
#Multiply S and L by 100 for plotting
fec2[which(fec2$variable=="S_adult"),"value"]=fec2[which(fec2$variable=="S_adult"),"value"]*100
fec2$index= as.factor(fec2$index)

#make component label
fec2$component<- "survival*100"
fec2$component[fec2$variable=="Long"]<- "longevity"
fec2$component[fec2$variable=="F"]<- "fecundity"

#FECUNDITY PLOT
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Whiteley2019/figures/")
pdf("PP_FecundPlot.pdf",height = 6, width = 8)
ggplot(data=fec2[fec2$Species %in% c("Plutella xylostella","Platyptilia carduidactyla"),], aes(x=Temp, y=value, color=component, lty=index))+facet_wrap(~Species)+geom_line(lwd=0.5)+
  guides(lty = FALSE)
dev.off()

#DEVELOPMENT PLOT
dat.sub1= dat.pp[which(!is.na(dat.pp$ET)&!is.na(dat.pp$LT)&!is.na(dat.pp$TP)&!is.na(dat.pp$DE)&!is.na(dat.pp$DLT)&!is.na(dat.pp$DP)),]
#mean across index
dat.sub1$lab= paste(dat.sub1$Species, round(dat.sub1$lat,2), sep=" ")
dat.sub1= aggregate(dat.sub1, list(dat.sub1$lab,dat.sub1$Species,dat.sub1$index,dat.sub1$lat,dat.sub1$Order), FUN="mean")
names(dat.sub1)[1:5]= c("lab","Species","index","lat","Order")

#PLOT DEV TIME
temps=10:40

for(k in 1:nrow(dat.sub1)){
  dat.sub= dat.sub1[k,]
  
  dt.e1=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
  dt.e= as.data.frame(cbind(temps, dt.e1))
  dt.e$lab= dat.sub$lab
  dt.e$index= dat.sub$index
  dt.e$species= dat.sub$Species
  dt.e$stage= "larvae"
  colnames(dt.e)[2]<-"devtime"
  
  dt.l1=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
  dt.l= as.data.frame(cbind(temps, dt.l1+dt.e1))
  dt.l$lab= dat.sub$lab
  dt.l$index= dat.sub$index
  dt.l$species= dat.sub$Species
  dt.l$stage= "pupae"
  colnames(dt.l)[2]<-"devtime"
  
  dt.p1=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
  dt.p= as.data.frame(cbind(temps, dt.p1+dt.l1+dt.e1))
  dt.p$lab= dat.sub$lab
  dt.p$index= dat.sub$index
  dt.p$species= dat.sub$Species
  dt.p$stage= "adult"
  colnames(dt.p)[2]<-"devtime"
  
  if(k==1) dt.dat= rbind(dt.e, dt.l, dt.p)
  if(k>1) dt.dat= rbind(dt.dat, rbind(dt.e, dt.l, dt.p) )
}
dt.dat$stage= factor(dt.dat$stage, levels=c("larvae","pupae","adult") )

#make labels
dat.sub1$label= paste("eggs: T0=", dat.sub1$ET, ", G=", dat.sub1$DE, ", larvae: T0=", dat.sub1$LT, ", G=", dat.sub1$DLT, ", pupae: T0=", dat.sub1$TP, ", G=", dat.sub1$DP, sep="")
#dt.dat$index= factor(dt.dat$index)

#plot
pdf("PP_DevTime.pdf",height = 8, width = 8)
ggplot(data=dt.dat, aes(x=temps, y=devtime, lty=stage))+geom_line(lwd=0.5)+facet_wrap(~lab) + xlim(15,35)+ylim(0,60)+
  ylab("Development time(days)")+ theme(legend.position="bottom")
#p<- p+  geom_text(data = dat.sub1,label=dat.sub1$label)
dev.off()

#Continue component plot
fec3= fec2[(fec2$Species %in% c("Plutella xylostella","Platyptilia carduidactyla")),]
comp= fec3[,c("Species", "Temp", "index", "value", "component")]

dev= dt.dat[dt.dat$stage=="adult",c("species", "temps", "index", "devtime")]
dev= dev[dev$index %in% c(1.1,924),]
names(dev)<- c("Species", "Temp", "index", "value")
dev$component= "devtime"
comp= rbind(comp,dev)
#restrict to two populations
comp=comp[comp$index %in% c(1.1,924),]

#estimate r
#cast
comp1 <- dcast(comp, Species+Temp~component, mean)
comp1$surv= comp1[,6]/100
#restrict to available fecundity
comp1 <- comp1[!is.na(comp1$fecundity),]
#estimate plutella survival ##FIX
inds= which(comp1$Species=="Platyptilia carduidactyla")
comp1[is.na(comp1$surv),"surv"]= 1-(1/comp1[inds, "longevity"])/10*comp1[inds, "devtime"]

#calculate r=ln(sf)/G
comp1$r= log(comp1$surv*comp1$fecundity)/comp1$devtime
#calculate R0
comp1$R0= comp1$surv*comp1$fecundity

#melt
comp2= melt(comp1, id.vars=c("Species","Temp") , measure.vars=c("devtime","fecundity","surv","r","R0"))
comp2= comp2[comp2$Species=="Platyptilia carduidactyla",]

# #make label
# comp2$label= "G: generation length (days)"
# comp2$label[comp2$variable=="fecundity"]="F: fecundity"
# comp2$label[comp2$variable=="surv"]="S: survival"
# comp2$label[comp2$variable=="r"]="r= ln(SF)/G"
# comp2$label[comp2$variable=="R0"]= "R0= SF"
# comp2$label= factor(comp2$label, levels=c("S: survival","F: fecundity","G: generation length (days)","r= ln(SF)/G","R0= SF"))

#FIGURE 4
#2x2

comp2$label= "G: generation length (days)"
comp2$label[comp2$variable=="fecundity" | comp2$variable=="R0"]="F: fecundity (solid), R0= SF (dashed)"
comp2$label[comp2$variable=="surv"]="S: survival"
comp2$label[comp2$variable=="r"]="r= ln(SF)/G"
comp2$label= factor(comp2$label, levels=c("S: survival","F: fecundity (solid), R0= SF (dashed)","G: generation length (days)","r= ln(SF)/G"))

comp2$metric=" "
comp2$metric[comp2$variable=="R0"]= "R0= SF"

pdf("Fig3__FitComponentPlot.pdf",height = 8, width = 6)
ggplot(data=comp2, aes(x=Temp, y=value, lty=metric))+facet_wrap(~label, scales="free_y")+geom_line(lwd=0.5)+
  ylab("fitness component value")+xlab("Temperature (C)")+theme_bw()+theme(legend.position = "none")
dev.off()

#----------------
#try to put r and R0 on same axis
#make label
comp2$label= "G: generation length (days)"
comp2$label[comp2$variable=="fecundity"]="F: fecundity"
comp2$label[comp2$variable=="surv"]="S: survival"
comp2$label[comp2$variable=="r"]="fitness"
comp2$label[comp2$variable=="R0"]= "fitness"
comp2$label= factor(comp2$label, levels=c("S: survival","F: fecundity","G: generation length (days)","fitness"))

comp2$metric=" "
comp2$metric[comp2$variable=="r"]="r= ln(SF)/G"
comp2$metric[comp2$variable=="R0"]= "R0= SF"

ggplot(data=comp2, aes(x=Temp, y=value, lty=metric))+facet_wrap(~label, scales="free_y")+geom_line(lwd=0.5)+
  ylab("fitness component value")+xlab("Temperature (C)")
#can't figure out two axes
