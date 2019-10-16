library(ggplot2)
library(reshape2)
library(nlme)

#DEVELOPMENT DATA
fdir= "/Volumes/GoogleDrive/My Drive/Seasonality/"
setwd(paste(fdir,"out/",sep="") )
dddat= read.csv("SeasonalityDatabase_MASTER.csv")
##Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]
#remove phys outliers
dddat$omit=NA
dddat[which(dddat$BDT.C< (-7)),"omit"]="y" #drops 3
dddat[which(dddat$EADDC>2000),"omit"]="y"  #drops 9

#----------------------
##ACROSS STAGES
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

#subset of orders
dddat.o= subset(dddat, dddat$Order %in% c("Coleoptera", "Lepidoptera","Diptera","Hymenoptera") )

#PLOTTING BY STAGE

#find species with all data
dddat1= dddat.o[which(!is.na(dddat.o$ET)&!is.na(dddat.o$LT)&!is.na(dddat.o$TP)),]
#just coleoptera and lepidoptera
#dddat1= dddat1[which(dddat1$Order %in% c("Coleoptera","Lepidoptera")),]
dddat1$ind=1:nrow(dddat1)
#mean across index
dddat1= aggregate(dddat1, list(dddat1$index,dddat1$Species,dddat1$lat,dddat1$Order), FUN="mean")
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
  ylab("lower developmental temperature, T0 (C)")+ theme(legend.position="bottom")#+
  scale_colour_gradientn(colours = heat.colors(10))
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

#------------------------
#DEVELOPMENT TIME
#plot out develop time in different stages across temperature
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

temps=10:40

dt= function(temp, T0, G){
dev=NA
if(temp>T0)  dev= G/(temp-T0)
return(dev)
}

#FIGURE 2
#PLOT DT BY STAGE
dat.sub=dddat.o[dddat.o$Genus=="Liriomyza",]
#dat.sub=dddat.o[dddat.o$Species=="Pieris brassicae",]

dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]
dat.sub1= dat.sub1[c(5,8,10),]

for(k in 1:nrow(dat.sub1)){
dat.sub= dat.sub1[k,]

dt.e1=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
dt.e= as.data.frame(cbind(temps, dt.e1))
dt.e$species= dat.sub$Species
dt.e$index= dat.sub$index
dt.e$stage= "larvae"
colnames(dt.e)[2]<-"devtime"

dt.l1=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
dt.l= as.data.frame(cbind(temps, dt.l1+dt.e1))
dt.l$species= dat.sub$Species
dt.l$index= dat.sub$index
dt.l$stage= "pupae"
colnames(dt.l)[2]<-"devtime"

dt.p1=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
dt.p= as.data.frame(cbind(temps, dt.p1+dt.l1+dt.e1))
dt.p$species= dat.sub$Species
dt.p$index= dat.sub$index
dt.p$stage= "adult"
colnames(dt.p)[2]<-"devtime"

if(k==1) dt.dat= rbind(dt.e, dt.l, dt.p)
if(k>1) dt.dat= rbind(dt.dat, rbind(dt.e, dt.l, dt.p) )
}
dt.dat$stage= factor(dt.dat$stage, levels=c("larvae","pupae","adult") )

#make labels
dat.sub1$label= paste("eggs: T0=", round(dat.sub1$ET,1), ", G=", round(dat.sub1$DE,1), ", larvae: T0=", round(dat.sub1$LT,1), ", G=", round(dat.sub1$DLT,1), ", pupae: T0=", round(dat.sub1$TP,1), ", G=", round(dat.sub1$DP,1), sep="")
dt.dat$lab=NA
dt.dat$lab= dat.sub1$label[match(dt.dat$index, dat.sub1$index)]

#plot
pdf("Fig2_DevTime.pdf",height = 6, width = 10)
ggplot(data=dt.dat, aes(x=temps, y=devtime, lty=stage))+geom_line(lwd=0.5)+facet_wrap(~lab, labeller = label_wrap_gen()) +xlim(15,35)+ylim(0,60)+
  ylab("Development time(days)") + theme(legend.position="bottom") #, strip.background = element_blank(), strip.text.x = element_blank()
dev.off()
#---------------
#FIGURE 3
#PLOT G BY T
dat.sub=dddat.o[dddat.o$Order=="Lepidoptera",]
dat.sub1= dat.sub[which(!is.na(dat.sub$ET)&!is.na(dat.sub$LT)&!is.na(dat.sub$TP)&!is.na(dat.sub$DE)&!is.na(dat.sub$DLT)&!is.na(dat.sub$DP)),]

temps=seq(10,25,0.3)

pdf("Fig3_PropLarvalDevTime.pdf",height = 6, width = 6)

#plot larval as proportion of total
for(k in 1:nrow(dat.sub1)){
  dat.sub= dat.sub1[k,]
  dt.e=unlist(lapply(temps,FUN=dt, T0=dat.sub$ET, G=dat.sub$DE))
  dt.l=unlist(lapply(temps,FUN=dt, T0=dat.sub$LT, G=dat.sub$DLT))
  dt.p=unlist(lapply(temps,FUN=dt, T0=dat.sub$TP, G=dat.sub$DP))
  if(k==1) plot(temps, dt.l/(dt.e+dt.l+dt.p), type="l", col=heat.colors(nrow(dat.sub1))[k], ylim=c(0,1),xlim=c(15,25), ylab="Proportion development time as larvae", xlab="Temperature (C)")
  if(k>1) points(temps, dt.l/(dt.e+dt.l+dt.p), type="l", col=heat.colors(nrow(dat.sub1))[k])
}
dev.off()

#==============================
#LOOK FOR FECUNDITY DATA

dat.sub=dddat[dddat$Species=="Plutella xylostella" | dddat$Species=="Platyptilia carduidactyla",]
#write out
#write.csv(dat.sub, "SeasonalityDatabase_sub.csv" )
#add development data for Iran population
#dddat.add= read.csv("SeasonalityDatabase_MASTER_add2019.csv")
#read back in
dat.pp= read.csv("SeasonalityDatabase_sub.csv")

#Look for papers with fecundity in title
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
#melt
comp2= melt(comp1, id.vars=c("Species","Temp") , measure.vars=c("devtime","fecundity","surv","r"))

#make label
comp2$label= "G: generation length (days)"
comp2$label[comp2$variable=="fecundity"]="f: fecundity"
comp2$label[comp2$variable=="surv"]="s: survival"
comp2$label[comp2$variable=="r"]="r= ln(sf)/G"
comp2$label= factor(comp2$label, levels=c("s: survival","f: fecundity","G: generation length (days)","r= ln(sf)/G") )

#FIGURE 4
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Whiteley2019/figures/")
pdf("Fig4__FitComponentPlot.pdf",height = 8, width = 6)
ggplot(data=comp2, aes(x=Temp, y=value))+facet_grid(label~Species, scales="free_y")+geom_line(lwd=0.5)+
  ylab("fitness component value")+xlab("Temperature (C)")
dev.off()

