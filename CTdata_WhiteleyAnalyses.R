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
fdir= "/Volumes/GoogleDrive/My Drive/Seasonality/"
setwd(paste(fdir,"out/",sep="") )
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
plot(dddat$LLT50, dddat$BDT.C, type="p", ylim=c(-10,20))
points(dddat$LLT100, dddat$BDT.C, type="p", pch="*")
plot(dddat$LLT50, dddat$EADDC, type="p", ylim=c(0,1800))
points(dddat$LLT100, dddat$EADDC, type="p", pch="*")

#----------------------
##ACROSS STAGES
#BDT.C, ET, LT, 
#EADDC, DE, DLT, DP

dddat$ET[which(dddat$ET=="-")]<-NA
dddat$ET= as.numeric(as.character(dddat$ET))
plot(dddat$ET, dddat$LT, xlim=c(-5,25), ylim=c(-5,25))
abline(a=0,b=1)
mod1=lm(dddat$LT~dddat$ET)
abline(mod1, lty="dashed")

dddat$DE[which(dddat$DE=="-")]<-NA
dddat$DE= as.numeric(as.character(dddat$DE))
plot(dddat$DE, dddat$DLT, xlim=c(0,300),ylim=c(0,1000))
plot(dddat$DLT, dddat$DP,xlim=c(0,1000))







