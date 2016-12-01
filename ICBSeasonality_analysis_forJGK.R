#SEASONALITY ANALYSIS

fdir= "C:\\Users\\Buckley\\Google Drive\\Seasonality\\"

#LOAD LIBRARIES
library(ggplot2)
library(tidyr)
library(grid)
library(plyr); library(dplyr)
library(nominatim) #for osm geocode
library(gridExtra)

#LOAD FUNCTIONS
#calculate slope
calcm= function(x, ys=ys){
  yd= which( !is.na(x) )
  if(length(yd)>0) mod1= lm( as.numeric(x[yd])~ys[yd] )   
  return(tryCatch(coef(summary(mod1))[2, ], error=function(e) c(NA,NA,NA,NA)))
}

#plot with shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

#--------------------------------
#LOAD DATA
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

#============================================
#Fig. 3 LATITUDINAL GRADIENTS IN DATA

#ANALYZE DEVELOPMENTAL CONSTRAINTS
dat=dddat

#plots by order

#restrict to orders with data
dat2=dat[dat$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]
#Neuroptera   Thysanoptera
dat2$pupal= as.factor(dat2$pupal)

p<- ggplot(data=dat2, aes(x=T0, y = DDD, shape=pupal, color=Family ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(0,1600) +xlim(-5,24) +scale_colour_discrete(guide = FALSE)+theme_bw()  #abs(lat)
p1= p + geom_point() 

#by latitude
p<- ggplot(data=dat2, aes(x=abs(lat), y = T0, shape=pupal, color=log(DDD) ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(-5,24) +xlab("Absolute latitude (°)")+theme_bw() #+ scale_color_gradientn(colours=matlab.like(10)) #+scale_shape_discrete(guide = FALSE)
p2= p + geom_point() +geom_smooth(method=lm, se=FALSE)

p<- ggplot(data=dat2, aes(x=abs(lat), y = DDD, shape=pupal, color=log(T0) ))+facet_grid(.~Order) + scale_shape_manual(values = c(1,19)) +ylim(0,1600)+xlab("Absolute latitude (°)") +theme_bw() #+ scale_color_gradientn(colours=matlab.like(10)) #+scale_shape_discrete(guide = FALSE)
p3= p + geom_point() 

#--------------------
setwd(paste(fdir,"figures\\",sep="") )
pdf("T0_DDD.pdf", height = 8, width = 8)

grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size="last"))

dev.off()

#FIG 6. FITNESS SURFACES
#number generations
ngens.ave= rowMeans(ngens, na.rm=TRUE)
#phen.dat2= as.data.frame( cbind(dddat[,c(1:9)],ngens.ave))
phen.dat2= as.data.frame( cbind(dddat[,c(1:6,7,9)],ngens.ave))
phen.dat2= phen.dat2[phen.dat2$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]

#----------------------------
#model

mod1= lm(dat$DDD~ poly(dat$T0) *dat$pupal * abs(dat$lat) *dat$Order)
mod1= lm(dat$DDD~ dat$Order * abs(dat$lat))

dat1= dat[!is.na(dat$lat),]
mod1= lm(dat$T0~ dat$pupal + abs(dat$lat) + abs(dat$lat)^2)
mod1= lm(dat$T0~ dat$Order * abs(dat$lat) )

#===============================
#Plot fitness surface

ngens[is.nan(ngens)] = NA 

phen.dat2= as.data.frame( cbind(1:nrow(dddat), dddat[,c("Species","Order","lon","lat","T0")],ngens))
phen.dat3= as.data.frame( cbind(1:nrow(dddat), dddat[,c("Species","Order","lon","lat","T0")],rowMeans(ngens, na.rm=T) ))
colnames(phen.dat3)[7]= "Ngen"

#restrict to orders with data
phen.dat3=phen.dat3[phen.dat3$Order %in% c("Coleoptera","Diptera","Hemiptera","Homoptera","Hymenoptera","Lepidoptera") ,]

#ave gens across years
setwd(paste(fdir,"figures\\",sep="") )
pdf("Ngen_byLDT.pdf", height = 4, width = 12)

p<- ggplot(data=phen.dat3, aes(x=T0, y = Ngen, color=abs(lat) ))+facet_grid(.~Order)  +ylim(0,20) +xlim(-5,25) +geom_smooth(data = phen.dat3, formula=y~x, aes(x=T0, y = Ngen), method=loess, se=TRUE)
p + geom_point()
dev.off()

#plot Ngen by latitude
setwd(paste(fdir,"figures\\",sep="") )
pdf("Ngen_byLat.pdf", height = 4, width = 12)
p<- ggplot(data=phen.dat3, aes(x=abs(lat), y = Ngen, color=T0 ))+facet_grid(.~Order)  +ylim(0,20) +xlim(0,60) +geom_smooth(data = phen.dat3, formula=y~x, aes(x=abs(lat), y = Ngen), method=loess, se=TRUE)
p + geom_point()
dev.off() 

#subset to Ngen data
phen.dat3= phen.dat3[which(!is.nan(phen.dat3$Ngen) ),]

#aggregate into latitudinal bins
phen.dat3$latbin= cut(phen.dat3$lat,breaks= c(-45,-20, 0,10, 20,30,40,50,61) ) #c(-45,0,20,40,61)

phen.w= phen.dat3 %>%
  group_by( latbin ) %>% 
  do(mod = lm(Ngen ~ T0,na.action=na.omit,  data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2]) %>%
  select(-mod)

#==========================================================
#FITNESS SURFACE

fld <- with(phen.dat3, interp(x = abs(lat), y = T0, z = Ngen, duplicate=TRUE))

gdat <- interp2xyz(fld, data.frame=TRUE)

p3d= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  coord_equal() +
  geom_contour(color = "white", alpha = 0.5) + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Number\ngenerations") + 
  theme_bw()+xlab("latitude (°)")+ylab("T0 (°C)")

#plot
setwd(paste(fdir,"figures\\",sep="") )
pdf("FitnessSurface.pdf", height = 8, width = 8)
p3d
dev.off()

#============================================




