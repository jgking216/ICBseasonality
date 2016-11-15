#CALCULATE DD

#fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ICBseasonality\\"
fdir= "C:\\Users\\Buckley\\Google Drive\\Seasonality\\"

#--------------------------------
#LOAD AND CLEAN DATA

setwd(paste(fdir,"Database\\",sep="") )

dat= read.csv("ThermalDatabase_all_9Nov2016.csv", na.strings="")

ind_new= paste(dat$"Species",dat$"Author",dat$"Year",dat$"Location_new",sep="_")
ind_old= paste(dat$"Species",dat$"Author",dat$"Year",dat$"Location",sep="_")

#Use old locations where new is na
na_new= is.na(dat$"Location_new")
ind_new[na_new]=ind_old[na_new]

ind= unique(ind_new)

dat$index= match(ind_new, ind)

#change columns to numeric #CHECK CONVERSION
dat$lon_new= as.numeric(as.character(dat$lon_new))
dat$lat_new= as.numeric(as.character(dat$lat_new))
dat$quality= as.numeric(as.character(dat$quality))

#--------------------------

#Add foreign language locations
dat.lang= read.csv("ThermalDatabase_Fall2016_Languages_10Nov2016.csv", na.strings="")

#change columns to numeric #CHECK CONVERSION
dat.lang$lon_new= as.numeric(as.character(dat.lang$lon_new))
dat.lang$lat_new= as.numeric(as.character(dat.lang$lat_new))
dat.lang$quality= as.numeric(as.character(dat.lang$quality))
 
match1= match(dat.lang$Order, dat$Order.1)
dat[match1, c("Location_new","lon_new","lat_new","quality","comments","omit", "colony")]= dat.lang[, c("Location_new","lon_new","lat_new","quality","comments","omit", "colony")]


#--------------------------
#Fill in missing new locations with old locations and code quality as 4.

dat$Location_new= as.character(dat$Location_new)
dat$Location= as.character(dat$Location)

inds= which(is.na(dat$Location_new ))
dat[inds,"Location_new"]= dat[inds,"Location"]
dat[inds,"quality"]= 4

#Write out
write.csv(dat, "ThermalDatabase_updatelocs.csv")

#=======================================================
# INITIAL CODE TO ADD COLUMNS

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

