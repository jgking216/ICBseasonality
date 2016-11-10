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

#--------------------------
#Write out

write.csv(dat, "ThermalDatabase_all_10Nov2016.csv")