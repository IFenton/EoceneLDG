## Environment
## 14 / 10 / 2014
## Isabel Fenton

## loading in netcdf files for Eocene analysis, and creating a dataframe

library(ncdf)
source("C:/Documents/Science/PhD/Code/sp_mat_2_df.R")
source("C:/Documents/Science/PhD/Code/palettes.R")
source("C:/Documents/Science/PhD/Code/maps.R")
setwd("C:/Documents/Science/PhD/Project/Eocene/Data/Environment")

# 1. Load in temperature data ---------------------------------------------
# open the file
anncl.dl <- open.ncdf("tdluao.pfclann.nc")

# print the variables it contains
print(anncl.dl)

# get the relevant variable
ocean.temp <- get.var.ncdf(anncl.dl, temp_mm_uo)
str(ocean.temp)
summary(ocean.temp)

# get the coordinates
str(anncl.dl$var$temp_mm_uo)
ocean.long <- anncl.dl$var$temp_mm_uo$dim[[1]]$vals
ocean.long <- ifelse(ocean.long > 180, ocean.long - 360, ocean.long)
ocean.lat <- anncl.dl$var$temp_mm_uo$dim[[2]]$vals

# convert this to the right format
ocean.temp.dat <- sp.mat.2.df(ocean.long, ocean.lat, t(ocean.temp))

# plot this up
with(ocean.temp.dat[!is.na(ocean.temp.dat$val),], distrib.map(Long, Lat, val, col.land = "steelblue2", pch = 15))


