## Environment
## 14 / 10 / 2014
## Isabel Fenton

## loading in netcdf files for Eocene analysis, and creating a dataframe

library(ncdf)
source("C:/Documents/Science/PhD/Code/sp_mat_2_df.R")
source("C:/Documents/Science/PhD/Code/palettes.R")
source("C:/Documents/Science/PhD/Code/maps.R")
setwd("C:/Documents/Science/PhD/Project/Eocene/Data/Environment")

## 1. Load in a test data file ---------------------------------------------
# open the file
anncl.dl <- open.ncdf("1.1Danian/annual/tdluao.pfclann.nc")

# print the variables it contains
print(anncl.dl)

# get the relevant variable
ocean.temp <- get.var.ncdf(anncl.dl, "temp_mm_uo")
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

close.ncdf(anncl.dl)


## 2. Work out what is in the different data files -------------------------

## 2a. what do pfcl, pfsd and pgcl mean? -----------------------------------
# open the file
tmp.1 <- open.ncdf("1.1Danian/annual/tdluao.pfclann.nc")
tmp.2 <- open.ncdf("1.1Danian/annual/tdluao.pfsdann.nc")
tmp.3 <- open.ncdf("1.1Danian/annual/tdluao.pgclann.nc")

# print the variables it contains
print(tmp.1)
print(tmp.2)
print(tmp.3)

# check the values
v.tmp1 <- get.var.ncdf(tmp.1, "temp_mm_uo")
summary(as.numeric(v.tmp1))
# range -0.58 to 35.17, so temperature

v.tmp2 <- get.var.ncdf(tmp.2, "temp_mm_uo")
summary(as.numeric(v.tmp2))
# range 0.158 to 1.501, so likely to be sd

v.tmp3 <- get.var.ncdf(tmp.3, "temp_ym_dpth")
summary(v.tmp3)
# range -1.37 to 35.17, so temperature

## 2b. Check that pfsd is standard deviation -------------------------------
# check the logic that this is sd by loading in the data by season
tmp.s1 <- open.ncdf("1.1 Danian/quarterly/tdluao.pfcldjf.nc")
tmp.s2 <- open.ncdf("1.1 Danian/quarterly/tdluao.pfclmam.nc")
tmp.s3 <- open.ncdf("1.1 Danian/quarterly/tdluao.pfcljja.nc")
tmp.s4 <- open.ncdf("1.1 Danian/quarterly/tdluao.pfclson.nc")

# print the variables it contains
print(tmp.s1)
print(tmp.s2)
print(tmp.s3)
print(tmp.s4)

v.tmps1 <- get.var.ncdf(tmp.s1, "temp_mm_uo")
v.tmps2 <- get.var.ncdf(tmp.s2, "temp_mm_uo")
v.tmps3 <- get.var.ncdf(tmp.s3, "temp_mm_uo")
v.tmps4 <- get.var.ncdf(tmp.s4, "temp_mm_uo")

# pick a random site (can't pick 1,1 as this is NA)
v.tmps1[50, 50]
v.tmps2[50, 50]
v.tmps3[50, 50]
v.tmps4[50, 50]

mean(c(v.tmps1[50, 50], v.tmps2[50, 50], v.tmps3[50, 50], v.tmps4[50, 50]))
v.tmp1[50, 50]

sd(c(v.tmps1[50, 50], v.tmps2[50, 50], v.tmps3[50, 50], v.tmps4[50, 50]))
v.tmp2[50, 50]

# this is not the sd of seasons. Is it the sd of months
tmp.m1 <- open.ncdf("1.1 Danian/monthly/tdluao.pfcljan.nc")
tmp.m2 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclfeb.nc")
tmp.m3 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclmar.nc")
tmp.m4 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclapr.nc")
tmp.m5 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclmay.nc")
tmp.m6 <- open.ncdf("1.1 Danian/monthly/tdluao.pfcljun.nc")
tmp.m7 <- open.ncdf("1.1 Danian/monthly/tdluao.pfcljul.nc")
tmp.m8 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclaug.nc")
tmp.m9 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclsep.nc")
tmp.m10 <- open.ncdf("1.1 Danian/monthly/tdluao.pfcloct.nc")
tmp.m11 <- open.ncdf("1.1 Danian/monthly/tdluao.pfclnov.nc")
tmp.m12 <- open.ncdf("1.1 Danian/monthly/tdluao.pfcldec.nc")

# print the variables it contains
print(tmp.m1)
print(tmp.m2)
print(tmp.m3)
print(tmp.m4)
print(tmp.m5)
print(tmp.m6)
print(tmp.m7)
print(tmp.m8)
print(tmp.m9)
print(tmp.m10)
print(tmp.m11)
print(tmp.m12)

v.tmpm1 <- get.var.ncdf(tmp.m1, "temp_mm_uo")
v.tmpm2 <- get.var.ncdf(tmp.m2, "temp_mm_uo")
v.tmpm3 <- get.var.ncdf(tmp.m3, "temp_mm_uo")
v.tmpm4 <- get.var.ncdf(tmp.m4, "temp_mm_uo")
v.tmpm5 <- get.var.ncdf(tmp.m5, "temp_mm_uo")
v.tmpm6 <- get.var.ncdf(tmp.m6, "temp_mm_uo")
v.tmpm7 <- get.var.ncdf(tmp.m7, "temp_mm_uo")
v.tmpm8 <- get.var.ncdf(tmp.m8, "temp_mm_uo")
v.tmpm9 <- get.var.ncdf(tmp.m9, "temp_mm_uo")
v.tmpm10 <- get.var.ncdf(tmp.m10, "temp_mm_uo")
v.tmpm11 <- get.var.ncdf(tmp.m11, "temp_mm_uo")
v.tmpm12 <- get.var.ncdf(tmp.m12, "temp_mm_uo")

# pick a random site (can't pick 1,1 as this is NA)
mean(c(v.tmpm1[50, 50], v.tmpm2[50, 50], v.tmpm3[50, 50], v.tmpm4[50, 50], v.tmpm5[50, 50], v.tmpm6[50, 50], v.tmpm7[50, 50], v.tmpm8[50, 50], v.tmpm9[50, 50], v.tmpm10[50, 50], v.tmpm11[50, 50], v.tmpm12[50, 50]))
v.tmp1[50, 50]

sd(c(v.tmpm1[50, 50], v.tmpm2[50, 50], v.tmpm3[50, 50], v.tmpm4[50, 50], v.tmpm5[50, 50], v.tmpm6[50, 50], v.tmpm7[50, 50], v.tmpm8[50, 50], v.tmpm9[50, 50], v.tmpm10[50, 50], v.tmpm11[50, 50], v.tmpm12[50, 50]))
v.tmp2[50, 50]

# doesn't appear to be that simple. So not quite sure what the sd values are, but probably don't need them.

# close the files
close.ncdf(tmp.1)
close.ncdf(tmp.2)
close.ncdf(tmp.3)
close.ncdf(tmp.s1)
close.ncdf(tmp.s2)
close.ncdf(tmp.s3)
close.ncdf(tmp.s4)
close.ncdf(tmp.m1)
close.ncdf(tmp.m2)
close.ncdf(tmp.m3)
close.ncdf(tmp.m4)
close.ncdf(tmp.m5)
close.ncdf(tmp.m6)
close.ncdf(tmp.m7)
close.ncdf(tmp.m8)
close.ncdf(tmp.m9)
close.ncdf(tmp.m10)
close.ncdf(tmp.m11)
close.ncdf(tmp.m12)


## 3. Set up three dataframes for the environmental variables -----------------------------------
# want values for mean SST, Mixed layer depth, 10deg isotherm, sd SST, sd Salinity, mean salinity
# so most of this data comes from pfclann. sd data comes from monthly pfcl. isotherm data comes from pgcl.

# get lat / long extent
ypr.env <- ocean.temp.dat[, 1:2]
lut.env <- ocean.temp.dat[, 1:2]
bar.env <- ocean.temp.dat[, 1:2]
pri.env <- ocean.temp.dat[, 1:2]

folders <- c("2.1Ypresian", "2.2Lutetian", "2.3Bartonian", "2.4Priabonian")
eoc.df <- c("ypr.env", "lut.env", "bar.env", "pri.env")


# 4. Create a function to add the values to a dataframe -------------------
ncdf.env <- function(filepath, filename, varname) {
  # input:  filepath (string)- set the wd to this; paste(folders[i], "/annual/", sep = "")
  #         filename (string) - ncdf file; grep("pfclann", dir(), value = T)
  #         varname (string) - from the ncdf file
  # output: dataframe with the lat, long, values 
  # save old wd
  oldwd <- getwd()
  # set the wd
  setwd(filepath)
  # open the relevant file
  tmp.ncdf <- open.ncdf(filename)
  on.exit(close.ncdf(tmp.ncdf)) 
  setwd(oldwd)
  # extract the relevant variable
  tmp.var <- get.var.ncdf(tmp.ncdf, varname)
  # extract the lat / long
  long <- get.var.ncdf(tmp.ncdf, "longitude")
  long <- ifelse(long > 180, long - 360, long)
  lat <- get.var.ncdf(tmp.ncdf, "latitude")
  # convert this to a dataframe
  if (length(dim(tmp.var)) == 2) {
    return(sp.mat.2.df(long, lat, t(tmp.var)))
  } else {
    depths <- get.var.ncdf(tmp.ncdf, "depth_1")
    df <- sp.mat.2.df(long, lat, t(tmp.var[, , 1]))
    names(df)[ncol(df)] <- round(depths[1])
    for(i in 2:dim(tmp.var)[3]) {
      df <- cbind(df, sp.mat.2.df(long, lat, t(tmp.var[, , i]))[3])
      names(df)[ncol(df)] <- round(depths[i])
    }
    return(df)
  }
}


## 5. mean SST ------------------------------------------------------------------

## 5a. work out the variable name for mean SST ----------------------------------
# open one of the files
anncl.dl <- open.ncdf("1.1 Danian/annual/tdluao.pfclann.nc")
print(anncl.dl)
# two possible variables
ocean.temp <- get.var.ncdf(anncl.dl, "temp_mm_uo") # ocean top-level temperature (K) - apparently though appears to be in deg C
ocean.temp2 <- get.var.ncdf(anncl.dl, "temp_mm_dpth") # potential temperature (ocean) (deg C)

# lots of columns of data, look at one
summary(ocean.temp)[,20]
summary(ocean.temp2)[,20]

# these appear to be identical, so check
sum(ocean.temp != ocean.temp2, na.rm = T)
# 3843, suggesting they are not always identical. Question is why?
tmp <- which(ocean.temp[, 20] != ocean.temp2[, 20])
cbind(ocean.temp[tmp, 20], ocean.temp2[tmp, 20])

# differences appear to only be in the last value of the decimals, so basically identical. Therefore pick one: ocean top-level temperature

close(anncl.dl)

## 5b. add mean SST as columns to the dataframes ---------------------------
tmp <- ncdf.env(paste(folders[1], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "temp_mm_uo")
ypr.env$mnSST <- tmp[, 3]

tmp <- ncdf.env(paste(folders[2], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "temp_mm_uo")
lut.env$mnSST <- tmp[, 3]

tmp <- ncdf.env(paste(folders[3], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "temp_mm_uo")
bar.env$mnSST <- tmp[, 3]

tmp <- ncdf.env(paste(folders[4], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "temp_mm_uo")
pri.env$mnSST <- tmp[, 3]

## 5c. check these look sensible on maps -----------------------------------
with(ypr.env[!is.na(ypr.env$mnSST),], distrib.map(Long, Lat, mnSST, col.land = "steelblue2", pch = 15))
with(lut.env[!is.na(lut.env$mnSST),], distrib.map(Long, Lat, mnSST, col.land = "steelblue2", pch = 15))
with(bar.env[!is.na(bar.env$mnSST),], distrib.map(Long, Lat, mnSST, col.land = "steelblue2", pch = 15))
with(pri.env[!is.na(pri.env$mnSST),], distrib.map(Long, Lat, mnSST, col.land = "steelblue2", pch = 15))

with(ypr.env, plot(mnSST ~ Lat, pch = "."))
with(lut.env, points(mnSST ~ Lat, pch = ".", col = 2))
with(bar.env, points(mnSST ~ Lat, pch = ".", col = 3))
with(pri.env, points(mnSST ~ Lat, pch = ".", col = 4))

with(ypr.env, points(-90:90, predict(gam(mnSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.env, points(-90:90, predict(gam(mnSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.env, points(-90:90, predict(gam(mnSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.env, points(-90:90, predict(gam(mnSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

with(ldg.m.data, points(meanSST.1deg ~ Lat, pch = ".", col = 5))
with(ldg.m.data, points(-90:90, predict(gam(meanSST.1deg ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))


## 6. Data for MLD -------------------------------------------------------

## 6a. Load in the data ----------------------------------------------------
print(anncl.dl)

tmp <- ncdf.env(paste(folders[1], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "mixLyrDpth_mm_uo")
ypr.env$MLD <- tmp[, 3]

tmp <- ncdf.env(paste(folders[2], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "mixLyrDpth_mm_uo")
lut.env$MLD <- tmp[, 3]

tmp <- ncdf.env(paste(folders[3], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "mixLyrDpth_mm_uo")
bar.env$MLD <- tmp[, 3]

tmp <- ncdf.env(paste(folders[4], "/annual/", sep = ""), grep("pfclann", dir(), value = T), "mixLyrDpth_mm_uo")
pri.env$MLD <- tmp[, 3]

## 6b. check these look sensible -------------------------------------------
with(ypr.env[!is.na(ypr.env$MLD),], distrib.map(Long, Lat, MLD, col.land = "steelblue2", pch = 15))
with(lut.env[!is.na(lut.env$MLD),], distrib.map(Long, Lat, MLD, col.land = "steelblue2", pch = 15))
with(bar.env[!is.na(bar.env$MLD),], distrib.map(Long, Lat, MLD, col.land = "steelblue2", pch = 15))
with(pri.env[!is.na(pri.env$MLD),], distrib.map(Long, Lat, MLD, col.land = "steelblue2", pch = 15))

with(ypr.env, plot(MLD ~ Lat, pch = "."))
with(lut.env, points(MLD ~ Lat, pch = ".", col = 2))
with(bar.env, points(MLD ~ Lat, pch = ".", col = 3))
with(pri.env, points(MLD ~ Lat, pch = ".", col = 4))

with(ypr.env, points(-90:90, predict(gam(MLD ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.env, points(-90:90, predict(gam(MLD ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.env, points(-90:90, predict(gam(MLD ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.env, points(-90:90, predict(gam(MLD ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

with(ldg.m.data, points(mean.pt ~ Lat, pch = ".", col = 5))
with(ldg.m.data, points(-90:90, predict(gam(mean.pt ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))


## 7. Data for 10deg -------------------------------------------------------

## 7a. Work out how to get the data ----------------------------------------
print(tmp.3)
str(v.tmp3)
tmp.depths <- get.var.ncdf(tmp.3, "depth_1")

# write a function to get the 10deg depth from a row
eoc.depth10 <- function (x) {
  x <- x[names(x) %in% round(tmp.depths)]
  if (sum(!is.na(x)) == 0) {
    return(NA)
  } else {
    tmp.depths[which(abs(x - 10) == min(abs(x - 10), na.rm = T))[1]]
  }
}

## 7b. Load in the data ----------------------------------------------------
ypr.depth <- ncdf.env(paste(folders[1], "/annual/", sep = ""), grep("pgclann", dir(), value = T), "insitu_T_ym_dpth")
ypr.depth$depth10 <- apply(ypr.depth, 1, eoc.depth10)
ypr.env$depth10 <- ypr.depth$depth10

lut.depth <- ncdf.env(paste(folders[2], "/annual/", sep = ""), grep("pgclann", dir(), value = T), "insitu_T_ym_dpth")
lut.depth$depth10 <- apply(lut.depth, 1, eoc.depth10)
lut.env$depth10 <- lut.depth$depth10

bar.depth <- ncdf.env(paste(folders[3], "/annual/", sep = ""), grep("pgclann", dir(), value = T), "insitu_T_ym_dpth")
bar.depth$depth10 <- apply(bar.depth, 1, eoc.depth10)
bar.env$depth10 <- bar.depth$depth10

pri.depth <- ncdf.env(paste(folders[4], "/annual/", sep = ""), grep("pgclann", dir(), value = T), "insitu_T_ym_dpth")
pri.depth$depth10 <- apply(pri.depth, 1, eoc.depth10)
pri.env$depth10 <- pri.depth$depth10

## 7c. check these look sensible -------------------------------------------
with(ypr.env[!is.na(ypr.env$depth10),], distrib.map(Long, Lat, depth10, col.land = "steelblue2", pch = 15))
with(lut.env[!is.na(lut.env$depth10),], distrib.map(Long, Lat, depth10, col.land = "steelblue2", pch = 15))
with(bar.env[!is.na(bar.env$depth10),], distrib.map(Long, Lat, depth10, col.land = "steelblue2", pch = 15))
with(pri.env[!is.na(pri.env$depth10),], distrib.map(Long, Lat, depth10, col.land = "steelblue2", pch = 15))

with(ypr.env, plot(-depth10 ~ Lat, pch = "."))
with(lut.env, points(-depth10 ~ Lat, pch = ".", col = 2))
with(bar.env, points(-depth10 ~ Lat, pch = ".", col = 3))
with(pri.env, points(-depth10 ~ Lat, pch = ".", col = 4))

with(ypr.env, points(-90:90, predict(gam(-depth10 ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.env, points(-90:90, predict(gam(-depth10 ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.env, points(-90:90, predict(gam(-depth10 ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.env, points(-90:90, predict(gam(-depth10 ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

with(ldg.m.data, points(-depth10deg ~ Lat, pch = ".", col = 5))
with(ldg.m.data, points(-90:90, predict(gam(-depth10deg ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))


## 8. sd SST ---------------------------------------------------------------

## 8a. check the data ------------------------------------------------------
dan.aug <- open.ncdf("1.1Danian/monthly/tdluao.pfclaug.nc")

# print the variables it contains
print(dan.aug)
# same as annual but monthly

## 8b. Get the data --------------------------------------------------------
ypr.SST.mon <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "temp_mm_uo")
names(ypr.SST.mon)[3] <- "jan"
ypr.SST.mon$feb <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$mar <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$apr <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$may <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$jun <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$jul <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$aug <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$sep <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$oct <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$nov <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$dec <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "temp_mm_uo")[, 3]
ypr.SST.mon$sd <- apply(ypr.SST.mon[,3:14], 1, sd, na.rm = T)
ypr.env$sdSST <- ypr.SST.mon$sd

lut.SST.mon <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "temp_mm_uo")
names(lut.SST.mon)[3] <- "jan"
lut.SST.mon$feb <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$mar <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$apr <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$may <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$jun <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$jul <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$aug <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$sep <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$oct <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$nov <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$dec <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "temp_mm_uo")[, 3]
lut.SST.mon$sd <- apply(lut.SST.mon[,3:14], 1, sd, na.rm = T)
lut.env$sdSST <- lut.SST.mon$sd

bar.SST.mon <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "temp_mm_uo")
names(bar.SST.mon)[3] <- "jan"
bar.SST.mon$feb <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$mar <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$apr <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$may <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$jun <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$jul <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$aug <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$sep <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$oct <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$nov <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$dec <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "temp_mm_uo")[, 3]
bar.SST.mon$sd <- apply(bar.SST.mon[,3:14], 1, sd, na.rm = T)
bar.env$sdSST <- bar.SST.mon$sd

pri.SST.mon <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "temp_mm_uo")
names(pri.SST.mon)[3] <- "jan"
pri.SST.mon$feb <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$mar <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$apr <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$may <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$jun <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$jul <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$aug <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$sep <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$oct <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$nov <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$dec <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "temp_mm_uo")[, 3]
pri.SST.mon$sd <- apply(pri.SST.mon[,3:14], 1, sd, na.rm = T)
pri.env$sdSST <- pri.SST.mon$sd

# 8c. Check the data ------------------------------------------------------
with(ypr.env[!is.na(ypr.env$sdSST),], distrib.map(Long, Lat, sdSST, col.land = "steelblue2", pch = 15))
with(lut.env[!is.na(lut.env$sdSST),], distrib.map(Long, Lat, sdSST, col.land = "steelblue2", pch = 15))
with(bar.env[!is.na(bar.env$sdSST),], distrib.map(Long, Lat, sdSST, col.land = "steelblue2", pch = 15))
with(pri.env[!is.na(pri.env$sdSST),], distrib.map(Long, Lat, sdSST, col.land = "steelblue2", pch = 15))

with(ypr.env, plot(-sdSST ~ Lat, pch = "."))
with(lut.env, points(-sdSST ~ Lat, pch = ".", col = 2))
with(bar.env, points(-sdSST ~ Lat, pch = ".", col = 3))
with(pri.env, points(-sdSST ~ Lat, pch = ".", col = 4))

with(ypr.env, points(-90:90, predict(gam(-sdSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.env, points(-90:90, predict(gam(-sdSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.env, points(-90:90, predict(gam(-sdSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.env, points(-90:90, predict(gam(-sdSST ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

with(ldg.m.data, points(-sdSST.1deg ~ Lat, pch = ".", col = 5))
with(ldg.m.data, points(-90:90, predict(gam(-sdSST.1deg ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))


# 9. SD salinity ----------------------------------------------------------
dan.aug <- open.ncdf("1.1Danian/monthly/tdluao.pfclaug.nc")

# print the variables it contains
print(dan.aug)
# same as annual but monthly

## 9b. Get the data --------------------------------------------------------
ypr.sal.mon <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "salinity_mm_dpth")
names(ypr.sal.mon)[3] <- "jan"
ypr.sal.mon$feb <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$mar <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$apr <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$may <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$jun <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$jul <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$aug <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$sep <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$oct <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$nov <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$dec <- ncdf.env(paste(folders[1], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "salinity_mm_dpth")[, 3]
ypr.sal.mon$sd <- apply(ypr.sal.mon[,3:14], 1, sd, na.rm = T)
ypr.env$sdsal <- ypr.sal.mon$sd

lut.sal.mon <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "salinity_mm_dpth")
names(lut.sal.mon)[3] <- "jan"
lut.sal.mon$feb <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$mar <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$apr <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$may <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$jun <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$jul <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$aug <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$sep <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$oct <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$nov <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$dec <- ncdf.env(paste(folders[2], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "salinity_mm_dpth")[, 3]
lut.sal.mon$sd <- apply(lut.sal.mon[,3:14], 1, sd, na.rm = T)
lut.env$sdsal <- lut.sal.mon$sd

bar.sal.mon <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "salinity_mm_dpth")
names(bar.sal.mon)[3] <- "jan"
bar.sal.mon$feb <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$mar <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$apr <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$may <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$jun <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$jul <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$aug <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$sep <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$oct <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$nov <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$dec <- ncdf.env(paste(folders[3], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "salinity_mm_dpth")[, 3]
bar.sal.mon$sd <- apply(bar.sal.mon[,3:14], 1, sd, na.rm = T)
bar.env$sdsal <- bar.sal.mon$sd

pri.sal.mon <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljan", dir(), value = T), "salinity_mm_dpth")
names(pri.sal.mon)[3] <- "jan"
pri.sal.mon$feb <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclfeb", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$mar <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclmar", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$apr <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclapr", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$may <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclmay", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$jun <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljun", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$jul <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcljul", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$aug <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclaug", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$sep <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclsep", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$oct <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcloct", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$nov <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfclnov", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$dec <- ncdf.env(paste(folders[4], "/monthly/", sep = ""), grep("pfcldec", dir(), value = T), "salinity_mm_dpth")[, 3]
pri.sal.mon$sd <- apply(pri.sal.mon[,3:14], 1, sd, na.rm = T)
pri.env$sdsal <- pri.sal.mon$sd

# 9c. Check the data ------------------------------------------------------
with(ypr.env[!is.na(ypr.env$sdsal),], distrib.map(Long, Lat, sdsal, col.land = "steelblue2", pch = 15))
with(lut.env[!is.na(lut.env$sdsal),], distrib.map(Long, Lat, sdsal, col.land = "steelblue2", pch = 15))
with(bar.env[!is.na(bar.env$sdsal),], distrib.map(Long, Lat, sdsal, col.land = "steelblue2", pch = 15))
with(pri.env[!is.na(pri.env$sdsal),], distrib.map(Long, Lat, sdsal, col.land = "steelblue2", pch = 15))

with(ypr.env, plot(sdsal ~ Lat, pch = "."))
with(lut.env, points(sdsal ~ Lat, pch = ".", col = 2))
with(bar.env, points(sdsal ~ Lat, pch = ".", col = 3))
with(pri.env, points(sdsal ~ Lat, pch = ".", col = 4))

with(ypr.env, points(-90:90, predict(gam(sdsal ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.env, points(-90:90, predict(gam(sdsal ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.env, points(-90:90, predict(gam(sdsal ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.env, points(-90:90, predict(gam(sdsal ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

with(ldg.m.data, points(sdSal.0m ~ Lat, pch = ".", col = 5))
with(ldg.m.data, points(-90:90, predict(gam(sdSal.0m ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))

