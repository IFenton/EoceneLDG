## Created 13 / 3 / 2014
## Isabel Fenton
## 
## Code to calculate diversity for Paleogene data
## 
## Input files:
## 
##
## Output files:
## 
## To work on:
## 

setwd("C:/Documents/Science/PhD/Project/Eocene/")

# 1. Load in the data -----------------------------------------------------
load("C:/Documents/Science/PhD/Project/Eocene/Outputs/OlComb.RData")
load("C:/Documents/Science/PhD/Project/Eocene/Outputs/EocComb.RData")

# dataset for sites with all species
EO_abun <- merge(EocComb, OlComb, all = TRUE)

nrow(EO_abun) == nrow(EocComb) + nrow(OlComb)
head(EO_abun)

# remove sites which don't have whole species counts
EO_abun <- EO_abun[EO_abun$All != "No", ]
dim(EO_abun)


# 2. Calculate species richness -------------------------------------------

# 2a. Generate a dataframe with each row as a site at a given age ---------
# create a dataframe containing rows for each unique site / age
# add column of site ID / age for sorting
EO_abun$site.age <- NA
for (i in 1:nrow(EO_abun)) {
  # split the ID
  tmp.split <- unlist(strsplit(EO_abun$ID[i], "_"))
  # paste the ID and the age
  # if it is not neptune
  if (length(tmp.split) == 3) {
    EO_abun$site.age[i] <- paste(tmp.split[2], "_", EO_abun$Age[i], sep = "") 
  } else {
    EO_abun$site.age[i] <- paste(tmp.split[2], ".", tmp.split[3], "_", EO_abun$Age[i], sep = "")
  }
}

# one row for each unique combination
EO_abun_div <- EO_abun[!duplicated(EO_abun$site.age), ] 
# remove rows for species and abundance
EO_abun_div <- EO_abun_div[, -c(which(names(EO_abun_div) == "Species"), which(names(EO_abun_div) == "Abundance"))]

# 2b. Calculate the number of species at each site for each age -----------
EO_abun_div$SR <- NA
for(i in 1:nrow(EO_abun_div)) {
  EO_abun_div$SR[i] <- sum(EO_abun$site.age == EO_abun_div$site.age[i])
}

head(EO_abun_div)

with(EO_abun_div, plot(Latitude, SR, pch = 16, col = as.factor(Period)))
par(ask = TRUE)
for (i in unique(EO_abun_div$Period)) {
  with(EO_abun_div[EO_abun_div$Period == i, ], plot(Latitude, SR, pch = 16, col = as.factor(Period), main = i))
}

# 2c. Generate a dataframe with sites by period ---------------------------
# create a dataframe of average richness by period
EO_abun_div$site.period <- NA
for (i in 1:nrow(EO_abun_div)) {
  tmp.split <- unlist(strsplit(EO_abun_div$site.age[i], "_"))
  EO_abun_div$site.period[i] <- paste(tmp.split[1], "_", EO_abun_div$Period[i], sep = "") 
}

EO_abun_age <- EO_abun_div[!duplicated(EO_abun_div$site.period), ] 
# remove sprich column
EO_abun_age <- EO_abun_age[, -21]

# 2d. Calculate average species richness for each site --------------------
EO_abun_age$SR <- NA
for(i in 1:nrow(EO_abun_age)) {
  EO_abun_age$SR[i] <- mean(EO_abun_div$SR[which(EO_abun_div$site.period == EO_abun_age$site.period[i])])
}

# 2e. Plot this up --------------------------------------------------------

with(EO_abun_age, plot(Latitude, SR, pch = 16, col = as.factor(Period)))
par(ask = TRUE)
for (i in unique(EO_abun_age$Period)) {
  with(EO_abun_age[EO_abun_age$Period == i, ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = i))
}

par(ask = FALSE)


# species richness by latitude split by period -------------
par(mfrow = c(2, 3))
with(EO_abun_age[EO_abun_age$Period == "Lower Eocene", ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = "Lower Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene", ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = "Middle Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Upper Eocene", ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = "Upper Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Lower Oligocene", ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = "Lower Oligocene"))
with(EO_abun_age[EO_abun_age$Period == "Upper Oligocene", ], plot(pal.lat.m, SR, pch = 16, col = as.factor(Period), main = "Upper Oligocene"))


# species richness by latitude split by hemisphere and period -------------
par(mfrow = c(2, 3))
with(EO_abun_age[EO_abun_age$Period == "Lower Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (abs(pal.lat.m) == pal.lat.m) + 1, main = "Lower Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (abs(pal.lat.m) == pal.lat.m) + 1, main = "Middle Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Upper Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (abs(pal.lat.m) == pal.lat.m) + 1, main = "Upper Eocene"))
with(EO_abun_age[EO_abun_age$Period == "Lower Oligocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (abs(pal.lat.m) == pal.lat.m) + 1, main = "Lower Oligocene"))
with(EO_abun_age[EO_abun_age$Period == "Upper Oligocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (abs(pal.lat.m) == pal.lat.m) + 1, main = "Upper Oligocene"))


# species richness with regression fitted by period -----------------------
for (i in unique(EO_abun_age$Period)) {
  with(EO_abun_age[EO_abun_age$Period == i, ], plot(abs(pal.lat.m), SR, pch = 16, col = as.factor(Period), main = i))
  abline(lm(SR ~ pal.lat.m, data = EO_abun_age[EO_abun_age$Period == i, ]))
}

head(EO_abun_age$site.period)


# species richness by latitude plotted by period --------------------------
par(mfrow = c(2, 3))
with(EO_abun_age[EO_abun_age$Period == "Lower Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))
with(EO_abun_age[EO_abun_age$Period == "Upper Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))
with(EO_abun_age[EO_abun_age$Period == "Lower Oligocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))
with(EO_abun_age[EO_abun_age$Period == "Upper Oligocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))


# species richness against latitude by period, coloured by hemisphere --------------
par(mfrow = c(2, 3))
with(EO_abun_age[EO_abun_age$Period == "Lower Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
with(EO_abun_age[EO_abun_age$Period == "Upper Eocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
with(EO_abun_age[EO_abun_age$Period == "Lower Oligocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
with(EO_abun_age[EO_abun_age$Period == "Upper Oligocene" & EO_abun_age$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))


# distribution maps of species richness by period -------------------------
with(EO_abun_age[EO_abun_age$Period == "Lower Eocene", ], distrib.map(mod_long51, mod_lat51, SR))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene", ], distrib.map(mod_long42, mod_lat42, SR))
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene", ], distrib.map(mod_long39, mod_lat39, SR))
with(EO_abun_age[EO_abun_adge$Period == "Upper Eocene", ], distrib.map(mod_long35, mod_lat35, SR))


# plot map of points at 39Ma ----------------------------------------------
PTR.39Mya <- readPNG("Output\\39Mya.png")
PTR.39Mya.img <- pixmapGrey(PTR.39Mya)
plot(PTR.39Mya.img)
# use locator() to work out these coordinates
points((-180 + 180) + 135.8, (-90 + 90) + 202, pch = 16)
points((-180 + 180) + 135.8, (90 + 90) /180 * 903.1  + 202, pch = 16)
points((180 + 180) / 360 * 1802.6 + 135.8, (-90 + 90) /180 * 903.1  + 202, pch = 16)
points((180 + 180) / 360 * 1802.6 + 135.8, (90 + 90) /180 * 903.1  + 202, pch = 16)
with(EO_abun_age[EO_abun_age$Period == "Middle Eocene",], points((mod_long39 + 180) / 360 * 1802.6 + 135.8, (mod_lat39 + 90) / 180 * 903.1  + 202, pch = 16, col = "black"))


#  create a smaller dataset -----------------------------------------------
# remove all almost and those with poor preservation
EO_abun_sub <- EO_abun_age[EO_abun_age$All != "Yes (almost)" &EO_abun_age$Preservation != 3 & EO_abun_age$Preservation != 4, ]
table(EO_abun_age$Preservation) # remove 3 and 4


# species richness in eocene ages with average maximum SR --------------------------------

par(mfrow = c(1, 3))
with(EO_abun_sub[EO_abun_sub$Period == "Lower Eocene" & EO_abun_sub$All != "Yes (almost)" , ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
mn.low <- rep(NA, 7)
for (i in 1:7) {
  mn.low[i] <- max(EO_abun_sub$SR[abs(EO_abun_sub$pal.lat.m) < (i * 10) & abs(EO_abun_sub$pal.lat.m) > ((i - 1) * 10) & EO_abun_sub$Period == "Lower Eocene"])
}
points(seq(5, 65, 10), mn.low, type = "l")

with(EO_abun_sub[EO_abun_sub$Period == "Middle Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
mn.mid <- rep(NA, 7)
for (i in 1:7) {
  mn.mid[i] <- max(EO_abun_sub$SR[abs(EO_abun_sub$pal.lat.m) < (i * 10) & abs(EO_abun_sub$pal.lat.m) > ((i - 1) * 10) & EO_abun_sub$Period == "Middle Eocene"])
}
points(seq(5, 65, 10), mn.mid, type = "l")


with(EO_abun_sub[EO_abun_sub$Period == "Upper Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
mn.up <- rep(NA, 7)
for (i in 1:7) {
  mn.up[i] <- max(EO_abun_sub$SR[abs(EO_abun_sub$pal.lat.m) < (i * 10) & abs(EO_abun_sub$pal.lat.m) > ((i - 1) * 10) & EO_abun_sub$Period == "Upper Eocene"])
}
points(seq(5, 65, 10), mn.up, type = "l")


# species richness in eocene ages with regression fitted ----------------------------

par(mfrow = c(1, 3))
with(EO_abun_sub[EO_abun_sub$Period == "Lower Eocene" & EO_abun_sub$All != "Yes (almost)" , ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modL <- with(EO_abun_sub[EO_abun_sub$Period == "Lower Eocene", ], lm(SR ~ poly(abs(pal.lat.m), 2)))
points(0:70, predict(modL, data.frame(pal.lat.m = 0:70)), type = "l")


with(EO_abun_sub[EO_abun_sub$Period == "Middle Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modM <- with(EO_abun_sub[EO_abun_sub$Period == "Middle Eocene", ], lm(SR ~ poly(abs(pal.lat.m), 2)))
points(0:70, predict(modM, data.frame(pal.lat.m = 0:70)), type = "l")

with(EO_abun_sub[EO_abun_sub$Period == "Upper Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modU <- with(EO_abun_sub[EO_abun_sub$Period == "Upper Eocene", ], lm(SR ~ poly(abs(pal.lat.m), 2)))
points(0:70, predict(modU, data.frame(pal.lat.m = 0:70)), type = "l")


# CBEP --------------------------------------------------------------------

# plot only those that are yes and preservation of 0 or n/a



# species richness by lat coloured by All? --------------------------------
par(mfrow = c(1, 3))
with(EO_abun_sub[EO_abun_sub$Period == "Lower Eocene" & EO_abun_sub$All != "Yes (almost)" , ], plot(abs(pal.lat.m), SR, pch = 16, col = All))

with(EO_abun_sub[EO_abun_sub$Period == "Middle Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))

with(EO_abun_sub[EO_abun_sub$Period == "Upper Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = All))

# species richness by lat coloured by Preservation? --------------------------------
par(mfrow = c(1, 3))
with(EO_abun_sub[EO_abun_sub$Period == "Lower Eocene" & EO_abun_sub$All != "Yes (almost)" , ], plot(abs(pal.lat.m), SR, pch = 16, col = Preservation))

with(EO_abun_sub[EO_abun_sub$Period == "Middle Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = Preservation))

with(EO_abun_sub[EO_abun_sub$Period == "Upper Eocene" & EO_abun_sub$All != "Yes (almost)", ], plot(abs(pal.lat.m), SR, pch = 16, col = Preservation))


# genus level richness ----------------------------------------------------
summary(EO_abun_div$Preservation)

# remove all poor preservation
EO_abun_CBEP <- EO_abun_div[EO_abun_div$Preservation != 2 & EO_abun_div$Preservation != 3 & EO_abun_div$Preservation != 4 & EO_abun_div$Preservation != 5, ]
summary(EO_abun_CBEP$Preservation)
dim(EO_abun_div)
dim(EO_abun_CBEP)

# remove all incomplete data
summary(EO_abun_CBEP$All)
EO_abun_CBEP <- EO_abun_CBEP[EO_abun_CBEP$All != "Yes (almost)" & EO_abun_CBEP$All != "Yes (maybe)", ]
summary(EO_abun_CBEP$All)

# calculate SR
EO_abun_CBEP$SR <- NA
for(i in 1:nrow(EO_abun_CBEP)) {
  EO_abun_CBEP$SR[i] <- length(unique(EO_abun$Corr.species[EO_abun$site.age == EO_abun_CBEP$site.age[i]]))
}
# coloured by hemisphere
par(mfrow = c(1, 3))
with(EO_abun_CBEP[EO_abun_CBEP$Period == "Lower Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modL <- with(EO_abun_CBEP[EO_abun_CBEP$Period == "Lower Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modL, data.frame(pal.lat.m = 0:70)), type = "l")
with(EO_abun_CBEP[EO_abun_CBEP$Period == "Lower Eocene", ], anova(modL, lm(SR ~ 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Middle Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modM <- with(EO_abun_CBEP[EO_abun_CBEP$Period == "Middle Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modM, data.frame(pal.lat.m = 0:70)), type = "l")
with(EO_abun_CBEP[EO_abun_CBEP$Period == "Middle Eocene", ], anova(modM, lm(SR ~ 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Upper Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modU <- with(EO_abun_CBEP[EO_abun_CBEP$Period == "Upper Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modU, data.frame(pal.lat.m = 0:70)), type = "l")
with(EO_abun_CBEP[EO_abun_CBEP$Period == "Upper Eocene", ], anova(modU, lm(SR ~ 1)))

# coloured by source
par(mfrow = c(1, 3))
with(EO_abun_CBEP[EO_abun_CBEP$Period == "Lower Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Middle Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Upper Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))


# add genus level richness
EO_abun_CBEP$GR <- NA
for(i in 1:nrow(EO_abun_CBEP)) {
  EO_abun_CBEP$GR[i] <- length(unique(grep("^[A-Z]", unlist(strsplit(EO_abun$Corr.species[EO_abun$site.age == EO_abun_CBEP$site.age[i]], " ")), value = TRUE)))
}

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Lower Eocene", ], plot(abs(pal.lat.m), GR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Middle Eocene", ], plot(abs(pal.lat.m), GR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))

with(EO_abun_CBEP[EO_abun_CBEP$Period == "Upper Eocene", ], plot(abs(pal.lat.m), GR, pch = 16, col = sapply(ID, function(x) length(grep("^M", x, value = TRUE)) + 1)))



# calculating maximum diversity for a given time period -------------------
# get a dataset which is the unique values of age & lat & long
names(EO_abun_CBEP)
max.abun <- EO_abun_CBEP[!duplicated(EO_abun_CBEP[, c(3:4, 7)]), ] 

# calculate max SR within these
for (i in 1:nrow(max.abun)) {
  max.abun$SR[i] <- max(EO_abun_CBEP$SR[which(EO_abun_CBEP$Latitude == max.abun$Latitude[i] & EO_abun_CBEP$Longitude == max.abun$Longitude[i]  & EO_abun_CBEP$Period == max.abun$Period[i])])
}
head(max.abun$SR)
EO_abun_CBEP[1:10,]

par(mfrow = c(1, 3))
with(max.abun[max.abun$Period == "Lower Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modL <- with(max.abun[max.abun$Period == "Lower Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modL, data.frame(pal.lat.m = 0:70)), type = "l")
with(max.abun[max.abun$Period == "Lower Eocene", ], anova(modL, lm(SR ~ 1)))

with(max.abun[max.abun$Period == "Middle Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modM <- with(max.abun[max.abun$Period == "Middle Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modM, data.frame(pal.lat.m = 0:70)), type = "l")
with(max.abun[max.abun$Period == "Middle Eocene", ], anova(modM, lm(SR ~ 1)))

with(max.abun[max.abun$Period == "Upper Eocene", ], plot(abs(pal.lat.m), SR, pch = 16, col = (pal.lat.m > 0) + 1))
modU <- with(max.abun[max.abun$Period == "Upper Eocene", ], lm(SR ~ abs(pal.lat.m)))
points(0:70, predict(modU, data.frame(pal.lat.m = 0:70)), type = "l")
with(max.abun[max.abun$Period == "Upper Eocene", ], anova(modU, lm(SR ~ 1)))


# 3. Calculate evenness ---------------------------------------------------


# 4. Calculate lineage ages --------------------------------------------------




