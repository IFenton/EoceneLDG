# modelling

library(spdep)
library(ncf)
source("../../../../Work/1311 LDGPaper/Code/140420SARerrOptimising_NC.R")

# 1. run an ols richness model -----------------------------------------------------
reoc.rsr.l0 <- lm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.rsr.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.rsr.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.rsr.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.rsr.l0.sac)

# run sar model
ldg.coords <- cbind(ldg.m.data$Long,ldg.m.data$Lat)
ldg.coords <- as.matrix(ldg.coords)
reoc.rsr.sarW <- with(ldg.m.data, sar.optimised(reoc.rsr.l0.sac$real$x.intercept, rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(reoc.rsr.sarW$obj, Nagelkerke = TRUE) # 
AIC(reoc.rsr.sarW$obj) # 

# predict
ypr.env$rarefy.sr <- sar.predict(reoc.rsr.sarW, newdata = ypr.env, olddata = ldg.m.data)
summary(ypr.env$rarefy.sr)
with(ypr.env, distrib.map(Long, Lat, rarefy.sr))

lut.env$rarefy.sr <- sar.predict(reoc.rsr.sarW, newdata = lut.env, olddata = ldg.m.data)
summary(lut.env$rarefy.sr)
with(lut.env, distrib.map(Long, Lat, rarefy.sr))

bar.env$rarefy.sr <- sar.predict(reoc.rsr.sarW, newdata = bar.env, olddata = ldg.m.data)
summary(bar.env$rarefy.sr)
with(bar.env, distrib.map(Long, Lat, rarefy.sr))

pri.env$rarefy.sr <- sar.predict(reoc.rsr.sarW, newdata = pri.env, olddata = ldg.m.data)
summary(pri.env$rarefy.sr)
with(pri.env, distrib.map(Long, Lat, rarefy.sr))


# 2. evenness -------------------------------------------------------------
reoc.eve.l0 <- lm(simpsonEve ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.eve.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.eve.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.eve.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.eve.l0.sac)

# run sar model
reoc.eve.sarW <- with(ldg.m.data, sar.optimised(reoc.eve.l0.sac$real$x.intercept, rarefysr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(reoc.eve.sarW$obj, Nagelkerke = TRUE) # 
AIC(reoc.eve.sarW$obj) # 

# predict
ypr.env$simpsonEve <- sar.predict(reoc.eve.sarW, newdata = ypr.env, olddata = ldg.m.data)
summary(ypr.env$simpsonEve)
with(ypr.env, distrib.map(Long, Lat, simpsonEve))

lut.env$simpsonEve <- sar.predict(reoc.eve.sarW, newdata = lut.env, olddata = ldg.m.data)
summary(lut.env$simpsonEve)
with(lut.env, distrib.map(Long, Lat, simpsonEve))

bar.env$simpsonEve <- sar.predict(reoc.eve.sarW, newdata = bar.env, olddata = ldg.m.data)
summary(bar.env$simpsonEve)
with(bar.env, distrib.map(Long, Lat, simpsonEve))

pri.env$simpsonEve <- sar.predict(reoc.eve.sarW, newdata = pri.env, olddata = ldg.m.data)
summary(pri.env$simpsonEve)
with(pri.env, distrib.map(Long, Lat, simpsonEve))

# 3. Average community age ------------------------------------------------
reoc.lna.l0 <- lm(MorphoAgeAbun ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.lna.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.lna.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.lna.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.lna.l0.sac)

# run sar model
reoc.lna.sarW <- with(ldg.m.data, sar.optimised(reoc.lna.l0.sac$real$x.intercept, rarefysr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
summary(reoc.lna.sarW$obj, Nagelkerke = TRUE) # 
AIC(reoc.lna.sarW$obj) # 

# predict
ypr.env$MorphoAgeAbun <- sar.predict(reoc.lna.sarW, newdata = ypr.env, olddata = ldg.m.data)
summary(ypr.env$MorphoAgeAbun)
with(ypr.env, distrib.map(Long, Lat, MorphoAgeAbun))

lut.env$MorphoAgeAbun <- sar.predict(reoc.lna.sarW, newdata = lut.env, olddata = ldg.m.data)
summary(lut.env$MorphoAgeAbun)
with(lut.env, distrib.map(Long, Lat, MorphoAgeAbun))

bar.env$MorphoAgeAbun <- sar.predict(reoc.lna.sarW, newdata = bar.env, olddata = ldg.m.data)
summary(bar.env$MorphoAgeAbun)
with(bar.env, distrib.map(Long, Lat, MorphoAgeAbun))

pri.env$MorphoAgeAbun <- sar.predict(reoc.lna.sarW, newdata = pri.env, olddata = ldg.m.data)
summary(pri.env$MorphoAgeAbun)
with(pri.env, distrib.map(Long, Lat, MorphoAgeAbun))
