# modelling

library(spdep)
library(ncf)
source("C:/Documents/Science/Work/1311 LDGPaper/Code/140420SARerrOptimising_NC.R")
source("C:/Documents/Science/PhD/Code/sar_predict.R")

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
reoc.op.nb <- dnearneigh(ldg.coords, 0, reoc.rsr.sarW$dist, longlat = TRUE)
reoc.op.w <- nb2listw(reoc.op.nb, glist = NULL, style = "W", zero.policy = TRUE)
reoc.rsr.op0 <- errorsarlm(reoc.rsr.sarW$mod, listw = reoc.op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(reoc.rsr.op0, Nagelkerke = TRUE) # 0.82286
AIC(reoc.rsr.op0) # 2738.962

# check predictions with recent data
recent.rsr <- sar.predict(reoc.rsr.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1]
with(ldg.p.data, distrib.map(Long, Lat, recent.rsr))
# these look reasonable

# predict
ypr.pred <- ypr.env
names(ypr.pred) <- c("Long", "Lat", "meanSST.1deg", "mean.pt", "depth10deg", "sdSST.1deg", "sdSal.0m", "meanSal.0m")
ypr.pred$dissolution <- 0
ypr.pred <- na.omit(ypr.pred)
ypr.pred$rarefy.sr <- sar.predict(reoc.rsr.op0, newdata = ypr.pred, olddata = ldg.m.data)[, 1]
summary(ypr.pred$rarefy.sr)
with(ypr.pred, distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))

lut.pred <- lut.env
names(lut.pred) <- c("Long", "Lat", "meanSST.1deg", "mean.pt", "depth10deg", "sdSST.1deg", "sdSal.0m", "meanSal.0m")
lut.pred$dissolution <- 0
lut.pred <- na.omit(lut.pred)
lut.pred$rarefy.sr <- sar.predict(reoc.rsr.op0, newdata = lut.pred, olddata = ldg.m.data)[, 1]
summary(lut.pred$rarefy.sr)
with(lut.pred, distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))

bar.pred <- bar.env
names(bar.pred) <- c("Long", "Lat", "meanSST.1deg", "mean.pt", "depth10deg", "sdSST.1deg", "sdSal.0m", "meanSal.0m")
bar.pred$dissolution <- 0
bar.pred <- na.omit(bar.pred)
bar.pred$rarefy.sr <- sar.predict(reoc.rsr.op0, newdata = bar.pred, olddata = ldg.m.data)[, 1]
summary(bar.pred$rarefy.sr)
with(bar.pred, distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))

pri.pred <- pri.env
names(pri.pred) <- c("Long", "Lat", "meanSST.1deg", "mean.pt", "depth10deg", "sdSST.1deg", "sdSal.0m", "meanSal.0m")
pri.pred$dissolution <- 0
pri.pred <- na.omit(pri.pred)
pri.pred$rarefy.sr <- sar.predict(reoc.rsr.op0, newdata = pri.pred, olddata = ldg.m.data)[, 1]
summary(pri.pred$rarefy.sr)
with(pri.pred, distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Lat, rarefy.sr, pch = "."))
with(lut.pred, points(Lat, rarefy.sr, pch = ".", col = 2))
with(bar.pred, points(Lat, rarefy.sr, pch = ".", col = 3))
with(pri.pred, points(Lat, rarefy.sr, pch = ".", col = 4))


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


# 4. predicting with only temperature -----------------------------------------------------
temp.mod <- lm(rarefy.sr ~ poly(meanSST.1deg, 3), data = ldg.m.data)
summary(temp.mod)

y.tmp.rsr <- predict(temp.mod, newdata = ypr.pred)
summary(y.tmp.rsr)
with(ypr.pred, distrib.map(Long, Lat, y.tmp.rsr, pch = 15, col.land = "steelblue2"))

l.tmp.rsr <- predict(temp.mod, newdata = lut.pred)
summary(l.tmp.rsr)
with(lut.pred, distrib.map(Long, Lat, l.tmp.rsr, pch = 15, col.land = "steelblue2"))

b.tmp.rsr <- predict(temp.mod, newdata = bar.pred)
summary(b.tmp.rsr)
with(bar.pred, distrib.map(Long, Lat, b.tmp.rsr, pch = 15, col.land = "steelblue2"))

p.tmp.rsr <- predict(temp.mod, newdata = pri.pred)
summary(p.tmp.rsr)
with(pri.pred, distrib.map(Long, Lat, p.tmp.rsr, pch = 15, col.land = "steelblue2"))

r.tmp.rsr <- predict(temp.mod, newdata = ldg.p.data)
summary(r.tmp.rsr)
with(ldg.p.data, distrib.map(Long, Lat, r.tmp.rsr, pch = 15, col.land = "steelblue2"))

with(ldg.p.data, plot(Lat, r.tmp.rsr, pch = ".", col = 5, ylim = c(0, 20)))
with(ypr.pred, points(Lat, y.tmp.rsr, pch = "."))
with(lut.pred, points(Lat, l.tmp.rsr, pch = ".", col = 2))
with(bar.pred, points(Lat, b.tmp.rsr, pch = ".", col = 3))
with(pri.pred, points(Lat, p.tmp.rsr, pch = ".", col = 4))

with(ldg.p.data, points(-90:90, predict(gam(r.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))
with(ypr.pred, points(-90:90, predict(gam(y.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.pred, points(-90:90, predict(gam(l.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.pred, points(-90:90, predict(gam(b.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.pred, points(-90:90, predict(gam(p.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

