# modelling

library(spdep)
library(ncf)
library(mgcv)
source("C:/Documents/Science/PhD/Work/1311 LDGPaper/Code/140420SARerrOptimising_NC.R")
source("C:/Documents/Science/PhD/Code/sar_predict.R")
load("../1311 LDGPaper/Reanalysis/Outputs/ldg_p_margo_pred.RData")
load("../1311 LDGPaper/Reanalysis/Outputs/margo_mod.RData")


# 1. try the model predictions -----------------------------------------------------
#run an ols richness model
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
# very high points, might be driven by salinity


## 2. Predicting with logs -------------------------------------------------
# the data is giving very strange predictions, so try with logged depths and richness (to remove negatives)

reoc.rsr.log.l0 <- lm(log(rarefy.sr + 1) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + log(mean.pt + 1) + log(depth10deg + 1) + meanSal.0m + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.rsr.log.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.rsr.log.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.rsr.log.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.rsr.log.l0.sac)

# run sar model
ldg.coords <- cbind(ldg.m.data$Long,ldg.m.data$Lat)
ldg.coords <- as.matrix(ldg.coords)
reoc.rsr.log.sarW <- with(ldg.m.data, sar.optimised(reoc.rsr.log.l0.sac$real$x.intercept,log(rarefy.sr + 1) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + log(mean.pt + 1) + log(depth10deg + 1) + meanSal.0m + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
reoc.log.op.nb <- dnearneigh(ldg.coords, 0, reoc.rsr.log.sarW$dist, longlat = TRUE)
reoc.log.op.w <- nb2listw(reoc.log.op.nb, glist = NULL, style = "W", zero.policy = TRUE)
reoc.rsr.log.op0 <- errorsarlm(reoc.rsr.log.sarW$mod, listw = reoc.log.op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(reoc.rsr.log.op0, Nagelkerke = TRUE) # 0.89483
AIC(reoc.rsr.log.op0) # -697.957

# check predictions with recent data
recent.log.rsr <- exp(sar.predict(reoc.rsr.log.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1] - 1)
with(ldg.p.data, distrib.map(Long, Lat, recent.log.rsr))
with(ldg.p.data[which(recent.log.rsr > 50), ], distrib.map(Long, Lat, recent.log.rsr[recent.log.rsr > 50]))
with(ldg.p.data[which(recent.log.rsr < 50), ], distrib.map(Long, Lat, recent.log.rsr[recent.log.rsr < 50]))
# this fails at extreme values, so get high richness predicted in northern high latitudes

# predict
ypr.pred$log.rsr <- exp(sar.predict(reoc.rsr.log.op0, newdata = ypr.pred, olddata = ldg.m.data)[, 1] - 1)
summary(ypr.pred$log.rsr)
with(ypr.pred, distrib.map(Long, Lat, log.rsr, pch = 15, col.land = "steelblue2"))

lut.pred$log.rsr <- exp(sar.predict(reoc.rsr.log.op0, newdata = lut.pred, olddata = ldg.m.data)[, 1] - 1)
summary(lut.pred$log.rsr)
with(lut.pred, distrib.map(Long, Lat, log.rsr, pch = 15, col.land = "steelblue2"))

bar.pred$log.rsr <- exp(sar.predict(reoc.rsr.log.op0, newdata = bar.pred, olddata = ldg.m.data)[, 1] - 1)
summary(bar.pred$log.rsr)
with(bar.pred, distrib.map(Long, Lat, log.rsr, pch = 15, col.land = "steelblue2"))

pri.pred$log.rsr <- exp(sar.predict(reoc.rsr.log.op0, newdata = pri.pred, olddata = ldg.m.data)[, 1] - 1)
summary(pri.pred$log.rsr)
with(pri.pred, distrib.map(Long, Lat, log.rsr, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Lat, log.rsr, pch = "."))
with(lut.pred, points(Lat, log.rsr, pch = ".", col = 2))
with(bar.pred, points(Lat, log.rsr, pch = ".", col = 3))
with(pri.pred, points(Lat, log.rsr, pch = ".", col = 4))
# very unconvincing


## 3. Only logging richness ------------------------------------------------
# the data is giving very strange predictions, so try with logged depths and richness (to remove negatives)
# linear model
reoc.rsr.log2.l0 <- lm(log(rarefy.sr + 1) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.rsr.log2.l0)

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.rsr.log2.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.rsr.log2.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.rsr.log2.l0.sac)

# run sar model
ldg.coords <- cbind(ldg.m.data$Long,ldg.m.data$Lat)
ldg.coords <- as.matrix(ldg.coords)
reoc.rsr.log2.sarW <- with(ldg.m.data, sar.optimised(reoc.rsr.log2.l0.sac$real$x.intercept,log(rarefy.sr + 1) ~ (poly(meanSST.1deg, 3) + sdSST.1deg + mean.pt + depth10deg + meanSal.0m + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
reoc.log2.op.nb <- dnearneigh(ldg.coords, 0, reoc.rsr.log2.sarW$dist, longlat = TRUE)
reoc.log2.op.w <- nb2listw(reoc.log2.op.nb, glist = NULL, style = "W", zero.policy = TRUE)
reoc.rsr.log2.op0 <- errorsarlm(reoc.rsr.log2.sarW$mod, listw = reoc.log2.op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(reoc.rsr.log2.op0, Nagelkerke = TRUE) # 0.89462 
AIC(reoc.rsr.log2.op0) # -696.6558

# check predictions with recent data
recent.log2.rsr <- exp(sar.predict(reoc.rsr.log2.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1] - 1)
with(ldg.p.data, distrib.map(Long, Lat, recent.log2.rsr))
with(ldg.p.data[which(recent.log2.rsr > 30), ], distrib.map(Long, Lat, recent.log2.rsr[recent.log2.rsr > 50]))
with(ldg.p.data[which(recent.log2.rsr < 30), ], distrib.map(Long, Lat, recent.log2.rsr[recent.log2.rsr < 50]))
# still giving some slightly strange predictions
# don't think sarlm can handle poisson and the data isn't integer

# predict
ypr.pred$log2.rsr <- exp(sar.predict(reoc.rsr.log2.op0, newdata = ypr.pred, olddata = ldg.m.data)[, 1] - 1)
summary(ypr.pred$log2.rsr)
with(ypr.pred, distrib.map(Long, Lat, log2.rsr, pch = 15, col.land = "steelblue2"))

lut.pred$log2.rsr <- exp(sar.predict(reoc.rsr.log2.op0, newdata = lut.pred, olddata = ldg.m.data)[, 1] - 1)
summary(lut.pred$log2.rsr)
with(lut.pred, distrib.map(Long, Lat, log2.rsr, pch = 15, col.land = "steelblue2"))

bar.pred$log2.rsr <- exp(sar.predict(reoc.rsr.log2.op0, newdata = bar.pred, olddata = ldg.m.data)[, 1] - 1)
summary(bar.pred$log2.rsr)
with(bar.pred, distrib.map(Long, Lat, log2.rsr, pch = 15, col.land = "steelblue2"))

pri.pred$log2.rsr <- exp(sar.predict(reoc.rsr.log2.op0, newdata = pri.pred, olddata = ldg.m.data)[, 1] - 1)
summary(pri.pred$log2.rsr)
with(pri.pred, distrib.map(Long, Lat, log2.rsr, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Lat, log2.rsr, pch = "."))
with(lut.pred, points(Lat, log2.rsr, pch = ".", col = 2))
with(bar.pred, points(Lat, log2.rsr, pch = ".", col = 3))
with(pri.pred, points(Lat, log2.rsr, pch = ".", col = 4))

# basically ridiculous


## 4. Try cutting out the extra salinity -----------------------------------
# use the same model, but exclude points where salinity is < 25psu
with(ypr.pred[ypr.pred$meanSal.0m > 25,], distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))
with(lut.pred[ypr.pred$meanSal.0m > 25,], distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))
with(bar.pred[ypr.pred$meanSal.0m > 25,], distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))
with(pri.pred[ypr.pred$meanSal.0m > 25,], distrib.map(Long, Lat, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(ypr.pred[ypr.pred$meanSal.0m > 25,], plot(Lat, rarefy.sr, pch = "."))
with(lut.pred[ypr.pred$meanSal.0m > 25,], points(Lat, rarefy.sr, pch = ".", col = 2))
with(bar.pred[ypr.pred$meanSal.0m > 25,], points(Lat, rarefy.sr, pch = ".", col = 3))
with(pri.pred[ypr.pred$meanSal.0m > 25,], points(Lat, rarefy.sr, pch = ".", col = 4))


## 5. predicting with only temperature -----------------------------------------------------
load("../1311 LDGPaper/Reanalysis/Outputs/ldg_margo_mod.RData")
load("../1311 LDGPaper/Reanalysis/Outputs/margo_mod.RData")

# this is what I did have
temp.mod.orig <- lm(rarefy.sr ~ poly(meanSST.1deg, 3), data = rsr.margo.mod)

temp.mod <- with(rsr.margo.mod, gam(rarefy.sr ~ s(Longitude, Latitude, k = 80) + s(meanSST.1deg), gamma = 1.4))
summary(temp.mod)

# look for spatial autocorrelation in the residuals
# using spline.correlog
temp.mod.sac <- with(rsr.margo.mod, spline.correlog(Long, Lat, temp.mod$residuals, latlon = TRUE, resamp = 10))
summary(temp.mod.sac)

# run sar model
ldg.coords <- cbind(rsr.margo.mod$Longitude, rsr.margo.mod$Latitude)
ldg.coords <- as.matrix(ldg.coords)
temp.mod.sarW <- with(rsr.margo.mod, sar.optimised(temp.mod.sac$real$x.intercept, rarefy.sr ~ poly(meanSST.1deg, 3), ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
temp.mod.nb <- dnearneigh(ldg.coords, 0, temp.mod.sarW$dist, longlat = TRUE)
temp.mod.w <- nb2listw(temp.mod.nb, glist = NULL, style = "W", zero.policy = TRUE)
temp.mod.op0 <- errorsarlm(temp.mod.sarW$mod, listw = temp.mod.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(temp.mod.op0, Nagelkerke = TRUE) # 0.86296 
AIC(temp.mod.op0) # 4496.966

# check predictions with recent data
r.tmp.rsr <- sar.predict(temp.mod.op0, newdata = ldg.p.data, olddata = rsr.margo.mod)[, 1]
summary(r.tmp.rsr)
with(ldg.p.data, distrib.map(Long, Lat, r.tmp.rsr))
# these look reasonable

names(eEoc.env)[3] <- "meanSST.1deg"
eEoc.env <- na.omit(eEoc.env)
eEoc.tmp.rsr <- sar.predict(temp.mod.op0, newdata = eEoc.env, olddata = rsr.margo.mod)[, 1]
summary(eEoc.tmp.rsr)
with(eEoc.env, distrib.map(Long, Lat, eEoc.tmp.rsr, pch = 15, col.land = "steelblue2"))

names(mEoc.env)[3] <- "meanSST.1deg"
mEoc.env <- na.omit(mEoc.env)
mEoc.tmp.rsr <- sar.predict(temp.mod.op0, newdata = mEoc.env, olddata = rsr.margo.mod)[, 1]
summary(mEoc.tmp.rsr)
with(mEoc.env, distrib.map(Long, Lat, mEoc.tmp.rsr, pch = 15, col.land = "steelblue2"))

names(lEoc.env)[3] <- "meanSST.1deg"
lEoc.env <- na.omit(lEoc.env)
lEoc.tmp.rsr <- sar.predict(temp.mod.op0, newdata = lEoc.env, olddata = rsr.margo.mod)[, 1]
summary(lEoc.tmp.rsr)
with(lEoc.env, distrib.map(Long, Lat, lEoc.tmp.rsr, pch = 15, col.land = "steelblue2"))

with(ldg.p.data, plot(Lat, r.tmp.rsr, pch = ".", col = 5, ylim = c(0, 20)))
with(eEoc.env, points(Lat, eEoc.tmp.rsr, pch = "."))
with(mEoc.env, points(Lat, mEoc.tmp.rsr, pch = ".", col = 2))
with(lEoc.env, points(Lat, lEoc.tmp.rsr, pch = ".", col = 3))

with(ldg.p.data, points(-90:90, predict(gam(r.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))
with(eEoc.env, points(-90:90, predict(gam(eEoc.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(mEoc.env, points(-90:90, predict(gam(mEoc.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(lEoc.env, points(-90:90, predict(gam(lEoc.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))



# compare these predictions with the ols model
r.tmp.rsr.ols <- predict(temp.mod, newdata = ldg.p.margo)
summary(r.tmp.rsr.ols)
with(ldg.p.margo, distrib.map(Longitude, Latitude, r.tmp.rsr.ols))

y.tmp.rsr.ols <- ypr.pred$temp_sr <- predict(temp.mod, newdata = ypr.pred)
summary(y.tmp.rsr.ols)
with(ypr.pred, distrib.map(Long, Lat, y.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

l.tmp.rsr.ols <- lut.pred$temp_sr <-  predict(temp.mod, newdata = lut.pred)
summary(l.tmp.rsr.ols)
with(lut.pred, distrib.map(Long, Lat, l.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

b.tmp.rsr.ols <- bar.pred$temp_sr <- predict(temp.mod, newdata = bar.pred)
summary(b.tmp.rsr.ols)
with(bar.pred, distrib.map(Long, Lat, b.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

p.tmp.rsr.ols <- pri.pred$temp_sr <- predict(temp.mod, newdata = pri.pred)
summary(p.tmp.rsr.ols)
with(pri.pred, distrib.map(Long, Lat, p.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

#proxy
e.tmp.rsr.ols <- predict(temp.mod, newdata = temp.bi.early)
u.tmp.rsr.ols <- predict(temp.mod, newdata = temp.bi.late)

# GCM2
eEoc.tmp.rsr.ols <- predict(temp.mod, newdata = eEoc.env)
summary(eEoc.tmp.rsr.ols)
with(eEoc.env, distrib.map(Long, Lat, eEoc.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

mEoc.tmp.rsr.ols <- predict(temp.mod, newdata = mEoc.env)
summary(mEoc.tmp.rsr.ols)
with(mEoc.env, distrib.map(Long, Lat, mEoc.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

lEoc.tmp.rsr.ols <- predict(temp.mod, newdata = lEoc.env)
summary(lEoc.tmp.rsr.ols)
with(lEoc.env, distrib.map(Long, Lat, lEoc.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

par(mfrow = c(1,2))
with(ldg.p.margo, plot(Latitude, r.tmp.rsr.ols, pch = ".", col = 5, ylim = c(0, 20), ylab = "predicted richness"))
with(ypr.pred, points(Lat, y.tmp.rsr.ols, pch = "."))
with(lut.pred, points(Lat, l.tmp.rsr.ols, pch = ".", col = 2))
with(bar.pred, points(Lat, b.tmp.rsr.ols, pch = ".", col = 3))
with(pri.pred, points(Lat, p.tmp.rsr.ols, pch = ".", col = 4))
with(temp.bi.early, points(Lat, e.tmp.rsr.ols, pch = 16, col = 6))
with(temp.bi.late, points(Lat, u.tmp.rsr.ols, pch = 16, col = 7))

with(ldg.p.margo, plot(Latitude, r.tmp.rsr.ols, pch = ".", col = 5, ylim = c(0, 20), ylab = "predicted richness"))
with(eEoc.env, points(Lat, eEoc.tmp.rsr.ols, pch = "."))
with(mEoc.env, points(Lat, mEoc.tmp.rsr.ols, pch = ".", col = 2))
with(lEoc.env, points(Lat, lEoc.tmp.rsr.ols, pch = ".", col = 3))


with(ldg.p.data, points(-65:70, predict(gam(r.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -65:70)), type = "l", col = 5))
with(ypr.pred, points(-90:90, predict(gam(y.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.pred, points(-90:90, predict(gam(l.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.pred, points(-90:90, predict(gam(b.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.pred, points(-90:90, predict(gam(p.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))
with(temp.bi.early, points(-70:80, predict(gam(e.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 6))
with(temp.bi.late, points(-70:80, predict(gam(u.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 7))


# basically the model predicts diversity drops rapidly with higher T
plot(-5:50, predict(temp.mod, newdata = data.frame(meanSST.1deg = -5:50)))


# 6. Andy's model ---------------------------------------------------------
reoc.rsr.log3.l0 <- lm(rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + log(mean.pt + 1) + log(depth10deg + 1) + poly(meanSal.0m,3) + sdSal.0m + dissolution)^2, data = ldg.m.data)
summary(reoc.rsr.log3.l0)

predict()

# look for spatial autocorrelation in the residuals
# using spline.correlog
reoc.rsr.log3.l0.sac <- with(ldg.m.data, spline.correlog(Long, Lat, reoc.rsr.log3.l0$residuals, latlon = TRUE, resamp = 10))
summary(reoc.rsr.log3.l0.sac)

# run sar model
ldg.coords <- cbind(ldg.m.data$Long,ldg.m.data$Lat)
ldg.coords <- as.matrix(ldg.coords)
reoc.rsr.log3.sarW <- with(ldg.m.data, sar.optimised(reoc.rsr.log3.l0.sac$real$x.intercept,rarefy.sr ~ (poly(meanSST.1deg, 3) + sdSST.1deg + log(mean.pt + 1) + log(depth10deg + 1) + poly(meanSal.0m,3) + sdSal.0m + dissolution)^2, ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
reoc.log3.op.nb <- dnearneigh(ldg.coords, 0, reoc.rsr.log3.sarW$dist, longlat = TRUE)
reoc.log3.op.w <- nb2listw(reoc.log3.op.nb, glist = NULL, style = "W", zero.policy = TRUE)


reoc.rsr.log3.op0 <- errorsarlm(reoc.rsr.log3.sarW$mod, listw = reoc.log3.op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(reoc.rsr.log3.op0, Nagelkerke = TRUE) # 0.89483

reoc.rsr.log3.op0 <- errorsarlm(update(reoc.rsr.log3.sarW$mod, ~, listw = reoc.log3.op.w, zero.policy = TRUE, tol.solve = 1e-18))
summary(reoc.rsr.log3.op0, Nagelkerke = TRUE) # 0.89483
AIC(reoc.rsr.log3.op0) # -697.957

# check predictions with recent data
recent.log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1]
with(ldg.p.data, distrib.map(Long, Lat, recent.log3.rsr))
with(ldg.p.data[recent.log3.rsr < 0 | recent.log3.rsr > 30, ], distrib.map(Long, Lat, recent.log3.rsr))
with(ldg.p.data[recent.log3.rsr > 0 & recent.log3.rsr < 30, ], distrib.map(Long, Lat, recent.log3.rsr[recent.log3.rsr > 0 & recent.log3.rsr < 30]))

# predict

ALE.prox.env, olddata = ldg.m.data)[, 1]
                                
summary(LE.prox.env$pred.rsr)

with(LE.prox.env, points(Lat, pred.rsr, pch = 16, col = 6))

UE.prox.env$pred.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = UE.prox.env, olddata = ldg.m.data)[, 1]

summary(UE.prox.env$pred.rsr)

with(UE.prox.env, points(Lat, pred.rsr, pch = 16, col = 8))

with(LE.prox.env, plot(pred.rsr, SR, xlim = c(0, 40), pch = 16))
with(UE.prox.env, points(pred.rsr, SR, col = 4, pch = 16))
legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene"), pch = 16, col = c(1, 2, 4))


# 7. Matching quantiles to proxy ------------------------------------------

# lower eocene
rq.LE <- rq(LE.SR ~ poly(abs(LE.lat),3), tau = 0.75)
rq.LE <- predict(rq.LE, newdata = data.frame(LE.lat = 0:70))

plot(abs(LE.lat), LE.SR, pch = 16, ylim = c(0, 35))
points(0:70, rq.LE, type = "l", lwd = 2)

with(temp.bi.early, points(-70:80, predict(gam(e.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 6, lwd = 2))

for(i in seq(0.01, 1, by = 0.01)) {
  tmp <- temp.bi.early
  tmp$meanSST.1deg <- i*tmp$meanSST.1deg
  tmp$tmp.pred <- predict(temp.mod, newdata = tmp)
  with(tmp, points(-70:80, predict(gam(tmp.pred ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = rainbow(120)[i*100]))
}

# adding a constant
plot(abs(LE.lat), LE.SR, pch = 16, ylim = c(0, 35))
points(0:70, rq.LE, type = "l", lwd = 2)

with(temp.bi.early, points(-70:80, predict(gam(e.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 6, lwd = 2))

for(i in 0:15) {
  tmp <- temp.bi.early
  tmp$meanSST.1deg <- tmp$meanSST.1deg - i
  tmp$tmp.pred <- predict(temp.mod, newdata = tmp)
  with(tmp, points(-70:80, predict(gam(tmp.pred ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = rainbow(20)[i+1]))
}

# upper eocene
rq.UE <- rq(UE.SR ~ poly(abs(UE.lat),3), tau = 0.75)
rq.UE <- predict(rq.UE, newdata = data.frame(UE.lat = 0:70))

plot(abs(UE.lat), UE.SR, pch = 16, ylim = c(0, 35))
points(0:70, rq.UE, type = "l", lwd = 2, col = 4)

with(temp.bi.late, points(-70:80, predict(gam(u.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 7, lwd = 2))

for(i in seq(0.01, 1, by = 0.01)) {
  tmp <- temp.bi.late
  tmp$meanSST.1deg <- i*tmp$meanSST.1deg
  tmp$tmp.pred <- predict(temp.mod, newdata = tmp)
  with(tmp, points(-70:80, predict(gam(tmp.pred ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = rainbow(120)[i*100], lwd = 2))
}

# adding a constant
with(temp.bi.late, plot(-70:80, predict(gam(u.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = 7, lwd = 2, ylab = "SR", xlab = "Lat"))
points(0:70, pred.UE.75, type = "l", lwd = 2, col = 4)
points(c(-70:0, 0:70), c(rev(pred.UE.75), pred.UE.75), type = "l", lwd = 2, col = 4)
                              
for(i in 0:15) {
  tmp <- temp.bi.late
  tmp$meanSST.1deg <- tmp$meanSST.1deg - i
  tmp$tmp.pred <- predict(temp.mod, newdata = tmp)
  with(tmp, points(-70:80, predict(gam(tmp.pred ~ s(Lat)), data.frame(Lat = -70:80)), type = "l", col = rainbow(20)[i+1]))
}


# 8. GAMs -----------------------------------------------------------------

# 8i. Richness ------------------------------------------------------------
reoc.rsr.gam <- with(rsr.margo.mod, gam(rarefy.sr ~ s(Longitude, Latitude, k = 80) + s(meanSST.1deg) + sdSST.1deg + mean.mld.t + depth10deg + absMnSal.0m + sdSal.0m + delta_carb_ion, gamma = 1.4))

# check predictions with recent data
recent.rsr <- predict.gam(reoc.rsr.gam, newdata = ldg.p.margo)

with(ldg.p.margo, distrib.map(Longitude, Latitude, recent.rsr))

ypr.pred <- ypr.env
names(ypr.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
ypr.pred$absMnSal.0m <- abs(ypr.pred$absMnSal.0m - 35.1)
ypr.pred$delta_carb_ion <- 0
ypr.pred <- na.omit(ypr.pred)
ypr.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = ypr.pred)
summary(ypr.pred$rarefy.sr)
with(ypr.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Latitude, rarefy.sr))

bar.pred <- bar.env
names(bar.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
bar.pred$absMnSal.0m <- abs(bar.pred$absMnSal.0m - 35.1)
bar.pred$delta_carb_ion <- 0
bar.pred <- na.omit(bar.pred)
bar.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = bar.pred)
summary(bar.pred$rarefy.sr)
with(bar.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(bar.pred, plot(Latitude, rarefy.sr))

pri.pred <- pri.env
names(pri.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
pri.pred$absMnSal.0m <- abs(pri.pred$absMnSal.0m - 35.1)
pri.pred$delta_carb_ion <- 0
pri.pred <- na.omit(pri.pred)
pri.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = pri.pred)
summary(pri.pred$rarefy.sr)
with(pri.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(pri.pred, plot(Latitude, rarefy.sr))

lut.pred <- lut.env
names(lut.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
lut.pred$absMnSal.0m <- abs(lut.pred$absMnSal.0m - 35.1)
lut.pred$delta_carb_ion <- 0
lut.pred <- na.omit(lut.pred)
lut.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = lut.pred)
summary(lut.pred$rarefy.sr)
with(lut.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(lut.pred, plot(Latitude, rarefy.sr))

eEoc.pred <- eEoc.env
names(eEoc.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
eEoc.pred$absMnSal.0m <- abs(eEoc.pred$absMnSal.0m - 35.1)
eEoc.pred$delta_carb_ion <- 0
eEoc.pred <- na.omit(eEoc.pred)
eEoc.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = eEoc.pred)
summary(eEoc.pred$rarefy.sr)
with(eEoc.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(eEoc.pred, plot(Latitude, rarefy.sr))

mEoc.pred <- mEoc.env
names(mEoc.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
mEoc.pred$absMnSal.0m <- abs(mEoc.pred$absMnSal.0m - 35.1)
mEoc.pred$delta_carb_ion <- 0
mEoc.pred <- na.omit(mEoc.pred)
mEoc.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = mEoc.pred)
summary(mEoc.pred$rarefy.sr)
with(mEoc.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(mEoc.pred, plot(Latitude, rarefy.sr))

lEoc.pred <- lEoc.env
names(lEoc.pred) <- c("Longitude", "Latitude", "meanSST.1deg", "mean.mld.t", "depth10deg", "sdSST.1deg", "sdSal.0m", "absMnSal.0m")
lEoc.pred$absMnSal.0m <- abs(lEoc.pred$absMnSal.0m - 35.1)
lEoc.pred$delta_carb_ion <- 0
lEoc.pred <- na.omit(lEoc.pred)
lEoc.pred$rarefy.sr <- predict.gam(reoc.rsr.gam, newdata = lEoc.pred)
summary(lEoc.pred$rarefy.sr)
with(lEoc.pred, distrib.map(Longitude, Latitude, rarefy.sr, pch = 15, col.land = "steelblue2"))

with(lEoc.pred, plot(Latitude, rarefy.sr))

eEoc.pred$temp_sr <- predict(temp.mod, newdata = eEoc.pred)
mEoc.pred$temp_sr <- predict(temp.mod, newdata = mEoc.pred)
lEoc.pred$temp_sr <- predict(temp.mod, newdata = lEoc.pred)

# 8ii. Evenness -----------------------------------------------------------
reoc.eve.gam <- with(eve.margo.mod, gam(simpsonEve ~ s(Longitude, Latitude, k = 80) + s(meanSST.1deg) + sdSST.1deg + mean.mld.t + depth10deg + absMnSal.0m + sdSal.0m + delta_carb_ion, gamma = 1.4))

# check predictions with recent data
recent.eve <- predict.gam(reoc.eve.gam, newdata = ldg.p.margo)

with(ldg.p.margo, distrib.map(Longitude, Latitude, recent.eve))

ypr.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = ypr.pred)
summary(ypr.pred$simpsonEve)
with(ypr.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Latitude, simpsonEve))

bar.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = bar.pred)
summary(bar.pred$simpsonEve)
with(bar.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(bar.pred, plot(Latitude, simpsonEve))

pri.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = pri.pred)
summary(pri.pred$simpsonEve)
with(pri.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(pri.pred, plot(Latitude, simpsonEve))

lut.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = lut.pred)
summary(lut.pred$simpsonEve)
with(lut.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(lut.pred, plot(Latitude, simpsonEve))

eEoc.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = eEoc.pred)
summary(eEoc.pred$simpsonEve)
with(eEoc.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(eEoc.pred, plot(Latitude, simpsonEve))

mEoc.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = mEoc.pred)
summary(mEoc.pred$simpsonEve)
with(mEoc.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(mEoc.pred, plot(Latitude, simpsonEve))

lEoc.pred$simpsonEve <- predict.gam(reoc.eve.gam, newdata = lEoc.pred)
summary(lEoc.pred$simpsonEve)
with(lEoc.pred, distrib.map(Longitude, Latitude, simpsonEve, pch = 15, col.land = "steelblue2"))

with(lEoc.pred, plot(Latitude, simpsonEve))

# 8iii. Lineage age -------------------------------------------------------
reoc.lna.gam <- with(eve.margo.mod, gam(LinAgeAbun ~ s(Longitude, Latitude, k = 80) + s(meanSST.1deg) + sdSST.1deg + mean.mld.t + depth10deg + absMnSal.0m + sdSal.0m + delta_carb_ion, gamma = 1.4))

# check predictions with recent data
recent.lna <- predict.gam(reoc.lna.gam, newdata = ldg.p.margo)

with(ldg.p.margo, distrib.map(Longitude, Latitude, recent.lna))

ypr.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = ypr.pred)
summary(ypr.pred$lin.age)
with(ypr.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Latitude, lin.age))

bar.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = bar.pred)
summary(bar.pred$lin.age)
with(bar.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(bar.pred, plot(Latitude, lin.age))

pri.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = pri.pred)
summary(pri.pred$lin.age)
with(pri.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(pri.pred, plot(Latitude, lin.age))

lut.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = lut.pred)
summary(lut.pred$lin.age)
with(lut.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(lut.pred, plot(Latitude, lin.age))

eEoc.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = eEoc.pred)
summary(eEoc.pred$lin.age)
with(eEoc.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(eEoc.pred, plot(Latitude, lin.age))

mEoc.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = mEoc.pred)
summary(mEoc.pred$lin.age)
with(mEoc.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(mEoc.pred, plot(Latitude, lin.age))

lEoc.pred$lin.age <- predict.gam(reoc.lna.gam, newdata = lEoc.pred)
summary(lEoc.pred$lin.age)
with(lEoc.pred, distrib.map(Longitude, Latitude, lin.age, pch = 15, col.land = "steelblue2"))

with(lEoc.pred, plot(Latitude, lin.age))

# 8iv. Functional richness ------------------------------------------------
reoc.FRic.gam <- with(fric.margo.mod, gam(FRic ~ s(Longitude, Latitude, k = 80) + s(meanSST.1deg) + sdSST.1deg + mean.mld.t + depth10deg + absMnSal.0m + sdSal.0m + delta_carb_ion, gamma = 1.4))

# check predictions with recent data
recent.FRic <- predict.gam(reoc.FRic.gam, newdata = ldg.p.margo)

with(ldg.p.margo, distrib.map(Longitude, Latitude, recent.FRic))

ypr.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = ypr.pred)
summary(ypr.pred$FRic)
with(ypr.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(ypr.pred, plot(Latitude, FRic))

bar.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = bar.pred)
summary(bar.pred$FRic)
with(bar.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(bar.pred, plot(Latitude, FRic))

pri.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = pri.pred)
summary(pri.pred$FRic)
with(pri.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(pri.pred, plot(Latitude, FRic))

lut.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = lut.pred)
summary(lut.pred$FRic)
with(lut.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(lut.pred, plot(Latitude, FRic))

eEoc.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = eEoc.pred)
summary(eEoc.pred$FRic)
with(eEoc.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(eEoc.pred, plot(Latitude, FRic))

mEoc.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = mEoc.pred)
summary(mEoc.pred$FRic)
with(mEoc.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(mEoc.pred, plot(Latitude, FRic))

lEoc.pred$FRic <- predict.gam(reoc.FRic.gam, newdata = lEoc.pred)
summary(lEoc.pred$FRic)
with(lEoc.pred, distrib.map(Longitude, Latitude, FRic, pch = 15, col.land = "steelblue2"))

with(lEoc.pred, plot(Latitude, FRic))


# 9. GAM predictions for diversity locations ------------------------------
match.2.5deg <- function(site, coords) {
  # site - single coordinate for the site (i.e. either lat or long)
  # coords - the list of all possible coordinates
  # create a blank list to hold the values
  pos <- which(abs(site - coords) == min(abs(site - coords)))
  output <- coords[pos]
  return(output)
}

# create a column for the tectonic data
EocComb_age$pred.sr.tect <- NA
EocComb_age$pred.temp.tect <- NA
EocComb_age$pred.eve.tect <- NA
EocComb_age$pred.lna.tect <- NA
EocComb_age$pred.fric.tect <- NA

# calculate the predicted values for the tectonic data
# Early Eocene (Ypresian)
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(ypr.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(ypr.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Early Eocene") {
    if (length(which(ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i])) > 0) {
    EocComb_age$pred.sr.tect[i] <- ypr.pred$rarefy.sr[ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i]]
    EocComb_age$pred.temp.tect[i] <- ypr.pred$temp_sr[ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i]]
    EocComb_age$pred.eve.tect[i] <- ypr.pred$simpsonEve[ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i]]
    EocComb_age$pred.lna.tect[i] <- ypr.pred$lin.age[ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i]]
    EocComb_age$pred.fric.tect[i] <- ypr.pred$FRic[ypr.pred$Longitude == eoc.long[i] & ypr.pred$Latitude == eoc.lat[i]]
    } 
  }
}
  
# middle Eocene (bartonian)
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(bar.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(bar.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Middle Eocene") {
    if (length(which(bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i])) > 0) {
      EocComb_age$pred.sr.tect[i] <- bar.pred$rarefy.sr[bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.temp.tect[i] <- bar.pred$temp_sr[bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.eve.tect[i] <- bar.pred$simpsonEve[bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.lna.tect[i] <- bar.pred$lin.age[bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.fric.tect[i] <- bar.pred$FRic[bar.pred$Longitude == eoc.long[i] & bar.pred$Latitude == eoc.lat[i]]
    } 
  }
}

# Late Eocene (Priabonian)
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(pri.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(pri.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Late Eocene") {
    if (length(which(pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i])) > 0) {
      EocComb_age$pred.sr.tect[i] <- pri.pred$rarefy.sr[pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.temp.tect[i] <- pri.pred$temp_sr[pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.eve.tect[i] <- pri.pred$simpsonEve[pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.lna.tect[i] <- pri.pred$lin.age[pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.fric.tect[i] <- pri.pred$FRic[pri.pred$Longitude == eoc.long[i] & pri.pred$Latitude == eoc.lat[i]]
    } 
  }
}

## for the CO2 data
EocComb_age$pred.sr.co2 <- NA
EocComb_age$pred.temp.co2 <- NA
EocComb_age$pred.eve.co2 <- NA
EocComb_age$pred.lna.co2 <- NA
EocComb_age$pred.fric.co2 <- NA

# Early Eocene
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(eEoc.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(eEoc.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Early Eocene") {
    if (length(which(eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i])) > 0) {
      EocComb_age$pred.sr.co2[i] <- eEoc.pred$rarefy.sr[eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.temp.co2[i] <- eEoc.pred$temp_sr[eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.eve.co2[i] <- eEoc.pred$simpsonEve[eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.lna.co2[i] <- eEoc.pred$lin.age[eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.fric.co2[i] <- eEoc.pred$FRic[eEoc.pred$Longitude == eoc.long[i] & eEoc.pred$Latitude == eoc.lat[i]]
    } 
  }
}

# Middle Eocene
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(mEoc.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(mEoc.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Middle Eocene") {
    if (length(which(mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i])) > 0) {
      EocComb_age$pred.sr.co2[i] <- mEoc.pred$rarefy.sr[mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.temp.co2[i] <- mEoc.pred$temp_sr[mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.eve.co2[i] <- mEoc.pred$simpsonEve[mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.lna.co2[i] <- mEoc.pred$lin.age[mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.fric.co2[i] <- mEoc.pred$FRic[mEoc.pred$Longitude == eoc.long[i] & mEoc.pred$Latitude == eoc.lat[i]]
    } 
  }
}

# Late Eocene
eoc.long <- sapply(EocComb_age$pal.long, match.2.5deg, unique(lEoc.pred$Longitude))
eoc.lat <- sapply(EocComb_age$pal.lat, match.2.5deg, unique(lEoc.pred$Latitude))

for (i in 1:nrow(EocComb_age)) {
  if (EocComb_age$Period[i] == "Late Eocene") {
    if (length(which(lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i])) > 0) {
      EocComb_age$pred.sr.co2[i] <- lEoc.pred$rarefy.sr[lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.temp.co2[i] <- lEoc.pred$temp_sr[lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.eve.co2[i] <- lEoc.pred$simpsonEve[lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.lna.co2[i] <- lEoc.pred$lin.age[lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i]]
      EocComb_age$pred.fric.co2[i] <- lEoc.pred$FRic[lEoc.pred$Longitude == eoc.long[i] & lEoc.pred$Latitude == eoc.lat[i]]
    } 
  }
}

# testing goodness of fit
# mean SR
with(EocComb_age, plot(pred.sr.tect, mean.SR, pch = 16, col = as.factor(Period)))
abline(0,1)
with(EocComb_age, plot(pred.sr.co2, mean.SR, pch = 16, col = as.factor(Period)))
abline(0,1)

# max SR
with(EocComb_age, plot(pred.sr.tect, max.SR, pch = 16, col = as.factor(Period)))
abline(0,1)
with(EocComb_age, plot(pred.sr.co2, max.SR, pch = 16, col = as.factor(Period)))
abline(0,1)

# rarefied SR
with(EocComb_age, plot(pred.sr.tect, rarefy.sr, pch = 16, col = as.factor(Period)))
abline(0,1)
with(EocComb_age, plot(pred.sr.co2, rarefy.sr, pch = 16, col = as.factor(Period)))
abline(0,1)

# calculate RMSE
rmse <- function(obs, pred)
{
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

## for entire dataset
# max.SR
with(EocComb_age, rmse(pred.sr.tect, max.SR)) # 8.182827
with(EocComb_age, rmse(pred.sr.co2, max.SR)) # 7.973676
# mean.SR
with(EocComb_age, rmse(pred.sr.tect, mean.SR)) # 6.487511
with(EocComb_age, rmse(pred.sr.co2, mean.SR)) # 6.863412
# rarefied SR
with(EocComb_age, rmse(pred.sr.tect, rarefy.sr)) # 8.367944
with(EocComb_age, rmse(pred.sr.co2, rarefy.sr)) # 7.794836

## split by period
## Early Eocene
# max.SR
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.tect, max.SR)) # 9.192648
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.co2, max.SR)) # 8.754297
# mean.SR
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.tect, mean.SR)) # 7.039098
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.co2, mean.SR)) # 7.645259
# rarefied SR
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.tect, rarefy.sr)) # 11.38624
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.sr.co2, rarefy.sr)) # 10.10226

# Middle Eocene
# max.SR
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.tect, max.SR)) # 7.605059
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.co2, max.SR)) # 7.683994
# mean.SR
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.tect, mean.SR)) # 5.933478
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.co2, mean.SR)) # 6.723788
# rarefied SR
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.tect, rarefy.sr)) # 4.928517
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.sr.co2, rarefy.sr)) # 4.520846

# Late Eocene
# max.SR
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.tect, max.SR)) # 7.088453
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.co2, max.SR)) # 6.801072
# mean.SR
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.tect, mean.SR)) # 6.265335
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.co2, mean.SR)) # 5.406247
# rarefied SR
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.tect, rarefy.sr)) # 4.77973
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.sr.co2, rarefy.sr)) # 4.184962


## for temperature only model
# tectonic
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.temp.tect, mean.SR)) # 5.761854
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.temp.tect, mean.SR)) # 5.390061
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.temp.tect, mean.SR)) # 5.252173

# co2 model
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.temp.co2, mean.SR)) # 6.179765
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.temp.co2, mean.SR)) #  6.616324
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.temp.co2, mean.SR)) # 5.008741

# for evenness model
# tectonic
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.eve.tect, simpsonEve)) # 0.2631576
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.eve.tect, simpsonEve)) # 0.2572157
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.eve.tect, simpsonEve)) # 0.2203721

# co2 model
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.eve.co2, simpsonEve)) # 0.231757
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.eve.co2, simpsonEve)) #  0.2174826
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.eve.co2, simpsonEve)) # 0.2119109

# for lineage age model
# tectonic
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.lna.tect, lin.age)) # 7.386006
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.lna.tect, lin.age)) # 4.127436
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.lna.tect, lin.age)) # 3.414303

# co2 model
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.lna.co2, lin.age)) # 5.980174
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.lna.co2,  lin.age)) # 3.843624
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.lna.co2, lin.age)) # 5.220653

# for lineage age model
# tectonic
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.fric.tect, FRic)) # 0.4456317
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.fric.tect, FRic)) # 0.4982162
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.fric.tect, FRic)) # 0.5690583

# co2 model
with(EocComb_age[EocComb_age$Period == "Early Eocene",], rmse(pred.fric.co2, FRic)) # 0.575143
with(EocComb_age[EocComb_age$Period == "Middle Eocene",], rmse(pred.fric.co2, FRic)) # 0.5667388
with(EocComb_age[EocComb_age$Period == "Late Eocene",], rmse(pred.fric.co2, FRic)) # 0.5704184


# 10. Save out -------------------------------------------------------------
save(ypr.pred, lut.pred, bar.pred, pri.pred, temp.bi.early, temp.bi.late, eEoc.pred, mEoc.pred, lEoc.pred, eEoc.env, mEoc.env, lEoc.env, y.tmp.rsr.ols, l.tmp.rsr.ols, b.tmp.rsr.ols, p.tmp.rsr.ols, e.tmp.rsr.ols, u.tmp.rsr.ols, eEoc.tmp.rsr.ols, mEoc.tmp.rsr.ols, lEoc.tmp.rsr.ols, file = "../../../../../Work/1404 Eocene diversity/Output/eoc_pred.RData")
save(temp.mod, file = "Output/temp_mod.RData")
