# modelling

library(spdep)
library(ncf)
source("C:/Documents/Science/PhD/Work/1311 LDGPaper/Code/140420SARerrOptimising_NC.R")
source("C:/Documents/Science/PhD/Code/sar_predict.R")

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
temp.mod <- lm(rarefy.sr ~ poly(meanSST.1deg, 3), data = ldg.m.data)
summary(temp.mod)

# look for spatial autocorrelation in the residuals
# using spline.correlog
temp.mod.sac <- with(ldg.m.data, spline.correlog(Long, Lat, temp.mod$residuals, latlon = TRUE, resamp = 10))
summary(temp.mod.sac)

# run sar model
ldg.coords <- cbind(ldg.m.data$Long,ldg.m.data$Lat)
ldg.coords <- as.matrix(ldg.coords)
temp.mod.sarW <- with(ldg.m.data, sar.optimised(temp.mod.sac$real$x.intercept, rarefy.sr ~ poly(meanSST.1deg, 3), ldg.coords, style = "W", tol = 4, longlat = TRUE, zero.policy = TRUE))
temp.mod.nb <- dnearneigh(ldg.coords, 0, temp.mod.sarW$dist, longlat = TRUE)
temp.mod.w <- nb2listw(temp.mod.nb, glist = NULL, style = "W", zero.policy = TRUE)
temp.mod.op0 <- errorsarlm(temp.mod.sarW$mod, listw = temp.mod.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(temp.mod.op0, Nagelkerke = TRUE) # 0.78294 
AIC(temp.mod.op0) # 2791.414

# check predictions with recent data
r.tmp.rsr <- sar.predict(temp.mod.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1]
summary(r.tmp.rsr)
with(ldg.p.data, distrib.map(Long, Lat, r.tmp.rsr))
# these look reasonable

y.tmp.rsr <- sar.predict(temp.mod.op0, newdata = ypr.pred, olddata = ldg.m.data)[, 1]
summary(y.tmp.rsr)
with(ypr.pred, distrib.map(Long, Lat, y.tmp.rsr, pch = 15, col.land = "steelblue2"))

l.tmp.rsr <- sar.predict(temp.mod.op0, newdata = lut.pred, olddata = ldg.m.data)[, 1]
summary(l.tmp.rsr)
with(lut.pred, distrib.map(Long, Lat, l.tmp.rsr, pch = 15, col.land = "steelblue2"))

b.tmp.rsr <- sar.predict(temp.mod.op0, newdata = bar.pred, olddata = ldg.m.data)[, 1]
summary(b.tmp.rsr)
with(bar.pred, distrib.map(Long, Lat, b.tmp.rsr, pch = 15, col.land = "steelblue2"))

p.tmp.rsr <- sar.predict(temp.mod.op0, newdata = pri.pred, olddata = ldg.m.data)[, 1]
summary(p.tmp.rsr)
with(pri.pred, distrib.map(Long, Lat, p.tmp.rsr, pch = 15, col.land = "steelblue2"))

with(ldg.p.data, plot(Lat, r.tmp.rsr, pch = ".", col = 5, ylim = c(0, 20)))
with(ypr.pred, points(Lat, y.tmp.rsr, pch = "."))
with(lut.pred, points(Lat, l.tmp.rsr, pch = ".", col = 2))
with(bar.pred, points(Lat, b.tmp.rsr, pch = ".", col = 3))
with(pri.pred, points(Lat, p.tmp.rsr, pch = ".", col = 4))

with(ldg.p.data, points(-90:90, predict(gam(r.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5, lwd = 2))
with(ypr.pred, points(-90:90, predict(gam(y.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", lwd = 2))
with(lut.pred, points(-90:90, predict(gam(l.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2, lwd = 2))
with(bar.pred, points(-90:90, predict(gam(b.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3, lwd = 2))
with(pri.pred, points(-90:90, predict(gam(p.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4, lwd = 2))

legend("topright", c("Ypresian", "Lutetian", "Bartonian", "Priabonian", "Recent"), col = 1:5, lwd = 2)

with(ldg.p.data, plot(Lat, r.tmp.rsr, pch = ".", col = 5, ylim = c(0, 20)))
with(ypr.pred, points(Lat, y.tmp.rsr, pch = "."))
with(lut.pred, points(Lat, l.tmp.rsr, pch = ".", col = 2))
with(pri.pred, points(Lat, p.tmp.rsr, pch = ".", col = 4))

with(ldg.p.data, points(-60:75, predict(gam(r.tmp.rsr ~ s(Lat)), data.frame(Lat = -60:75)), type = "l", col = 5, lwd = 2))
with(ypr.pred, points(-90:90, predict(gam(y.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", lwd = 2))
with(lut.pred, points(-90:90, predict(gam(l.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2, lwd = 2))
with(pri.pred, points(-90:90, predict(gam(p.tmp.rsr ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4, lwd = 2))

legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene", "Recent"), col = c(1:2,4:5), lwd = 2)




# compare these predictions with the ols model
r.tmp.rsr.ols <- predict(temp.mod, newdata = ldg.p.data)
summary(r.tmp.rsr.ols)
with(ldg.p.data, distrib.map(Long, Lat, r.tmp.rsr.ols))

y.tmp.rsr.ols <- predict(temp.mod, newdata = ypr.pred)
summary(y.tmp.rsr.ols)
with(ypr.pred, distrib.map(Long, Lat, y.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

l.tmp.rsr.ols <- predict(temp.mod, newdata = lut.pred)
summary(l.tmp.rsr.ols)
with(lut.pred, distrib.map(Long, Lat, l.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

b.tmp.rsr.ols <- predict(temp.mod, newdata = bar.pred)
summary(b.tmp.rsr.ols)
with(bar.pred, distrib.map(Long, Lat, b.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

p.tmp.rsr.ols <- predict(temp.mod, newdata = pri.pred)
summary(p.tmp.rsr.ols)
with(pri.pred, distrib.map(Long, Lat, p.tmp.rsr.ols, pch = 15, col.land = "steelblue2"))

with(ldg.p.data, plot(Lat, r.tmp.rsr.ols, pch = ".", col = 5, ylim = c(0, 20), ylab = "predicted richness"))
with(ypr.pred, points(Lat, y.tmp.rsr.ols, pch = "."))
with(lut.pred, points(Lat, l.tmp.rsr.ols, pch = ".", col = 2))
with(bar.pred, points(Lat, b.tmp.rsr.ols, pch = ".", col = 3))
with(pri.pred, points(Lat, p.tmp.rsr.ols, pch = ".", col = 4))

with(ldg.p.data, points(-90:90, predict(gam(r.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 5))
with(ypr.pred, points(-90:90, predict(gam(y.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l"))
with(lut.pred, points(-90:90, predict(gam(l.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.pred, points(-90:90, predict(gam(b.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.pred, points(-90:90, predict(gam(p.tmp.rsr.ols ~ s(Lat)), data.frame(Lat = -90:90)), type = "l", col = 4))

# basically the model predicts diversity drops rapidly with higher T
plot(-5:50, predict(temp.mod, newdata = data.frame(meanSST.1deg = -5:50)))


# 5. Andy's model ---------------------------------------------------------
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

reoc.rsr.log3.op0 <- errorsarlm(update(reoc.rsr.log3.sarW$mod, ~, listw = reoc.log3.op.w, zero.policy = TRUE, tol.solve = 1e-18)
summary(reoc.rsr.log3.op0, Nagelkerke = TRUE) # 0.89483



AIC(reoc.rsr.log3.op0) # -697.957

# check predictions with recent data
recent.log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = ldg.p.data, olddata = ldg.m.data)[, 1]
with(ldg.p.data, distrib.map(Long, Lat, recent.log3.rsr))
with(ldg.p.data[recent.log3.rsr < 0 | recent.log3.rsr > 30, ], distrib.map(Long, Lat, recent.log3.rsr))
with(ldg.p.data[recent.log3.rsr > 0 & recent.log3.rsr < 30, ], distrib.map(Long, Lat, recent.log3.rsr[recent.log3.rsr > 0 & recent.log3.rsr < 30]))

# predict
ypr.pred$log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = ypr.pred, olddata = ldg.m.data)[, 1]
summary(ypr.pred$log3.rsr)
with(ypr.pred[ypr.pred$meanSal.0m > 30 & ypr.pred$log3.rsr > 0 & ypr.pred$log3.rsr < 100, ], distrib.map(Long, Lat, log3.rsr, pch = 15, col.land = "steelblue2"))

lut.pred$log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = lut.pred, olddata = ldg.m.data)[, 1]
summary(lut.pred$log3.rsr)
with(lut.pred, distrib.map(Long, Lat, log3.rsr, pch = 15, col.land = "steelblue2"))

bar.pred$log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = bar.pred, olddata = ldg.m.data)[, 1]
summary(bar.pred$log3.rsr)
with(bar.pred, distrib.map(Long, Lat, log3.rsr, pch = 15, col.land = "steelblue2"))

pri.pred$log3.rsr <- sar.predict(reoc.rsr.log3.op0, newdata = pri.pred, olddata = ldg.m.data)[, 1]
summary(pri.pred$log3.rsr)
with(pri.pred, distrib.map(Long, Lat, log3.rsr, pch = 15, col.land = "steelblue2"))

with(ldg.p.data[recent.log3.rsr > 0 & recent.log3.rsr < 30, ], plot(Lat, recent.log3.rsr[recent.log3.rsr > 0 & recent.log3.rsr < 30], pch = ".", col = 5))
#with(ldg.m.data, points(Lat, rarefy.sr, pch = 16, col = 6))
with(ypr.pred[ypr.pred$meanSal.0m > 30 & ypr.pred$log3.rsr > 0 & ypr.pred$log3.rsr < 60, ], points(Lat, log3.rsr, pch = "."))
with(lut.pred[lut.pred$meanSal.0m > 30 & lut.pred$log3.rsr > 0 & lut.pred$log3.rsr < 60, ], points(Lat, log3.rsr, pch = ".", col = 2))
with(bar.pred[bar.pred$meanSal.0m > 30 & bar.pred$log3.rsr > 0 & bar.pred$log3.rsr < 60, ], points(Lat, log3.rsr, pch = ".", col = 3))
with(pri.pred[pri.pred$meanSal.0m > 30 & pri.pred$log3.rsr > 0 & pri.pred$log3.rsr < 60, ], points(Lat, log3.rsr, pch = ".", col = 4))

with(ldg.p.data[recent.log3.rsr > 0 & recent.log3.rsr < 30, ], points(-90:90, predict(gam(recent.log3.rsr[recent.log3.rsr > 0 & recent.log3.rsr < 30] ~ s(Lat, k = 20)), data.frame(Lat = -90:90)), type = "l", col = 5))
with(ypr.pred[ypr.pred$meanSal.0m > 30 & ypr.pred$log3.rsr > 0 & ypr.pred$log3.rsr < 60, ], points(-90:90, predict(gam(log3.rsr ~ s(Lat, k = 20)), data.frame(Lat = -90:90)), type = "l"))
with(lut.pred[lut.pred$meanSal.0m > 30 & lut.pred$log3.rsr > 0 & lut.pred$log3.rsr < 60, ], points(-90:90, predict(gam(log3.rsr ~ s(Lat, k = 20)), data.frame(Lat = -90:90)), type = "l", col = 2))
with(bar.pred[bar.pred$meanSal.0m > 30 & bar.pred$log3.rsr > 0 & bar.pred$log3.rsr < 60, ], points(-90:90, predict(gam(log3.rsr ~ s(Lat, k = 20)), data.frame(Lat = -90:90)), type = "l", col = 3))
with(pri.pred[pri.pred$meanSal.0m > 30 & pri.pred$log3.rsr > 0 & pri.pred$log3.rsr < 60, ], points(-90:90, predict(gam(log3.rsr ~ s(Lat, k = 20)), data.frame(Lat = -90:90)), type = "l", col = 4))

