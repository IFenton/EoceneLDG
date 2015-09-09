## Phil Trans Paper
## Date created: 6 / 9 / 2015
## Last modified: 9 / 9 / 2015
##
## Code used in the Phil Trans paper to produce the analyses.
##
## Project: 1404 Eocene diversity
## Previous file: eocene_div.R (Eocene.RProj)
## Next file:
##

## Inputs: -----------------------------------------------------------------
## EocComb_Abun, EocComb_age, EocComb_age.mn

## Outputs: ----------------------------------------------------------------

## Libraries ---------------------------------------------------------------
load("../../Project/Eocene/Outputs/EocComb.RData")
load("../1408 Diversity measures/Outputs/pca.RData")


# 1. Diversity graphs -----------------------------------------------------
# N.b. some of the data in the EocComb_age is NA, as these sites don't have abundance data. 

## 1i. Species richness ----------------------------------------------------
png("Figures/1_SR_Recent.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(ldg.margo.mod, plot(abs(Latitude), sp.rich, pch = 16, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", ylim = c(0, 32), cex.main = 2))
lm.Rec <- with(ldg.margo.mod, lm(sp.rich ~ poly(abs(Latitude), 2)))
p.Rec <- predict(lm.Rec, newdata = data.frame(Latitude = 0:70))
points(0:70, p.Rec, type = "l", lwd = 2)
dev.off()

lm.EE <- with(EocComb_age[EocComb_age$Period == "Early Eocene", ], lm(max.SR ~ poly(abs(pal.lat), 2)))
lm.ME <- with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], lm(max.SR ~ poly(abs(pal.lat), 2)))
lm.LE <- with(EocComb_age[EocComb_age$Period == "Late Eocene", ], lm(max.SR ~ poly(abs(pal.lat), 2)))
p.EE <- predict(lm.EE, newdata = data.frame(pal.lat = 0:65))
p.ME <- predict(lm.ME, newdata = data.frame(pal.lat = 0:65))
p.LE <- predict(lm.LE, newdata = data.frame(pal.lat = 0:65))

png("Figures/1_SR_Eocene.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), max.SR, pch = 15, ylab = "Species richness", ylim = c(0, 40), bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2))
with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], points(abs(pal.lat), max.SR, pch = 16, col = "grey50"))
with(EocComb_age[EocComb_age$Period == "Late Eocene", ], points(abs(pal.lat), max.SR, pch = 17, col = "grey80"))
legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene", "Recent"), col = c(1, "grey50", "grey80", 1), lty = c(1,1,1,2), pch = c(15:17, NA))
points(0:65, p.EE, type = "l", lwd = 2)     
points(0:65, p.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, p.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, lty = 2 )
dev.off()

png("Figures/1_SR_Eocene2.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), max.SR, pch = 15, ylab = "Species richness", ylim = c(0, 30), bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2, type = "n"))
legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene", "Recent"), col = c(1, "grey50", "grey80", 1), lty = c(1,1,1,2), lwd = 2)
points(0:65, p.EE, type = "l", lwd = 2)     
points(0:65, p.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, p.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, lty = 2 )
dev.off()

## 1ii. Species evenness ---------------------------------------------------
png("Figures/1_eve_recent.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(ldg.margo.mod, plot(abs(Latitude), simpsonEve, pch = 16, ylab = "Simpsons Evenness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", ylim = c(0, 1), cex.main = 2))
lm.Rec <- with(ldg.margo.mod, lm(simpsonEve ~ abs(Latitude)))
p.Rec <- predict(lm.Rec, newdata = data.frame(Latitude = 0:70))
points(0:70, p.Rec, type = "l", lwd = 2)
dev.off()

lm.EE.eve <- with(EocComb_age[EocComb_age$Period == "Early Eocene" & EocComb_age$simpson != Inf & !is.na(EocComb_age$simpson), ], lm(simpson/max.SR ~ abs(pal.lat)))
lm.ME.eve <- with(EocComb_age[EocComb_age$Period == "Middle Eocene" & EocComb_age$simpson != Inf & !is.na(EocComb_age$simpson), ], lm(simpson/max.SR ~ abs(pal.lat)))
lm.LE.eve <- with(EocComb_age[EocComb_age$Period == "Late Eocene" & EocComb_age$simpson != Inf & !is.na(EocComb_age$simpson), ], lm(simpson/max.SR ~ abs(pal.lat)))
p.EE.eve <- predict(lm.EE.eve, newdata = data.frame(pal.lat = 0:65))
p.ME.eve <- predict(lm.ME.eve, newdata = data.frame(pal.lat = 0:65))
p.LE.eve <- predict(lm.LE.eve, newdata = data.frame(pal.lat = 0:65))

png("Figures/1_eve_eocene.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), simpson/max.SR, pch = 15, col = 1, ylab = "Simpsons Evenness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", ylim = c(0, 1), cex.main = 2))
with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], points(abs(pal.lat), simpson/max.SR, pch = 16, col = "grey50"))
with(EocComb_age[EocComb_age$Period == "Late Eocene", ], points(abs(pal.lat), simpson/max.SR, pch = 17, col = "grey80"))
points(0:65, p.EE.eve, type = "l", lwd = 2)     
points(0:65, p.ME.eve, type = "l", lwd = 2, col = "grey50")
points(0:65, p.LE.eve, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, lty = 2 )
dev.off()

png("Figures/1_eve_eocene2.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), simpson/max.SR, pch = 15, col = 1, ylab = "Simpsons Evenness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", ylim = c(0, 1), cex.main = 2, type = "n"))
points(0:65, p.EE.eve, type = "l", lwd = 2)     
points(0:65, p.ME.eve, type = "l", lwd = 2, col = "grey50")
points(0:65, p.LE.eve, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, lty = 2 )
dev.off()

# 1iii. Lineage age -------------------------------------------------------
png("Figures/1_lna_recent.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(ldg.margo.mod, plot(abs(Latitude), LinAgeAbun, pch = 16, ylab = "Average age", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2, ylim = c(0, 35)))
lm.Rec <- with(ldg.margo.mod, lm(LinAgeAbun ~ poly(abs(Latitude), 2)))
p.Rec <- predict(lm.Rec, newdata = data.frame(Latitude = 0:70))
points(0:70, p.Rec, type = "l", lwd = 2, col = 1)
dev.off()

lna.EE <- with(EocComb_age[EocComb_age$Period == "Early Eocene", ], lm(lin.age ~ abs(pal.lat)))
lna.ME <- with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], lm(lin.age ~ abs(pal.lat)))
lna.LE <- with(EocComb_age[EocComb_age$Period == "Late Eocene", ], lm(lin.age ~ abs(pal.lat)))
plna.EE <- predict(lna.EE, newdata = data.frame(pal.lat = 0:65))
plna.ME <- predict(lna.ME, newdata = data.frame(pal.lat = 0:65))
plna.LE <- predict(lna.LE, newdata = data.frame(pal.lat = 0:65))

png("Figures/1_lna_eocene.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), lin.age, pch = 15, col = 1, ylab = "Average age", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2, ylim = c(0, 25)))
with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], points(abs(pal.lat), lin.age, pch = 16, col = "grey50"))
with(EocComb_age[EocComb_age$Period == "Late Eocene", ], points(abs(pal.lat), lin.age, pch = 17, col = "grey80"))
points(0:65, plna.EE, type = "l", lwd = 2)     
points(0:65, plna.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, plna.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, col = 1, lty = 2)
dev.off()

png("Figures/1_lna_eocene2.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), lin.age, pch = 15, col = 1, ylab = "Average age", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2, ylim = c(0, 25), type = "n"))
points(0:65, plna.EE, type = "l", lwd = 2)     
points(0:65, plna.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, plna.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, col = 1, lty = 2)
dev.off()

## 1iv. Functional richness ------------------------------------------------
png("Figures/1_FRic_recent.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(ldg.margo.mod, plot(abs(Latitude), FRic, pch = 16, ylab = "Functional richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2))
lm.Rec <- with(ldg.margo.mod, lm(FRic ~ poly(abs(Latitude), 2)))
p.Rec <- predict(lm.Rec, newdata = data.frame(Latitude = 0:70))
points(0:70, p.Rec, type = "l", lwd = 2, col = 1)
dev.off()

fric.EE <- with(EocComb_age[EocComb_age$Period == "Early Eocene", ], lm(FRic ~ abs(pal.lat)))
fric.ME <- with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], lm(FRic ~ abs(pal.lat)))
fric.LE <- with(EocComb_age[EocComb_age$Period == "Late Eocene", ], lm(FRic ~ abs(pal.lat)))
pfric.EE <- predict(fric.EE, newdata = data.frame(pal.lat = 0:65))
pfric.ME <- predict(fric.ME, newdata = data.frame(pal.lat = 0:65))
pfric.LE <- predict(fric.LE, newdata = data.frame(pal.lat = 0:65))

png("Figures/1_FRic_eocene.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), FRic, pch = 15, col = 1, ylab = "Functional richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2))
with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], points(abs(pal.lat), FRic, pch = 16, col = "grey50"))
with(EocComb_age[EocComb_age$Period == "Late Eocene", ], points(abs(pal.lat), FRic, pch = 17, col = "grey80"))
points(0:65, pfric.EE, type = "l", lwd = 2)     
points(0:65, pfric.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, pfric.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, col = 1, lty = 2)
dev.off()

png("Figures/1_FRic_eocene2.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), FRic, pch = 15, col = 1, ylab = "Functional richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2, type = "n"))
points(0:65, pfric.EE, type = "l", lwd = 2)     
points(0:65, pfric.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, pfric.LE, type = "l", lwd = 2, col = "grey80")
points(0:70, p.Rec, type = "l", lwd = 2, col = 1, lty = 2)
dev.off()


# 2. PCA constructions ----------------------------------------------------

# 2i. PC1-PC2 -------------------------------------------------------------
png("Figures/2_pca.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(robPCA2$scores[, 1], robPCA2$scores[, 2], pch = 16, col = as.numeric(eoc.pca.margo.plot$eco) + 2, xlab = "PCA 1", ylab = "PCA2", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, ylim = c(-15, 10), xlim = c(-13, 7))
legend("bottomleft", legend = levels(eoc.pca.margo.plot$eco), pch = 16, col = 3:7)
par(new = TRUE)
plot(robPCA2$loadings[, 1], robPCA2$loadings[, 2], pch = 16, col = pca.margo.plot$Ocean2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", xlim = c(-1.04, 0.56), ylim = c(-1.2, 0.8), bty = "l")
for (i in 1:nrow(robPCA2$loadings)){
  points(c(0, robPCA2$loadings[i, 1]), c(0, robPCA2$loadings[i, 2]), type = "l", lwd = 2)
  text(robPCA2$loadings[i, 1] - .12, robPCA2$loadings[i, 2] + .05, rownames(robPCA2$loadings)[i], font = 2)
}
dev.off()

png("Figures/2_PCA_eoc.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(robPCA3$scores[, 1], robPCA3$scores[, 2], pch = c(15, 17, 16)[as.numeric(as.factor(eoc.pca.col$Period))], xlab = "PCA 1", ylab = "PCA2", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, col = as.numeric(eoc.pca.col$eco) + 2, xlim = c(15, -15), ylim = c(4, -4))
legend("topright", legend = c("Early Eocene", "Middle Eocene", "Late Eocene"), pch = 15:17)
par(new = TRUE)
plot(robPCA3$loadings[, 1], robPCA3$loadings[, 2], pch = 16, col = pca.margo.plot$Ocean2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", xlim = c(3, -3), ylim = c(0.8, -0.8), bty = "l")
for (i in 1:nrow(robPCA3$loadings)){
  points(c(0, robPCA3$loadings[i, 1]), c(0, robPCA3$loadings[i, 2]), type = "l", lwd = 2)
  text(robPCA3$loadings[i, 1] - 0.35, robPCA3$loadings[i, 2] - 0.03, rownames(robPCA3$loadings)[i], font = 2)
}
dev.off()

# 2ii. PC1-PC3 -----------------------------------------------------------------
png("Figures/2_pca13.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(robPCA2$scores[, 1], robPCA2$scores[, 3], pch = 16, col = as.numeric(eoc.pca.margo.plot$eco) + 2, xlab = "PCA 1", ylab = "PCA3", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, ylim = c(-6, 8), xlim = c(-13, 7))
legend("bottomleft", legend = levels(eoc.pca.margo.plot$eco), pch = 16, col = 3:7)
par(new = TRUE)
plot(robPCA2$loadings[, 1], robPCA2$loadings[, 3], pch = 16, col = pca.margo.plot$Ocean2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", xlim = c(-1.04, 0.56), ylim = c(-0.48, 0.64), bty = "l")
for (i in 1:nrow(robPCA2$loadings)){
  points(c(0, robPCA2$loadings[i, 1]), c(0, robPCA2$loadings[i, 3]), type = "l", lwd = 2)
  text(robPCA2$loadings[i, 1] - .12, robPCA2$loadings[i, 3] - .025, rownames(robPCA2$loadings)[i], font = 2)
}
dev.off()

png("Figures/2_PCA_eoc13.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(robPCA3$scores[, 1], robPCA3$scores[, 3], pch = c(15, 17, 16)[as.numeric(as.factor(eoc.pca.col$Period))], xlab = "PCA 1", ylab = "PCA3", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, col = as.numeric(eoc.pca.col$eco) + 2, xlim = c(15, -15), ylim = c(12, -12))
legend("topright", legend = c("Early Eocene", "Middle Eocene", "Late Eocene"), pch = 15:17)
par(new = TRUE)
plot(robPCA3$loadings[, 1], robPCA3$loadings[, 3], pch = 16, col = pca.margo.plot$Ocean2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", xlim = c(.75, -.75), ylim = c(.6, -.6), bty = "l")
for (i in 1:nrow(robPCA3$loadings)){
  points(c(0, robPCA3$loadings[i, 1]), c(0, robPCA3$loadings[i, 3]), type = "l", lwd = 2)
  text(robPCA3$loadings[i, 1] - 0.1, robPCA3$loadings[i, 3] - 0.03, rownames(robPCA3$loadings)[i], font = 2)
}
dev.off()


# 3. GAM modelling --------------------------------------------------------
load("../../1404 Eocene diversity/Output/eoc_pred.RData")
setwd("../../../../../Work/1404 Eocene diversity/")
png("Figures/3_gam_eEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))

with(ypr.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "purple"))
with(eEoc.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "orange"))
with(ypr.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(eEoc.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.EE, type = "l", lwd = 2)   
dev.off()

png("Figures/3_gam_mEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))

with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))

with(bar.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "purple"))
with(mEoc.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "orange"))
with(bar.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(mEoc.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.ME, type = "l", lwd = 2) 
dev.off()

png("Figures/3_gam_lEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))

with(EocComb_age[EocComb_age$Period == "Late Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))
with(pri.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "purple"))
with(lEoc.pred, points(abs(Latitude), rarefy.sr, pch = ".", col = "orange"))
with(pri.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(lEoc.pred, points(abs(c(-90:90)), predict(gam(rarefy.sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.LE, type = "l", lwd = 2) 
dev.off()


# 4. Coccolithophores -----------------------------------------------------
load("../1408 Coccoliths/Outputs/Coccolithophores.RData")

lm.coc.EE <- with(age.eo.cocco[age.eo.cocco$age == "Early Eocene", ], lm(sp.rich ~ abs(pal.lat)))
lm.coc.ME <- with(age.eo.cocco[age.eo.cocco$age == "Middle Eocene", ], lm(sp.rich ~ abs(pal.lat)))
lm.coc.LE <- with(age.eo.cocco[age.eo.cocco$age == "Late Eocene", ], lm(sp.rich ~ poly(abs(pal.lat), 2)))
p.coc.EE <- predict(lm.coc.EE, newdata = data.frame(pal.lat = 0:65))
p.coc.ME <- predict(lm.coc.ME, newdata = data.frame(pal.lat = 0:65))
p.coc.LE <- predict(lm.coc.LE, newdata = data.frame(pal.lat = 0:65))

png("Figures/4_SR_Cocco.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(age.eo.cocco[age.eo.cocco$age == "Early Eocene", ], plot(abs(pal.lat), sp.rich, pch = 15, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(Latitude)", cex.main = 2))
with(age.eo.cocco[age.eo.cocco$age == "Middle Eocene", ], points(abs(pal.lat), sp.rich, pch = 16, col = "grey50"))
with(age.eo.cocco[age.eo.cocco$age == "Late Eocene", ], points(abs(pal.lat), sp.rich, pch = 17, col = "grey80"))
legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene"), col = c(1, "grey50", "grey80"), pch = c(15:17))
points(0:65, p.coc.EE, type = "l", lwd = 2)     
points(0:65, p.coc.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, p.coc.LE, type = "l", lwd = 2, col = "grey80")
dev.off()

png("Figures/4_SR_Cocco2.png", width = 400, height = 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(age.eo.cocco[age.eo.cocco$age == "Early Eocene", ], plot(abs(pal.lat), sp.rich, pch = 15, ylab = "Species richness", ylim = c(0, 30), bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "abs(pal.lat)", cex.main = 2, type = "n"))
legend("topright", c("Early Eocene", "Middle Eocene", "Late Eocene"), col = c(1, "grey50", "grey80"), lwd = 2)
points(0:65, p.coc.EE, type = "l", lwd = 2)     
points(0:65, p.coc.ME, type = "l", lwd = 2, col = "grey50")
points(0:65, p.coc.LE, type = "l", lwd = 2, col = "grey80")
dev.off()


# 5. Using only temperature -----------------------------------------------
png("Figures/5_temp_eEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))
with(EocComb_age[EocComb_age$Period == "Early Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))

with(ypr.pred, points(abs(Latitude), temp_sr, pch = ".", col = "purple"))
with(eEoc.pred, points(abs(Latitude), temp_sr, pch = ".", col = "orange"))
with(ypr.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(eEoc.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.EE, type = "l", lwd = 2)   
dev.off()

png("Figures/5_temp_mEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))

with(EocComb_age[EocComb_age$Period == "Middle Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))

with(lut.pred, points(abs(Latitude), temp_sr, pch = ".", col = "purple"))
with(mEoc.pred, points(abs(Latitude), temp_sr, pch = ".", col = "orange"))
with(lut.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(mEoc.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.ME, type = "l", lwd = 2) 
dev.off()

png("Figures/5_temp_lEoc.png", 400, 400)
par(mar = c(5.1, 5.1, 3.1, 2.1))

with(EocComb_age[EocComb_age$Period == "Late Eocene", ], plot(abs(pal.lat), max.SR, pch = 16, col = 1, ylab = "Species richness", bty = "l", las = 1, cex.axis = 1.2, cex.lab = 1.8, xlab = "Latitude", ylim = c(0, 32), cex.main = 2))

with(pri.pred, points(abs(Latitude), temp_sr, pch = ".", col = "purple"))
with(lEoc.pred, points(abs(Latitude), temp_sr, pch = ".", col = "orange"))
with(pri.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "purple"))
with(lEoc.pred, points(abs(c(-90:90)), predict(gam(temp_sr ~ s(Latitude)), data.frame(Latitude = -90:90)), type = "l", lwd = 2, col = "orange"))
points(0:65, p.LE, type = "l", lwd = 2) 
dev.off()


# 6. PF extinctions -------------------------------------------------------
load("../../Project/Foraminifera/Outputs/150318_pf_traits.RData")
head(pf.traits)

# create a dataset of the species that went extinct during the Eocene
eoc.ext.traits <- pf.traits[pf.traits$aL.end < 56 & pf.traits$aL.end > 33.8, ]
# add a column specifying which period they went extinct in
eoc.ext.traits$period <- NA
for (i in 4:6) {
  eoc.ext.traits$period[eoc.ext.traits$aL.end < dates$Start.Date[i] & eoc.ext.traits$aL.end > dates$End.Date[i]] <- as.character(dates$Epoch[i])
}
eoc.ext.traits$period <- factor(eoc.ext.traits$period, levels = c("Early Eocene", "Middle Eocene", "Late Eocene"))

# how many went extinct at different periods
table(eoc.ext.traits$period)

# how were these extinctions split by ecogroup
table(eoc.ext.traits$period, eoc.ext.traits$eco)

eEoc.traits <- pf.traits[pf.traits$aL.end < 56 & pf.traits$aL.start > 47.80, ]
mEoc.traits <- pf.traits[pf.traits$aL.end < 47.80 & pf.traits$aL.start > 38, ]
lEoc.traits <- pf.traits[pf.traits$aL.end < 38 & pf.traits$aL.start > 33.9, ]

# changes in eco group
table(eEoc.traits$eco)
table(mEoc.traits$eco)
table(lEoc.traits$eco)
table(bfd.traits$eco)
table(eoc.ext.traits$period, eoc.ext.traits$eco)
