## Proxies for productivity
## 14 / 10 / 2014
## Isabel Fenton

## Moore et al (2014) suggests a couple of possible proxies for productivity
## Check these in the recent before trying to use them in the past

# ldg.data
library(fields)
library(sp)

## 1. Benthic forams -------------------------------------------------------
# extract benthic foram column from bfd data
# need to remove extra rows
GCdeci <- function (x, y) ifelse(x >= 0 & y < 0, -as.numeric(paste(x, abs(round(y / 60 * 10000, 0)), sep= ".")), as.numeric(paste(x, abs(round(y / 60 * 10000, 0)), sep = ".")))

tmp <- which(point.in.polygon(GCdeci(ldg.bfd$Longd, ldg.bfd$Longm), GCdeci(ldg.bfd$Latd, ldg.bfd$Latm), world.dat$x, world.dat$y) != 0)
tmp <- c(tmp, which(is.na(ldg.bfd$Water.Depth)))

bfs <- ldg.bfd$Benthics[-tmp]

# plot that against chla
plot(bfs, ldg.env$meanChla.1deg)

# no pattern, but lots of reasons why not. How about plotting BF/PF ratio
plot(bfs / ldg.data$Total.Planktics, ldg.env$meanChla.1deg)

# try to separate out the accumulation at the bottom by logging
plot(log(bfs / ldg.data$Total.Planktics), log(ldg.env$meanChla.1deg))

# slight pattern, but not clear. Check this r2 of a model
summary(lm(log(ldg.env$meanChla.1deg) ~ log(bfs / ldg.data$Total.Planktics)))
# r2 = 0.109, so not a very good model

# check the global map
distrib.map(ldg.data$Long, ldg.data$Lat, log(bfs / ldg.data$Total.Planktics))

## 2. Radiolaria --------------------------------------------------------------
# extract radiolarian column from bfd data
# need to remove extra rows
rdl <- ldg.bfd$Radiolarians[-tmp]

# plot that against chla
plot(rdl, ldg.env$meanChla.1deg)

# no pattern, but lots of reasons why not. How about plotting BF/PF ratio
plot(rdl / ldg.data$Total.Planktics, ldg.env$meanChla.1deg)

# try to separate out the accumulation at the bottom by logging
plot(log(rdl / ldg.data$Total.Planktics), log(ldg.env$meanChla.1deg))

# less clear pattern. Check this r2 of a model
summary(lm(log(ldg.env$meanChla.1deg) ~ log(rdl / ldg.data$Total.Planktics)))
# r2 = 0.0129, so not really worth it

