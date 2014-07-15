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

# 2a. Generate a dataframe with each row as a site at a given age ----------
EO_abun <- 

# 2b. Calculate the number of species at each site for each age -----------

# 2c. Generate a dataframe with sites by period ---------------------------

# 2d. Calculate average species richness for each site --------------------

# 2e. Plot this up --------------------------------------------------------


# 3. Calculate evenness ---------------------------------------------------



# 4. Calculate lineage ages --------------------------------------------------


