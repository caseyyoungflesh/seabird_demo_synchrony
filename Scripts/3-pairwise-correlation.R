#################
# Youngflesh et al. - 3 - Pairwise correlations between species
#################

# Lack of synchronized breeding success in a seabird community: extreme events, niche separation, and environmental variability

# Casey Youngflesh, Yun Li, Heather J. Lynch, Karine Delord, Christophe Barbraud, Rubao Ji, Stephanie Jenouvrier


# Load packages -----------------------------------------------------------

library(dplyr)
library(purrr)


# Data process ------------------------------------------------------------------

setwd('Data')

data <- read.csv('data_pub.csv', header = TRUE)

#separate by species
ADPE <- dplyr::filter(data, SPECIES == 'ADPE')
SOFU <- dplyr::filter(data, SPECIES == 'SOFU')
CAPE <- dplyr::filter(data, SPECIES == 'CAPE')
SNPE <- dplyr::filter(data, SPECIES == 'SNPE')
SPSK <- dplyr::filter(data, SPECIES == 'SPSK')

#left_join all data.frames
mrg <- list(ADPE[,c('YEAR', 'BS')], SOFU[,c('YEAR', 'BS')], 
            CAPE[,c('YEAR', 'BS')], SNPE[,c('YEAR', 'BS')], 
            SPSK[,c('YEAR', 'BS')]) %>% 
  purrr::reduce(dplyr::full_join, 
                by = 'YEAR')
colnames(mrg) <- c('YEAR', 'ADPE', 'SOFU', 
                   'CAPE', 'SNPE', 'SPSK')

rna <- which(is.na(mrg), arr.ind = TRUE)[,1]
d_idx <- rna[which(!duplicated(rna))]

mrg2 <- mrg[-d_idx,]


#detrend - output residuals
ADPE_res <- scale(residuals(lm(mrg2$ADPE ~ mrg2$YEAR)))[,1]
SOFU_res <- scale(residuals(lm(mrg2$SOFU ~ mrg2$YEAR)))[,1]
CAPE_res <- scale(residuals(lm(mrg2$CAPE ~ mrg2$YEAR)))[,1]
SNPE_res <- scale(residuals(lm(mrg2$SNPE ~ mrg2$YEAR)))[,1]
SPSK_res <- scale(residuals(lm(mrg2$SPSK ~ mrg2$YEAR)))[,1]

#merge residuals
mrg3 <- data.frame(YEAR = mrg2$YEAR, ADPE_res, SOFU_res, 
                   CAPE_res, SNPE_res, SPSK_res)


# Calculate pairwise correlation and bootstrap CI  --------------------------------

#remove YEAR
mrg4 <- mrg3[,-1]
narr <- array(NA, dim = c(5, 5, 10000))
for (i in 1:10000)
{
  rs1 <- base::sample(1:NROW(mrg4), size = NROW(mrg4), replace = TRUE)
  tt <- mrg4[rs1,]
  narr[,,i] <- cor(tt)
}

corr_mn <- cor(mrg4)
corr_LCI <- apply(narr, c(1,2), function(x) quantile(x, probs = 0.025))
corr_UCI <- apply(narr, c(1,2), function(x) quantile(x, probs = 0.975))

#get species names
cn <- substr(colnames(corr_mn), 1, 4)
sp1 <- c(rep(cn[1], 4), rep(cn[2], 3), 
         rep(cn[3], 2), rep(cn[4], 1))
sp2 <- c(cn[2:5], cn[3:5], cn[4:5], cn[5])
#indices to extract corr values from mat
mat_idx <- which(upper.tri(corr_mn), arr.ind = TRUE)
mat_idx2 <- mat_idx[order(mat_idx[,1]),]
#fill df with correlations
corr_df <- data.frame(sp1 = sp1, sp2 = sp2, 
                      corr_mn = corr_mn[mat_idx2], 
                      corr_LCI = corr_LCI[mat_idx2], 
                      corr_UCI = corr_UCI[mat_idx2])
