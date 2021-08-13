#################
# Youngflesh et al. - 4 - Synchrony comparisons
#################

# Lack of synchronized breeding success in a seabird community: extreme events, niche separation, and environmental variability

# Casey Youngflesh, Yun Li, Heather J. Lynch, Karine Delord, Christophe Barbraud, Rubao Ji, Stephanie Jenouvrier


# load packages -----------------------------------------------------------

library(dplyr)
library(purrr)
library(ggplot2)
library(rstanarm)
library(MCMCvis)


# calculate synchrony for each population ---------------------------------

# Pointe Géologie
setwd('Data')

PG_bs <- read.csv('data_pub.csv', header = TRUE)

#only years where all species have BS data
ADPE <- dplyr::filter(PG_bs, SPECIES == 'ADPE')
SOFU <- dplyr::filter(PG_bs, SPECIES == 'SOFU')
CAPE <- dplyr::filter(PG_bs, SPECIES == 'CAPE')
SNPE <- dplyr::filter(PG_bs, SPECIES == 'SNPE')
SPSK <- dplyr::filter(PG_bs, SPECIES == 'SPSK')


#left_join all data.frames
PG_bs1 <- list(ADPE[,c('YEAR', 'BS')], SOFU[,c('YEAR', 'BS')], 
               CAPE[,c('YEAR', 'BS')], SNPE[,c('YEAR', 'BS')], 
               SPSK[,c('YEAR', 'BS')]) %>% 
  purrr::reduce(dplyr::full_join, 
                by = 'YEAR')
colnames(PG_bs1) <- c('YEAR', 'ADPE', 'SOFU', 'CAPE', 'SNPE', 'SPSK')

rna <- which(is.na(PG_bs1), arr.ind = TRUE)[,1]
d_idx <- rna[which(!duplicated(rna))]

PG_bs2 <- PG_bs1[-d_idx,]

ADPE_res <- scale(residuals(lm(ADPE ~ YEAR, data = PG_bs2)))[,1]
SOFU_res <- scale(residuals(lm(SOFU ~ YEAR, data = PG_bs2)))[,1]
CAPE_res <- scale(residuals(lm(CAPE ~ YEAR, data = PG_bs2)))[,1]
SNPE_res <- scale(residuals(lm(SNPE ~ YEAR, data = PG_bs2)))[,1]
SPSK_res <- scale(residuals(lm(SPSK ~ YEAR, data = PG_bs2)))[,1]

PG_bs3 <- data.frame(YEAR = PG_bs2$YEAR, ADPE_res, SOFU_res, 
                     CAPE_res, SNPE_res, SPSK_res)

#correlation matrix
PG_cc <- cor(PG_bs3[,-1])

#no extreme year
ey_idx <- which(PG_bs2$YEAR == 2013)
PG_bs2_ne <- PG_bs2[-ey_idx,]
ADPE_res_ne <- scale(residuals(lm(ADPE ~ YEAR, data = PG_bs2_ne)))[,1]
SOFU_res_ne <- scale(residuals(lm(SOFU ~ YEAR, data = PG_bs2_ne)))[,1]
CAPE_res_ne <- scale(residuals(lm(CAPE ~ YEAR, data = PG_bs2_ne)))[,1]
SNPE_res_ne <- scale(residuals(lm(SNPE ~ YEAR, data = PG_bs2_ne)))[,1]
SPSK_res_ne <- scale(residuals(lm(SPSK ~ YEAR, data = PG_bs2_ne)))[,1]

PG_bs3_ne <- data.frame(YEAR = PG_bs2_ne$YEAR, ADPE_res_ne, SOFU_res_ne, 
                        CAPE_res_ne, SNPE_res_ne, SPSK_res_ne)

#correlation matrix
PG_cc_ne <- cor(PG_bs3_ne[,-1])


# Isle of May
setwd('Data')
# source: Lahoz-Monfort et al. 2013
IM_bs <- read.csv('IM_data.csv')
IM_bs$BS <- IM_bs$CHICKS / IM_bs$ABUN

usp <- sort(unique(IM_bs$SPECIES))

ATPU <- dplyr::filter(IM_bs, SPECIES == usp[1])
BLKI <- dplyr::filter(IM_bs, SPECIES == usp[2])
COMU <- dplyr::filter(IM_bs, SPECIES == usp[3])
EUSH <- dplyr::filter(IM_bs, SPECIES == usp[4])
RAZO <- dplyr::filter(IM_bs, SPECIES == usp[5])

#detrend
ATPU_res <- scale(residuals(lm(BS ~ YEAR, data = ATPU)))[,1]
BLKI_res <- scale(residuals(lm(BS ~ YEAR, data = BLKI)))[,1]
COMU_res <- scale(residuals(lm(BS ~ YEAR, data = COMU)))[,1]
EUSH_res <- scale(residuals(lm(BS ~ YEAR, data = EUSH)))[,1]
RAZO_res <- scale(residuals(lm(BS ~ YEAR, data = RAZO)))[,1]

IM_bs1 <- data.frame(YEAR = 1:length(ATPU_res), ATPU_res, 
                     BLKI_res, COMU_res, EUSH_res, RAZO_res)

#correlation matrix
IM_cc <- cor(IM_bs1[,-1])


# Southeast Farallon Island
# source: https://data.prbo.org/cadc2/index.php?page=colony-data
SF_bs <- read.csv('SF_data.csv')
SF_bs$ABUN <- as.numeric(SF_bs$ABUN)

#only years where all species have BS data
usp <- unique(SF_bs$SPECIES)
ASSP <- dplyr::filter(SF_bs, SPECIES == usp[1])
BRCO <- dplyr::filter(SF_bs, SPECIES == usp[2])
CAAU <- dplyr::filter(SF_bs, SPECIES == usp[3])
COMU2 <- dplyr::filter(SF_bs, SPECIES == usp[4])
PECO <- dplyr::filter(SF_bs, SPECIES == usp[5])
PIGU <- dplyr::filter(SF_bs, SPECIES == usp[6])
RHAU <- dplyr::filter(SF_bs, SPECIES == usp[7])
WEGU <- dplyr::filter(SF_bs, SPECIES == usp[8])


#left_join all data.frames
SF_bs1 <- list(ASSP[,c('YEAR', 'BS')], BRCO[,c('YEAR', 'BS')], 
               CAAU[,c('YEAR', 'BS')], COMU2[,c('YEAR', 'BS')], 
               PECO[,c('YEAR', 'BS')], PIGU[,c('YEAR', 'BS')],
               RHAU[,c('YEAR', 'BS')], WEGU[,c('YEAR', 'BS')]) %>% 
  purrr::reduce(dplyr::full_join, 
                by = 'YEAR')
colnames(SF_bs1) <- c('YEAR', 'ASSP', 'BRCO', 'CAAU', 
                      'COMU2', 'PECO', 'PIGU', 'RHAU',
                      'WEGU')

rna <- which(is.na(SF_bs1), arr.ind = TRUE)[,1]
d_idx <- rna[which(!duplicated(rna))]

SF_bs2 <- SF_bs1[-d_idx,]

ASSP_res <- scale(residuals(lm(ASSP ~ YEAR, data = SF_bs2)))[,1]
BRCO_res <- scale(residuals(lm(BRCO ~ YEAR, data = SF_bs2)))[,1]
CAAU_res <- scale(residuals(lm(CAAU ~ YEAR, data = SF_bs2)))[,1]
COMU2_res <- scale(residuals(lm(COMU2 ~ YEAR, data = SF_bs2)))[,1]
PECO_res <- scale(residuals(lm(PECO ~ YEAR, data = SF_bs2)))[,1]
PIGU_res <- scale(residuals(lm(PIGU ~ YEAR, data = SF_bs2)))[,1]
RHAU_res <- scale(residuals(lm(RHAU ~ YEAR, data = SF_bs2)))[,1]
WEGU_res <- scale(residuals(lm(WEGU ~ YEAR, data = SF_bs2)))[,1]

SF_bs3 <- data.frame(YEAR = SF_bs2$YEAR, ASSP_res, BRCO_res, 
                     CAAU_res, COMU2_res, PECO_res, PIGU_res,
                     RHAU_res, WEGU_res)

#correlation matrix
SF_cc <- cor(SF_bs3[,-1])


# Tern Island
# source: Dearborn et al. 2001
TI_bs <- read.csv('TI_data.csv')

#only years where all species have BS data
usp <- unique(TI_bs$SPECIES)
BFAL <- dplyr::filter(TI_bs, SPECIES == usp[1])
BFBO <- dplyr::filter(TI_bs, SPECIES == usp[2])
BLNO <- dplyr::filter(TI_bs, SPECIES == usp[3])
LAAL <- dplyr::filter(TI_bs, SPECIES == usp[4])
RTTR <- dplyr::filter(TI_bs, SPECIES == usp[5])
WHTE <- dplyr::filter(TI_bs, SPECIES == usp[6])

#left_join all data.frames
TI_bs1 <- list(BFAL[,c('YEAR', 'BS')], BFBO[,c('YEAR', 'BS')], 
               BLNO[,c('YEAR', 'BS')], LAAL[,c('YEAR', 'BS')], 
               RTTR[,c('YEAR', 'BS')], WHTE[,c('YEAR', 'BS')]) %>% 
  purrr::reduce(dplyr::full_join, 
                by = 'YEAR')
colnames(TI_bs1) <- c('YEAR', 'BFAL', 'BFBO', 'BLNO', 'LAAL', 'RTTR', 'WHTE')

rna <- which(is.na(TI_bs1), arr.ind = TRUE)[,1]
d_idx <- rna[which(!duplicated(rna))]

TI_bs2 <- TI_bs1[-d_idx,]

BFAL_res <- scale(residuals(lm(BFAL ~ YEAR, data = TI_bs2)))[,1]
BFBO_res <- scale(residuals(lm(BFBO ~ YEAR, data = TI_bs2)))[,1]
BLNO_res <- scale(residuals(lm(BLNO ~ YEAR, data = TI_bs2)))[,1]
LAAL_res <- scale(residuals(lm(LAAL ~ YEAR, data = TI_bs2)))[,1]
RTTR_res <- scale(residuals(lm(RTTR ~ YEAR, data = TI_bs2)))[,1]
WHTE_res <- scale(residuals(lm(WHTE ~ YEAR, data = TI_bs2)))[,1]

TI_bs3 <- data.frame(YEAR = TI_bs2$YEAR, BFAL_res, BFBO_res, 
                     BLNO_res, LAAL_res, RTTR_res, WHTE_res)

#correlation matrix
TI_cc <- cor(TI_bs3[,-1])


#mean correlation
mn_PG <- mean(PG_cc[upper.tri(PG_cc)])
mn_PG_ne <- mean(PG_cc_ne[upper.tri(PG_cc_ne)])
mn_IM <- mean(IM_cc[upper.tri(IM_cc)])
mn_SF <- mean(SF_cc[upper.tri(SF_cc)])
mn_TI <- mean(TI_cc[upper.tri(TI_cc)])


# compare env variability and synchrony -----------------------------------

#processed from globcolour product - GSM method algo
chl_data <- read.csv('chl.csv')

#average Chl over specified months for each site/year
#PG: Dec - March
#IM: April - Sep
#SF: March - July
#TI: Dec - July

tf1_PG <- dplyr::filter(chl_data, SITE == 'PG', 
                        MONTH %in% 12)
#add one to year for Dec
tf1_PG$YEAR <- as.numeric(tf1_PG$YEAR) + 1
tf2_PG <- dplyr::filter(chl_data, SITE == 'PG', 
                        MONTH %in% 1:3)
tf3_PG <- rbind(tf1_PG, tf2_PG)

chl_PG <- group_by(tf3_PG, YEAR) %>% 
  summarize(mean_CHL = mean(CHL, na.rm = TRUE))

tf_IM <- dplyr::filter(chl_data, SITE == 'IM', 
                       MONTH %in% 4:9)

chl_IM <- group_by(tf_IM, YEAR) %>% 
  summarize(mean_CHL = mean(CHL, na.rm = TRUE))

tf_SF <- dplyr::filter(chl_data, SITE == 'SF',
                       MONTH %in% 3:7)
chl_SF <- group_by(tf_SF, YEAR) %>% 
  summarize(mean_CHL = mean(CHL, na.rm = TRUE))

tf1_TI <- dplyr::filter(chl_data, SITE == 'TI', 
                        MONTH %in% 12)
#add one to year for Dec
tf1_TI$YEAR <- as.numeric(tf1_TI$YEAR) + 1
tf2_TI <- dplyr::filter(chl_data, SITE == 'TI',
                        MONTH %in% 1:7)
tf3_TI <- rbind(tf1_TI, tf2_TI)
chl_TI <- group_by(tf3_TI, YEAR) %>% 
  summarize(mean_CHL = mean(CHL, na.rm = TRUE))

pl_df <- data.frame(SITE = c('PG', 'IM', 'SF', 'TI'),
                    CHL_CV = c(sd(chl_PG$mean_CHL) / mean(chl_PG$mean_CHL),
                               sd(chl_IM$mean_CHL) / mean(chl_IM$mean_CHL),
                               sd(chl_SF$mean_CHL) / mean(chl_SF$mean_CHL),
                               sd(chl_TI$mean_CHL) / mean(chl_TI$mean_CHL)),
                    SYN = c(mn_PG, mn_IM, 
                            mn_SF, mn_TI))

#fit model
fit_syn <- rstanarm::stan_glm(SYN ~ CHL_CV, 
                              data = pl_df,
                              iter = 10000)

chl_cv_ch <- MCMCvis::MCMCchains(fit_syn, params = 'CHL_CV')


#prob slope < 0
ch_l0 <- sum(chl_cv_ch[,1] < 0) / NROW(chl_cv_ch)

#summary
fsummary <- MCMCvis::MCMCsummary(fit_syn, round = 3)


#posterior for int and slope
ab_ch <- MCMCvis::MCMCchains(fit_syn, param = c('(Intercept)', 'CHL_CV'))

#simulated x to plot best fit line and uncertainty
x_sim <- seq(range(pl_df$CHL_CV)[1] - 0.01, range(pl_df$CHL_CV)[2] + 0.01, length.out = 100)
mf <- matrix(NA, nrow = NROW(ab_ch), ncol = length(x_sim))
for (i in 1:NROW(ab_ch))
{
  mf[i,] <- ab_ch[i,1] + ab_ch[i,2] * x_sim
}

fit_mn <- apply(mf, 2, mean)
fit_LCI <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
fit_UCI <- apply(mf, 2, function(x) quantile(x, probs = 0.975))

f_pl_df <- data.frame(fit_mn,
                      fit_LCI,
                      fit_UCI,
                      x_sim)

#plot
fig_B <- ggplot(data = pl_df, aes(CHL_CV, SYN)) +
  geom_point(size = 3) +
  geom_line(data = f_pl_df, aes(x_sim, fit_mn),
            col = 'red', alpha = 0.5, size = 2,
            inherit.aes = FALSE, lty = 2) +
  xlab(expression(paste('Coefficient of variation chl-a (', 
                        CV[chl], ')'))) +
  ylab(expression(paste('Community synchrony (', bar(rho), ')'))) +
  annotate('text', x = 0.430, 0.208, label = 'Pointe Géologie') +
  annotate('text', x = 0.208, 0.252, label = 'Isle of May') +
  annotate('text', x = 0.270, 0.317, label = 'Southeast Farallon Island') +
  annotate('text', x = 0.117, 0.325, label = 'Tern Island') +
  theme(
    panel.grid.major = element_line(color = 'grey80'),
    panel.grid.minor = element_line(color = 'grey90'),
    panel.background = element_blank(),
    panel.border = element_rect(fill= NA, color = 'black', size = 1.5),
    axis.ticks.length = unit(0.2, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18))


# evenness ----------------------------------------------------------------

#range body mass from Billerman et al. 2020 - use mean
#3800 - 8200g
mADPE <- median(ADPE$ABUN, na.rm = TRUE) * mean(c(3800, 8200))
#720 - 1180g
mSOFU <- median(SOFU$ABUN, na.rm = TRUE) * mean(c(720, 1180))
#340 - 528g
mCAPE <- median(CAPE$ABUN, na.rm = TRUE) * mean(c(340, 528))
#202 - 570g
mSNPE <- median(SNPE$ABUN, na.rm = TRUE) * mean(c(202, 570))
#900 - 1600g
mSPSK <- median(SPSK$ABUN, na.rm = TRUE) * mean(c(900, 1600))

s_all_PG <- sum(sort(c(mADPE, mSOFU, mCAPE, mSNPE, mSPSK), decreasing = TRUE))

temp_PG <- c((mADPE/s_all_PG), (mSOFU/s_all_PG), (mCAPE/s_all_PG),
             (mSNPE/s_all_PG), (mSPSK/s_all_PG))

evenness_PG <- (s_all_PG/(s_all_PG-1))*(1-(sum(temp_PG^2)))


#Isle of May
#310 - 550g
mATPU <- median(ATPU$ABUN, na.rm = TRUE) * mean(c(310, 550))
#365 - 400g
mBLKI <- median(BLKI$ABUN, na.rm = TRUE) * mean(c(365, 400))
#800 – 1125g
mCOMU <- median(COMU$ABUN, na.rm = TRUE) * mean(c(800, 1125))
#1407 – 2154g
mEUSH <- median(EUSH$ABUN, na.rm = TRUE) * mean(c(1407, 2154))
#505 – 890g
mRAZO <- median(RAZO$ABUN, na.rm = TRUE) * mean(c(505, 890))

s_all_IM <- sum(sort(c(mATPU, mBLKI, mCOMU, mEUSH, mRAZO), decreasing = TRUE))

temp_IM <- c((mATPU/s_all_IM), (mBLKI/s_all_IM), (mCOMU/s_all_IM),
             (mEUSH/s_all_IM), (mRAZO/s_all_IM))

evenness_IM <- (s_all_IM/(s_all_IM-1))*(1-(sum(temp_IM^2)))


# Southeast Farallon Island
#NO ABUNDANCE DATA FOR ASSP OR RHAU - EXCLUDE FROM EVENNESS CALCULATION
#35 - 40g
#mASSP <- median(ASSP$ABUN, na.rm = TRUE) * mean(c(35, 40))
#1930 - 2570g
mBRCO <- median(BRCO$ABUN, na.rm = TRUE) * mean(c(1930, 2570))
#150 – 200g
mCAAU <- median(CAAU$ABUN, na.rm = TRUE) * mean(c(150, 200))
#800 - 1125g
mCOMU2 <- median(COMU2$ABUN, na.rm = TRUE) * mean(c(800, 1125))
#1531 - 2034
mPECO <- median(PECO$ABUN, na.rm = TRUE) * mean(c(1531, 2034))
#450 - 550g
mPIGU <- median(PIGU$ABUN, na.rm = TRUE) * mean(c(450, 550))
#389 – 569g
#mRHAU <- median(RHAU$ABUN, na.rm = TRUE) * mean(c(389, 569))
#1050 - 1250g
mWEGU <- median(WEGU$ABUN, na.rm = TRUE) * mean(c(1050, 1250))

s_all_SF <- sum(sort(c(mBRCO, mCAAU, mCOMU2, mPECO, mPIGU, mWEGU), decreasing = TRUE), na.rm = TRUE)

temp_SF <- c((mBRCO/s_all_SF), (mCAAU/s_all_SF), (mCOMU2/s_all_SF), (mPECO/s_all_SF), 
             (mPIGU/s_all_SF), (mWEGU/s_all_SF))

evenness_SF <- (s_all_SF/(s_all_SF-1))*(1-(sum(temp_SF^2, na.rm = TRUE)))


# Tern Island
#2990 - 3400g
mBFAL <- median(BFAL$ABUN, na.rm = TRUE) * mean(c(2990, 3400))
#1348 - 1702g
mBFBO <- median(BFBO$ABUN, na.rm = TRUE) * mean(c(1348, 1702))
#85 – 140g
mBLNO <- median(BLNO$ABUN, na.rm = TRUE) * mean(c(85, 140))
#1900 - 3100
mLAAL <- median(LAAL$ABUN, na.rm = TRUE) * mean(c(1900, 3100))
#600 - 835
mRTTR <- median(RTTR$ABUN, na.rm = TRUE) * mean(c(600, 835))
#77 - 157g
mWHTE <- median(WHTE$ABUN, na.rm = TRUE) * mean(c(77, 157))

s_all_TI <- sum(sort(c(mBFAL, mBFBO, mBLNO, mLAAL, mRTTR, mWHTE), decreasing = TRUE))

temp_TI <- c((mBFAL/s_all_TI), (mBFBO/s_all_TI), (mBLNO/s_all_TI),
             (mLAAL/s_all_TI), (mRTTR/s_all_TI), (mWHTE/s_all_TI))

evenness_TI <- (s_all_TI/(s_all_TI-1))*(1-(sum(temp_TI^2)))


#combine evenness estimates
e_df <- data.frame(SITE = c('PG', 'IM', 'SF', 'TI'),
                   evenness = c(evenness_PG, evenness_IM, evenness_SF, evenness_TI))

#Sugihara et al. 2003 Fig. 4 data
sea2003 <- read.csv('Sugihara_et_al_2003_fig_4.csv')

cols <- colorspace::qualitative_hcl(n = 4)
ggplot(sea2003, aes(Evenness)) +
  geom_vline(data = e_df, aes(xintercept = evenness, color = factor(SITE)), 
             #col = 'red', 
             size = 2, alpha = 0.8, lty = 2) +
  geom_histogram(binwidth = 0.1, colour = "black", alpha = 0.4) + 
  theme_bw() +
  theme(
    #axis.text = element_text(size = 14), #axis label size
    panel.grid.major = element_line(color = 'grey80'), #lower # is darker
    panel.grid.minor = element_line(color = 'grey90'),
    panel.background = element_blank(),
    panel.border = element_rect(fill= NA, color = 'black', size = 1.5),
    axis.ticks.length = unit(0.2, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)) +
  xlim(0, 1) +
  ylab('Count') +
  xlab('Community evenness') +
  labs(color = 'Site') +
  scale_color_manual(labels = c('Isle of May', 'Pointe Géologie', 'Southeast Farallon Island', 'Tern Island'), values = cols)


# effect of extreme years on sync estimates -------------------------------

ES_fun <- function(INPUT, NSP)
{
  #get years with no missing BS
  nyr <- dplyr::filter(INPUT, !is.na(BS)) %>%
    dplyr::count(YEAR)
  
  #filter for those years
  d3 <- INPUT %>%
    dplyr::filter(YEAR %in% nyr$YEAR[nyr$n == NSP])  
  
  ###########
  #fit model - intercept for species and intercept for year
  fit2 <- rstanarm::stan_glmer(BS ~ 1 + (1 | SPECIES) + (1 | YEAR), 
                               data = d3,
                               iter = 10000,
                               adapt_delta = 0.99)
  #extract chains for intercepts
  b_ch <- MCMCvis::MCMCchains(fit2, params = 'b')
  
  #year columns
  yr_cols <- grep('YEAR', colnames(b_ch))
  
  #exclude final year col (new year)
  nb_ch <- b_ch[,yr_cols[-length(yr_cols)]]
  
  #calculate median absolute deviation
  madv <- rep(NA, NROW(nb_ch))
  for (i in 1:NROW(nb_ch))
  {
    #i <- 1
    madv[i] <- mad(nb_ch[i,])
  }
  
  #posterior means for delta
  nb_mn <- apply(nb_ch, 2, mean)
  
  #extreme years
  EYR <- unique(d3$YEAR)[which((abs(nb_mn - median(nb_mn)) / mean(madv)) > 2)]
  ###########
  
  
  #if there are extreme years
  if (length(EYR) > 0) 
  {
    #remove extreme year
    '%ni%' <- Negate('%in%')
    
    d4 <- d3 %>%
      dplyr::filter(YEAR %ni% EYR)
    
    usp <- unique(d4$SPECIES)
    t_res_ne <- data.frame(YEAR = unique(d4$YEAR))
    for (k in 1:length(usp))
    {
      #k <- 1
      tsp <- dplyr::filter(d4, SPECIES == usp[k])
      tres <- residuals(lm(BS ~ YEAR, data = tsp))
      t_res_ne <- cbind(t_res_ne, tres)
    }
    
    colnames(t_res_ne) <- c('YEAR', usp)
    
    #correlation matrix
    cc_ne <- cor(t_res_ne[,-1])
    mn_ne <- mean(cc_ne[upper.tri(cc_ne)])
    
  } else {
    EYR <- NA
    mn_ne <- NA
  }
  
  usp <- unique(d3$SPECIES)
  t_res_e <- data.frame(YEAR = unique(d3$YEAR))
  for (k in 1:length(usp))
  {
    #k <- 1
    tsp_e <- dplyr::filter(d3, SPECIES == usp[k])
    tres_e <- residuals(lm(BS ~ YEAR, data = tsp_e))
    t_res_e <- cbind(t_res_e, tres_e)
  }
  
  colnames(t_res_e) <- c('YEAR', usp)
  
  #correlation matrix
  cc_e <- cor(t_res_e[,-1])
  mn_e <- mean(cc_e[upper.tri(cc_e)])
  
  out <- list(EYR = EYR,
              cor = mn_e,
              ne_cor = mn_ne)
 
  return(out)
}

#PG
PG_c <- ES_fun(INPUT = PG_bs, NSP = 5)

#IM
IM_c <- ES_fun(INPUT = IM_bs, NSP = 5)

#SF
SF_c <- ES_fun(INPUT = SF_bs, NSP = 8)

#TI
TI_c <- ES_fun(INPUT = TI_bs, NSP = 6)
