#################
# Youngflesh et al. - 5 - Simulation study
#################

# Lack of synchronized breeding success in a seabird community: extreme events, niche separation, and environmental variability

# Casey Youngflesh, Yun Li, Heather J. Lynch, Karine Delord, Christophe Barbraud, Rubao Ji, Stephanie Jenouvrier


# load packages -----------------------------------------------------------

library(boot)
library(dplyr)
library(reshape2)
library(ggplot2)


# run simulation -------------------------------------------------------

#number of species
NS <- 6
#number of years
YRS <- 100

#function to simulate data
sim_fun <- function(seed, beta1_sd, beta2_sd, d_sd, s_sd, eps_sd)
{
  set.seed(seed)
  
  #simulate relationship between env and pro
  beta1 <- rnorm(NS, 0, beta1_sd)
  beta2 <- rnorm(1, 0, beta2_sd) 
  
  #simulate env
  env_d <- rnorm(YRS, 0, d_sd)
  env_s <- rnorm(YRS, 0, s_sd)
  
  #simulate pro
  psim <- matrix(nrow = YRS, ncol = NS)
  for (i in 1:NS)
  {
    #i <- 1
    psim[,i] <- boot::inv.logit((beta1[i] * env_d) + (beta2 * env_s) + rnorm(YRS, 0, eps_sd))
  }
  
  out <- list(seed = seed,
              beta1 = beta1,
              beta2 = beta2,
              env_d = env_d,
              env_s = env_s,
              psim = psim)
  
  return(out)
}

#Run each set of scenarios 100 times
mn_rho_hd_ls <- rep(NA, 100)
mn_rho_ld_hs <- rep(NA, 100)
for (i in 1:100)
{
  #i <- 1
  #high var env_d, low var env_s
  hd_ls <- sim_fun(seed = i, d_sd = 5, s_sd = 1,
                   beta1_sd = 0.1, beta2_sd = 0.1, eps_sd = 0.3)
  #low var env_d, high var env_s
  ld_hs <- sim_fun(seed = i, d_sd = 1, s_sd = 5,
                   beta1_sd = 0.1, beta2_sd = 0.1, eps_sd = 0.3)
  
  #mean correlation for high var env_d - low var env_s
  rho_hd_ls <- cor(hd_ls$psim)
  mn_rho_hd_ls[i] <- mean(rho_hd_ls[upper.tri(rho_hd_ls)])
  
  #mean correlation for low var env_d - high var env_s
  rho_ld_hs <- cor(ld_hs$psim)
  mn_rho_ld_hs[i] <- mean(rho_ld_hs[upper.tri(rho_ld_hs)])
}

#specify colors
cols <- colorspace::qualitative_hcl(n = 6)

#simulate data with seed = 51 to plot
#high var env_d, low var env_s
SEED <- 51
hd_ls <- sim_fun(seed = SEED, d_sd = 5, s_sd = 1,
                 beta1_sd = 0.1, beta2_sd = 0.1, eps_sd = 0.3)
#low var env_d, high var env_s
ld_hs <- sim_fun(seed = SEED, d_sd = 1, s_sd = 5,
                 beta1_sd = 0.1, beta2_sd = 0.1, eps_sd = 0.3)

#calculate rho bar for seed
rho_hd_ls <- cor(hd_ls$psim)
mn_rho_hd_ls <- mean(rho_hd_ls[upper.tri(rho_hd_ls)])

rho_ld_hs <- cor(ld_hs$psim)
mn_rho_ld_hs <- mean(rho_ld_hs[upper.tri(rho_ld_hs)])

#simulation 1
reshape2::melt(data.frame(time = 1:100, hd_ls$psim), 'time') %>%
  ggplot(aes(time, value, color = variable)) +
  geom_line(alpha = 0.4, size = 1) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'grey80'), #lower # is darker
    panel.grid.minor = element_line(color = 'grey90'),
    panel.background = element_blank(),
    panel.border = element_rect(fill= NA, color = 'black', size = 1.5),
    axis.ticks.length = unit(0.2, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)) +
  labs(color = 'Hypothetical Species') +
  scale_color_manual(labels = paste0(1:6), values = cols) +
  ggtitle(latex2exp::TeX(paste0('Simulation 1 - d = 5, s = 1; $\\bar{\\rho} = ', round(mn_rho_hd_ls, 2)))) +
  #ggtitle(latex2exp::TeX('$\\bar{\\rho}')) +
  xlab('Time Step') +
  ylab(latex2exp::TeX('$p_{sim}')) +
  ylim(c(0,1))


#simulation 2
reshape2::melt(data.frame(time = 1:100, ld_hs$psim), 'time') %>%
  ggplot(aes(time, value, color = variable)) +
  geom_line(alpha = 0.4, size = 1) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'grey80'), #lower # is darker
    panel.grid.minor = element_line(color = 'grey90'),
    panel.background = element_blank(),
    panel.border = element_rect(fill= NA, color = 'black', size = 1.5),
    axis.ticks.length = unit(0.2, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18)) +
  labs(color = 'Hypothetical Species') +
  scale_color_manual(labels = paste0(1:6), values = cols) +
  ggtitle('Simulation 2 - d = 1, s = 5') +
  ggtitle(latex2exp::TeX(paste0('Simulation 2 - d = 1, s = 5; $\\bar{\\rho} = ', round(mn_rho_ld_hs, 2)))) +
  xlab('Time Step') +
  ylab(latex2exp::TeX('$p_{sim}')) +
  ylim(c(0,1))
