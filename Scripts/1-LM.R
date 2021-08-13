#################
# Youngflesh et al. - 1 - LM model - synchrony in community
#################

# Lack of synchronized breeding success in a seabird community: extreme events, niche separation, and environmental variability

# Casey Youngflesh, Yun Li, Heather J. Lynch, Karine Delord, Christophe Barbraud, Rubao Ji, Stephanie Jenouvrier


# Load packages -----------------------------------------------------------

library(dplyr)
library(rjags)
library(parallel)
library(MCMCvis)


# Data process ------------------------------------------------------------------

setwd('Data')

data <- read.csv('data_pub.csv', header = TRUE)

#separate by species
ADPE <- dplyr::filter(data, SPECIES == 'ADPE')
SOFU <- dplyr::filter(data, SPECIES == 'SOFU')
CAPE <- dplyr::filter(data, SPECIES == 'CAPE')
SNPE <- dplyr::filter(data, SPECIES == 'SNPE')
SPSK <- dplyr::filter(data, SPECIES == 'SPSK')


#sort for complete years only
BS_mrg <- data.frame(YEAR = ADPE$YEAR, 
                     ADPE_BS = ADPE$BS, 
                     SOFU_BS = SOFU$BS, 
                     CAPE_BS = CAPE$BS, 
                     SNPE_BS = SNPE$BS, 
                     SPSK_BS = SPSK$BS)

#which rows have no NAs (data for every species)
BS_ind <- which(table(which(!is.na(BS_mrg), arr.ind = TRUE)[,1]) == 6)

ADPE_f <- ADPE[BS_ind, ]
SOFU_f <- SOFU[BS_ind, ]
CAPE_f <- CAPE[BS_ind, ]
SNPE_f <- SNPE[BS_ind, ]
SPSK_f <- SPSK[BS_ind, ]


#Data for JAGS model

LM_DATA <- list(
  E1 = ADPE_f$ABUN * 2, #two eggs
  E2 = SOFU_f$ABUN, #one egg
  E3 = CAPE_f$ABUN, #one egg
  E4 = SNPE_f$ABUN, #one egg
  E5 = SPSK_f$ABUN * 2, #two eggs
  F1 = ADPE_f$CHICKS,
  F2 = SOFU_f$CHICKS,
  F3 = CAPE_f$CHICKS,
  F4 = SNPE_f$CHICKS,
  F5 = SPSK_f$CHICKS,
  N = NROW(ADPE_f))

DATA <- LM_DATA
{
  sink("LM_model.jags")
  
  cat("
      
      model {
      
      # Likelihood ------------------------------------------------------------
      
      for (i in 1:N)
      {
      F1[i] ~ dbin(rho[1,i], E1[i])
      F2[i] ~ dbin(rho[2,i], E2[i])
      F3[i] ~ dbin(rho[3,i], E3[i])
      F4[i] ~ dbin(rho[4,i], E4[i])  
      F5[i] ~ dbin(rho[5,i], E5[i])
      # where Es[i] = # eggs laid for that species * number of breeding pairs in year i
      }
      
      
      # Priors and regressions -----------------------------------------------
      for (i in 1:N)
      {
      # Logistic regressions on productivity (with random terms)
      #bs are offsets for each species
      
      logit(rho[1,i]) <- b0[1] + d[i] + e[1,i]
      logit(rho[2,i]) <- b0[2] + d[i] + e[2,i]
      logit(rho[3,i]) <- b0[3] + d[i] + e[3,i]
      logit(rho[4,i]) <- b0[4] + d[i] + e[4,i]
      logit(rho[5,i]) <- b0[5] + d[i] + e[5,i]
      
      #priors    
      d[i] ~ dnorm(0, tau.d)
      
      for (j in 1:5)
      {
      e[j, i] ~ dnorm(0, tau.e[j])
      }
      }
      
      tau.d  <- 1/var.d
      var.d <- sd.d * sd.d
      sd.d ~ dunif(0, 3)
      
      for (j in 1:5)
      {
      #priors for species effects
      b0[j] ~ dunif(-5,5)
      
      tau.e[j] <- 1/var.e[j]
      var.e[j] <- sd.e[j] * sd.e[j]
      sd.e[j] ~ dunif(0,3)
      
      # Synchrony indices
      I[j] <- var.d/(var.e[j] + var.d)
      
      }
      
      # mean synchrony index
      I_bar <- mean(I)

      #for mean productivity
      mn_pro[1] <- mean(rho[1,]) * 2
      mn_pro[2] <- mean(rho[2,])
      mn_pro[3] <- mean(rho[3,])
      mn_pro[4] <- mean(rho[4,])
      mn_pro[5] <- mean(rho[5,]) * 2
      
      }",fill = TRUE)

  sink()
}


# Starting values ---------------------------------------------------------


Inits_1 <- list(b0 = rep(0, 5),
                d = rnorm(DATA$N),
                e = matrix(rnorm(DATA$N*5), nrow = 5),
                sd.e = runif(5, 0, 2),
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(b0 = rep(0, 5),
                d = rnorm(DATA$N),
                e = matrix(rnorm(DATA$N*5), nrow = 5),
                sd.e = runif(5, 0, 2),
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(b0 = rep(0, 5),
                d = rnorm(DATA$N),
                e = matrix(rnorm(DATA$N*5), nrow = 5),
                sd.e = runif(5, 0, 2),
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('b0',
          'd',
          'var.d',
          'var.e',
          'mn_pro',
          'I',
          'I_bar')


# Inputs for MCMC ---------------------------------------------------------

JAGS_FILE <- 'LM_model.jags'
n_adapt <- 8000  # number for initial adapt
n_burn <- 30000000 # number burnin
n_draw <- 100000000  # number of final draws to make
n_thin <- 500    # thinning rate
n_chain <- 3  # number of chains


# Run model (parallel) ---------------------------------------------------------------

#number of chains
cl <- parallel::makeCluster(n_chain)

pid <- NA
for(i in 1:n_chain)
{
  pidNum <- capture.output(cl[[i]])
  start <- regexpr("pid", pidNum)[[1]]
  end <- nchar(pidNum)
  pid[i] <- substr(pidNum, (start + 4), end)
}

parallel::clusterExport(cl,
                        c('DATA',
                          'n_adapt',
                          'n_burn',
                          'n_draw',
                          'n_thin',
                          'Pars',
                          'pid',
                          'F_Inits',
                          'JAGS_FILE'
                        ))


ptm <- proc.time()
out.1 <- parallel::clusterEvalQ(cl,
                                {
                                  require(rjags)
                                  processNum <- which(pid==Sys.getpid())
                                  m.inits <- F_Inits[[processNum]]
                                  
                                  jm = jags.model(data = DATA,
                                                  file = paste0(JAGS_FILE),
                                                  inits = m.inits,
                                                  n.chains = 1,
                                                  n.adapt = n_adapt)
                                  
                                  update(jm,
                                         n.iter = n_burn)
                                  
                                  samples = coda.samples(jm,
                                                         n.iter = n_draw,
                                                         variable.names = Pars,
                                                         thin = n_thin)
                                  return(samples)
                                })


out <- coda::mcmc.list(out.1[[1]][[1]],
                       out.1[[2]][[1]],
                       out.1[[3]][[1]])

n_total <- n_burn + n_draw
n_extra <- 0

stopCluster(cl)
n_final <- floor((n_draw + n_extra)/n_thin)
tt <- (proc.time() - ptm)[3]/60 #minutes


# Analyze results ---------------------------------------------------------

MCMCvis::MCMCsummary(out, round = 4)


# identify extreme years --------------------------------------------------

d_ch <- MCMCvis::MCMCchains(out, params = 'd')

#calculate median absolute deviation
#https://www.wikiwand.com/en/Median_absolute_deviation
mad <- rep(NA, NROW(d_ch))

for (i in 1:NROW(d_ch))
{
  #i <- 1
  mad[i] <- median(abs(d_ch[i,] - median(d_ch[i,])))
}

#posterior means for delta
d_mn <- MCMCvis::MCMCpstr(out, params = 'd')[[1]]

#extreme years - sensu Palmer et al. 2017 - Eq 2.1
(abs((d_mn - median(d_mn))) / mean(mad)) > 2
