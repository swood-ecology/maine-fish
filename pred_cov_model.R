##################################
#Ex 6 State Space - JAGS
##################################

#add hierarchy to get at relationship with env covariates
#add density dependence



# set wd ------------------------------------------------------------------

setwd("~/Google_Drive/R/NEFI_course/EF_Activities/data")



# load packages -----------------------------------------------------------

library(rjags)
library(jagsRun)
library(MCMCvis)



# load data ---------------------------------------------------------------

setwd('~/Google_Drive/R/NEFI_course/Fish_project/Data/')
data <- readRDS('Allyn_EcoForecastingProjectData.rds')


#filter by season and species
species <- unique(data$CommonName)

#species that move less:
#lobseter
#black seabass
#probs squid

SP_NUM <- 8
d_spr <- dplyr::filter(data, 
                         CommonName == species[SP_NUM],
                         Season == 'SPRING')

d_fall <- dplyr::filter(data, 
                       CommonName == species[SP_NUM],
                       Season == 'FALL')

yrs <- unique(d_spr$Year)

#sort data into df - group by year over space
out <- data.frame()
for (i in 1:length(yrs))
{
  #i <- 1
  temp <- dplyr::filter(d_spr, Year == yrs[i])
  BM <- sum(temp$Biomass)
  sst <- mean(temp$SeasonalSST)
  depth <- mean(temp$Depth)
  tt <- data.frame(YEAR = yrs[i],
                   SPECIES = temp$CommonName[1],
                   BM = BM,
                   SST = sst,
                   DEPTH = depth)
  out <- rbind(out, tt)
}

N <- NROW(out)

BM2 <- out$BM
BM2[(N - 2):N] <- NA

DATA <- list(
  y = log(BM2),
  year = out$YEAR,
  depth = out$DEPTH,
  sst = out$SST,
  N = N,
  x_ic = mean(log(out$BM)),
  tau_ic = 10,
  a_obs = 1,
  r_obs = 1,
  a_add = 1,
  r_add = 1)


# Model -------------------------------------------------------------------

setwd("~/Google_Drive/R/NEFI_course/Results")

{
  sink('RW.jags')
  
  cat("
  
  model {
  
  #observation model
  for (t in 1:N)
  {
    y[t] ~ dnorm(x[t], tau_obs)
  }
  
  #data model
  for (t in 2:N)
  {
    x[t] ~ dnorm(mu[t], tau_add)
    mu[t] <- x[t-1] + beta_0 + beta_x * x[t-1] + 
              beta_1 * depth[t] + beta_2 * sst[t]
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  
  beta_0 ~ dnorm(0, 0.001)
  beta_x ~ dnorm(0, 0.001)
  beta_1 ~ dnorm(0, 0.001)
  beta_2 ~ dnorm(0, 0.001)
  tau_obs ~ dgamma(a_obs, r_obs)
  tau_add ~ dgamma(a_add, r_add)

  }",fill = TRUE)
  
  sink()
}



# Starting values ---------------------------------------------------------

Inits_1 <- list(
                .RNG.name = "base::Mersenne-Twister", 
                .RNG.seed = 1)

Inits_2 <- list(
                .RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = 2)

Inits_3 <- list(
                .RNG.name = "base::Marsaglia-Multicarry", 
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('tau_obs',
          'tau_add',
          'beta_0',
          'beta_x',
          'beta_1',
          'beta_2',
          'x')


# Run model ---------------------------------------------------------------


#make sure model compiles
# jagsRun(jagsData = DATA,
#         jagsModel = 'pwatch_surv.jags',
#         jagsInits = F_Inits,
#         DEBUG = TRUE)

fit <- jagsRun(jagsData = DATA, 
               jagsModel = 'RW.jags',
               jagsInits = F_Inits,
               params = Pars,
               jagsID = 'BLSE_1',
               jagsDsc = 'random walk with cov',
               obj_out = TRUE,
               n_chain = 3,
               n_adapt = 2000,
               n_burn = 40000,
               n_draw = 40000,
               n_thin = 1,
               EXTRA = FALSE,
               Rhat_max = 1.1,
               n_max = 100000)


#summaries
MCMCvis::MCMCsummary(fit, params = c('tau_obs', 'tau_add',
                                     'beta_x', 'beta_0',  
                                     'beta_1', 'beta_2', 'x'), round = 2)
MCMCvis::MCMCplot(fit, params = 'x')
MCMCvis::MCMCplot(fit, params = c('beta_1', 'beta_2'))

#MCMCvis::MCMCtrace(fit, params = c('tau_obs', 'tau_add', 'x'), ind = TRUE)



#extract posteriors
x_med <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) median(exp(x)))[[1]]
x_LCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.025))[[1]]
x_UCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.975))[[1]]

com <- c(x_med, x_LCI, x_UCI)

#plot
plot(DATA$year, x_med, type = 'n', ylim = range(com, na.rm = TRUE), 
     ylab = "Biomass", xlim = range(DATA$year))
#mean
lines(DATA$year, x_med, col = 'red', lwd = 2)
#LCI
lines(DATA$year, x_LCI, col = 'red', lwd = 2, lty = 2)
#UCI
lines(DATA$year, x_UCI, col = 'red', lwd = 2, lty = 2)
#observed data
lines(DATA$year, exp(DATA$y), lty = 2)

#plot held out data
na.val <- which(is.na(BM2))
points(DATA$year[na.val], out$BM[na.val], pch = 19)


