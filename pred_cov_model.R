###############
#Predict lobster abundance with state space and covariates
#
###############


# set wd ------------------------------------------------------------------

setwd("~/Google_Drive/R/NEFI_course/EF_Activities/data")



# load packages -----------------------------------------------------------

library(rjags)
library(MCMCvis)



# load data ---------------------------------------------------------------

setwd('~/Google_Drive/R/NEFI_course/Fish_project/Data/')
data <- readRDS('LobsterEcoForecastingProjectData.rds')

d_spr <- dplyr::filter(data,
                       Season == 'FALL')

yrs <- unique(d_spr$Year)
zones <- unique(d_spr$ZONEID)

#sort data into df - group by year over space
out <- matrix(NA, nrow = length(yrs), ncol =  length(zones))
SST <- matrix(NA, nrow = length(yrs), ncol =  length(zones))
for (i in 1:length(yrs))
{
  #i <- 1
  temp <- dplyr::filter(d_spr, Year == yrs[i], )
  for (j in 1:length(zones))
  {
    #j <- 1
    temp_z <- dplyr::filter(temp, ZONEID == zones[j])
    out[i,j] <- mean(temp_z$Biomass)
    SST[i,j] <- mean(temp_z$SeasonalSST)
  }
}

#replace 0 with small value (log does like 0)
z_ind <- which(out == 0, arr.ind = TRUE)
out[z_ind] <- 0.01

#number of years
NY <- NROW(out)

#number of zones
NZ <- NCOL(out)

#i = year
#j = zone
out2 <- out
to.na <- (NY-3):NY
out2[to.na,] <- NA

DATA <- list(
  y = log(out2),
  sst = SST,
  NY = NY,
  NZ = NZ,
  x_ic = mean(log(out)), #mean of biomass for initial value
  tau_ic = 1/(var(as.vector(log(out)))*100)) #variance biomass times 100



# Model -------------------------------------------------------------------

setwd("~/Google_Drive/R/NEFI_course/Results")

{
  sink('mat.jags')
  
  cat("
  
  model {
  
  #observation model
  for (i in 1:NY)
  {
    for (j in 1:NZ)
    {
      y[i,j] ~ dnorm(x[i,j], tau_obs)
      ysim[i,j] ~ dnorm(x[i,j], tau_obs)
    }
    #sum ysim across all sites
    ysim_all[i] <- mean(ysim[i,])
  }
  
  
  #process model
  for (i in 2:NY)
  {
    for (j in 1:NZ)
    {
      x[i,j] ~ dnorm(mu[i,j], tau_add)
      mu[i,j] <- x[i-1,j] + eta + gamma[j] + #alpha[i] +
                  beta_x * x[i-1,j] + kappa * sst[i,j]
    }
    x_all[i] <- mean(x[i,])
  }
  
  #### Priors
  #intial state
  for (j in 1:NZ)
  {
    x[1,j] ~ dnorm(x_ic, tau_ic)
  }
  #first time step sum
  x_all[1] <- mean(x[1,])
  
  #eta - grand mean
  eta ~ dnorm(0, 0.01)
  
  # #alpha - year effect
  # for (i in 1:NY)
  # {
  #   alpha[i] ~ dnorm(0, tau_alpha)
  # }
 
  #gamma - zone effect
  for (j in 1:NZ)
  {
    gamma[j] ~ dnorm(0, tau_gamma)
  }
 
  #beta_x - temporal autocorrelation effect
  beta_x ~ dnorm(0, 0.01)
  
  #k - carrying capacity
  # for (i in 2:NY)
  # {
  #   for (j in 1:NZ)
  #   {
  #     # k[i,j] ~ dnorm(mu_k[i,j], tau_k)
  #     # mu_k[i,j] = alpha_k + beta_k * sst[i,j]
  #     #k[i,j] = alpha_k + beta_k * sst[i,j]
  #   }
  # }
  # k ~ dnorm(0, 0.01)
  
  #kappa - SST effect
  kappa ~ dnorm(0, 0.01)
  
  #alpha_k - intercept SST effect on k
  # for (j in 1:NZ)
  # {
  #   alpha_k[j] ~ dnorm(mu_alpha_k, tau_alpha_k)
  # }
  # mu_alpha_k ~ dnorm(0, 0.01)
  # tau_alpha_k ~ dgamma(0.01, 0.01)
  # alpha_k ~ dnorm(0, 0.01)
  # 
  # #beta_k - slope SST effect on k
  # beta_k ~ dnorm(0, 0.01)
  # 
  # #tau_k - precision k
  # tau_k ~ dgamma(0.01, 0.01)
 
  #observation error
  tau_obs ~ dgamma(0.1, 0.1)
 
  #process eror
  tau_add ~ dgamma(0.1, 0.1)
 
  #prec alpha
  tau_alpha ~ dgamma(0.1, 0.1)
 
  #prec gamma
  tau_gamma ~ dgamma(0.1, 0.1)

  }",fill = TRUE)
  
  sink()
}



# Starting values ---------------------------------------------------------

x_init <- matrix(nrow = DATA$NY, ncol = DATA$NZ)
x_init[1,] <- rep(1, DATA$NZ)

Inits_1 <- list(x = x_init,
  .RNG.name = "base::Mersenne-Twister",
  .RNG.seed = 1)

Inits_2 <- list(x = x_init,
  .RNG.name = "base::Wichmann-Hill",
  .RNG.seed = 2)

Inits_3 <- list(x = x_init,
  .RNG.name = "base::Marsaglia-Multicarry",
  .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('tau_obs',
          'tau_add',
          'eta',
          'alpha',
          'gamma',
          'tau_alpha',
          'tau_gamma',
          'beta_x',
          # 'k',
          # 'alpha_k',
          # 'mu_alpha_k',
          # 'tau_alpha_k',
          # 'beta_k',
          # 'tau_k',
          'kappa',
          'ysim',
          'ysim_all',
          'x_all',
          'x')


# Run model ---------------------------------------------------------------


#make sure model compiles
# jagsRun(jagsData = DATA,
#         jagsModel = 'pwatch_surv.jags',
#         jagsInits = F_Inits,
#         DEBUG = TRUE)


rjags::load.module("glm")

#compile model
jm <- rjags::jags.model(data = DATA, 
                        file = 'mat.jags', 
                        inits = F_Inits, 
                        n.chains = 3, 
                        n.adapt = 5000)

#burn-in
stats::update(jm, 
              n.iter = 50000)

#draw samples
fit <- rjags::coda.samples(jm, 
                               n.iter = 50000, 
                               variable.names = Pars, 
                               thin = 1)


#summarize output
MCMCvis::MCMCsummary(fit, excl = c('x', 'x_all', 'ysim', 'ysim_all'), round = 2)
#MCMCvis::MCMCplot(fit, params = 'alpha')
MCMCvis::MCMCplot(fit, params = 'gamma')

#extract median and CI for x
x_med <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) median(exp(x)))[[1]]
x_LCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.025))[[1]]
x_UCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.975))[[1]]

#extract median and CI for ysim (for prediction interval)
ys_med <- MCMCvis::MCMCpstr(fit, params = 'ysim', 
                           func = function(x) median(exp(x)))[[1]]
ys_LCI <- MCMCvis::MCMCpstr(fit, params = 'ysim', 
                           func = function(x) quantile(exp(x), probs = 0.025))[[1]]
ys_UCI <- MCMCvis::MCMCpstr(fit, params = 'ysim', 
                           func = function(x) quantile(exp(x), probs = 0.975))[[1]]


#combine to get range for plot
#com <- c(x_med, x_LCI, x_UCI)
com <- c(ys_med, ys_LCI, ys_UCI)


plt_df <- data.frame(zone = rep(zones, each = NROW(x_med)),
                     year = rep(yrs, NCOL(x_med)),
                     x_med = as.vector(x_med),
                     x_LCI = as.vector(x_LCI),
                     x_UCI = as.vector(x_UCI),
                     ys_med = as.vector(ys_med),
                     ys_LCI = as.vector(ys_LCI),
                     ys_UCI = as.vector(ys_UCI))

#data witheld
out_wh <- log(out[to.na,])

#plot
time <- 1:DATA$NY
# plot(time, x_med[,1], type = 'n', ylim = range(com, na.rm = TRUE),
#      ylab = "Biomass")
par(mfrow = c(3, 3))
for (i in 1:NCOL(x_med))
{
  plot(time, x_med[,1], type = 'n', ylim = range(com, na.rm = TRUE),
       ylab = "Biomass", main = paste0('Zone ', zones[i]))
  #i <- 2
  #pred
  polygon(cbind(c(time, rev(time), time[1]), 
                c(ys_LCI[,i], rev(ys_UCI[,i]), ys_LCI[1,i])), 
          border = FALSE, col = rgb(0,0,1,0.3))
  
  #CI
  polygon(cbind(c(time, rev(time), time[1]), 
                c(x_LCI[,i], rev(x_UCI[,i]), x_LCI[1,i])), 
          border = FALSE, col = rgb(1,0,0,0.3))
  
  lines(time, x_med[,i]) #model mean
  
  #actual data
  points(time, exp(DATA$y[,i]), pch = '+', col = 'black')
  #witheld data
  points(to.na, exp(out_wh[,i]), pch = '+', col = 'red')
}


#probability integral transform plot
pin_vec <- rep(NA, 19)
ch_vec <- seq(0.05, 0.95, length = 19)
for (i in 1:length(ch_vec))
{
  #i <- 1
  rn <- ch_vec[i] / 2
  temp_LCI <- MCMCvis::MCMCpstr(fit, params = 'ysim', 
                                func = function(x) quantile(exp(x), probs = 0.5 - rn))[[1]]
  temp_UCI <- MCMCvis::MCMCpstr(fit, params = 'ysim', 
                                func = function(x) quantile(exp(x), probs = 0.5 + rn))[[1]]
  pin_vec[i] <- sum(exp(out_wh) > temp_LCI[to.na,] &
        exp(out_wh) < temp_UCI[to.na,]) / length(out_wh)
}

#coverage plot
par(mfrow = c(1,1))
plot(ch_vec, pin_vec, ylab = 'Prop witheld data in prediction interval',
     xlab = 'Prediction interval')
abline(a = 0, b = 1, col = 'red', lty = 2)

#RMSE = 30.6
# x_ch <- MCMCvis::MCMCchains(fit, params = 'x')
# (RMSE <- sqrt(mean((x_med[to.na,] - exp(out_wh))^2)))
# RMSE / diff(range(out))

#Predictive loss
#PL = var(resid) + var(pred)
PL <- var()


#one step ahead prediction

#extract median and CI for x
x_med_a <- MCMCvis::MCMCpstr(fit, params = 'x_all',
                           func = function(x) median(exp(x)))[[1]]
x_LCI_a <- MCMCvis::MCMCpstr(fit, params = 'x_all',
                           func = function(x) quantile(exp(x), probs = 0.025))[[1]]
x_UCI_a <- MCMCvis::MCMCpstr(fit, params = 'x_all',
                           func = function(x) quantile(exp(x), probs = 0.975))[[1]]

#extract median and CI for ysim (for prediction interval)
ys_med_a <- MCMCvis::MCMCpstr(fit, params = 'ysim_all',
                            func = function(x) median(exp(x)))[[1]]
ys_LCI_a <- MCMCvis::MCMCpstr(fit, params = 'ysim_all',
                            func = function(x) quantile(exp(x), probs = 0.025))[[1]]
ys_UCI_a <- MCMCvis::MCMCpstr(fit, params = 'ysim_all',
                            func = function(x) quantile(exp(x), probs = 0.975))[[1]]


par(mfrow = c(1, 1))
#all sites mean
plot(time, x_med_a, type = 'n', ylim = range(c(ys_LCI_a, ys_UCI_a), na.rm = TRUE),
     ylab = "Biomass")
#PI
polygon(cbind(c(time, rev(time), time[1]),
              c(ys_LCI_a, rev(ys_UCI_a), ys_LCI_a[1])),
        border = FALSE, col = rgb(0,0,1,0.3))
#CI
polygon(cbind(c(time, rev(time), time[1]),
              c(x_LCI_a, rev(x_UCI_a), x_LCI_a[1])),
        border = FALSE, col = rgb(1,0,0,0.3))

lines(1:DATA$NY, x_med_a) #model mean



#GP
#check PL
#check PL and coverage on random walk
