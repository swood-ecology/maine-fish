


# set wd ------------------------------------------------------------------

setwd("~/Google_Drive/R/NEFI_course/EF_Activities/data")



# load packages -----------------------------------------------------------

library(rjags)
library(jagsRun)
library(MCMCvis)



# load data ---------------------------------------------------------------

setwd('~/Google_Drive/R/NEFI_course/Fish_project/Data/')
data <- readRDS('LobsterEcoForecastingProjectData.rds')

d_spr <- dplyr::filter(data,
                       Season == 'FALL')

head(d_spr)
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

DATA <- list(
  y = log(out),
  #sst = sst,
  NY = NY,
  NZ = NZ,
  x_ic = mean(log(out)), #mean of biomass for initial value
  tau_ic = 1/(var(as.vector(log(out)))*10)) #variance biomass times 10



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
    }
  }
  
  #data model
  for (i in 2:NY)
  {
    for (j in 1:NZ)
    {
      x[i,j] ~ dnorm(mu[i,j], tau_add)
      mu[i,j] <- x[i-1,j] #+ COVARIATES AND OTHER THINGS GO HERE
    }
  }
  
  #### Priors
  #intial state
  for (j in 1:NZ)
  {
    x[1,j] ~ dnorm(x_ic, tau_ic)
  }
  #observation error
  tau_obs ~ dgamma(0.1, 0.1)
  #process eror
  tau_add ~ dgamma(0.1, 0.1)

  }",fill = TRUE)
  
  sink()
}



# Starting values ---------------------------------------------------------
# 
# Inits_1 <- list(
#   .RNG.name = "base::Mersenne-Twister", 
#   .RNG.seed = 1)
# 
# Inits_2 <- list(
#   .RNG.name = "base::Wichmann-Hill", 
#   .RNG.seed = 2)
# 
# Inits_3 <- list(
#   .RNG.name = "base::Marsaglia-Multicarry", 
#   .RNG.seed = 3)
# 
# F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('tau_obs',
          'tau_add',
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
                        #inits = start, 
                        n.chains = 3, 
                        n.adapt = 2000)
#burn-in
stats::update(jm, 
              n.iter = 50000)

#draw samples
fit <- rjags::coda.samples(jm, 
                               n.iter = 50000, 
                               variable.names = Pars, 
                               thin = 1)


#summarize output
MCMCvis::MCMCsummary(fit, excl = 'x')


#extract median and CI for x
x_med <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) median(exp(x)))[[1]]
x_LCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.025))[[1]]
x_UCI <- MCMCvis::MCMCpstr(fit, params = 'x', 
                           func = function(x) quantile(exp(x), probs = 0.975))[[1]]

#combine to get range for plot
com <- c(x_med, x_LCI, x_UCI)




#plot - NEEDS WORK
plot(1:DATA$NY, x_med[,1], type = 'n', ylim = range(com, na.rm = TRUE), 
     ylab = "Biomass")
for (i in 1:NCOL(x_med))
{
  #i <- 1
  lines(1:DATA$NY, x_med[,i]) #model mean
  # lines(1:DATA$NY, x_LCI[,i], col = 'red', lwd = 2, lty = 2) #model LCI
  # lines(1:DATA$NY, x_UCI[,i], col = 'red', lwd = 2, lty = 2) #model UCI
}
