library(tidyverse)
library(rjags)

# Read in data
lobstah.data <- 
  readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/LobsterEcoForecastingProjectData.rds")
lobstah.data <- lobstah.data %>% select(Year:Depth, ZONEID)

# Visualize data
hist(lobstah.data$Biomass)
hist(log(lobstah.data$Biomass+0.000001))

# Take log and drop outliers
lobstah.data$logBiomass <- log(lobstah.data$Biomass + 0.000001) 
lobstah <- lobstah.data %>%
  filter(logBiomass > mean(lobstah.data$logBiomass) - 2*sd(lobstah.data$logBiomass))

ggplot(lobstah, aes(y = logBiomass, x = Year)) + 
  geom_point(alpha = 0.5) + 
  theme_bw()

# Plot all data averaged together
ggplot(aggregate(cbind(logBiomass, Depth, SeasonalSST) ~ Year, lobstah, mean), 
       aes(y = logBiomass, x = Year)) + 
  geom_line() + 
  theme_bw()

# Plot seasons separately
ggplot(aggregate(cbind(logBiomass, Depth, SeasonalSST) ~ Year + Season, lobstah, 
                 mean), 
       aes(y = logBiomass, x = Year, group=Season, color=Season)) + 
  geom_line() + 
  theme_bw()

# Plot zones separately
ggplot(aggregate(cbind(logBiomass, Depth, SeasonalSST) ~ Year + ZONEID, 
                 lobstah %>% filter(Season == 'SPRING'), 
                 mean), 
       aes(y = logBiomass, x = Year, group=ZONEID, color=ZONEID)) + 
  geom_line() + 
  theme_bw()

# Aggregate data to year
year.lobstah.data <- aggregate(logBiomass ~ Year, 
                               lobstah %>% filter(Season == 'SPRING'), 
                               mean)


######################################
# Implement simple random walk model #
######################################

attach(year.lobstah.data)

# Define Random Walk model
RandomWalk = "
model{
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }

  #### Process Model
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1],tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

# Define data and parameters
data <- list(y=logBiomass,
             n=length(logBiomass),
             x_ic=2.5,
             tau_ic=1,a_obs=1.5,r_obs=1,a_add=1.5,r_add=1)

# Define initial values
nchain = 4
init <- list()
for(i in 1:nchain){
  y.samp = sample(logBiomass,length(logBiomass),
                  replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=1/var(log(y.samp)))
}

# Run model
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 4)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
gelman.plot(jags.out)

burnin = 1000                                ## determine convergence
jags.burn <- window(jags.out, start=burnin)  ## remove burn-in
summary(jags.burn) ## check diagnostics post burn-in

# CI plot
rw.out <- as.matrix(jags.burn)

x.cols <- grep("^x",colnames(rw.out)) ## grab all columns that start with the letter x
ci <- apply(rw.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(Year,ci[2,],type='n',ylim=c(0,7))
ecoforecastR::ciEnvelope(Year,ci[1,],ci[3,],
                         col=ecoforecastR::col.alpha("lightBlue",0.75))
points(Year,logBiomass,pch="+",cex=0.5)


#####################################
# Forecast simple random walk model #
#####################################

# Define function to print
plot.run <- function(){
  x.cols <- grep("^x",colnames(rw.out)) ## grab all columns that start with the letter x
  ci <- apply(rw.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
  plot(Year,ci[2,],type='n',ylim=c(0,7),xlim=c(time[1], time[length(time)]))
  ecoforecastR::ciEnvelope(Year,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
  points(Year,logBiomass,pch="+",cex=0.5)
}

# Define global parameters: NT = number of years to forecast; Nmc = number of MC sims
NT <- 20
Nmc <- 1000
time1 = Year       ## calibration period
time2 = seq(from=Year[length(Year)],to=Year[length(Year)]+(NT-1), by=1)   ## forecast period
time = c(time1, time2)
N.cols <- c("grey50","red","green","blue","orange") ## set colors
trans <- 0.8       ## set transparancy

# Define forecasting function
##` @param IC   Initial Conditions
##` @param n    Size of Monte Carlo ensemble
##` @param Q    Variance of prediction
forecastRW <- function(IC,Q,n=Nmc){
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    N[,t] <- rnorm(n, Nprev, Q)     ## predict next step
    Nprev <- N[,t]                  ## update IC
  }
  return(N)
}

# Deterministic forecast
det <- forecastRW(IC=mean(rw.out[,ncol(rw.out)]),
          Q=0,
          n=1)
plot.run()
lines(time2, det, col="purple", lwd=3)

# Depends on Q
det <- forecastRW(IC=mean(rw.out[,ncol(rw.out)]),
                  Q=sd(rw.out[,ncol(rw.out)]),
                  n=1)
plot.run()
lines(time2, det, col="purple", lwd=3)

# Depends on IC
prow = sample.int(nrow(rw.out),Nmc,replace=TRUE)

IC.model <- forecastRW(IC=rw.out[prow,ncol(rw.out)],
                  Q=0,
                  n=Nmc)
plot.run()
IC.ci = apply(IC.model,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,IC.ci[1,],IC.ci[3,],col=ecoforecastR::col.alpha(N.cols[1],trans))
lines(time2,IC.ci[2,],lwd=0.5)

# Depends on IC and Q
prow = sample.int(nrow(rw.out),Nmc,replace=TRUE)

IC.Q.model <- forecastRW(IC=rw.out[prow,ncol(rw.out)],
                       Q=sd(rw.out[,ncol(rw.out)]),
                       n=Nmc)
plot.run()
IC.Q.ci = apply(IC.Q.model,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,IC.Q.ci[1,],IC.Q.ci[3,],col=ecoforecastR::col.alpha(N.cols[1],trans))
lines(time2,IC.Q.ci[2,],lwd=0.5)

detach(year.lobstah.data)


############################################
# Implement hierarchical random walk model #
############################################

zone.lobstah.data <- aggregate(cbind(logBiomass, Depth, SeasonalSST) ~ Year + ZONEID, 
                                lobstah %>% filter(Season == "SPRING"), 
                               mean)

attach(zone.lobstah.data)

# Define the model
HigherRW = "
model{
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }

  #### Random Effect
  for(g in 1:ng){
    alpha_site[g] ~ dnorm(0,tau_site)
  }
  
  #### Process Model
  for(t in 2:n){
    mu[t] <- x[t-1] + alpha_site[zoneid[t]]
    x[t] ~ dnorm(mu[t], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  tau_site ~ dgamma(0.1,0.1)  ## site random effect precision
}
"

# Define data
hrw.data <- list(y=logBiomass, 
                 n = logBiomass %>% length(), 
                 ng = ZONEID %>% unique() %>% length(), 
                 zoneid=ZONEID %>% as.numeric(),
                 x_ic=2.5,
                 tau_ic=1,a_obs=1.5,r_obs=1,a_add=1.5,r_add=1)

# Define initial values
nchain = 4
init <- list()
for(i in 1:nchain){
  y.samp = sample(logBiomass,length(logBiomass),
                  replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=1/var(log(y.samp)))
}

# Run models
hrw.model <- jags.model(data = hrw.data, 
                        file = textConnection(HigherRW), 
                        n.chains = 4,
                        n.adapt = 1000)

hrw.samp <- coda.samples (model = hrw.model,
                           variable.names = c("x","alpha_site","tau_add","tau_obs","tau_site"),
                           n.iter = 10000)

burnin = 1000
hrw.burn <- window(hrw.samp, burnin)


###########################################
# Forecast hierarchical random walk model #
###########################################

# Store model runs as matrices
hrw.params <- as.matrix(hrw.burn)
hrw.x <- hrw.params[ ,11:ncol(hrw.params)]
hrw.params <- hrw.params[, 1:10]

# Define global parameters: NT = number of years to forecast; Nmc = number of MC sims
NT = 20
NG = ZONEID %>% unique() %>% length()
Nmc = 1000
time1 = Year       ## calibration period
time2 = seq(from=Year[length(Year)]+1,to=Year[length(Year)]+NT, by=1)   ## forecast period
time = c(time1, time2)
IC = hrw.x[sample.int(nrow(hrw.x),Nmc,replace=TRUE),]
hrw.alpha = hrw.params[,1:NG]
alpha = hrw.alpha[sample.int(nrow(hrw.alpha),Nmc,replace=TRUE),]

# Define forecasting function
##` @param IC     Initial Conditions
##` @param Nmc    Size of Monte Carlo ensemble
##` @param Q      Variance of prediction
##` @param alpha  Vector of group parameters
##` @param NT     Number of years to forecast
##` @param groups List of groups
forecastHRW <- function(IC,Q,Nmc,NT,alpha,groups){
  N <- matrix(NA,Nmc,NT)                    ## storage
  mu <- matrix(NA,Nmc,NT)                   ## storage
  Nprev <- IC[ , ncol(IC)]                  ## initialize
  
  for(t in 1:NT){
    for(n in 1:ncol(IC)){
      mu[,t] <- Nprev + alpha[ ,groups[n]]  ## Define mean as depending on group
      N[,t] <- rnorm(Nmc, mu[ ,t], Q)       ## predict next step
      Nprev <- N[,t]                        ## update IC
    }
  }
  return(N)
}

# Depend on global SD
hrw.fore <- forecastHRW(IC = IC, Q = sd(hrw.x[,ncol(hrw.x)]), Nmc=Nmc, 
                        alpha = alpha, NT = NT, groups = ZONEID %>% as.numeric())

detach(zone.lobstah.data)
