library(tidyverse)
library(rjags)


# Read in data
fish.data <- readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/Allyn_EcoForecastingProjectData.rds")

# Filter out lobster in Maine
maine.lobstah <- fish.data %>% 
  select(-Survey, -X) %>%
  filter(CommonName == "AMERICAN LOBSTER") %>%
  filter(Y > 42.874)

# Plot all data
ggplot(maine.lobstah, aes(y = Biomass, x = Year)) + 
  geom_point(alpha = 0.5) + 
  ylim(0,500) +
  theme_bw()

# Take annual average
mean.maine.lobstah <- aggregate(cbind(Biomass, Depth, SeasonalSST) ~ Year, 
                                maine.lobstah, mean)

# Plot average biomass trend
ggplot(mean.maine.lobstah, aes(x = Year, y = Biomass)) + 
  geom_line() + 
  theme_bw()


######################################
# Implement simple random walk model #
######################################

attach(mean.maine.lobstah)

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
data <- list(y=Biomass,
             n=length(Biomass),
             x_ic=20,
             tau_ic=0.5,a_obs=1,r_obs=1,a_add=1,r_add=1)



# Define initial values
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(Biomass,length(Biomass),
                  replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
}

# Run model
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
gelman.plot(jags.out)

burnin = 1000                                ## determine convergence
jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
summary(jags.burn) ## check diagnostics post burn-in

# CI plot
rw.out <- as.matrix(jags.burn)

x.cols <- grep("^x",colnames(rw.out)) ## grab all columns that start with the letter x
ci <- apply(rw.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(mean.maine.lobstah$Year,ci[2,],type='n',ylim=c(0,50))
ecoforecastR::ciEnvelope(mean.maine.lobstah$Year,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(mean.maine.lobstah$Year,mean.maine.lobstah$Biomass,pch="+",cex=0.5)


#####################################
# Forecast simple random walk model #
#####################################

# Define function to print
plot.run <- function(){
  x.cols <- grep("^x",colnames(rw.out)) ## grab all columns that start with the letter x
  ci <- apply(rw.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
  plot(Year,ci[2,],type='n',ylim=c(0,75),xlim=c(time[1], time[length(time)]))
  ecoforecastR::ciEnvelope(Year,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
  points(Year,Biomass,pch="+",cex=0.5)
}

# Define global parameters: NT = number of years to forecast; Nmc = number of MC sims
NT <- 20
Nmc <- 1000
time1 = Year       ## calibration period
time2 = seq(from=Year[length(Year)],to=Year[length(Year)]+(NT-1), by=1)   ## forecast period
time = c(time1, time2)

# Define forecasting function
##` @param IC   Initial Conditions
##` @param n    Size of Monte Carlo ensemble
##` @param Q    Variance of prediction
forecastRW <- function(IC,Q,n=Nmc){
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    N[,t] <- rnorm(n,Nprev,Q)       ## predict next step
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
ecoforecastR::ciEnvelope(time2,IC.ci[1,],IC.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,IC.ci[2,],lwd=0.5)

# Depends on IC and Q
prow = sample.int(nrow(rw.out),Nmc,replace=TRUE)

IC.Q.model <- forecastRW(IC=rw.out[prow,ncol(rw.out)],
                       Q=sd(rw.out[,ncol(rw.out)]),
                       n=Nmc)
plot.run()
IC.Q.ci = apply(IC.Q.model,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,IC.Q.ci[1,],IC.Q.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,IC.Q.ci[2,],lwd=0.5)

detach(mean.maine.lobstah)

############################################
# Implement hierarchical random walk model #
############################################

# Read data
lobstah.data <- 
  readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/LobsterEcoForecastingProjectData.rds")

lobstah.data <- lobstah.data %>% select(Year:Depth, ZONEID)

mean.lobstah.data <- aggregate(cbind(Biomass, Depth, SeasonalSST) ~ Year + ZONEID, 
                                lobstah.data, mean)

attach(mean.lobstah.data)

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

hrw.data <- list(y = Biomass, 
                 n = Biomass %>% length(), 
                 ng = ZONEID %>% unique() %>% length(), 
                 zoneid=ZONEID %>% as.numeric(),
                 x_ic = 20, tau_ic = 0.5, a_obs = 1, r_obs = 1, a_add = 1, r_add = 1)

hrw.model <- jags.model(data = hrw.data, 
                          file = textConnection(HigherRW), 
                          n.chains = 3, 
                          n.adapt = 2000)

hrw.out <- coda.samples (model = hrw.model,
                           variable.names = c("tau_add","tau_obs","tau_site","alpha_site"),
                           n.iter = 10000)

hrw.out.x <- coda.samples (model = hrw.model,
                         variable.names = c("x"),
                         n.iter = 10000)


