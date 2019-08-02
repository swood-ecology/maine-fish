library(tidyverse)
library(rjags)


################################################
# Implement hierarchical, dynamic linear model #
################################################

# Read data
lobstah.data <- 
  readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/LobsterEcoForecastingProjectData.rds")

lobstah.data <- lobstah.data %>% select(Year:Depth, ZONEID)
lobstah.data$logBiomass <- log(lobstah.data$Biomass + 0.000001) 

lobstah <- lobstah.data %>%
  filter(logBiomass > mean(lobstah.data$logBiomass) - 
           2*sd(lobstah.data$logBiomass))
lobstah$numZone <- lobstah$ZONEID %>% as.numeric()

mean.lobstah.data <- aggregate(cbind(logBiomass, Depth, SeasonalSST) ~ Year + ZONEID, 
                               lobstah %>% filter(Season == "SPRING"), 
                               mean)

attach(mean.lobstah.data)

# Define the model
HigherDLM = "
  model{
  #### Data Model
  for(i in 1:n){
    y[i] ~ dnorm(x[i],tau_obs)
  }

  #### Random Effect
  for(g in 1:ng){
    alpha_site[g] ~ dnorm(0,tau_site)
  }

  #### Process Model
  for(i in 2:n){
    mu[i] <- x[i-1] + alpha_site[zoneid[i]] + beta*sst[i]
    x[i] ~ dnorm(mu[i], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  beta ~ dnorm(0,1)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  tau_site ~ dgamma(0.1,0.1)  ## site random effect precision
}
"

h.dlm.data <- list(y = logBiomass, sst = SeasonalSST,
                 n = logBiomass %>% length(), 
                 ng = ZONEID %>% unique() %>% length(), 
                 zoneid = ZONEID %>% as.numeric(), 
                 x_ic = 2.5, tau_ic = 0.5, 
                 a_obs = 1, r_obs = 1, a_add = 1, r_add = 1)

h.dlm.model <- jags.model(data = h.dlm.data, 
                          file = textConnection(HigherDLM),
                          n.chains = 4)

h.dlm.out <- coda.samples (model = h.dlm.model,
                           variable.names = c("alpha_site","beta", 
                                              "tau_add","tau_obs","tau_site",
                                              "x"),
                           n.iter = 10000)

gelman.plot(h.dlm.out)

burnin = 5000
hrw.burn <- window(h.dlm.out, burnin)
summary(hrw.burn)


# Fit DLM model with ecoforecastR
ef.out <- ecoforecastR::fit_dlm(
  model = list(obs = "logBiomass", fixed = "~ 1 + X + SeasonalSST"),
  mean.lobstah.data
)

## parameter diagnostics
params <- window(ef.out$params, start = 1000) ## remove burn-in
summary(params)
cor(as.matrix(params))

## confidence interval
dlm.out <- as.matrix(ef.out$predict)
dlm.ci <- apply(dlm.out, 2, quantile, c(0.025, 0.5, 0.975))
plot(mean.lobstah.data$Year,
  dlm.ci[2, ],
  type = "n"
)
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(mean.lobstah.data$Year, dlm.ci[1, ], dlm.ci[3, ], col = ecoforecastR::col.alpha("lightBlue", 0.75))
points(mean.lobstah.data$Year, mean.lobstah.data$logBiomass, pch = "+", cex = 0.5)


########################
# Fit random slope DLM #
########################

# Define the model
RSDLM = "
model{
  #### Data Model
  for(i in 1:n){
    y[i] ~ dnorm(x[i],tau_obs)
  }

  #### Random Effect
  for(g in 1:ng){
    beta_site[g] ~ dnorm(0.2, 0.1)
  }

  #### Process Model
  for(i in 2:n){
    mu[i] <- x[i-1] + alpha + beta_site[zoneid[i]]*sst[i]
    x[i] ~ dnorm(mu[i], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(0.1,0.1)
  tau_add ~ dgamma(0.1,0.1)
  alpha ~ dnorm(2, 0.1)
}
"

h.dlm.data <- list(y = logBiomass, sst = SeasonalSST,
                   n = logBiomass %>% length(), 
                   ng = ZONEID %>% unique() %>% length(),
                   zoneid = ZONEID %>% as.numeric(), 
                   x_ic = 2.5, tau_ic = 0.1)

h.dlm.model <- jags.model(data = h.dlm.data, 
                          file = textConnection(RSDLM),
                          n.chains = 4)

h.dlm.out <- coda.samples (model = h.dlm.model,
                           variable.names = c("alpha","beta_site",
                                              "tau_add","tau_obs","x"),
                           n.iter = 10000)
gelman.plot(h.dlm.out)

burnin = 6000
hrw.burn <- window(h.dlm.out, burnin)
summary(hrw.burn)

detach(mean.lobstah.data)


##################################
# Run for all data, not averages #
##################################
RSDLM.AD = "
model{
  #### Data Model
  for(i in 1:n){
    y[i] ~ dnorm(x[i],tau_obs)
  }

  #### Random slope
  for(g in 1:ng){
    beta_site[g] ~ dnorm(0.2, 0.1)
  }

  #### Year-based random walk value
  for(t in year1:(year1+ny)){
    x_year[t] ~ dnorm(0.2, 0.1)
  }

  #### Process Model
  for(i in 2:n){
    x[i] ~ dnorm(mu[i], tau_add)  
    mu[i] <- x_year[year[i-1]] + alpha + beta_site[zoneid[i]]*sst[i]
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(0.1,0.1)
  tau_add ~ dgamma(0.1,0.1)
  alpha ~ dnorm(2, 0.1)
}
"

attach(lobstah)
h.dlm.data <- list(y = logBiomass, sst = SeasonalSST,
                   n = logBiomass %>% length(), 
                   ng = ZONEID %>% unique() %>% length(),
                   ny = Year %>% unique() %>% length(),
                   year1 = Year %>% as.numeric() %>% min(),
                   year = Year %>% as.numeric(),
                   zoneid = numZone, 
                   x_ic = 2.5, tau_ic = 0.1)

h.dlm.model <- jags.model(data = h.dlm.data, 
                          file = textConnection(RSDLM.AD),
                          n.chains = 4)

h.dlm.out <- coda.samples (model = h.dlm.model,
                           variable.names = c("alpha","beta_site",
                                              "tau_add","tau_obs","x_year"),
                           n.iter = 10000)
gelman.plot(h.dlm.out)

burnin = 6000
hrw.burn <- window(h.dlm.out, burnin)
summary(hrw.burn)

out <- as.matrix(hrw.burn)[,2:8]

hist <- data.frame(out[,1],
           out[,2],
           out[,3],
           out[,4],
           out[,5],
           out[,6],
           out[,7]
)

ggplot(reshape::melt(hist)) + geom_freqpoly(aes(x = value,
                                 y = ..density.., colour = variable)) +
  theme_bw()

detach(lobstah)

