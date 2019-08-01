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
  filter(logBiomass > mean(lobstah.data$logBiomass) - 2*sd(lobstah.data$logBiomass))

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
    mu[i] <- x[i-1] + alpha_site[zoneid[i]] + beta1*x1[i] + beta2*x2[i]
    x[i] ~ dnorm(mu[i], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  beta1 ~ dnorm(0,1)
  beta2 ~ dnorm(0,1)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  tau_site ~ dgamma(0.1,0.1)  ## site random effect precision
}
"

h.dlm.data <- list(y = logBiomass, x1 = Depth, x2 = SeasonalSST,
                 n = logBiomass %>% length(), 
                 ng = ZONEID %>% unique() %>% length(), 
                 zoneid = ZONEID %>% as.numeric(), 
                 x_ic = 2.5, tau_ic = 0.5, 
                 a_obs = 1, r_obs = 1, a_add = 1, r_add = 1)

h.dlm.model <- jags.model(data = h.dlm.data, 
                          file = textConnection(HigherDLM),
                          n.chains = 4)

h.dlm.out <- coda.samples (model = h.dlm.model,
                           variable.names = c("alpha_site","beta1", "beta2",
                                              "tau_add","tau_obs","tau_site",
                                              "x"),
                           n.iter = 10000)

gelman.plot(h.dlm.out)

burnin = 1000
hrw.burn <- window(h.dlm.out, burnin)


# Fit DLM model with ecoforecastR
ef.out <- ecoforecastR::fit_dlm(
  model = list(obs = "logBiomass", fixed = "~ 1 + X + Depth + SeasonalSST"),
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

detach(mean.lobstah.data)
