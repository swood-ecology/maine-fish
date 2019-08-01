library(tidyverse)
library(rjags)


# Read in data
fish.data <- readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/LobsterEcoForecastingProjectData.rds")

# Filter out lobster in Maine
maine.lobstah <- fish.data %>%
  select(Year:Depth, ZONEID)

# Plot all data
ggplot(maine.lobstah, aes(y = Biomass, x = Year)) +
  geom_point(alpha = 0.5) +
  ylim(0, 500) +
  theme_bw()

# Take annual average
mean.maine.lobstah <- aggregate(
  cbind(Biomass, Depth, SeasonalSST) ~ Year+ZONEID+Season,
  maine.lobstah, mean
)

# Fit DLM model
## fit the model
ef.out <- ecoforecastR::fit_dlm(
  model = list(obs = "Biomass", fixed = "~ 1 + X + Depth + SeasonalSST"),
  mean.maine.lobstah
)

## parameter diagnostics
params <- window(ef.out$params, start = 1000) ## remove burn-in
summary(params)
cor(as.matrix(params))

## confidence interval
dlm.out <- as.matrix(ef.out$predict)
dlm.ci <- apply(dlm.out, 2, quantile, c(0.025, 0.5, 0.975))
plot(mean.maine.lobstah$Year,
  dlm.ci[2, ],
  type = "n"
)
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(mean.maine.lobstah$Year, dlm.ci[1, ], dlm.ci[3, ], col = ecoforecastR::col.alpha("lightBlue", 0.75))
points(mean.maine.lobstah$Year, mean.maine.lobstah$Biomass, pch = "+", cex = 0.5)


# Take seasonal average
zone.lobstah <- aggregate(
  cbind(Biomass) ~ Year + ZONEID,
  maine.lobstah, mean
)

# Plot average biomass trend
ggplot(aggregate(
  cbind(Biomass) ~ Year + ZONEID,
  maine.lobstah, mean
), 
aes(x = Year, y = Biomass, facet = ZONEID)) +
  geom_line(aes(color = ZONEID)) +
  theme_bw()



lobstah_rand <- "
model{                                                                                                 
                                                                                                       
  #### Priors                                                                                            
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)

  #### Random Effects                                                                                    
  tau_alpha~dgamma(0.1,0.1)                                                                     
  for(i in 1:nseasons){                                                                             
    alpha[i]~dnorm(0,tau_alpha)                                                                  
  }

  #### Fixed Effects                                                                                     
  betaX ~dnorm(0,0.001)                                                                           
  betaIntercept~dnorm(0,0.001)                                                                       
  betaDepth~dnorm(0,0.001)                                                                            
  betaSeasonalSST~dnorm(0,0.001)                                                                      
  for(j in 1: 3 ){                                                                                      
    muXf[j] ~ dnorm(0,0.001)                                                                             
    tauXf[j] ~ dgamma(0.01,0.01)                                                                         
  }                                                                                                       

  #### Data Model                                                                                        
  for(t in 1:n){                                                                                         
    OBS[t] ~ dnorm(x[t],tau_obs)                                                                         
    Xf[t,1] ~ dnorm(muXf[1],tauXf[1])                                                                   
    Xf[t,2] ~ dnorm(muXf[2],tauXf[2])                                                                        
    Xf[t,3] ~ dnorm(muXf[3],tauXf[3])                                                                        
  }                                                                                                      

  #### Process Model                                                                                    
  for(t in 2:n){                                                                                         
    mu[t] <- x[t-1]  + betaX*x[t-1] + betaIntercept*Xf[t,1] + betaDepth*Xf[t,2] + betaSeasonalSST*Xf[t,3]
    x[t]~dnorm(mu[t],tau_add)                                                                            
  }                                                                                                      
}   
"

# Define data and parameters
season.data <- list(
  y = season.maine.lobstah$Biomass,
  nrep = unique(season.maine.lobstah$Season),
  n = length(season.maine.lobstah$Biomass),
  x_ic = 20,
  tau_ic = 0.5, a_obs = 1, r_obs = 1, a_add = 1, r_add = 1
)

# Define initial values
nchain <- 3
init <- list()
for (i in 1:nchain) {
  y.samp <- sample(mean.maine.lobstah$Biomass, length(mean.maine.lobstah$Biomass),
    replace = TRUE
  )
  init[[i]] <- list(tau_add = 1 / var(diff(log(y.samp))), tau_obs = 5 / var(log(y.samp)))
}

# Run model
j.model <- jags.model(
  file = textConnection(RandomWalk),
  data = data,
  inits = init,
  n.chains = 3
)


#########################################################

# Read in lobster data with management groupings
lobstah.spatial <- readRDS("~/Box Sync/Courses/NEFI/maine-fish/data/LobsterEcoForecastingProjectData.rds")

# Select only needed columns
lobstah.spatial %>%
  select(Year:Depth, ZONEID, Long, Lat) -> lobstah.spatial
