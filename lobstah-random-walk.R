library(tidyverse)
library(rjags)


# Read in data
fish.data <- readRDS("~/Box Sync/Courses/NEFI/Allyn_EcoForecastingProjectData.rds")

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
data <- list(y=mean.maine.lobstah$Biomass,
             n=length(mean.maine.lobstah$Biomass),
             x_ic=20,
             tau_ic=0.5,a_obs=1,r_obs=1,a_add=1,r_add=1)

# Define initial values
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(mean.maine.lobstah$Biomass,length(mean.maine.lobstah$Biomass),
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
out <- as.matrix(jags.burn)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot(mean.maine.lobstah$Year,ci[2,],type='n',ylim=c(0,50))
ecoforecastR::ciEnvelope(mean.maine.lobstah$Year,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(mean.maine.lobstah$Year,mean.maine.lobstah$Biomass,pch="+",cex=0.5)
