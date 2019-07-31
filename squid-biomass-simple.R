library(readr)
library(dplyr)
library(ggplot2)
library(rjags)

all_fish <- readRDS("Allyn_EcoForecastingProjectData.rds")

#surface temperature, lobster, abundance across data type
#group by year, season, species

squid <- all_fish %>% 
  filter(CommonName == "LONGFIN SQUID") %>%
  group_by(Year) %>%
  summarize(mean_biomass = mean(Biomass))

ggplot(squid, aes(x = Year, y = mean_biomass)) +
  geom_line()
hist(squid$mean_biomass)

#prior has gamma distribution?
squid_model <- "
model {
  mu ~ dnorm(mu0, T) 
  for(i in 1:n){
    x[i] ~ dnorm(mu, s)
  }
}
"

data <- list(x = squid$mean_biomass, n = length(squid), mu0 = 5, 
             s = 0.01, T = 0.00001)

inits <- list()
inits[[1]] <- list(mu = 10)
inits[[2]] <- list(mu = 15)
inits[[3]] <- list(mu = 8)

model <- jags.model(file = textConnection(squid_model), 
                    data = data, 
                    inits = inits, 
                    n.chains = 3)

out <- coda.samples(model = model, 
                    variable.names = "mu", 
                    n.iter = 10000)

gelman.diag(out)
gelman.plot(out)

plot(out)
summary(out)
