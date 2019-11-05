########## ########## Vergleich der Simulationen in WinBUGS und STAN bei einfachen Modellen ########## ########## 
########## Test mit Stan und Normalverteilung ##########



##### Clear data
rm(list=ls())



#### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Ergebnisse/zusätzliche Experimente/Vergleich der Simulationen bei einfachen Modellen")


 
#### Requiering stan
library("rstan")
library("rstantools")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(1))
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

inits1 <- list(a=0, b=0, c=0, c=0, d=0, e=0, f=0, g=0, h=0)
inits2 <- list(a=1, b=1, c=1, c=1, d=1, e=1, f=1, g=1, h=1)
inits3 <- list(a=-1, b=-1, c=-1, c=-1, d=-1, e=-1, f=-1, g=-1, h=-1)

all.inits <- list(inits1, inits2, inits3)

#data_list <- list()
m <- stan_model('Test Stan NV.stan')



# Simulation 
stan_samples <- sampling(m, iter=20000, init=all.inits, verbose=T, chain=3, warmup= 10000, control = list(max_treedepth = 12, adapt_delta = 0.90)) # !! iter nachher erhöhen




Stan_summary <- summary(stan_samples, pars = c("a", "b", "c", "d", "e", "f", "g", "h"))$summary 
Stan_summary





########## ########## Simulation beendet ########## ########## 