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



#Data
#l=1
#m=1
#n=1
o=1
p=1
q=1
data_list= list( o=o, p=p, q=q)

# Initialwerte
inits1 <- list(a=0, b=0, c=0, c=0, d=2, e=2, f=2, g=2, h=2, i=2, j=2, k=2, o=2, p=2, q=2 )
inits2 <- list(a=1, b=1, c=1, c=1, d=1, e=1, f=1, g=1, h=1, i=1, j=1, k=1, o=1, p=1, q=1 )
inits3 <- list(a=1, b=1, c=1, c=1, d=1, e=1, f=1, g=1, h=1, i=1, j=1, k=1, o=1, p=1, q=1 )
#, l=1, m=2, n=2

all.inits <- list(inits1, inits2, inits3)
all.inits <- list(inits3)


#data_list <- list()
m <- stan_model('weitere Verteilungen.stan')


# Simulation 
stan_samples <- sampling(m, data=data_list, init=all.inits,iter=30000,  verbose=TRUE, chain=1, warmup= 15000, control = list(max_treedepth = 12, adapt_delta = 0.90)) 


#pairs(stan_samples, pars = c("a", "b", "c", "d", "e", "f", "g",  "lp__"), las = 1)
# bringt kaum was

Stan_summary <- summary(stan_samples, pars = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k" ))$summary 
Stan_summary
# Error in check_pars(allpars, pars) : no parameter o, p, q




########## ########## Simulation beendet ########## ########## 