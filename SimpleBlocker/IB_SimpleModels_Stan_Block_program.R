########## ########## Vergleich des einfachen Blocker Beispiels in WinBUGS und STAN  ########## ########## 




##### Clear data
rm(list=ls())



#### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Ergebnisse/zusätzliche Experimente/einfacher Blocker Vergleich WinBUGS und STAN")



#### Requiering stan
library("rstan")
library("rstantools")
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores(1))
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



#### Assignment data to stan
rt = c(3, 7, 5, 102, 28, 4, 98, 60, 25, 138, 64, 45, 9, 57, 25, 33, 
    28, 8, 6, 32, 27, 22)
nt =  c(38, 114, 69, 1533, 355, 59, 945, 632, 278, 1916, 873, 263, 
    291, 858, 154, 207, 251, 151, 174, 209, 391, 680)
rc = c(3, 14, 11, 127, 27, 6, 152, 48, 37, 188, 52, 47, 16, 45, 31, 
    38, 12, 6, 3, 40, 43, 39)
nc =c(39, 116, 93, 1520, 365, 52, 939, 471, 282, 1921, 583, 266, 
    293, 883, 147, 213, 122, 154, 134, 218, 364, 674)
N = 22


data_list <- list(N=N, rt=rt, nt=nt, rc=rc, nc=nc)



#### Read in inits

inits1 <- function(chain_id = 1) {
  list(d <- 0,
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               sigmasq_delta <- 1,
               mu <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta_new <- 0)
}

inits2 <- function(chain_id = 2) {
  list(d <- 0,
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               sigmasq_delta <- 1,
               mu <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta_new <- 0)
}

inits3 <- function(chain_id = 3) {
  list(d <- 0,
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               sigmasq_delta <- 1,
               mu <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               delta_new <- 0)
}

all.inits <- list(inits1, inits2, inits3)



#### Compiling
m <- stan_model('IB_SimpleModels_Stan_Block_Stan.stan')



#### Simulation 
stan_samples <- sampling(m, data=data_list, iter=20000, init="all.inits", verbose=T, chain=3, warmup= 10000, control = list(max_treedepth = 10, adapt_delta = 0.85)) # !! iter nachher erhöhen


#### Summary
Stan_summary <- summary(stan_samples, pars = c("d", "delta_new", "sigma_delta", "sigmasq_delta" ))$summary 
Stan_summary





########## ########## Simulation beendet ########## ########## 