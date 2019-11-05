########## ########## MCMC Comparison of the different samplers and their packages ########## ########## 
########## using the Blocker example with the Random Effects Model



########## testing running time
#### to test the running time, the Samplers are running separately and the package "tictoc" are used


########## testing other diagnostica
#### because of development of compareMCMCs, similar functions from other packages e.g. coda-package are used




#### Clear data
rm(list=ls())



#### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Ergebnisse/zusätzliche Experimente/MCMC Vergleiche")



### loading libraries
library(R2WinBUGS)
library(rjags) 
library(jagsUI)
library(R2jags)
library(runjags)
library(nimble)
#library("rstan")
#library("rstantools")

library(coda)  
library(tictoc)
library(lattice) 
library(random)
library(matrixStats) # zus?tzl Paket, berechnet Median
library(car)
library(fitR)
library(mcmcr)
#library(BayesianTools) 

# load RNGs
load.module("glm") 
load.module("lecuyer")
#list.factories(type="rng")
#parallel.seeds("lecuyer::RngStream", 5);
load.module("dic")



##### Read the data into R.
#data = read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F)
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
head(data) # Shows the first six entries
data2  = read.table("Data_Blocker_Rest.txt")
head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
# ns # check
nt <-  ncol(data[,5:6])
#nt # check
na <- data[,7]
# na  # check
r <- data[,1:2]
# r # Check
n <-  data[,3:4]
# n # Check
t <-  data[,5:6]
# t # Check

dat <- list("ns", "nt", "na", "r", "n", "t")  # names list of numbers 
dat2 <- list(ns=ns, nt=nt, na=na, r=r, n=n, t=t) 



##### Parameter to monitor/save
# analog https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/some_BUGS_comparisons/blocker/nimble_blocker_comparisons.html,
# all known parameters were monitored
params <- c("mu", "d", "sd", "A" )



##### read in inits with chains
inits1 <- list(d=c( NA, 0), 
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=2)

inits3 <- list(d=c( NA, 2), 
               sd=2, 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=3 ) 

all.inits <- list(inits1, inits2, inits2)



##### model Code
# the model code is named as "Blocker_Model_Random.txt"





# Part 1: Testing running time of R code of WinBUGS ----------------------------------




tic("Total WinBUGS")
tic("Compiling and Simulation/Sampling")
test_WinBUGS <- bugs(data=dat, inits=all.inits, parameters=params, model.file="Blocker_Model_Random.bug",
                    n.chains=3, n.iter=20000, n.burnin =10000,
                    bugs.directory="C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/WinBUGS/WinBUGS - Programm und Skript/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14/"
)
toc()
tic("generating summary")
print(test_WinBUGS)
#plot(test_WinBUGS)
toc()
toc()





# Part 2: Testing running time of R code of JAGS ----------------------------------




##### ##### jagsUI package ##### ##### 


##### Inits for jagsUI
inits1 <- list(d=c( NA, 0), 
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               A=0,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               A=1,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=2)

inits3 <- list(d=c( NA, 2), 
               sd=2, 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
               A=2,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=3 ) 

# Note: "lecuyer::RngStream" does not work with jagsUI

all.inits.jagsUI <- list(inits1, inits2, inits2)



tic("Total jagsUI")
tic("Compiling and Simulation/Sampling")
jags.m.jagsUI <- jags(data=dat, 
               inits=all.inits,
               parameters.to.save=params,
               model.file="Blocker_Model_Random.txt", 
               n.chains=3,
               n.iter=20000, n.burnin=10000) # store.data=TRUE
toc()
tic("generating summary")
jags.View(jags.m.jagsUI)
jags.m.jagsUI
toc()
toc()




##### ##### R2jags package ##### ##### 


tic("Total R2jags")
tic("Compiling and Simulation/Sampling")
jags.m.R2jags <- jags(data=dat,  inits=all.inits, parameters.to.save=params,  n.chains = 3, n.iter = 30000, n.burnin = 10000,
               model.file="Blocker_Model_Random.txt", DIC=TRUE,  jags.module = c("glm","dic") )  
toc()
tic("generating summary")
print(jags.m.R2jags)
jags.m.R2jags[["BUGSoutput"]][["median"]]
toc()
toc()




##### ##### rjags package ##### ##### 


##### read in inits with chains for rjags
inits.rjags <- function(){
  #chain 1
  list(d=c( NA, 0), 
       sd=1, 
       mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
       .RNG.name="lecuyer::RngStream", .RNG.seed=1
  )  
  #chain 2
  list(d=c( NA, -1), 
       sd=4, 
       mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
       .RNG.name="lecuyer::RngStream", .RNG.seed=2 
  )
  #chain 3
  list(d=c( NA, 2), 
       sd=2, 
       mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
       .RNG.name="lecuyer::RngStream", .RNG.seed=3
  ) 
}


tic("Total rjags")
tic("Compiling")
jags.m.rjags <- jags.model( file = "Blocker_Model_Random.txt", data=dat2, inits=inits.rjags, n.chains=3, n.adapt=1500)
update(jags.m.rjags, 10000) # burn in 
toc()
tic("Simulation/Sampling")
samps_coda <- coda.samples( jags.m.rjags, variable.names=c("mu", "d", "sd", "A"),   n.iter=30000, DIC=T  ) 
toc()
tic("generating summary")
summary(window(samps_coda, start=10001)) # burnin = 10000
#### Median Calculation
# Median for T1
median(as.matrix(samps_coda[,1]))
# Median for T2
median(as.matrix(samps_coda[,2]))
# Median for d[2]
median(as.matrix(samps_coda[,3]))
# Median for sd
median(as.matrix(samps_coda[,4]))
toc()
toc()




##### ##### runjags package ##### ##### 


tic("Total runjags")
tic("Compiling and Simulation/Sampling")
jags.m.runjags <- run.jags(model="Blocker_Model_Random.txt", monitor=c("mu", "d", "sd", "A"  ),
                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , n.chains=3, inits=all.inits, burnin = 10000, sample = 30000, adapt = 1500)
toc()
tic("generating summary")
print(jags.m.runjags)
summary(jags.m.runjags$mcmc) # f?r 2.5 - 97.5  CrI
toc()
toc()





# Part 3: Testing running time of R code of NIMBLE ----------------------------------




##### Adaptions for NIMBLE

Nimble_constants = list(ns=ns, nt=nt, na=na, t=t)
Nimble_data = list(r=r, n=n)
Nimble_inits = list(d = c( NA, 0),
                    mu = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    # zus?tzlich noch einzuf?gen:
                    A=0)
Code_Model<- nimbleCode( {
  for(i in 1:ns){								 # LOOP THROUGH STUDIES 
    w[i,1] <- 0 								# adjustment for multi-arm trials is zero for control arm 
    delta[i,1] <- 0 								# treatment effect is zero for control arm 
    mu[i] ~ dnorm(0,.0001) 				# vague priors for all trial baselines 
    for (k in 1:na[i]) { 								# LOOP THROUGH ARMS 
      r[i,k] ~ dbin(p[i,k],n[i,k]) 				# binomial likelihood 
      logit(p[i,k]) <- mu[i] + delta[i,k]				 # model for linear predictor 
      rhat[i,k] <- p[i,k] * n[i,k] 				# expected value of the numerators 
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))+  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))  #Deviance contribution 
    } 
    resdev[i] <- sum(dev[i,1:na[i]]) 				# summed residual deviance contribution for this trial 
    for (k in 2:na[i]) { 								# LOOP THROUGH ARMS 
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) 				# trial-specific LOR distributions 
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] 				# mean of LOR distributions (with multi-arm trial correction) 
      taud[i,k] <- tau *2*(k-1)/k 								# precision of LOR distributions (with multi-arm trial correction) 
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) 				# adjustment for multi-arm RCTs 
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) 				# cumulative adjustment for multi-arm trials 
    } 
  }    
  totresdev <- sum(resdev[1:22]) 								#Total Residual Deviance 
  d[1] <- 0												 # treatment effect is zero for reference treatment 
  sd ~ dunif(0,5) 												# vague prior for between-trial SD.
  tau <- pow(sd,-2) 								# between-trial precision = (1/between-trial variance) 
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale  
  # Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 
  
  for (k in 2:nt){ d[k] ~ dnorm(0,0.0001) }
  for (k in 1:nt) { logit(T[k]) <- A + d[k] }
  A ~ dnorm(-2.2, 3.3) 
})




##### ##### nimbleModel package ##### ##### 


tic("Total nimbleModel")
tic("Compiling and Simulation/Sampling")
Model_Nimble_nimbleModel <- nimbleModel(code = Code_Model, name = "ProcessedModel", constants = Nimble_constants,
                            data = Nimble_data, inits = Nimble_inits) 
mcmc.out.nimbleModel <- nimbleMCMC(code = Code_Model, constants = Nimble_constants,
                       data = Nimble_data, inits = Nimble_inits,
                       nchains = 3, niter = 20000, nburnin = 10000, 
                       summary = TRUE, WAIC = F,
                       monitors = c("mu", "d", "sd", "A"))
toc()
tic("generating summary")
mcmc.out.nimbleModel[["summary"]][["all.chains"]]
toc()
toc()
mcmc.out.nimbleModel_II <- nimbleMCMC(code = Code_Model, constants = Nimble_constants,
                                   data = Nimble_data, inits = Nimble_inits,
                                   nchains = 3, niter = 20000, nburnin = 10000, 
                                   summary = F, WAIC = F,  samplesAsCodaMCMC = T,
                                   monitors = c("mu", "d", "sd", "A"))




##### ##### readBUGSModel package ##### ##### 


tic("Total readBUGSModel")
tic("Compiling and Simulation/Sampling")
readBUGS_Model <- readBUGSmodel(model='Blocker_Random_Model_Nimble.bug', data = 'Blocker_Data_Nimble.R',
                                inits = 'Blocker_Inits_Nimble.R' )  
mcmc.out.readBUGSModel <- nimbleMCMC(code = readBUGS_Model,
                       nchains = 3, niter = 20000, nburnin = 10000,
                       summary = TRUE, WAIC = F,
                       monitors = c("mu", "d", "sd", "A"))
toc()
tic("generating summary")
mcmc.out.readBUGSModel[["summary"]][["all.chains"]]
toc()
toc()
mcmc.out.readBUGSModel_II <- nimbleMCMC(code = readBUGS_Model,
                                     nchains = 3, niter = 20000, nburnin = 10000,
                                     summary = F, samplesAsCodaMCMC = T, WAIC = F,
                                     monitors = c("mu", "d", "sd", "A"))





# Part 4: Testing further diagnostica with the coda package ---------------




##### ##### Testing of the effective Size ##### ##### 


## of WinBUGS
effectiveSize(test_WinBUGS$sims.matrix) 


## of JAGS
# jagsUI
effectiveSize(jags.m.jagsUI)
# R2jags
effectiveSize(jags.m.R2jags)
# rjags
effectiveSize(samps_coda)
# runjags
effectiveSize(jags.m.runjags)


## of NIMBLE
# nimbleModel
EES1 <- effectiveSize(mcmc.out.nimbleModel[["samples"]][["chain1"]])
EES2 <- effectiveSize(mcmc.out.nimbleModel[["samples"]][["chain2"]])
EES3 <- effectiveSize(mcmc.out.nimbleModel[["samples"]][["chain3"]])
e <- (EES1+EES2+EES3)/3
e
# readBUGSModel
EES1_II <- effectiveSize(mcmc.out.readBUGSModel[["samples"]][["chain1"]])
EES2_II <- effectiveSize(mcmc.out.readBUGSModel[["samples"]][["chain2"]])
EES3_II <- effectiveSize(mcmc.out.readBUGSModel[["samples"]][["chain3"]])
f <- (EES1_II+EES2_II+EES3_II)/3
f





# Part 5: Testing trace- and density plot with the plot function        ---------------------------------------------------------
## additional tests
## the orders have to be activated for seeing single plots

## the plots are generated with the plot()-function if possible
## else they are generated by a separate function

## the files are saved in external pdf files


#par(ask=F)

## of WinBUGS
#pdf(file="Traceplots_WinBUGS.pdf")
#plot(test_WinBUGS, x=max)
#R2jags::traceplot(test_WinBUGS, ask = F) 
#plot(density(test_WinBUGS$sims.matrix)) # better only direct via programm
#dev.off()
#traceplot(test_WinBUGS)


## of JAGS
# jagsUI
#pdf(file="Traceplots_jagsUI.pdf")
#plot(jags.m.jagsUI)
#R2jags::traceplot(jags.m.jagsUI, ask = F) 
#densplot()
#dev.off()

# R2jags
#pdf(file="Traceplots_R2jags.pdf")
#plot(jags.m.R2jags)
#R2jags::traceplot(jags.m.R2jags, ask = F) 
#dev.off()

# rjags
#pdf(file="Traceplots_rjags.pdf")
#plot(samps_coda)
#dev.off()

# runjags
#pdf(file="Traceplots_runjags.pdf")
#plot(jags.m.runjags)
#R2jags::traceplot(jags.m.runjags$mcmc, ask = F) 
#dev.off()


## of nimble
# nimbleModel
#pdf(file="Traceplots_nimbleModel.pdf")
#plot(mcmc.out.nimbleModel_II)
#dev.off()
#R2jags::traceplot(mcmc.out.nimbleModel_II, ask = F) 

# readBUGSModel
#pdf(file="Traceplots_readBUGSModel.pdf")
#plot(mcmc.out.readBUGSModel_II)
#dev.off()









# Part 6: Testing autocorrelation plot with the plot function        ---------------------------------------------------------

## additional tests
## the orders have to be activated for seeing single plots


#par(ask=F)

## of WinBUGS
#pdf(file="autocorrplot_WinBUGS.pdf")
#autocorr.plot(test_WinBUGS$sims.matrix, ask = F)
#dev.off()


## of JAGS
# jagsUI
#pdf(file="autocorrplot_jagsUI.pdf")
#autocorr.plot(jags.m.jagsUI, ask = F)
#dev.off()

# R2jags
#pdf(file="autocorrplot_R2jags.pdf")
#autocorr.plot(jags.m.R2jags, ask = F)
#dev.off()

# rjags
#pdf(file="autocorrplot_rjags.pdf")
#autocorr.plot(samps_coda, ask = F)
#dev.off()

# runjags
#pdf(file="autocorrplot_runjags.pdf")
#autocorr.plot(jags.m.runjags, ask = F)
#dev.off()


## of nimble
# nimblemodel
#pdf(file="autocorrplot_nimblemodel.pdf")
#autocorr.plot(mcmc.out.nimbleModel_II, ask = F)
#dev.off()
# readBUGSmodel
#pdf(file="autocorrplot_nimblemodel.pdf")
#autocorr.plot(mcmc.out.readBUGSModel_II, ask = F)
#dev.off()





# Part 7: Comparison of Rhat ----------------------------------------------



### use of rhat function of mcmcr


## of WinBUGS
rhat(as.mcmc(test_WinBUGS$sims.matrix), by = "term")


## of JAGS
# jagsUI
rhat(as.mcmc(jags.m.jagsUI$BUGSoutput$sims.matrix), by = "term" )

# R2jags
rhat(as.mcmc(jags.m.R2jags$BUGSoutput$sims.matrix), by = "term"  )

# rjags
# no function available
rhat(samps_coda, by = "term", as_df = FALSE)

# runjags
# psrf = rhat
rhat(jags.m.runjags$mcmc, by = "term")


## of nimble
# nimblemodel
rhat(mcmc.out.nimbleModel_II, by = "term", as_df = FALSE)

# readBUGSmodel
rhat(mcmc.out.readBUGSModel_II, by = "term", as_df = FALSE)





# Part 8: Testing with functions of the sbfnk/fitR package  ---------------------------------------------------


## sbfnk/fitR: Tool box for fitting dynamic infectious disease models to time series. 




##### ##### burnAndThin()  function ##### ##### 

## additional tests
## the orders have to be activated for seeing single plots

## this test compares to "mixing" diagnostics and should show the effect of increasing of thinning,
## e.g. to reduce the amount of memory and storage space in long chains. 
## Afterwards, the autocorrelation plots will compared to the unthinned.


# BaT = Burn and Thin


#par(ask=F)


## of WinBUGS
#BaT_BUGS <- burnAndThin(test_WinBUGS$sims.matrix)
#pdf(file="BaT_autocorrplot_WinBUGS.pdf")
#autocorr.plot(BaT_BUGS, ask = F)
#dev.off()


## of JAGS
# jagsUI
#BaT_jagsUI <- burnAndThin(jags.m.jagsUI$BUGSoutput$sims.matrix) 
#pdf(file="BaT_autocorrplot_jagsUI.pdf")
#autocorr.plot(BaT_jagsUI, ask = F)
#dev.off()

# R2jags
#BaT_R2jags <- burnAndThin(jags.m.R2jags$BUGSoutput$sims.matrix) 
#pdf(file="BaT_autocorrplot_R2jags.pdf")
#autocorr.plot(BaT_R2jags, ask = F)
#dev.off()

# rjags
#BaT_rjags <- burnAndThin(samps_coda) 
#pdf(file="BaT_autocorrplot_rjags.pdf")
#autocorr.plot(BaT_rjags, ask = F)
#dev.off()

# runjags
#BaT_runjags <- burnAndThin(jags.m.runjags$mcmc) 
#pdf(file="BaT_autocorrplot_runjags.pdf")
#autocorr.plot(BaT_runjags, ask = F)
#dev.off()


## of nimble
# nimblemodel
#BaT_nimblemodel <- burnAndThin(mcmc.out.nimbleModel_II) 
#pdf(file="BaT_autocorrplot_nimblemodel.pdf")
#autocorr.plot(BaT_nimblemodel, ask = F)
#dev.off()


# readBUGSmodel
#BaT_readBUGSmodel <- burnAndThin(mcmc.out.readBUGSModel_II) 
#pdf(file="BaT_autocorrplot_nimblemodel.pdf")
#autocorr.plot(BaT_readBUGSmodel, ask = F)
#dev.off()





##### ##### effective sample size (ESS)  function ##### ##### 



## this is a way to estimate  for burn-in cut-off.
## EEs is the number of independent samples equivalent to our number of autocorrelated samples
## s.: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html#1_objectives
## Aim: which package achieve earlier the state for sampling
## -> reduction of expenses

## Note: plots need many space,
## only ESS from the mcmcr package will be better

#### ESS 
## of WinBUGS
ess(as.mcmc(test_WinBUGS$sims.matrix), by = "term")
help("ess")

## of JAGS
# jagsUI
ess(as.mcmc(jags.m.jagsUI), by = "term")

# R2jags
ess(as.mcmc(jags.m.R2jags), by = "term")

# rjags
ess(samps_coda, by = "term")

# runjags
ess(jags.m.runjags$mcmc, by = "term")


## of nimble
# nimblemodel
ess(mcmc.out.nimbleModel_II, by = "term")

# readBUGSmodel
ess(mcmc.out.readBUGSModel_II, by = "term")





# Part 9: Convergence diagnostic with the Gelman-Rubin Diagnostic    ------------------------------------------


## with the use of the BayesianTools package
## Note: gelmanDiagnostics() need know an object of mcmcSampler or mcmcSamplerList
## so the functions of the coda package are used
## Note: Rhat = 'potential scale reduction factor' from the gelman.diag()
## however the gelman.diag() generate an upper C.I. for Rhat and can be used for Verification 
## so this Part can be continued if necessary 





########## ########## ##########   Test finished   ########## ########## ########## 


