########## ########## Simulation of a possible application of NMA in biotechnology  ########## ##########
########## use of the runjags package
########## NMA in the wastewater treatment
##### Influence of A-Toxin to MO B




# Part 1:  Preparation ------------------------------------------------------------------




##### Clear data
rm(list=ls())



##### load libraries
library(rjags)
library(runjags)
library(coda)
library(random)
library(readxl)
library(matrixStats) # additional package, calculates Median
load.module("glm") 
load.module("lecuyer")
load.module("dic")



#####  Setting Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Ergebnisse/zusätzliche Experimente/Anwendung an die Biotec/Ansatz 1")



##### Generation of Data
## optimal conditions
set.seed(1)
Data_opt <- rnorm(12, mean = 6000, sd = 100)
#Data_opt
Data_opt[4:6] <- Data_opt[4:6]/10
Data_opt[7:9] <- Data_opt[7:9]/100
Data_opt[10:12] <- Data_opt[10:12]/1000
#Data_opt
Data_opt <- as.integer(  Data_opt)
#Data_opt

## Toxin A-free conditions
set.seed(2)
Data_tox_free <- rnorm(12, mean = 4000, sd = 100)
#Data_tox_free
Data_tox_free[4:6] <- Data_tox_free[4:6]/10
Data_tox_free[7:9] <- Data_tox_free[7:9]/100
Data_tox_free[10:12] <- Data_tox_free[10:12]/1000
#Data_tox_free
Data_tox_free <- as.integer(  Data_tox_free)
#Data_tox_free

## Toxin A containing conditions
set.seed(3)
Data_toxin <- rnorm(12, mean = 3800, sd = 100)
#Data_toxin
Data_toxin[4:6] <- Data_toxin[4:6]/(10*2)
Data_toxin[7:9] <- Data_toxin[7:9]/(100*3)
Data_toxin[10:12] <- Data_toxin[10:12]/(1000*4)
#Data_toxin
Data_toxin <- as.integer(  Data_toxin)
#Data_toxin

## all Data
Data_opt
Data_tox_free
Data_toxin
# -> Data seems plausible and are saved in the Daten.xlsx-file



##### Read the data into R.
#data = read_xlsx("Daten.xlsx")
data = as.matrix(read.table("Daten.txt", sep = "", header=F))
head(data) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <- nrow(data)
# ns # check
nt <- ncol(data[,5:6])
#nt # check

na <- data[,7]
# na  # check
r <- data[,1:2]
# r # Check
n <-  data[,3:4]
# n # Check
t <-  data[,5:6]
# t # Check



##### read in inits with chains
inits1 <- list(d=c( NA, 0), 
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=2)

inits3 <- list(d=c( NA, 2), 
               sd=2, 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=3 ) 

all.inits <- list(inits1, inits2, inits2)



##### define JAGS model within R
cat("model{								 # *** PROGRAM STARTS 
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
    totresdev <- sum(resdev[]) 								#Total Residual Deviance 
    d[1] <- 0												 # treatment effect is zero for reference treatment 
    sd ~ dunif(0,5) 												# vague prior for between-trial SD.
    tau <- pow(sd,-2) 								# between-trial precision = (1/between-trial variance) 

    # Provide estimates of treatment effects T[k] on the natural (probability) scale  
    # Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 

    for (k in 2:nt){ d[k] ~ dnorm(0,0.0001) }
    for (k in 1:nt) { logit(T[k]) <- A + d[k] }
    A ~ dnorm(-2.2, 3.3) 
    } ",
    file="Wastewater_Random.txt")





# Part 2: Simulation with rujags ----------------------------------------------------------------



##### Set up the JAGS model and settings
jags.m <- run.jags(model="Wastewater_Random.txt", monitor=c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , n.chains=3, inits=all.inits, burnin = 10000, sample = 20000, adapt = 1500)



#### optional, if not convergated:
#jags.m <- autorun.jags(model="Wastewater_Random.txt", monitor=c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
#                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , n.chains=3, 
#                   inits=all.inits, burnin = 10000, sample = 30000, adapt = 1500, , max.time="1hr")






# Part 3: Output posterior values and calculation median --------------------------




##### summarize posterior samples 
print(jags.m)
summary(jags.m$mcmc) # für 2.5 - 97.5  CrI



##### Generate plots and save to separate file
png(filename="Plot_Mo_B.png")
plot(jags.m)
dev.off()





##### ##### ##### simulation finished ##### ##### ##### 

