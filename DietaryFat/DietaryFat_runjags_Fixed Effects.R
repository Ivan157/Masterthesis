########## ########## Simulation Dietary Fat Beispiel mit Fixed Effects mit dem runjags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 





# Teil Simulation mit JAGS ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
library(rjags) 
library(runjags)
library(random)
library(coda)
load.module("glm") 
load.module("lecuyer")
load.module("dic")



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit runjags/DietaryFat")



##### Read the data into R.
#data = read.table("DietaryFat_Data.txt", sep = "", header=F)
data = as.matrix(read.table("DietaryFat_Data.txt", sep = "", header=T))
head(data) # Shows the first six entries
#data2  = as.data.frread.table("DietaryFat_Data_Rest.txt")
#head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
#ns # check
nt <-  ncol(data[,7:9])
#nt # check
na <- data[,10]
#na  # check
r <- data[,7:9]
#r # Check
E <-  data[,4:6]
#E # Check
t <-  data[,1:3]
#t # Check



##### read in inits with chains
inits1 <- list(d=c( NA, 0, 0), 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1, -1), 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=2)

inits3 <- list(d=c( NA, 2, 2), 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=3 ) 

all.inits <- list(inits1, inits2, inits2)



##### define JAGS model within R
cat("model{								# *** PROGRAM STARTS 
        for(i in 1:ns){                 # LOOP THROUGH STUDIES 
            w[i,1] <- 0                 # adjustment for multi-arm trials is zero for control arm 
            delta[i,1] <- 0             # treatment effect is zero for control arm 
            mu[i] ~ dnorm(0,.0001)      # vague priors for all trial baselines 
            
            for (k in 1:na[i]) {            # LOOP THROUGH ARMS 
                r[i,k] ~ dpois(theta[i,k])  # Poisson likelihood 
                theta[i,k] <- lambda[i,k]*E[i,k]            # failure rate * exposure 
                log(lambda[i,k]) <- mu[i] + delta[i,k]      # model for linear predictor 
                dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k]))     #Deviance contribution 
            } 
            resdev[i] <- sum(dev[i,1:na[i]])        # summed residual deviance contribution for this trial 
            for (k in 2:na[i]) {                    # LOOP THROUGH ARMS 
                delta[i,k] ~ dnorm(md[i,k],taud[i,k])           # trial-specific LOR distributions 
                md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k]      # mean of LOR distributions (with multi-arm trial correction
                taud[i,k] <- tau *2*(k-1)/k                     # precision of LOR distributions (with multi-arm trial correction
                w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) 				# adjustment for multi-arm RCTs 
                sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) 				# cumulative adjustment for multi-arm trials 
            } 
        }    
        totresdev <- sum(resdev[])          #Total Residual Deviance 
        d[1]<-0                             # treatment effect is zero for reference treatment 
        for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }    # vague priors for treatment effects 
        sd ~ dunif(0,5)                             # vague prior for between-trial SD 
        tau <- pow(sd,-2)                   # between-trial precision = (1/between-trial variance) 
    
        # zusätzlich eingefügt
        A ~ dnorm(-3,1.77) 
        for (k in 1:nt) { log(T[k]) <- A + d[k]  }
    } ",
    file="DietaryFat_Fixed.txt")



##### Set up the JAGS model and settings
jags.m <- run.jags(model="DietaryFat_Fixed.txt", monitor=c("d[2]", "T[1]", "T[2]",  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "E"=E, "t"=t) , n.chains=3, inits=all.inits, burnin = 20000, sample = 20000, adapt = 5000)



#### optional, falls nicht konvergiert:
#jags.m <- autorun.jags(model="DietaryFat_Random.txt", monitor=c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
#                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "E"=E, "t"=t) , n.chains=3, 
#                   inits=all.inits, burnin = 100000, sample = 100000, adapt = 2000, , max.time="1hr")





# Ausgabe posteriore Werte und Berechnung DIC --------------------------



print(jags.m)
summary(jags.m$mcmc) # für 2.5 - 97.5  CrI




########## ########## ########## Simulation beendet ########## ########## ##########

