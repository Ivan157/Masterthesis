########## ########## Simulation Dietary Fat Beispiel mit Fixed Effects mit dem rjags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 





# Teil Simulation mit JAGS ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
library(rjags)
library(coda)
library(random)
library(matrixStats) # zusätzl Paket, berechnet Median
load.module("glm") 
load.module("lecuyer")
load.module("dic")



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit rjags/DietaryFat")



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

dat <- list(ns=ns, nt=nt, na=na, r=r, E=E, t=t)  # names list of numbers 
#dat



##### read in inits with chains
inits1 <- list(d=c( NA, 0, 0), 
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=2)

inits3 <- list(d=c( NA, 2, 2), 
               sd=2, 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=3 ) 

all.inits <- list(inits1, inits2, inits2)



##### define JAGS model within R
cat("model{								# *** PROGRAM STARTS 
        for(i in 1:ns){ # LOOP THROUGH STUDIES 
          mu[i] ~ dnorm(0,.0001)    # vague priors for all trial baselines 
          for (k in 1:na[i]) {      # LOOP THROUGH ARMS 
             r[i,k] ~ dpois(theta[i,k])         # Poisson likelihood 
              theta[i,k] <- lambda[i,k]*E[i,k]  # event rate * exposure 
              log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]         # model for linear predictor 
              dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution 
          } 
          resdev[i] <- sum(dev[i,1:na[i]])      # summed residual deviance contribution for this trial 
        }    
        totresdev <- sum(resdev[])              #Total Residual Deviance 
       d[1]<-0          # treatment effect is zero for reference treatment 
      for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }  # vague priors for treatment effects 
      
      # zusätzlich eingefügt
      A ~ dnorm(-3,1.77) 
      for (k in 1:nt) { log(T[k]) <- A + d[k]  }
    } ",
    file="DietaryFat_Fixed.txt")



##### Set up the JAGS model and settings
jags.m <- jags.model( file = "DietaryFat_Random.txt", data=dat, inits=all.inits, n.chains=3, n.adapt=5000)



##### Initialization
update(jags.m, 20000) # burn in 



##### run JAGS and save posterior samples
samps_coda <- coda.samples( jags.m, variable.names=c("d[2]", "T[1]", "T[2]"),   n.iter=20000, DIC=T  ) 





# Ausgabe posteriore Werte, Berechnung Median und Berechnung DIC --------------------------




##### summarize posterior samples # funktioniert nur mit coda.samples
summary(window(samps_coda, start=20001)) # burnin = 20000
# Anmerkung: mean für totresdev ist D_res 



#### Trace Monitore
#plot(samps_coda) 



#### Median-Berechnung
# Median für T1
median(as.matrix(samps_coda[,1]))
# Median für T2
median(as.matrix(samps_coda[,2]))
# Median für d[2]
median(as.matrix(samps_coda[,3]))



##### run JAGS and save posterior samples
samps_jags <- jags.samples( jags.m, variable.names=c("dev",  "totresdev", "theta", "deviance"),   n.iter=100000, DIC=T ) 



##### Berechnung DIC und pD
D_res <-  summary(samps_jags$totresdev[])
D_res 
pD <- var(samps_jags$deviance)/2
pD 
DIC = D_res + pD
DIC




########## ########## ########## Simulation beendet ########## ########## ##########

