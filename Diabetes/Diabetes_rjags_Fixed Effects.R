########## ########## Simulation Diabetes Beispiel mit Random Effects mit dem rjags package
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
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit rjags/Ex3 Diabetes")



##### Read the data into R.
#data = read.table("Diabetes_Data_neu sortiert.txt", sep = "", header=F)
data = as.matrix(read.table("Diabetes_Data_neu sortiert.txt", sep = "", header=T))
head(data) # Shows the first six entries
#data2  = as.data.frread.table("Diabetes_Data_Rest.txt")
#head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
#ns # check
nt <-  6
#nt # check
na <- data[,11]
#na  # check
r <- data[,5:7]
#r # Check
time <-  data[,1]
#time # Check
t <-  data[,2:4]
#t # Check
n <-  data[,8:10]
#n # Check

dat <- list(ns=ns, nt=nt, na=na, r=r, time=time, t=t, n=n)  # names list of numbers 
dat



##### read in inits with chains
# Anmerkung: da cloglog als Link: wurden die Inits von JAGS generiert.
# Die manuellen Inits sind zwar drin, können aber rausgelassen werden

inits1 <- list(d=c(NA,0,0,0,0,0), 
               sd=1,  
               mu=c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0), 
               A=0 ,
               delta= structure(.Data= c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
                                         0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
                                         NA, 0, NA), .Dim=c(22, 3)))  
# .RNG.name="base::Wichmann-Hill", .RNG.seed=1

inits2 <- list(d=c(NA,-1,4,-1,2,3), 
               sd=3,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,0,0,   0,1,0,0,0,  1,1), 
               A=1 ,
               delta= structure(.Data= c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
                                         0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
                                         NA, 0, NA), .Dim=c(22, 3))) 

inits3 <- list(d=c(NA,1,4,-3,-2,3),  
               
               sd=4.5,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,-2,0,   0,1,0,-2,0,  1,1), 
               A=2 ,
               delta= structure(.Data= c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
                                         0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
                                         NA, 0, NA), .Dim=c(22, 3))) 
#,      .RNG.name="base::Wichmann-Hill", .RNG.seed=3

all.inits <- list(inits1, inits2, inits3)
# all.inits <- list(inits2)



##### define JAGS model within R
cat("model{             # *** PROGRAM STARTS 
        for(i in 1:ns){ # LOOP THROUGH STUDIES 
          mu[i] ~ dnorm(0,.0001) # vague priors for all trial baselines 
          for (k in 1:na[i]) { # LOOP THROUGH ARMS 
            r[i,k] ~ dbin(p[i,k],n[i,k]) # Binomial likelihood 
            cloglog(p[i,k]) <- log(time[i]) + mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor 
            rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
            dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
          } 
        resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
        }    
        totresdev <- sum(resdev[]) #Total Residual Deviance 
        d[1]<-0 # treatment effect is zero for reference treatment 
        for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 
        A ~ dnorm(-4.2,1.11) 
        for (k in 1:nt) { 
            cloglog(T[k]) <- log(3) + A + d[k]  
        }
      } ",
    file="Diabetes_Fixed.txt")



##### Set up the JAGS model and settings
jags.m <- jags.model( file = "Diabetes_Fixed.txt", data=dat,  n.chains=3, n.adapt=5000)



##### Initialization
update(jags.m, 10000) # burn in 



##### run JAGS and save posterior samples
samps_coda <- coda.samples( jags.m, variable.names=c("d[2]", "d[3]",  "d[4]", "d[5]", "d[6]", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]", "T[6]" ,  "totresdev"),   n.iter=20000, DIC=T  ) 





# Ausgabe posteriore Werte, Berechnung Median und Berechnung DIC --------------------------




##### summarize posterior samples # funktioniert nur mit coda.samples
summary(window(samps_coda, start=10001)) # burnin = 100000
# Anmerkung: mean für totresdev ist D_res 



#### Trace Monitore
#plot(samps_coda) 



#### Median-Berechnung

#as.matrix(samps_coda[])

# Median für T1
median(as.matrix(samps_coda[,1]))
# Median für T2
median(as.matrix(samps_coda[,2]))
# Median für T3
median(as.matrix(samps_coda[,3]))
# Median für T4
median(as.matrix(samps_coda[,4]))
# Median für T5
median(as.matrix(samps_coda[,5]))
# Median für T6
median(as.matrix(samps_coda[,6]))

# Median für d[2]
median(as.matrix(samps_coda[,7]))
# Median für d[3]
median(as.matrix(samps_coda[,8]))
# Median für d[4]
median(as.matrix(samps_coda[,9]))
# Median für d[5]
median(as.matrix(samps_coda[,10]))
# Median für d[6
median(as.matrix(samps_coda[,11]))



##### run JAGS and save posterior samples
samps_jags <- jags.samples( jags.m, variable.names=c("dev",  "totresdev", "rhat", "deviance"),   n.iter=20000, DIC=T ) 



##### Berechnung DIC und pD
D_res <-  summary(samps_jags$totresdev[])
D_res 
pD <- var(samps_jags$deviance)/2
pD 
DIC = D_res + pD
DIC





########## ########## ########## Simulation beendet ########## ########## ##########

