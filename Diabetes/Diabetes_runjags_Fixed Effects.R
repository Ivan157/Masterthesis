########## ########## Simulation Diabetes Beispiel mit Random Effects mit dem runjags package
########## Die Working Directory muss auf Ihre Bed�rfnisse angepasst werden 





# Teil Simulation mit JAGS ------------------------------------------------------------------




##### Clear data
rm(list=ls())



##### load libraries
library(rjags)
library(coda)
library(random)
library(matrixStats) # zus�tzl Paket, berechnet Median
load.module("glm") 
load.module("lecuyer")
load.module("dic")



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit runjags/Ex3 Diabetes")



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



##### read in inits with chains
# Anmerkung: da cloglog als Link: wurden die Inits von JAGS generiert.
# Die manuellen Inits sind zwar drin, k�nnen aber rausgelassen werden

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
jags.m <- run.jags(model="Diabetes_Fixed.txt", monitor=c("d[2]", "d[3]",  "d[4]", "d[5]", "d[6]", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]", "T[6]" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "time"=time, "t"=t, "n"=n) , n.chains=3,  burnin = 10000, sample = 20000, adapt = 5000)



#### optional, falls nicht konvergiert:
#jags.m <- autorun.jags(model="Diabetes_Fixed.txt", monitor=c("d[2]", "d[3]",  "d[4]", "d[5]", "d[6]", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]", "T[6]",  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
#                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "time"=time, "t"=t, "n"=n) , n.chains=3, inits=all.inits, burnin = 10000, sample = 20000, adapt = 5000, max.time="1hr")





# Ausgabe posteriore Werte und Berechnung DIC --------------------------




print(jags.m)
summary(jags.m$mcmc) # f�r 2.5 - 97.5  CrI





########## ########## ########## Simulation beendet ########## ########## ##########

