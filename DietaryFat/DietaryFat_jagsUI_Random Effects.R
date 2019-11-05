########## ########## Simulation Dietary Fat Beispiel mit Random Effects mit dem jagsUI package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 





# Teil Simulation mit JAGS ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
#library(rjags) # jagsUI benötigt dieses Paket
library(lattice) 
#library(coda)  
library(jagsUI)
#library(random)



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit jagsUI/Dietary fat")



##### Read the data into R.
#data = read.table("DietaryFat_Data.txt", sep = "", header=F)
data = as.matrix(read.table("DietaryFat_Data.txt", sep = "", header=T))
head(data) # Shows the first six entries
#data2  = as.data.frread.table("DietaryFat_Data_Rest.txt")
#head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
ns # check
nt <-  ncol(data[,7:9])
nt # check
na <- data[,10]
na  # check
r <- data[,7:9]
r # Check
E <-  data[,4:6]
E # Check
t <-  data[,1:3]
t # Check

dat <- list("ns", "nt", "na", "r", "E", "t")  # names list of numbers 
dat



##### Parameter to monitor/save
params <- c("d[2]", "T[1]", "T[2]", "sd" , "dev", "totresdev", "theta" )
params



##### read in inits with chains
inits1 <- list(d=c(NA, 0, 0), # Änderung Ausgangsdaten zu WinBUGS: Hinzufügen eines weiteren Initialwertes für d
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3) ,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=2)

inits3 <- list(d=c( NA, 2, 2), 
               sd=2, 
               mu= c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=3 ) 

# Achtung: "lecuyer::RngStream" funktioniert nicht mit jagsUI
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
    file="DietaryFat_Random.txt")



##### Set up the JAGS model and settings
jags.m <- jags(data=dat, 
               inits=all.inits,
               parameters.to.save=params,
               model.file="DietaryFat_Random.txt", 
               n.chains=3,
               n.iter=110000, n.burnin=90000,
               store.data=TRUE) 
    #bei n.iter=100k & burn-in=100k: 
        #Error in process.input(data, parameters.to.save, inits, n.chains, n.iter,  : 
        #Number of iterations must be larger than burn-in.
#Warning messages:
#1: In process.input(data, parameters.to.save, inits, n.chains, n.iter,  :
#Suppling a list of character strings to the data argument will be deprecated in the next version
#2: In process.input(data, parameters.to.save, inits, n.chains, n.iter,  :
#Suppling a character vector to the data argument will be deprecated in the next version
# => muss demnächst angepasst werden





# Anzeigen der posterioren Werte und des Medians --------------------------




#traceplot(jags.m) # zeigt Abbildungen einzeln nach Eingabe der Entertaste an
jags.View(jags.m)
jags.m
# mit jags-Funktion keine weiterführende Diagnostik möglich: das Objekt jags.m zeigt nur eine Liste von 24 Elementen an.
# ebenso wenig mit jags.basic



# Versuch mit jagsbasic
jagsbasic.m <- jags.basic(data=dat, 
                          inits=all.inits,
                          parameters.to.save=params,
                          model.file="DietaryFat_Random.txt", 
                          n.chains=3,
                          n.iter=110000, n.burnin=90000) 
#ja, funktioniert
jagsbasic.m_II <- do.call(rbind.data.frame, jagsbasic.m)



#### Berechnung Median
# Median von T1
median(jagsbasic.m_II$`T[1]`)

# Median von T2
median(jagsbasic.m_II$`T[2]`)

# Median von d
median(jagsbasic.m_II$`d[2]`)

# Median von sd
median(jagsbasic.m_II$`sd`)






########## ########## ########## Simulation beendet ########## ########## ##########

