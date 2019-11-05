########## ########## Simulation Dietary Fat Beispiel mit Fixed Effects mit dem jagsUI package
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
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit jagsUI/Ex3 Diabetes")



##### Read the data into R.
#data = read.table("Diabetes_Data_neu sortiert.txt", sep = "", header=F)
data = as.matrix(read.table("Diabetes_Data_neu sortiert.txt", sep = "", header=T))
head(data) # Shows the first six entries
#data2  = as.data.frread.table("Diabetes_Data_Rest.txt")
#head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
ns # check
nt <-  6
nt # check
na <- data[,11]
na  # check
r <- data[,5:7]
r # Check
time <-  data[,1]
time # Check
t <-  data[,2:4]
t # Check
n <-  data[,8:10]
n

dat <- list("ns", "nt", "na", "r", "time", "t", "n" )  # names list of numbers
dat



##### Parameter to monitor/save
params <- c("d[2]", "d[3]",  "d[4]", "d[5]", "d[6]", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]", "T[6]",  "totresdev"," dev", "rhat" )
params



##### read in inits with chains
# Anmerkung: da cloglog als Link: wurden die Inits von JAGS generiert.
# Die manuellen Inits sind zwar drin, können aber rausgelassen werden

inits1 <- list(d=c(NA,0,0,0,0,0), 
               sd=1,  
               mu=c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0), 
               A=0 )  
# .RNG.name="base::Wichmann-Hill", .RNG.seed=1

inits2 <- list(d=c(NA,-1,4,-1,2,3), 
               sd=3,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,0,0,   0,1,0,0,0,  1,1), 
               A=1 ) 

inits3 <- list(d=c(NA,1,4,-3,-2,3),  
               
               sd=4.5,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,-2,0,   0,1,0,-2,0,  1,1), 
               A=2 ) 
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
jags.m <- jags(data=dat, 
               inits=NULL,
               parameters.to.save=params,
               model.file="Diabetes_Fixed.txt", 
               n.chains=3,
               n.iter=20000, n.burnin=10000,
               store.data=TRUE) 

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
                          inits=NULL,
                          parameters.to.save=params,
                          model.file="Diabetes_Fixed.txt", 
                          n.chains=3,
                          n.iter=20000, n.burnin=10000) 

#jagsbasic.m
jagsbasic.m_II <- do.call(rbind.data.frame, jagsbasic.m)



#### Berechnung Median

# Median von T1
median(jagsbasic.m_II$`T[1]`)

# Median von T2
median(jagsbasic.m_II$`T[2]`)

# Median von T3
median(jagsbasic.m_II$`T[3]`)

# Median von T4
median(jagsbasic.m_II$`T[4]`)

# Median von T5
median(jagsbasic.m_II$`T[5]`)

# Median von T6
median(jagsbasic.m_II$`T[6]`)

# Median von d2
median(jagsbasic.m_II$`d[2]`)

# Median von d3
median(jagsbasic.m_II$`d[3]`)

# Median von d4
median(jagsbasic.m_II$`d[4]`)

# Median von d5
median(jagsbasic.m_II$`d[5]`)

# Median von d6
median(jagsbasic.m_II$`d[6]`)





########## ########## ########## Simulation beendet ########## ########## ##########

