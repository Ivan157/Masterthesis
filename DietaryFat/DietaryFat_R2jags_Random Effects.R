########## ########## Simulation Dietary Fat Beispiel mit Random Effects mit dem R2jags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 





# Teil Simulation mit JAGS ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
library(rjags) # R2jags benötigt rjags
library(R2jags)
library(random)
#load.module("glm") 
load.module("lecuyer")
#load.module("dic")



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit R2jags/DietaryFat")



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



dat <- list("ns", "nt", "na", "r", "E", "t")  # names list of numbers 
#dat



##### Parameter to monitor/save
params <- c("d[2]", "T[1]", "T[2]", "sd" , "dev", "totresdev", "theta" )
params



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
jags.m <- jags(data=dat,  inits=all.inits, parameters.to.save=params,  n.chains = 3, n.iter = 110000, n.burnin = 90000,
               model.file="DietaryFat_Random.txt", DIC=TRUE,  jags.module = c("glm","dic") )  
# zusätzlich noch mehrere Argumente standardmäßig dabei, v.a. interessant: DIC, jags.module



#### optional, falls nicht konvergiert:
#jags.m.upd <- autojags(jags.m)





# Ausgabe posteriore Werte, Berechnung Median und Berechnung DIC --------------------------



print(jags.m)



#### Median
jags.m[["BUGSoutput"]][["median"]]





# nachträgliche Berechnung von pD --------------------------------------



#jags.m_levPlot[["BUGSoutput"]][["summary"]]
out_lePlo <- capture.output( jags.m[["BUGSoutput"]][["summary"]])
cat("Hilf_pD", out_lePlo, file="Hilf.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf.txt", sep = "", header=F, skip=5, nrows=21)
#Hilf_data
Hilf_dev <- Hilf_data[1:21,2]
#Hilf_dev



# manuelle Berechnung von pD
#dev ist Std-Abweichung jedes einzelnen Werts
# insg 21 Werte
Var_manuell <-  sum(Hilf_dev)^2/21
pD_manuell <- Var_manuell/2
pD_manuell





# zusätzliche Diagnostik, bei Bedarf aktivieren --------------------------------------------------------------------



#pdf("DietaryFat_Random_trace.pdf")
#plot(jags.m)
#traceplot(jags.m)
#dev.off()



# Generate MCMC object for analysis
#Anm: es scheint, dass die Zeile "jags.m.mcmc <- as.mcmc(jags.m) " manuell ausgeführt werden muss
#jags.m.mcmc <- as.mcmc(jags.m) 
#jags.m.mcmc



#pdf("jags.m.mcmc.autocorr.pdf") # Autocorrelation plot
#autocorr.plot(jags.m.mcmc)
#dev.off()



# Other diagnostics using CODA:
#gelman.plot(jags.m.mcmc)
#geweke.diag(jags.m.mcmc)
#geweke.plot(jags.m.mcmc)
#raftery.diag(jags.m.mcmc)
#heidel.diag(jags.m.mcmc)





########## ########## ########## Simulation beendet ########## ########## ##########

