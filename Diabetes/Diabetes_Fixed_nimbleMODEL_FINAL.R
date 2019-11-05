########## ########## Simulation Diabetes Beispiel mit Fixed Effects mit NIMBLE ########## ########## 
########## Verwendung nimbleModel
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden.




# Teil 1 Creating a model and Simulation ------------------------------------------------------------------



##### Clear data
rm(list=ls())



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/Ex3 Diabetes")




##### load libraries
library(car)
library(nimble)
library(coda)
#library(igraph)



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



##### Zuordnen der Argumente für NIBMLE


### Zuordnen Konstanten
Nimble_constants = list(ns=ns, nt=nt, na=na, t=t, time=time, n=n)
Nimble_constants


### Zuordnen data
Nimble_data = list(r=r) 
Nimble_data


############################

### Zuordnen Inits
Nimble_inits = list(d=c(NA,0,0,0,0,0), 
                    mu=c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0), 
                    A=0 )  


### Create Model Code
Code_Model<- nimbleCode( {
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
  totresdev <- sum(resdev[1:22]) #Total Residual Deviance 
  d[1]<-0 # treatment effect is zero for reference treatment 
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 
  A ~ dnorm(-4.2,1.11) 
  for (k in 1:nt) { 
    cloglog(T[k]) <- log(3) + A + d[k]  
  }
})



# nimbleModel prozessiert BUGS-Modellcode und optionale Konstanten, Daten und Initialwerte. Liefert ein NIMBLE-Modell zurück.
# dieser Schritt ist bei den hier getesteten BUGS Modellen nicht nötig 
Model_Nimble <- nimbleModel(code = Code_Model, name = "ProcessedModel", constants = Nimble_constants,
                            data = Nimble_data, inits = Nimble_inits) 

Model_Nimble$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = Code_Model, constants = Nimble_constants,
                       data = Nimble_data, inits = Nimble_inits,
                       nchains = 3, niter = 100000, nburnin = 50000,
                       summary = TRUE, WAIC = F,
                       monitors = c("totresdev", "T", "d"))
#help(buildMCMC)





# Teil 2: Anzeigen Ergebnisse der Simulation ------------------------------




#### Zusammenfassung posterioreer Werte
mcmc.out[["summary"]][["all.chains"]]






#### Berechnung der CrI für T


## Berechnung CrI von T[1]
T1_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[1]"] , c(0.025, 0.975))
T1_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[1]"] , c(0.025, 0.975))
T1_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[1]"] , c(0.025, 0.975))

# CrI von T[1]
(T1_1 + T1_2 + T1_3)/3


## Berechnung CrI von T[2]
T2_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[2]"] , c(0.025, 0.975))
T2_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[2]"] , c(0.025, 0.975))
T2_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[2]"] , c(0.025, 0.975))

# CrI von T[2]
(T2_1 + T2_2 + T2_3)/3


## Berechnung CrI von T[3]
T3_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[3]"] , c(0.025, 0.975))
T3_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[3]"] , c(0.025, 0.975))
T3_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[3]"] , c(0.025, 0.975))

# CrI von T[3]
(T3_1 + T3_2 + T3_3)/3


## Berechnung CrI von T[4]
T4_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[4]"] , c(0.025, 0.975))
T4_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[4]"] , c(0.025, 0.975))
T4_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[4]"] , c(0.025, 0.975))

# CrI von T[4]
(T4_1 + T4_2 + T4_3)/3


## Berechnung CrI von T[5]
T5_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[5]"] , c(0.025, 0.975))
T5_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[5]"] , c(0.025, 0.975))
T5_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[5]"] , c(0.025, 0.975))

# CrI von T[5]
(T5_1 + T5_2 + T5_3)/3


## Berechnung CrI von T[6]
T6_1 <- quantile(mcmc.out$samples[["chain1"]][,"T[6]"] , c(0.025, 0.975))
T6_2 <- quantile(mcmc.out$samples[["chain2"]][,"T[6]"] , c(0.025, 0.975))
T6_3 <- quantile(mcmc.out$samples[["chain3"]][,"T[6]"] , c(0.025, 0.975))


# CrI von T[6]
(T6_1 + T6_2 + T6_3)/3



#### Berechnung der CrI für d


## Berechnung CrI von d[2]
d2_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[2]"] , c(0.025, 0.975))
d2_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[2]"] , c(0.025, 0.975))
d2_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[2]"] , c(0.025, 0.975))

# CrI von d[2]
(d2_1 + d2_2 + d2_3)/3


## Berechnung CrI von d[3]
d3_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[3]"] , c(0.025, 0.975))
d3_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[3]"] , c(0.025, 0.975))
d3_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[3]"] , c(0.025, 0.975))

# CrI von d[3]
(d3_1 + d3_2 + d3_3)/3


## Berechnung CrI von d[4]
d4_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[4]"] , c(0.025, 0.975))
d4_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[4]"] , c(0.025, 0.975))
d4_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[4]"] , c(0.025, 0.975))

# CrI von d[4]
(d4_1 + d4_2 + d4_3)/3


## Berechnung CrI von d[5]
d5_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[5]"] , c(0.025, 0.975))
d5_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[5]"] , c(0.025, 0.975))
d5_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[5]"] , c(0.025, 0.975))


# CrI von d[5]
(d5_1 + d5_2 + d5_3)/3


## Berechnung CrI von d[6]
d6_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[6]"] , c(0.025, 0.975))
d6_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[6]"] , c(0.025, 0.975))
d6_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[6]"] , c(0.025, 0.975))

# CrI von d[6]
(d6_1 + d6_2 + d6_3)/3




#Model_Nimble$dev


out_lePlo <- capture.output( Model_Nimble$dev)
cat("Hilf_pD", out_lePlo, file="Hilf2.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf2.txt", sep = "", header=F, skip=2)
Hilf_data


Hilf_dev <- c (Hilf_data[,2], Hilf_data[,3], Hilf_data[1,4], Hilf_data[5,4], Hilf_data[15,4], Hilf_data[16,4])
#Hilf_dev


# manuelle Berechnung von pD
# dev ist Std-Abweichung jedes einzelnen Werts
# insg 48 Werte
Var_manuell <-  sum(Hilf_dev)^2/48
pD_manuell <- Var_manuell/2
pD_manuell



#### Plot of model
#directed acyclic graph
#durch igraph
Model_Nimble$plotGraph() # Anweisung geht nicht bei nimbleMCMC



########## ########## ##########  Simulation beendet ########## ########## ##########

