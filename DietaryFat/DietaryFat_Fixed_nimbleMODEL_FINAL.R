########## ########## Simulation DietaryFat Beispiel mit Fixed Effects mit NIMBLE ########## ########## 
########## Verwendung nimbleModel
########## Die Working Directory muss auf Ihre Bed�rfnisse angepasst werden.




# Teil 1 Creating a model and Simulation ------------------------------------------------------------------



##### Clear data
rm(list=ls())



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/DietaryFat")



##### load libraries
library(car)
library(nimble)
library(coda)
#library(igraph)



##### Read the data into R.
data = as.matrix(read.table("DietaryFat_Data.txt", sep = "", header=T))
#head(data) # Shows the first six entries
data2  = read.table("DietaryFat_Data_Rest.txt")
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



##### Zuordnen der Argumente f�r NIBMLE

### Zuordnen Konstanten
Nimble_constants = list(ns=ns, nt=nt, na=na, t=t, E=E)
Nimble_constants


### Zuordnen data
Nimble_data = list(r=r) # , dev=dev
Nimble_data


### Zuordnen Inits
Nimble_inits = list(d=c( NA, 0, 0), 
                    mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    # zus�tzlich noch einzuf�gen:
                    A=0)


### Create Model Code
Code_Model<- nimbleCode( {
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
  totresdev <- sum(resdev[1:10])              #Total Residual Deviance 
  d[1]<-0          # treatment effect is zero for reference treatment 
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) }  # vague priors for treatment effects 
  
  # zus�tzlich eingef�gt
  A ~ dnorm(-3,1.77) 
  for (k in 1:nt) { log(T[k]) <- A + d[k]  }
})



# nimbleModel prozessiert BUGS-Modellcode und optionale Konstanten, Daten und Initialwerte. Liefert ein NIMBLE-Modell zur�ck.
# dieser Schritt ist bei den hier getesteten BUGS Modellen nicht n�tig 
Model_Nimble <- nimbleModel(code = Code_Model, name = "ProcessedModel", constants = Nimble_constants,
                            data = Nimble_data, inits = Nimble_inits) 

Model_Nimble$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = Code_Model, constants = Nimble_constants,
                       data = Nimble_data, inits = Nimble_inits,
                       nchains = 3, niter = 21000, nburnin = 19000,
                       summary = TRUE, WAIC = F,
                       monitors = c("totresdev", "T", "d"))
#help(buildMCMC)





# Teil 2: Anzeigen Ergebnisse der Simulation ------------------------------




#### Zusammenfassung posterioreer Werte
mcmc.out[["summary"]][["all.chains"]]



#### Berechnung der CrI

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


## Berechnung CrI von d[2]
d2_1 <- quantile(mcmc.out$samples[["chain1"]][,"d[2]"] , c(0.025, 0.975))
d2_2 <- quantile(mcmc.out$samples[["chain2"]][,"d[2]"] , c(0.025, 0.975))
d2_3 <- quantile(mcmc.out$samples[["chain3"]][,"d[2]"] , c(0.025, 0.975))

# CrI von d[2]
(d2_1 + d2_2 + d2_3)/3





# Teil 3: Nachtr�gliche Berechnung von pD und Erzeugung DAG des Modelles -----------------


#Model_Nimble$dev


out_lePlo <- capture.output( Model_Nimble$dev)
cat("Hilf_pD", out_lePlo, file="Hilf2.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf2.txt", sep = "", header=F, skip=2)
#Hilf_data


Hilf_dev <- c (Hilf_data[,2], Hilf_data[,3], Hilf_data[2,4])
#Hilf_dev


# manuelle Berechnung von pD
# dev ist Std-Abweichung jedes einzelnen Werts
# insg 21 Werte
Var_manuell <-  sum(Hilf_dev)^2/21
pD_manuell <- Var_manuell/2
pD_manuell



#### Plot of model
#directed acyclic graph
#durch igraph
Model_Nimble$plotGraph() # Anweisung geht nicht bei nimbleMCMC



########## ########## ##########  Simulation beendet ########## ########## ##########

