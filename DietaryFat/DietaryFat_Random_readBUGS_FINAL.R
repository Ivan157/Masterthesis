########## ########## Simulation DietaryFat Beispiel mit Random Effects mit NIMBLE
########## Verwendung readBUGS
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden.





# Teil 1 Creating a model ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
library(nimble)
library(car)
#library(igraph)
library(coda)



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/DietaryFat")



##### Read the data into R.
data = as.matrix(read.table("DietaryFat_Data.txt", sep = "", header=T))
# dieser Schritt ist für die Erzeugung des leverage Plotes und des pD notwendig



##### Definierung Model Code, seiner Konstanten, Daten, und initialen Werte für MCMC.
# help(readBUGSmodel) # additionelle Infos
readBUGS_Model <- readBUGSmodel(model='DietaryFat_Random_Model_Nimble.bug', data = 'DietaryFat_Data_Nimble.R',
                               inits = 'DietaryFat_Inits_Nimble.R' )  
readBUGS_Model$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = readBUGS_Model,
                       nchains = 3, niter = 110000, nburnin = 90000,
                       summary = TRUE, WAIC = F,
                       monitors = c("totresdev", "T", "d", "sd"))





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


## Berechnung CrI von sd
sd_1 <- quantile(mcmc.out$samples[["chain1"]][,"sd"] , c(0.025, 0.975))
sd_2 <- quantile(mcmc.out$samples[["chain2"]][,"sd"] , c(0.025, 0.975))
sd_3 <- quantile(mcmc.out$samples[["chain3"]][,"sd"] , c(0.025, 0.975))

# CrI von sd
(sd_1 + sd_2 + sd_3)/3





# Teil 3: Nachträgliche Berechnung von pD und Erzeugung DAG des Modelles -----------------


#Model_Nimble$dev



out_lePlo <- capture.output( readBUGS_Model$dev)
cat("Hilf_pD", out_lePlo, file="Hilf3.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf3.txt", sep = "", header=F, skip=2)
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
readBUGS_Model$plotGraph() # Anweisung geht nicht bei nimbleMCMC



########## ########## ##########  Simulation beendet ########## ########## ##########

