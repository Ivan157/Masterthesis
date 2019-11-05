########## ########## Simulation Diabetes Beispiel mit Fixed Effects mit NIMBLE
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
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/Ex3 Diabetes")



##### Read the data into R.
data = as.matrix(read.table("Diabetes_Data_neu sortiert.txt", sep = "", header=T))
# dieser Schritt ist für die Erzeugung des leverage Plotes und des pD notwendig



##### Definierung Model Code, seiner Konstanten, Daten, und initialen Werte für MCMC.
# help(readBUGSmodel) # additionelle Infos
readBUGS_Model <- readBUGSmodel(model='Diabetes_Fixed_Model_Nimble.bug', data = 'Diabetes_Data_Nimble.R',
                                inits = 'Diabetes_Inits_Fixed_Nimble.R' )  
readBUGS_Model$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = readBUGS_Model,
                       nchains = 3, niter = 100000, nburnin = 50000,
                       summary = TRUE, WAIC = F,
                       monitors = c("totresdev", "T", "d"))


readBUGS_Model$r


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





# Teil 3: Nachträgliche Berechnung von pD und Erzeugung DAG des Modelles -----------------



#Model_Nimble$dev


out_lePlo <- capture.output( readBUGS_Model$dev)
cat("Hilf_pD", out_lePlo, file="Hilf4.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf3.txt", sep = "", header=F, skip=2)
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
readBUGS_Model$plotGraph() # Anweisung geht nicht bei nimbleMCMC



########## ########## ##########  Simulation beendet ########## ########## ##########

