########## ########## Simulation Blocker Beispiel mit Fixed Effects mit NIMBLE
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
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/Blocker")



##### Read the data into R.
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
# dieser Schritt ist für die Erzeugung des leverage Plotes notwendig



##### Definierung Model Code, seiner Konstanten, Daten, und initialen Werte für MCMC.
# help(readBUGSmodel) # additionelle Infos
readBUGS_Model <- readBUGSmodel(model='Blocker_Fixed_Model_II.bug', data = 'Blocker_Data_Nimble.R',
                                inits = 'Blocker_Inits_Nimble_Fixed.R' )  
readBUGS_Model$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = readBUGS_Model,
                       nchains = 3, niter = 20000, nburnin = 10000,
                       summary = TRUE, WAIC = F,
                       monitors = c("rhat",'dev', "totresdev", "T", "d"))





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





# Teil 3: Leverage Plot und nachträgliche Berechnung von pD -----------------




out_lePlo <- capture.output( mcmc.out$summary)
cat("Hilf_lePlo", out_lePlo, file="Hilf4.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf4.txt", sep = "", header=F, skip=295, nrows=88)
#Hilf_data



#### Berechnung w_ik

Hilf_dev <- Hilf_data[1:44,3]
#Hilf_dev
Hilf_dev_k1 <- Hilf_dev[1:22]
#Hilf_dev_k1
Hilf_dev_k2 <- Hilf_dev[23:44]
#Hilf_dev_k2
Hilf_dev_II <- cbind(Hilf_dev_k1, Hilf_dev_k2)
#Hilf_dev_II
Hilf_dev_III <- cbind(Hilf_dev_II, total = rowMeans(Hilf_dev_II))
#Hilf_dev_III
w_ik <- sqrt(Hilf_dev_III[,3])
w_ik_neg <- -w_ik
fertige_Daten_w_ik <- cbind(Hilf_dev_III, w_ik_neg, w_ik)
fertige_Daten_w_ik



# manuelle Berechnung von pD
# dev ist Std-Abweichung jedes einzelnen Werts
# insg 22 Wertepaare
Var_manuell <-  sum(Hilf_dev)^2/44
pD_manuell <- Var_manuell/2
pD_manuell



#### Berechnung leverage_ik

Hilf_rhat <- Hilf_data[45:88,3]
#Hilf_rhat
Hilf_rhat_k1 <- Hilf_rhat[1:22]
#Hilf_rhat_k1
Hilf_rhat_k2 <- Hilf_rhat[23:44]
#Hilf_rhat_k2


dev_tilde_erst_k1 <- data[,1]*log(data[,1]/Hilf_rhat_k1)
#dev_tilde_erst_k1
dev_tilde_zweit_k1 <- (data[,3]-data[,1])*log((data[,3]-data[,1])/(data[,3]-Hilf_rhat_k1))
#dev_tilde_zweit_k1
dev_tilde_gesamt_k1 <- 2*(dev_tilde_erst_k1+dev_tilde_zweit_k1)
#dev_tilde_gesamt_k1
leverage_k1 <-  fertige_Daten_w_ik[,1] - dev_tilde_gesamt_k1

leverage_k1


dev_tilde_erst_k2 <- data[,2]*log(data[,2]/Hilf_rhat_k2)
#dev_tilde_erst_k2
dev_tilde_erst_k2 <- data[,2]*log(data[,2]/Hilf_rhat_k2)
#dev_tilde_erst_k2
dev_tilde_zweit_k2 <- (data[,4]-data[,2])*log((data[,4]-data[,2])/(data[,4]-Hilf_rhat_k2))
#dev_tilde_zweit_k2
dev_tilde_gesamt_k2 <- 2*(dev_tilde_erst_k2+dev_tilde_zweit_k2)
#dev_tilde_gesamt_k2
leverage_k2 <-  fertige_Daten_w_ik[,2] - dev_tilde_gesamt_k2

leverage_k2



#### Erzeugen leverage plot

scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_k1, leverage_k2),  
            main="Leverage plot for the fixed effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), 
            ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F, grid=TRUE )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 



#### Plot of model
#directed acyclic graph
#durch igraph
readBUGS_Model$plotGraph() # Anweisung geht nicht bei nimbleMCMC





########## ########## ########## Simulation beendet ########## ########## ##########

