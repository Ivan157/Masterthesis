########## ########## Simulation Blocker Beispiel mit Random Effects mit NIMBLE ########## ########## 
########## Verwendung nimbleModel
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden.




# Teil 1 Creating a model and Simulation ------------------------------------------------------------------



##### Clear data
rm(list=ls())



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/NIMBLE/Nachrechnen TSD2/Blocker")



##### load libraries
library(nimble)
library(car)
library(coda)
#library(igraph)



##### Read the data into R.
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
#head(data) # Shows the first six entries
data2  = read.table("Data_Blocker_Rest.txt")
#head(data2) # Shows the first six entries



##### Values for simulation, prepare dat for NIMBLE (allocation values from data)
ns <-  nrow(data)
#ns # check
nt <-  ncol(data[,5:6])
#nt # check
na <- data[,7]
#na  # check
r <- data[,1:2]
#r # Check
n <-  data[,3:4]
#n # Check
t <-  data[,5:6]
#t # Check



##### Zuordnen der Argumente für NIBMLE


### Zuordnen Konstanten
Nimble_constants = list(ns=ns, nt=nt, na=na, t=t)
#Nimble_constants


### Zuordnen data
Nimble_data = list(r=r, n=n)
#Nimble_data



### Zuordnen Inits
Nimble_inits = list(d = c( NA, 0),
                    sd = c( 1),
                    mu = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    # zusätzlich noch einzufügen:
                    A=0,
                    delta=matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol=2))



### Create Model Code
Code_Model<- nimbleCode( {
  for(i in 1:ns){								 # LOOP THROUGH STUDIES 
    w[i,1] <- 0 								# adjustment for multi-arm trials is zero for control arm 
    delta[i,1] <- 0 								# treatment effect is zero for control arm 
    mu[i] ~ dnorm(0,.0001) 				# vague priors for all trial baselines 
    for (k in 1:na[i]) { 								# LOOP THROUGH ARMS 
      r[i,k] ~ dbin(p[i,k],n[i,k]) 				# binomial likelihood 
      logit(p[i,k]) <- mu[i] + delta[i,k]				 # model for linear predictor 
      rhat[i,k] <- p[i,k] * n[i,k] 				# expected value of the numerators 
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution 
                       +  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))  
    } 
    resdev[i] <- sum(dev[i,1:na[i]]) 				# summed residual deviance contribution for this trial 
    for (k in 2:na[i]) { 								# LOOP THROUGH ARMS 
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) 				# trial-specific LOR distributions 
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] 				# mean of LOR distributions (with multi-arm trial correction) 
      taud[i,k] <- tau *2*(k-1)/k 								# precision of LOR distributions (with multi-arm trial correction) 
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) 				# adjustment for multi-arm RCTs 
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) 				# cumulative adjustment for multi-arm trials; Modifikation
    } 
  }    
  totresdev <- sum(resdev[1:22]) 								#Total Residual Deviance; Modifikation
  
  d[1] <- 0												 # treatment effect is zero for reference treatment 
  sd ~ dunif(0,5) 												# vague prior for between-trial SD.
  tau <- pow(sd,-2) 								# between-trial precision = (1/between-trial variance) 
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale  
  # Given a Mean Effect, meanA, for â€˜standardâ€™ treatment 1, with precision (1/variance) precA 
  
  
  for (k in 2:nt){ d[k] ~ dnorm(0,0.0001) }
  for (k in 1:nt) { logit(T[k]) <- A + d[k] }
  A ~ dnorm(-2.2, 3.3) 
})

# nimbleModel prozessiert BUGS-Modellcode und optionale Konstanten, Daten und Initialwerte. Liefert ein NIMBLE-Modell zurück.
# dieser Schritt ist bei den hier getesteten BUGS Modellen nicht nötig 
Model_Nimble <- nimbleModel(code = Code_Model, name = "ProcessedModel", constants = Nimble_constants,
                            data = Nimble_data, inits = Nimble_inits) 

Model_Nimble$initializeInfo()



##### Simulation
mcmc.out <- nimbleMCMC(code = Code_Model, constants = Nimble_constants,
                       data = Nimble_data, inits = Nimble_inits,
                       nchains = 3, niter = 20000, nburnin = 10000,
                       summary = TRUE, WAIC = F,
                       monitors = c("rhat",'dev', "totresdev", "T", "d", "sd"))
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


## Berechnung CrI von sd
sd_1 <- quantile(mcmc.out$samples[["chain1"]][,"sd"] , c(0.025, 0.975))
sd_2 <- quantile(mcmc.out$samples[["chain2"]][,"sd"] , c(0.025, 0.975))
sd_3 <- quantile(mcmc.out$samples[["chain3"]][,"sd"] , c(0.025, 0.975))

# CrI von sd
(sd_1 + sd_2 + sd_3)/3





# Teil 3: Leverage Plot und nachträgliche Berechnung von pD -----------------




out_lePlo <- capture.output( mcmc.out$summary)
cat("Hilf_lePlo", out_lePlo, file="Hilf.txt", sep="\n", append=TRUE)


Hilf_data = read.table("Hilf.txt", sep = "", header=F, skip=298, nrows=88)
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
            main="Leverage plot for the random effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), 
            ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F, grid=TRUE )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 



#### Plot of model
#directed acyclic graph
#durch igraph
Model_Nimble$plotGraph() # Anweisung geht nicht bei nimbleMCMC



########## ########## ##########  Simulation beendet ########## ########## ##########

