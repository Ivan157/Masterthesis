########## ########## Simulation Blocker Beispiel mit Random Effects mit dem runjags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden.


# Teil Simulation mit JAGS ------------------------------------------------------------------


##### Clear data
rm(list=ls())


##### load libraries
library(rjags) 
library(runjags)
library(random)
library(coda)
load.module("glm") 
load.module("lecuyer")
load.module("dic")


#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit runjags/Blocker")


##### Read the data into R.
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
#head(data) # Shows the first six entries
data2  = read.table("Data_Blocker_Rest.txt")
#head(data2) # Shows the first six entries


##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <-  nrow(data)
# ns # check
nt <-  ncol(data[,5:6])
#nt # check

na <- data[,7]
# na  # check
r <- data[,1:2]
# r # Check
n <-  data[,3:4]
# n # Check
t <-  data[,5:6]
# t # Check


##### read in inits with chains
inits1 <- list(d=c( NA, 0), 
               sd=1, 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="lecuyer::RngStream", .RNG.seed=1)  

inits2 <- list(d=c( NA, -1), 
               sd=4, 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=2)

inits3 <- list(d=c( NA, 2), 
               sd=2, 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
               .RNG.name="lecuyer::RngStream", .RNG.seed=3 ) 

all.inits <- list(inits1, inits2, inits2)


##### define JAGS model within R
model <- cat("model{								 # *** PROGRAM STARTS 
        for(i in 1:ns){								 # LOOP THROUGH STUDIES 
            w[i,1] <- 0 								# adjustment for multi-arm trials is zero for control arm 
            delta[i,1] <- 0 								# treatment effect is zero for control arm 
            mu[i] ~ dnorm(0,.0001) 				# vague priors for all trial baselines 
            for (k in 1:na[i]) { 								# LOOP THROUGH ARMS 
                r[i,k] ~ dbin(p[i,k],n[i,k]) 				# binomial likelihood 
                logit(p[i,k]) <- mu[i] + delta[i,k]				 # model for linear predictor 
                rhat[i,k] <- p[i,k] * n[i,k] 				# expected value of the numerators 
                dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))+  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))  #Deviance contribution 
            } 
            resdev[i] <- sum(dev[i,1:na[i]]) 				# summed residual deviance contribution for this trial 
            for (k in 2:na[i]) { 								# LOOP THROUGH ARMS 
                delta[i,k] ~ dnorm(md[i,k],taud[i,k]) 				# trial-specific LOR distributions 
                md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] 				# mean of LOR distributions (with multi-arm trial correction) 
                taud[i,k] <- tau *2*(k-1)/k 								# precision of LOR distributions (with multi-arm trial correction) 
                w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) 				# adjustment for multi-arm RCTs 
                sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) 				# cumulative adjustment for multi-arm trials 
            } 
        }    
    totresdev <- sum(resdev[]) 								#Total Residual Deviance 
    d[1] <- 0												 # treatment effect is zero for reference treatment 
    sd ~ dunif(0,5) 												# vague prior for between-trial SD.
    tau <- pow(sd,-2) 								# between-trial precision = (1/between-trial variance) 

    # Provide estimates of treatment effects T[k] on the natural (probability) scale  
    # Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 

    for (k in 2:nt){ d[k] ~ dnorm(0,0.0001) }
    for (k in 1:nt) { logit(T[k]) <- A + d[k] }
    A ~ dnorm(-2.2, 3.3) 
    } ",
    file="Blocker_Random.txt")



##### Set up the JAGS model and settings
jags.m <- run.jags(model="Blocker_Random.txt", monitor=c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , n.chains=3, inits=all.inits, burnin = 10000, sample = 30000, adapt = 1500)


#### optional, falls nicht konvergiert:
#jags.m <- autorun.jags(model="Blocker_Random.txt", monitor=c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev", "deviance", "pd", "dic" , "full.pd"  ),
#                   data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , n.chains=3, 
#                   inits=all.inits, burnin = 10000, sample = 30000, adapt = 1500, , max.time="1hr")


# Ausgabe posteriore Werte --------------------------------------


print(jags.m)
summary(jags.m$mcmc) # für 2.5 - 97.5  CrI

# Berechnung Daten für leverage plot --------------------------------------

#### erneute Simulation für die Erzeugung von dev und rhat
jags.m_levPlot <- run.jags(model="Blocker_Random.txt", monitor=c("dev", "rhat"), 
                           data=list("ns"=ns, "nt"=nt, "na"=na, "r"=r, "n"=n, "t"=t) , 
                           n.chains=3, inits=all.inits, burnin = 10000, sample = 20000, adapt = 1000)


jags.m_levPlot_summary <- summary(jags.m_levPlot$mcmc)
#jags.m_levPlot_summary


out_lePlo <- capture.output(   jags.m_levPlot_summary)
cat("Hilf_lePlo", out_lePlo, file="Hilf.txt", sep="\n", append=TRUE)
Hilf_data = read.table("Hilf.txt", sep = "", header=F, skip=11, nrows=88)
#Hilf_data


#### Berechnung w_ik

Hilf_dev <- Hilf_data[1:44,2]
Hilf_dev
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
#fertige_Daten_w_ik

# manuelle Berechnung von pD
#dev ist Std-Abweichung jedes einzelnen Werts
# insg 22 Wertepaare
Var_manuell <-  sum(Hilf_dev)^2/44
pD_manuell <- Var_manuell/2
pD_manuell

#### Berechnung leverage_ik

Hilf_rhat <- Hilf_data[45:88,2]
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

#leverage_k1

dev_tilde_erst_k2 <- data[,2]*log(data[,2]/Hilf_rhat_k2)
#dev_tilde_erst_k2
dev_tilde_erst_k2 <- data[,2]*log(data[,2]/Hilf_rhat_k2)
#dev_tilde_erst_k2
dev_tilde_zweit_k2 <- (data[,4]-data[,2])*log((data[,4]-data[,2])/(data[,4]-Hilf_rhat_k2))
#dev_tilde_zweit_k2
dev_tilde_gesamt_k2 <- 2*(dev_tilde_erst_k2+dev_tilde_zweit_k2)
#dev_tilde_gesamt_k2
leverage_k2 <-  fertige_Daten_w_ik[,2] - dev_tilde_gesamt_k2

#leverage_k2


#### Erzeugen leverage plot
library(car)
scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_k1, leverage_k2),  main="Leverage plot for the random effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 





########## ########## ########## Simulation beendet ########## ########## ##########
