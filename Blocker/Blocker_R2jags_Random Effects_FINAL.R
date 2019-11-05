########## ########## Simulation Blocker Beispiel mit Random Effects mit dem R2jags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden "samps_coda" aktiviert werden.


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
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit R2jags/Blocker")


##### Read the data into R.
#data = read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F)
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
head(data) # Shows the first six entries
data2  = read.table("Data_Blocker_Rest.txt")
head(data2) # Shows the first six entries


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

dat <- list("ns", "nt", "na", "r", "n", "t")  # names list of numbers 


##### Parameter to monitor/save
params <- c("d[2]", "T[1]", "T[2]", "sd" ,  "totresdev" )


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
cat("model{								 # *** PROGRAM STARTS 
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
jags.m <- jags(data=dat,  inits=all.inits, parameters.to.save=params,  n.chains = 3, n.iter = 30000, n.burnin = 10000,
               model.file="Blocker_Random.txt", DIC=TRUE,  jags.module = c("glm","dic") )  
# zusätzlich noch mehrere Argumente standardmäßig dabei, v.a. interessant: DIC, jags.module


#### optional, falls nicht konvergiert:
#jags.m.upd <- autojags(jags.m)



# Ausgabe posteriore Werte, Berechnung Median und Berechnung DIC --------------------------

print(jags.m)

#### Median
jags.m[["BUGSoutput"]][["median"]]



# Berechnung Daten für leverage plot --------------------------------------

#### erneute Simulation für die Erzeugung von dev und rhat
params_lev_plot <- c("dev", "rhat" )
jags.m_levPlot <- jags(data=dat,  inits=all.inits, parameters.to.save=params_lev_plot,  n.chains = 3, n.iter=10000, n.burnin=0,
                       model.file="Blocker_Random.txt", DIC=TRUE,  jags.module = c("glm","dic") )  


#jags.m_levPlot[["BUGSoutput"]][["summary"]]
out_lePlo <- capture.output( jags.m_levPlot[["BUGSoutput"]][["summary"]])
cat("Hilf_lePlo", out_lePlo, file="Hilf.txt", sep="\n", append=TRUE)

Hilf_data = read.table("Hilf.txt", sep = "", header=F, skip=2, nrows=89)
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

Hilf_rhat <- Hilf_data[46:89,2]
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
library(car)
scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_k1, leverage_k2),  main="Leverage plot for the random effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 



# zusätzliche Diagnostik, bei Bedarf aktivieren --------------------------------------------------------------------


#pdf("Blocker_Random_trace.pdf")
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

