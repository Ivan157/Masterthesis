########## ########## Simulation Blocker Beispiel mit Fixed Effects mit dem rjags package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden


# Teil Simulation mit JAGS ------------------------------------------------------------------


##### Clear data
rm(list=ls())


##### load libraries
library(rjags)
library(coda)
library(random)
library(matrixStats) # zusätzl Paket, berechnet Median
load.module("glm") 
load.module("lecuyer")
load.module("dic")

#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit rjags/Blocker")


##### Read the data into R.
data = as.matrix(read.table("Blocker_Data_neu sortiert.txt", sep = "", header=F))
head(data) # Shows the first six entries
data2  = read.table("Data_Blocker_Rest.txt")
head(data2) # Shows the first six entries


##### Values for simulation, prepare dat for JAGS (allocation values from data)
ns <- nrow(data)
# ns # check
nt <- ncol(data[,5:6])
#nt # check

na <- data[,7]
# na  # check
r <- data[,1:2]
# r # Check
n <-  data[,3:4]
# n # Check
t <-  data[,5:6]
# t # Check
dat <- list(ns=ns, nt=nt, na=na, r=r, n=n, t=t)  # names list of numbers 


##### read in inits with chains
inits <- function(){
  #chain 1
  list(d=c( NA, 0), 
       sd=1, 
       mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
       .RNG.name="lecuyer::RngStream", .RNG.seed=1
  )  
  #chain 2
  list(d=c( NA, -1), 
       sd=4, 
       mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
       .RNG.name="lecuyer::RngStream", .RNG.seed=2 
  )
  #chain 3
  list(d=c( NA, 2), 
       sd=2, 
       mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
       .RNG.name="lecuyer::RngStream", .RNG.seed=3
  ) 
}


##### define JAGS model within R
cat("model{ # *** PROGRAM STARTS 
      for(i in 1:ns){ # LOOP THROUGH STUDIES 
          mu[i] ~ dnorm(0,.0001) # vague priors for all trial baselines 
          for (k in 1:na[i]) { # LOOP THROUGH ARMS 
              r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood 
              logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor 
              rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
              dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))     #Deviance contribution 
          } 
          resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
      }    
      totresdev <- sum(resdev[]) #Total Residual Deviance 
      d[1]<-0 # treatment effect is zero for reference treatment 
      for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 

      # Provide estimates of treatment effects T[k] on the natural (probability) scale  
      # Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 

      for (k in 1:nt) { logit(T[k]) <- A + d[k] }
      A ~ dnorm(-2.2, 3.3) }",
    file="Blocker_Fixed.txt")


##### Set up the JAGS model and settings
jags.m <- jags.model( file = "Blocker_Random.txt", data=dat, inits=inits, n.chains=3, n.adapt=1500)


##### Initialization
update(jags.m, 10000) # burn in 


##### run JAGS and save posterior samples
samps_coda <- coda.samples( jags.m, variable.names=c("d[2]", "T[1]", "T[2]" ),   n.iter=30000, DIC=T  ) 
samps_jags <- jags.samples( jags.m, variable.names=c("rhat", "dev",  "totresdev", "deviance"),   n.iter=20000, DIC=T ) 



# Ausgabe posteriore Werte, Berechnung Median und Berechnung DIC --------------------------


##### summarize posterior samples # funktioniert nur mit coda.samples
summary(window(samps_coda, start=10001)) # burnin = 10000
# Anmerkung: mean für totresdev ist D_res 


#### Trace Monitore
plot(samps_coda) 


#### Median-Berechnung
# Median für T1
median(as.matrix(samps_coda[,1]))
# Median für T2
median(as.matrix(samps_coda[,2]))
# Median für d[2]
median(as.matrix(samps_coda[,3]))



##### Berechnung DIC und pD
D_res <-  summary(samps_jags$totresdev[])
D_res 
pD <- var(samps_jags$deviance)/2
pD 
DIC = D_res + pD
DIC





# Berechnung Daten für leverage plot --------------------------------------



#### Berechnung w_ik
samps_jags[["dev"]]
out_dev <- capture.output(samps_jags[["dev"]])
cat("Hilf_dev", out_dev, file="Hilf.txt", sep="\n", append=TRUE)
Hilf = read.table("Hilf.txt", sep = "", header=F, skip=3, nrows=22, colClasses=c("NULL",NA,NA))
Hilf
cbind(Hilf, total = rowMeans(Hilf))
Hilf2 <- cbind(Hilf, total = rowMeans(Hilf))
w_ik <- sqrt(Hilf2[,3])
w_ik_neg <- -w_ik
fertige_Daten_w_ik <- cbind(Hilf2, w_ik_neg, w_ik)
fertige_Daten_w_ik


#### Berechnung leverage_ik
samps_jags[["rhat"]]
out_rhat <- capture.output(samps_jags$rhat)
out_rhat
cat("Hilf_out_rhat", out_rhat, file="Hilf_rhat.txt", sep="\n", append=TRUE)
Hilf_rhat2 = read.table("Hilf_rhat.txt", sep = "", header=F, skip=3, nrows=22, colClasses=c("NULL",NA,NA))
Hilf_rhat2


#### Berechnung leverage
dev_tilde_erst_k1 <- data[,1]*log(data[,1]/Hilf_rhat2[,1])
#dev_tilde_erst_k1
dev_tilde_zweit_k1 <- (data[,3]-data[,1])*log((data[,3]-data[,1])/(data[,3]-Hilf_rhat2[,1]))
#dev_tilde_zweit_k1
dev_tilde_gesamt_k1 <- 2*(dev_tilde_erst_k1+dev_tilde_zweit_k1)
#dev_tilde_gesamt_k1
leverage_k1 <-  fertige_Daten_w_ik[,1] - dev_tilde_gesamt_k1
leverage_k1

dev_tilde_erst_k2 <- data[,2]*log(data[,2]/Hilf_rhat2[,2])
#dev_tilde_erst_k2
dev_tilde_zweit_k2 <- (data[,4]-data[,2])*log((data[,4]-data[,2])/(data[,4]-Hilf_rhat2[,2]))
#dev_tilde_zweit_k2
dev_tilde_gesamt_k2 <- 2*(dev_tilde_erst_k2+dev_tilde_zweit_k2)
#dev_tilde_gesamt_k2
leverage_k2 <-  fertige_Daten_w_ik[,2] - dev_tilde_gesamt_k2
leverage_k2


#### Erzeugen leverage plot
library(car)
scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_k1, leverage_k2),  main="Leverage plot for the fixed effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 




########## ########## ########## Simulation beendet ########## ########## ##########