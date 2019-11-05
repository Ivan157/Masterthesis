########## ########## Simulation Blocker Beispiel mit Fixed Effects mit dem jagsUI package
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 


# Teil Simulation mit JAGS ------------------------------------------------------------------



##### Clear data
rm(list=ls())



##### load libraries
#library(rjags) # jagsUI benötigt dieses Paket
library(lattice) 
#library(coda)  
library(jagsUI)
#library(random)



#####  Sichergehen richtiger Working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/Programm JAGS/Nachrechnen TSD2-Dokument/Nachrechnen mit jagsUI/Blocker")



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
params <- c("d[2]", "T[1]", "T[2]", "totresdev" )


##### read in inits with chains
inits1 <- list(d=c( NA, 0), 
               mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=1)  


inits2 <- list(d=c( NA, -1), 
               mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=2)


inits3 <- list(d=c( NA, 2), 
               mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3),
               .RNG.name="base::Wichmann-Hill", .RNG.seed=3 ) 

# Achtung: "lecuyer::RngStream" funktioniert nicht mit jagsUI

all.inits <- list(inits1, inits2, inits2)



##### define JAGS model within R
cat("model{ # *** PROGRAM STARTS 
      for(i in 1:ns){ # LOOP THROUGH STUDIES 
          mu[i] ~ dnorm(0,.0001) # vague priors for all trial baselines 
          for (k in 1:na[i]) { # LOOP THROUGH ARMS 
              r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood 
              logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor 
              rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
              dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) +  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
          } 
          resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
      }    
      totresdev <- sum(resdev[]) #Total Residual Deviance 
      d[1]<-0 # treatment effect is zero for reference treatment 
      for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 

      # Provide estimates of treatment effects T[k] on the natural (probability) scale  
      # Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 

      for (k in 1:nt) { logit(T[k]) <- A + d[k] }
      A ~ dnorm(-2.2, 3.3) 
      } ",
    file="Blocker_Fixed.txt")



##### Set up the JAGS model and settings

jags.m <- jags(data=dat, 
               inits=all.inits,
               parameters.to.save=params,
               model.file="Blocker_Fixed.txt", 
               n.chains=3,
               n.iter=20000, n.burnin=10000,
               store.data=TRUE) 

#Warning messages:
#1: In process.input(data, parameters.to.save, inits, n.chains, n.iter,  :
#Suppling a list of character strings to the data argument will be deprecated in the next version
#2: In process.input(data, parameters.to.save, inits, n.chains, n.iter,  :
#Suppling a character vector to the data argument will be deprecated in the next version
# => muss demnächst angepasst werden




# Anzeigen der posterioren Werte und des Medians --------------------------


#traceplot(jags.m) # zeigt Abbildungen einzeln nach Eingabe der Entertaste an
jags.View(jags.m)
jags.m
# mit jags-Funktion keine weiterführende Diagnostik möglich: das Objekt jags.m zeigt nur eine Liste von 24 Elementen an.
# ebenso wenig mit jags.basic


# Versuch mit jagsbasic
jagsbasic.m <- jags.basic(data=dat, 
                          inits=all.inits,
                          parameters.to.save=params,
                          model.file="Blocker_Fixed.txt", 
                          n.chains=3,
                          n.iter=20000, n.burnin=10000) 
#ja, funktioniert
#jagsbasic.m
jagsbasic.m_II <- do.call(rbind.data.frame, jagsbasic.m)


#### Berechnung Median
# Median von T1
median(jagsbasic.m_II$`T[1]`)

# Median von T2
median(jagsbasic.m_II$`T[2]`)

# Median von d
median(jagsbasic.m_II$`d[2]`)



# Leverage Plot -----------------------------------------------------------

## Erzeugen Daten für dev und rhat
params_levplot <- c("dev", "rhat")
jagsbasic.m_levplot <- jags.basic(data=dat, 
                                  inits=all.inits,
                                  parameters.to.save=params_levplot,
                                  model.file="Blocker_Fixed.txt", 
                                  n.chains=3,
                                  n.iter=10000, n.burnin=5000) 
jagsbasic.m_levplot <- do.call(rbind.data.frame, jagsbasic.m_levplot)


## Berechnung w_ik
dev_1.1 <- mean(jagsbasic.m_levplot$`dev[1,1]`)
dev_2.1 <- mean(jagsbasic.m_levplot$`dev[2,1]`)
dev_3.1 <- mean(jagsbasic.m_levplot$`dev[3,1]`)
dev_4.1 <- mean(jagsbasic.m_levplot$`dev[4,1]`)
dev_5.1 <- mean(jagsbasic.m_levplot$`dev[5,1]`)
dev_6.1 <- mean(jagsbasic.m_levplot$`dev[6,1]`)
dev_7.1 <- mean(jagsbasic.m_levplot$`dev[7,1]`)
dev_8.1 <- mean(jagsbasic.m_levplot$`dev[8,1]`)
dev_9.1 <- mean(jagsbasic.m_levplot$`dev[9,1]`)
dev_10.1 <- mean(jagsbasic.m_levplot$`dev[10,1]`)
dev_11.1 <- mean(jagsbasic.m_levplot$`dev[11,1]`)
dev_12.1 <- mean(jagsbasic.m_levplot$`dev[12,1]`)
dev_13.1 <- mean(jagsbasic.m_levplot$`dev[13,1]`)
dev_14.1 <- mean(jagsbasic.m_levplot$`dev[14,1]`)
dev_15.1 <- mean(jagsbasic.m_levplot$`dev[15,1]`)
dev_16.1 <- mean(jagsbasic.m_levplot$`dev[16,1]`)
dev_17.1 <- mean(jagsbasic.m_levplot$`dev[17,1]`)
dev_18.1 <- mean(jagsbasic.m_levplot$`dev[18,1]`)
dev_19.1 <- mean(jagsbasic.m_levplot$`dev[19,1]`)
dev_20.1 <- mean(jagsbasic.m_levplot$`dev[20,1]`)
dev_21.1 <- mean(jagsbasic.m_levplot$`dev[21,1]`)
dev_22.1 <- mean(jagsbasic.m_levplot$`dev[22,1]`)
dev_1.2 <- mean(jagsbasic.m_levplot$`dev[1,2]`)
dev_2.2 <- mean(jagsbasic.m_levplot$`dev[2,2]`)
dev_3.2 <- mean(jagsbasic.m_levplot$`dev[3,2]`)
dev_4.2 <- mean(jagsbasic.m_levplot$`dev[4,2]`)
dev_5.2 <- mean(jagsbasic.m_levplot$`dev[5,2]`)
dev_6.2 <- mean(jagsbasic.m_levplot$`dev[6,2]`)
dev_7.2 <- mean(jagsbasic.m_levplot$`dev[7,2]`)
dev_8.2 <- mean(jagsbasic.m_levplot$`dev[8,2]`)
dev_9.2 <- mean(jagsbasic.m_levplot$`dev[9,2]`)
dev_10.2 <- mean(jagsbasic.m_levplot$`dev[10,2]`)
dev_11.2 <- mean(jagsbasic.m_levplot$`dev[11,2]`)
dev_12.2 <- mean(jagsbasic.m_levplot$`dev[12,2]`)
dev_13.2 <- mean(jagsbasic.m_levplot$`dev[13,2]`)
dev_14.2 <- mean(jagsbasic.m_levplot$`dev[14,2]`)
dev_15.2 <- mean(jagsbasic.m_levplot$`dev[15,2]`)
dev_16.2 <- mean(jagsbasic.m_levplot$`dev[16,2]`)
dev_17.2 <- mean(jagsbasic.m_levplot$`dev[17,2]`)
dev_18.2 <- mean(jagsbasic.m_levplot$`dev[18,2]`)
dev_19.2 <- mean(jagsbasic.m_levplot$`dev[19,2]`)
dev_20.2 <- mean(jagsbasic.m_levplot$`dev[20,2]`)
dev_21.2 <- mean(jagsbasic.m_levplot$`dev[21,2]`)
dev_22.2 <- mean(jagsbasic.m_levplot$`dev[22,2]`)
dev_1 <- mean(c(dev_1.1,dev_1.2 ))
dev_2 <- mean(c(dev_2.1,dev_2.2 ))
dev_3 <- mean(c(dev_3.1,dev_3.2 ))
dev_4 <- mean(c(dev_4.1,dev_4.2 ))
dev_5 <- mean(c(dev_5.1,dev_5.2 ))
dev_6 <- mean(c(dev_6.1,dev_6.2 ))
dev_7 <- mean(c(dev_7.1,dev_7.2 ))
dev_8 <- mean(c(dev_8.1,dev_8.2 ))
dev_9 <- mean(c(dev_9.1,dev_9.2 ))
dev_10 <- mean(c(dev_10.1,dev_10.2 ))
dev_11 <- mean(c(dev_11.1,dev_11.2 ))
dev_12 <- mean(c(dev_12.1,dev_12.2 ))
dev_13 <- mean(c(dev_13.1,dev_13.2 ))
dev_14 <- mean(c(dev_14.1,dev_14.2 ))
dev_15 <- mean(c(dev_15.1,dev_15.2 ))
dev_16 <- mean(c(dev_16.1,dev_16.2 ))
dev_17 <- mean(c(dev_17.1,dev_17.2 ))
dev_18 <- mean(c(dev_18.1,dev_18.2 ))
dev_19 <- mean(c(dev_19.1,dev_19.2 ))
dev_20 <- mean(c(dev_20.1,dev_20.2 ))
dev_21 <- mean(c(dev_21.1,dev_21.2 ))
dev_22 <- mean(c(dev_22.1,dev_22.2 ))

dev_hilf <- rbind(dev_1, dev_2, dev_3, dev_4, dev_5, dev_6, dev_7, dev_8, dev_9, dev_10, dev_11, dev_12, dev_13, dev_14, dev_15, dev_16, dev_17, dev_18, dev_19, dev_20, dev_21, dev_22 )

w_ik <- sqrt(dev_hilf[,1])
w_ik_neg <- -w_ik

fertige_Daten_w_ik <- cbind( w_ik_neg, w_ik)
#fertige_Daten_w_ik



#### Berechnung leverage_ik
dev_bar_k1 <- rbind(dev_1.1, dev_2.1 , dev_3.1, dev_4.1, dev_5.1 , dev_6.1 , dev_7.1 , dev_8.1, dev_9.1 , dev_10.1 , dev_11.1, dev_12.1 , dev_13.1 , dev_14.1 , dev_15.1 , dev_16.1 , dev_17.1 , dev_18.1, dev_19.1, dev_20.1, dev_21.1, dev_22.1)
#dev_bar_k1
dev_bar_k2 <- rbind(dev_1.2, dev_2.2 , dev_3.2, dev_4.2, dev_5.2 , dev_6.2 , dev_7.2 , dev_8.2, dev_9.2 , dev_10.2 , dev_11.2, dev_12.2 , dev_13.2 , dev_14.2 , dev_15.2 , dev_16.2 , dev_17.2 , dev_18.2, dev_19.2, dev_20.2, dev_21.2, dev_22.2)
#dev_bar_k2

rhat_1.1 <- mean(jagsbasic.m_levplot$`rhat[1,1]`)
rhat_2.1 <- mean(jagsbasic.m_levplot$`rhat[2,1]`)
rhat_3.1 <- mean(jagsbasic.m_levplot$`rhat[3,1]`)
rhat_4.1 <- mean(jagsbasic.m_levplot$`rhat[4,1]`)
rhat_5.1 <- mean(jagsbasic.m_levplot$`rhat[5,1]`)
rhat_6.1 <- mean(jagsbasic.m_levplot$`rhat[6,1]`)
rhat_7.1 <- mean(jagsbasic.m_levplot$`rhat[7,1]`)
rhat_8.1 <- mean(jagsbasic.m_levplot$`rhat[8,1]`)
rhat_9.1 <- mean(jagsbasic.m_levplot$`rhat[9,1]`)
rhat_10.1 <- mean(jagsbasic.m_levplot$`rhat[10,1]`)
rhat_11.1 <- mean(jagsbasic.m_levplot$`rhat[11,1]`)
rhat_12.1 <- mean(jagsbasic.m_levplot$`rhat[12,1]`)
rhat_13.1 <- mean(jagsbasic.m_levplot$`rhat[13,1]`)
rhat_14.1 <- mean(jagsbasic.m_levplot$`rhat[14,1]`)
rhat_15.1 <- mean(jagsbasic.m_levplot$`rhat[15,1]`)
rhat_16.1 <- mean(jagsbasic.m_levplot$`rhat[16,1]`)
rhat_17.1 <- mean(jagsbasic.m_levplot$`rhat[17,1]`)
rhat_18.1 <- mean(jagsbasic.m_levplot$`rhat[18,1]`)
rhat_19.1 <- mean(jagsbasic.m_levplot$`rhat[19,1]`)
rhat_20.1 <- mean(jagsbasic.m_levplot$`rhat[20,1]`)
rhat_21.1 <- mean(jagsbasic.m_levplot$`rhat[21,1]`)
rhat_22.1 <- mean(jagsbasic.m_levplot$`rhat[22,1]`)
rhat_1.2 <- mean(jagsbasic.m_levplot$`rhat[1,2]`)
rhat_2.2 <- mean(jagsbasic.m_levplot$`rhat[2,2]`)
rhat_3.2 <- mean(jagsbasic.m_levplot$`rhat[3,2]`)
rhat_4.2 <- mean(jagsbasic.m_levplot$`rhat[4,2]`)
rhat_5.2 <- mean(jagsbasic.m_levplot$`rhat[5,2]`)
rhat_6.2 <- mean(jagsbasic.m_levplot$`rhat[6,2]`)
rhat_7.2 <- mean(jagsbasic.m_levplot$`rhat[7,2]`)
rhat_8.2 <- mean(jagsbasic.m_levplot$`rhat[8,2]`)
rhat_9.2 <- mean(jagsbasic.m_levplot$`rhat[9,2]`)
rhat_10.2 <- mean(jagsbasic.m_levplot$`rhat[10,2]`)
rhat_11.2 <- mean(jagsbasic.m_levplot$`rhat[11,2]`)
rhat_12.2 <- mean(jagsbasic.m_levplot$`rhat[12,2]`)
rhat_13.2 <- mean(jagsbasic.m_levplot$`rhat[13,2]`)
rhat_14.2 <- mean(jagsbasic.m_levplot$`rhat[14,2]`)
rhat_15.2 <- mean(jagsbasic.m_levplot$`rhat[15,2]`)
rhat_16.2 <- mean(jagsbasic.m_levplot$`rhat[16,2]`)
rhat_17.2 <- mean(jagsbasic.m_levplot$`rhat[17,2]`)
rhat_18.2 <- mean(jagsbasic.m_levplot$`rhat[18,2]`)
rhat_19.2 <- mean(jagsbasic.m_levplot$`rhat[19,2]`)
rhat_20.2 <- mean(jagsbasic.m_levplot$`rhat[20,2]`)
rhat_21.2 <- mean(jagsbasic.m_levplot$`rhat[21,2]`)
rhat_22.2 <- mean(jagsbasic.m_levplot$`rhat[22,2]`)

rhat_k1 <- rbind(rhat_1.1, rhat_2.1,  rhat_3.1,  rhat_4.1, rhat_5.1, rhat_6.1, rhat_7.1, rhat_8.1, rhat_9.1, rhat_10.1, rhat_11.1, rhat_12.1, rhat_13.1, rhat_14.1 , rhat_15.1 , rhat_16.1, rhat_17.1, rhat_18.1, rhat_19.1, rhat_20.1 , rhat_21.1 , rhat_22.1)
#rhat_k1

dev_tilde_erst_k1 <- data[,1]*log(data[,1]/rhat_k1[,1])
#dev_tilde_erst_k1
dev_tilde_zweit_k1 <- (data[,3]-data[,1])*log((data[,3]-data[,1])/(data[,3]-rhat_k1))
#dev_tilde_zweit_k1
dev_tilde_gesamt_k1 <- 2*(dev_tilde_erst_k1+dev_tilde_zweit_k1)

leverage_k1 <-  dev_bar_k1[,1] - dev_tilde_gesamt_k1

#leverage_k1


rhat_k2 <- rbind(rhat_1.2, rhat_2.2,  rhat_3.2,  rhat_4.2, rhat_5.2, rhat_6.2, rhat_7.2, rhat_8.2, rhat_9.2, rhat_10.2, rhat_11.2, rhat_12.2, rhat_13.2, rhat_14.2 , rhat_15.2 , rhat_16.2, rhat_17.2, rhat_18.2, rhat_19.2, rhat_20.2 , rhat_21.2 , rhat_22.2)
#rhat_k2

dev_tilde_erst_k2 <- data[,2]*log(data[,2]/rhat_k2[,1])
#dev_tilde_erst_k2
dev_tilde_zweit_k2 <- (data[,4]-data[,2])*log((data[,4]-data[,2])/(data[,4]-rhat_k2[,1]))
#dev_tilde_zweit_k2
dev_tilde_gesamt_k2 <- 2*(dev_tilde_erst_k2+dev_tilde_zweit_k2)

leverage_k2 <-  dev_bar_k2[,1] - dev_tilde_gesamt_k2

#leverage_k2




#### Erzeugen leverage plot
library(car)
scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_k1, leverage_k2),  main="Leverage plot for the Fixed effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 





########## ########## ########## Simulation beendet ########## ########## ##########

