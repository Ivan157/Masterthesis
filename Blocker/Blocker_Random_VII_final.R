########## ########## Simulation Blocker Beispiel mit random effects mit STAN ########## ########## 



##### Clear data
rm(list=ls())



#### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/STAN/Nachrechnen TSD2/Blocker")


#### Requiering stan
library("rstan")
library("rstantools")
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



#### Read in data
data = read.csv2("Blocker_Data_neu sortiert_VI.csv", header=TRUE, sep = ";", quote = "\"",  dec = ",", fill = TRUE, comment.char = "")
#data



#### Assignment data to stan
NO =nrow(data)
NT=max(data$ï..Treatment_t )
NS=max(data$Studie)
n=data$Gesamtanzahl_n
r=data$Erfolge_r
t=data$ï..Treatment_t
s=data$Studie
base=data$ï..Treatment_t

data_list <- list(NO=NO, NT=NT, NS=NS,  n=n, r=r, t=t,  s=s,base=base)



#### Read in inits
inits1 <- function(chain_id = 1) {
  list(d=c( NA, 0), 
       sd=1, 
       mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
}

inits2 <- function(chain_id = 2) {
  list(d=c( NA, -1), 
       sd=4, 
       mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3))
}

inits3 <- function(chain_id = 3) {
  list(d=c( NA, 2), 
       sd=2, 
       mu=c(-3, 5, -1, -3, 7, -3, -4, -3, -3, 0, -3, -3,0, 3, 5, -3, -3, -1, -3, -7, -3, -3)) 
}

all.inits <- list(inits1, inits2, inits3) 
#all.inits



# Compiling
m <- stan_model('Model_Random_VIII_final.stan')
m <- stan_model('Model_Random_VIII_final.stan')

# Simulation 
stan_samples <- sampling(m, data = data_list, iter=40000, verbose=T, chain=4) # !! iter nachher erhöhen



#### Plotting and summarizing the posterior distribution
stan_samples # = print(stan_samples)
plot(stan_samples)
Stan_summary <- summary(stan_samples, pars = c("d[2]","d_II[2]", "T[1]", "T[2]", "s_d", "totalresdev"), probs = c(0.025, 0.975))$summary 
Stan_summary



#### Additional Lines for Median
Median_d2 <- median(as.matrix(stan_samples, pars = c("d_II[2]")))
Median_d2
Median_T1 <- median(as.matrix(stan_samples, pars = c("T[1]")))
Median_T1
Median_T2 <- median(as.matrix(stan_samples, pars = c("T[2]")))
Median_T2
Median_sd <- median(as.matrix(stan_samples, pars = c("s_d")))
Median_sd





# Section for Convergence Diagnostic --------------------------------------------------


# konkreter Vergleich mit BUGS nicht möglich ->kein pD in dem Sinne
# DIC (und damit pD) ist veraltet
# => loo()-Fkt und WAIC
#(pD Code für Python)
#allerdings andere Konvergenz - Diagnostika:
# Diagnostik mir rstan Paket
sampler_params <- get_sampler_params(stan_samples, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
# each chain separately
lapply(sampler_params, summary, digits = 2)
Stan_summary_lp__ <- summary(stan_samples, pars = c("lp__"), probs = c(0.025, 0.975))$summary # sigmasq_delta entspricht sd, nachher ändern
Stan_summary_lp__
# weitere Möglichkeit: Package 'shinystan'





# Section for leverage plot -----------------------------------------------



#### Read in single values for dev and rhat

SingeValues_dev <- summary(stan_samples, pars = c("dev[1]", "dev[2]", "dev[3]", "dev[4]", "dev[5]", "dev[6]", "dev[7]", "dev[8]", "dev[9]", "dev[10]", "dev[11]", "dev[12]", "dev[13]", "dev[14]", "dev[15]", "dev[16]", "dev[17]", "dev[18]", "dev[19]", "dev[20]", "dev[21]", "dev[22]", "dev[23]", "dev[24]", "dev[25]", "dev[26]", "dev[27]", "dev[28]", "dev[29]", "dev[30]", "dev[31]", "dev[32]", "dev[33]", "dev[34]", "dev[35]", "dev[36]", "dev[37]", "dev[38]", "dev[39]", "dev[40]", "dev[41]", "dev[42]", "dev[43]", "dev[44]"
))$summary
#SingeValues_dev
SingeValues_rhat <- summary(stan_samples, pars = c("rhat[1]", "rhat[2]", "rhat[3]", "rhat[4]", "rhat[5]", "rhat[6]", "rhat[7]", "rhat[8]", "rhat[9]", "rhat[10]", "rhat[11]", "rhat[12]", "rhat[13]", "rhat[14]", "rhat[15]", "rhat[16]", "rhat[17]", "rhat[18]", "rhat[19]", "rhat[20]", "rhat[21]", "rhat[22]", "rhat[23]", "rhat[24]", "rhat[25]", "rhat[26]", "rhat[27]", "rhat[28]", "rhat[29]", "rhat[30]", "rhat[31]", "rhat[32]", "rhat[33]", "rhat[34]", "rhat[35]", "rhat[36]", "rhat[37]", "rhat[38]", "rhat[39]", "rhat[40]", "rhat[41]", "rhat[42]", "rhat[43]", "rhat[44]"
))$summary
#SingeValues_rhat

out_lePlo <- capture.output(   SingeValues_dev)
cat("Hilf_lePlo", out_lePlo, file="Hilf_lePlo.txt", sep="\n", append=TRUE)
out_lePlo <- capture.output(   SingeValues_rhat)
cat("Hilf_lePlo", out_lePlo, file="Hilf_lePlo.txt", sep="\n", append=TRUE)



Hilf_dev_I <-  read.table("Hilf_lePlo.txt", sep = "", header=F, skip=2, nrows=22)
#Hilf_dev_I
Hilf_dev_II <-  read.table("Hilf_lePlo.txt", sep = "", header=F, skip=24, nrows=22)
#Hilf_dev_II
Hilf_rhat_I <-  read.table("Hilf_lePlo.txt", sep = "", header=F, skip=48, nrows=22)
#Hilf_rhat_I
Hilf_rhat_II <-  read.table("Hilf_lePlo.txt", sep = "", header=F, skip=70, nrows=22)
#Hilf_rhat_II



#### Berechnung w_ik
Hilf_dev_1 <- cbind(Hilf_dev_I[,2], Hilf_dev_II[,2])
#Hilf_dev_1
Hilf_dev_2 <- cbind(Hilf_dev_1, total = rowMeans(Hilf_dev_1))
#Hilf_dev_2
w_ik <- sqrt(Hilf_dev_2[,3])
w_ik_neg <- -w_ik

fertige_Daten_w_ik <- cbind(Hilf_dev_2, w_ik_neg, w_ik)
#fertige_Daten_w_ik



#### Berechnung leverage_ik
dev_tilde_erst_I <- data[1:22,2]*log(data[1:22,2]/Hilf_rhat_I[,2])
#dev_tilde_erst_I
dev_tilde_zweit_I <- (data[1:22,3]-data[1:22,2])*log((data[1:22,3]-data[1:22,2])/(data[1:22,3]-Hilf_rhat_I[,2]))
#dev_tilde_zweit_I
dev_tilde_gesamt_I <- 2*(dev_tilde_erst_I+dev_tilde_zweit_I)
#dev_tilde_gesamt_I

dev_tilde_erst_II <- data[23:44,2]*log(data[23:44,2]/Hilf_rhat_II[,2])
#dev_tilde_erst_II
dev_tilde_zweit_II <- (data[23:44,3]-data[23:44,2])*log((data[23:44,3]-data[23:44,2])/(data[23:44,3]-Hilf_rhat_II[,2]))
#dev_tilde_zweit_II
dev_tilde_gesamt_II <- 2*(dev_tilde_erst_II+dev_tilde_zweit_II)
#dev_tilde_gesamt_II

leverage_I <-  fertige_Daten_w_ik[,1] - dev_tilde_gesamt_I
leverage_II <-  fertige_Daten_w_ik[,2] - dev_tilde_gesamt_II
#leverage_I
#leverage_II



#### Erzeugen leverage plot
library(car)
scatterplot(c(fertige_Daten_w_ik[,"w_ik_neg"], fertige_Daten_w_ik[,"w_ik"]), c(leverage_I, leverage_II),  main="Leverage plot for the random effects model", xlim=c(-3,3), ylim=c(0,4.5), xlab=expression('w'[ik]), ylab=expression('leverage'[ik]), regLine =F, smooth=F, boxplots=F )                
curve(-x^2+1, from=-3, to=3, col="red", lty="solid", add=T) 
curve(-x^2+2, from=-3, to=3, col="green", lty="dashed", add=T) 
curve(-x^2+3, from=-3, to=3, col="blueviolet", lty="dotted", add=T) 
curve(-x^2+4, from=-3, to=3, col="blue", lty="dotdash", add=T) 



# manuelle Berechnung von pD ----------------------------------------------


#dev ist Std-Abweichung jedes einzelnen Werts
# insg 22*2 Wertepaare, also 44 Werte
Hilf_dev_pD <- read.table("Hilf_lePlo.txt", sep = "", header=F, skip=2, nrows=44)
#Hilf_dev_pD
Var_manuell <-  sum(Hilf_dev_pD[,2])^2/44
pD_manuell <- Var_manuell/2
pD_manuell



########## ########## ########## Simulation beendet ########## ########## ##########
