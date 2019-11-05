########## ########## Simulation Blocker Beispiel mit fixed effects mit STAN ########## ########## 



##### Clear data
rm(list=ls())



#### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/STAN/Nachrechnen TSD2/DietaryFat")



#### Requiering stan
library("rstan")
library("rstantools")
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



#### Read in data
data = read.csv2("DietaryFat_Data_neu sortiert_II_ohneNA.csv", header=TRUE, sep = ";", quote = "\"",  dec = ",", fill = TRUE, comment.char = "")
#data



#### Assignment data to stan
NO =nrow(data)
NT=max(data$Treatment_t )
NS=max(data$Studie)
E=data$Explosure_time_E
r=data$Erfolge_r
t=data$Treatment_t
s=data$Studie
base=data$Treatment_t

data_list <- list(NO=NO, NT=NT, NS=NS,  E=E, r=r, t=t,  s=s,base=base)



#### Read in inits
inits1 <- function(chain_id = 1) { #*
  list(d=c( NA, 0, 0), 
       sd=1, 
       mu=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
}

inits2 <- function(chain_id = 2) {
  list(d=c( NA, -1, -1), 
       sd=4, 
       mu=c(-3, -3, -3, -3, -3, -3, -3, -3, -3, -3))
}

inits3 <- function(chain_id = 3) { #*
  list(d=c( NA, 2, 2), 
       sd=2, 
       mu=c(3, 5, 1, 3, 7, 3, 4, 3, 3, 0)) 
}

all.inits <- list(inits1, inits2, inits3) 
#all.inits

#* Hinweis:
# Chain 3: Rejecting initial value:
# Chain 3:   Log probability evaluates to log(0), i.e. negative infinity.
# Chain 3:   Stan can't start sampling from this initial value.
# => negative Werte zu positiv geändert
# inits für Kette 1+2+4 auch. Aber augenscheinlich keine Probleme



# Compiling
m <- stan_model('Model_Fixed.stan')
m <- stan_model('Model_Fixed.stan')



# Simulation 
#stan_samples <- sampling(m, data = data_list, iter=20000, verbose=T, chain=4, control = list(adapt_delta = 0.99)) # dauert zu lange
stan_samples <- sampling(m, data = data_list, iter=30000, verbose=T, chain=4, control = list(adapt_delta = 0.90)) # optimalere Lösung



#### Plotting and summarizing the posterior distribution
stan_samples # = print(stan_samples)
plot(stan_samples)
Stan_summary <- summary(stan_samples, pars = c("d[2]","d_II[2]", "T[1]", "T[2]", "totalresdev", "dev", "theta"), probs = c(0.025, 0.975))$summary 
Stan_summary



#### Additional Lines for Median
Median_d2 <- median(as.matrix(stan_samples, pars = c("d_II[2]")))
Median_d2
Median_T1 <- median(as.matrix(stan_samples, pars = c("T[1]")))
Median_T1
Median_T2 <- median(as.matrix(stan_samples, pars = c("T[2]")))
Median_T2






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





# manuelle Berechnung von pD ----------------------------------------------




SingeValues_dev <- summary(stan_samples, pars = c("dev[1]", "dev[2]", "dev[3]", "dev[4]", "dev[5]", "dev[6]", "dev[7]", "dev[8]", "dev[9]", "dev[10]", "dev[11]", "dev[12]", "dev[13]", "dev[14]", "dev[15]", "dev[16]", "dev[17]", "dev[18]", "dev[19]", "dev[20]", "dev[21]"))$summary
out_lePlo <- capture.output(   SingeValues_dev)
cat("Hilf_pD", out_lePlo, file="Hilf_pD.txt", sep="\n", append=TRUE)
Hilf_dev_pD <-  read.table("Hilf_pD.txt", sep = "", header=F, skip=2, nrows=21)
#dev ist Std-Abweichung jedes einzelnen Werts
# insg 21 Werte
Var_manuell <-  sum(Hilf_dev_pD[,2])^2/21
Var_manuell
pD_manuell <- Var_manuell/2
pD_manuell




########## ########## ########## Simulation beendet ########## ########## ##########

