########## ########## Simulation Dietary Fat Beispiel mit Random Effects mit STAN ########## ########## 
########## Die Working Directory muss auf Ihre Bedürfnisse angepasst werden 




##### Clear data
rm(list=ls())



##### Setting working directory
setwd("C:/Users/IvanB/Desktop/Masterarbeit/Statistische Programme und Gibbs Sampler/STAN/Nachrechnen TSD2/Ex3 Diabetes")



##### Requiering stan
library("rstan")
library("rstantools")
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')



##### Read in data
data = read.csv2("Diabetes_Data_neu sortiert_II_ohneNA.csv", header=TRUE, sep = ";", quote = "\"",  dec = ",", fill = TRUE, comment.char = "")
head(data) # Shows the first six entries
#data2  = read.table("Diabetes_Data_Rest.txt")
#head(data2) # Shows the first six entries



##### Assignment data to stan
NO =nrow(data)
#NO =48
NT=max(data$Treatment_t )
NS=max(data$Studie)
time=data$time
r=data$Erfolge_r
t=data$Treatment_t
s=data$Studie
base=data$Treatment_t
n=data$Gesamtanzahl_n

data_list <- list(NO=NO, NT=NT, NS=NS,  time=time, r=r, t=t,  s=s,base=base, n=n)



##### Read in inits
# Anmerkung: da cloglog als Link: wurden die Inits von STAN generiert.
# Die manuellen Inits sind zwar drin, können aber rausgelassen werden

inits1 <-  list(d=c(NA,0,0,0,0,0), 
                sd=1,  
                mu=c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0), 
               
                
                delta=  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0 ))
          # die anderen Inits müssen nach der Art von Inits1 angepasst werden
                
          #      delta= c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
          #               0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
          #               NA, 0, NA))  # structure(.Data= , .Dim=c(22, 3))
  # mu=c(0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0), 
#  A=c(0, 0, 0, 0, 0, 0) ,

inits2 <-   list(d=c(NA,-1,4,-1,2,3), 
               sd=3,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,0,0,   0,1,0,0,0,  1,1), 
               A=1 ,
               delta= c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
                                         0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
                                         NA, 0, NA)) # structure(.Data= , .Dim=c(22, 3))

inits3 <-   list(d=c(NA,1,4,-3,-2,3),  
               
               sd=4.5,  
               mu=c(1,1,0,1,0,    0,1,0,0,0,    1,1,0,-2,0,   0,1,0,-2,0,  1,1), 
               A=2 ,
               
              delta=  c(NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 
                                 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0,0, NA, 0, 0, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, 
                              NA, 0, NA)) # structure(.Data= , .Dim=c(22, 3))

all.inits <- list(inits1, inits2, inits3) 
#function(chain_id = 1) {



##### Compiling
m <- stan_model('Model_Random_final.stan')
m <- stan_model('Model_Random_final.stan')



stan_samples <- sampling(m, data = data_list, iter=20000, verbose=T, chain=3, init_r=0.1) 
#stan_samples <- sampling(m, data = data_list, iter=20000, verbose=T, chain=4, control = list(adapt_delta = 0.99)) # dauert zu lange
#stan_samples <- sampling(m, data = data_list, i





# Ausgabe posteriore Werte, Berechnung Median --------------------------




##### Plotting and summarizing the posterior distribution
stan_samples # = print(stan_samples)
plot(stan_samples)
Stan_summary <- summary(stan_samples, pars = c("d[2]", "d_II[2]", "d[3]",  "d[4]", "d[5]", "d[6]", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]", "T[6]",  "s_d" ,  "totalresdev"), probs = c(0.025, 0.975))$summary 
Stan_summary



##### Additional Lines for Median
Median_d2 <- median(as.matrix(stan_samples, pars = c("d_II[2]")))
Median_d2
Median_d3 <- median(as.matrix(stan_samples, pars = c("d_II[3]")))
Median_d3
Median_d4 <- median(as.matrix(stan_samples, pars = c("d_II[4]")))
Median_d4
Median_d5 <- median(as.matrix(stan_samples, pars = c("d_II[5]")))
Median_d5
Median_d6 <- median(as.matrix(stan_samples, pars = c("d_II[6]")))
Median_d6

Median_T1 <- median(as.matrix(stan_samples, pars = c("T[1]")))
Median_T1
Median_T2 <- median(as.matrix(stan_samples, pars = c("T[2]")))
Median_T2
Median_T3 <- median(as.matrix(stan_samples, pars = c("T[3]")))
Median_T3
Median_T4 <- median(as.matrix(stan_samples, pars = c("T[4]")))
Median_T4
Median_T5 <- median(as.matrix(stan_samples, pars = c("T[5]")))
Median_T5
Median_T6 <- median(as.matrix(stan_samples, pars = c("T[6]")))
Median_T6

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





# manuelle Berechnung von pD ----------------------------------------------




SingeValues_dev <- summary(stan_samples, pars = c("dev[1]", "dev[2]", "dev[3]", "dev[4]", "dev[5]", "dev[6]", "dev[7]", "dev[8]", "dev[9]", "dev[10]", "dev[11]", "dev[12]", "dev[13]", "dev[14]", "dev[15]", "dev[16]", "dev[17]", "dev[18]", "dev[19]", "dev[20]", "dev[21]" , "dev[22]", "dev[23]", "dev[24]", "dev[25]", "dev[26]", "dev[27]", "dev[28]", "dev[29]", "dev[30]", "dev[31]", "dev[32]", "dev[33]", "dev[34]", "dev[35]", "dev[36]", "dev[37]", "dev[38]", "dev[39]", "dev[40]", "dev[41]", "dev[42]", "dev[43]", "dev[44]", "dev[45]", "dev[46]", "dev[47]", "dev[48]"))$summary
out_lePlo <- capture.output(   SingeValues_dev)
cat("Hilf_pD", out_lePlo, file="Hilf_pD.txt", sep="\n", append=TRUE)
Hilf_dev_pD <-  read.table("Hilf_pD.txt", sep = "", header=F, skip=2, nrows=21)
#dev ist Std-Abweichung jedes einzelnen Werts
# insg 48 Werte
Var_manuell <-  sum(Hilf_dev_pD[,3])^2/48
Var_manuell
pD_manuell <- Var_manuell/2
pD_manuell

#Var_manuell <-  sum(Hilf_dev_pD[,4])/48
#Var_manuell


########## ########## ########## Simulation beendet ########## ########## ##########

