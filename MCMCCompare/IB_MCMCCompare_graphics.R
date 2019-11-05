##### Datei zum Erstellen Grafiken

setwd("C:/Users/IvanB/Desktop/Masterarbeit/Ergebnisse/zusätzliche Experimente/MCMC Vergleiche")

library(ggplot2)
library("readxl")
library("RColorBrewer")


# Part 1 ------------------------------------------------------------------


##### Grafik Verteilungen

Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

head(Daten)

Verteilungen <- Daten$Verteilung
#Verteilungen
Verteilung <- as.character(Verteilungen)
#Verteilungen
WinBUGS <- Daten$WinBUGS
#WinBUGS
jagsUI <- Daten$jagsUI
#jagsUI
R2jags <- Daten$R2jags
#R2jags
rjags <- Daten$rjags
#rjags
runjags <- Daten$runjags
#runjags
NimbleModel <- Daten$nimbleModel
#NimbelModel
readBUGSModel <- Daten$readBUGSModel
#readBUGSModel


# Data generation
x  <- seq(-2, 2, 0.05)
y1 <- pnorm(x)
y2 <- pnorm(x,1,1)
df <- data.frame(x,y1,y2)


# Data generation
x  <- Verteilung
#x
y1 <- WinBUGS
y2 <- jagsUI
y3 <- R2jags
y4 <- rjags
y5 <- runjags
y6 <- NimbleModel
y7 <-readBUGSModel
df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7)

# make x an ordered factor
x <- factor(x, levels =x)
#x
g <- ggplot(df, aes(x, compare, group = 1))
g <- g + geom_line(aes(y=y1), colour="green",stat="identity") 
g <- g + geom_line(aes(y=y2), colour="tomato",stat="identity")
g <- g + geom_line(aes(y=y3), colour="coral",stat="identity")
g <- g + geom_line(aes(y=y4), colour="sienna1",stat="identity")
g <- g + geom_line(aes(y=y5), colour="chocolate1",stat="identity")
g <- g + geom_line(aes(y=y6), colour="blue",stat="identity")
g <- g + geom_line(aes(y=y7), colour="deepskyblue",stat="identity")
g


g + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +scale_x_discrete(limits=x)


###
g <- ggplot(df, aes(x, compare, group = 1))
g <- g + geom_line(aes(y=y1), colour="green",stat="identity") 
g <- g + geom_line(aes(y=y2), colour="red",stat="identity")
g <- g + geom_line(aes(y=y3), colour="orange",stat="identity")
g <- g + geom_line(aes(y=y4), colour="yellow",stat="identity")
g <- g + geom_line(aes(y=y5), colour="chocolate1",stat="identity")
g <- g + geom_line(aes(y=y6), colour="blue",stat="identity")
g <- g + geom_line(aes(y=y7), colour="deepskyblue",stat="identity")
g


g + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +scale_x_discrete(limits=x)

###

g


###

geom_point(shape=18)

g <- ggplot(df, aes(x, compare, group = 1))
g <- g + geom_line(aes(y=y1), colour="green",stat="identity") + geom_point( aes( colour="green", shape=4))
g <- g + geom_line(aes(y=y2), colour="tomato",stat="identity")+ geom_point( aes(colour="tomato", shape=4))
g <- g + geom_line(aes(y=y3), colour="coral",stat="identity")+ geom_point( aes(colour="coral",shape=4))
g <- g + geom_line(aes(y=y4), colour="sienna1",stat="identity")+ geom_point( aes(colour="sienna1",shape=4))
g <- g + geom_line(aes(y=y5), colour="chocolate1",stat="identity")+ geom_point( aes(colour="chocolate1",shape=4))
g <- g + geom_line(aes(y=y6), colour="blue",stat="identity")+ geom_point( aes(colour="blue",shape=4))
g <- g + geom_line(aes(y=y7), colour="deepskyblue",stat="identity")+ geom_point( aes(colour="deepskyblue",shape=4))
g


g + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +scale_x_discrete(limits=x)
###


demo("colors")




# Part 1b: Grafik Verteilungen_II ------------------------------------------



Daten_1b <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen_II")
Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

Parameter <- Daten$Verteilung 
Parameter
  
head(Daten_1b)
Sampler <- Daten_1b$Sampler
Sampler
mu_1 <- Daten_1b$`mu[1]`
#mu_1
mu_2 <- Daten_1b$`mu[2]`
mu_3 <- Daten_1b$`mu[3]`
mu_4 <- Daten_1b$`mu[4]`
mu_5 <- Daten_1b$`mu[5]`
mu_6 <- Daten_1b$`mu[6]`
mu_7 <- Daten_1b$`mu[7]`
mu_8 <- Daten_1b$`mu[8]`
mu_9 <- Daten_1b$`mu[9]`
mu_10 <- Daten_1b$`mu[10]`
mu_11 <- Daten_1b$`mu[11]`
mu_12 <- Daten_1b$`mu[12]`
mu_13 <- Daten_1b$`mu[13]`
mu_14 <- Daten_1b$`mu[14]`
mu_15 <- Daten_1b$`mu[15]`
mu_16 <- Daten_1b$`mu[16]`
mu_17 <- Daten_1b$`mu[17]`
mu_18 <- Daten_1b$`mu[18]`
mu_19 <- Daten_1b$`mu[19]`
mu_20 <- Daten_1b$`mu[20]`
mu_21 <- Daten_1b$`mu[21]`
mu_22 <- Daten_1b$`mu[22]`
d_2 <- Daten_1b$`d[2]`
sd <- Daten_1b$sd
A <- Daten_1b$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Values of Parameters")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g
## => gespeichert unter "Bild_MCMC-Vergleich _ simErgeb.jpeg"





# Part 2: effective Size _ coda    ----------------------------------------------------------------



Daten_2 <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "effectiveSize")
#Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

#Parameter <- Daten$Verteilung 
#Parameter

head(Daten_2)


Sampler <- Daten_2$Sampler
Sampler
mu_1 <- Daten_2$`mu[1]`
#mu_1
mu_2 <- Daten_2$`mu[2]`
mu_3 <- Daten_2$`mu[3]`
mu_4 <- Daten_2$`mu[4]`
mu_5 <- Daten_2$`mu[5]`
mu_6 <- Daten_2$`mu[6]`
mu_7 <- Daten_2$`mu[7]`
mu_8 <- Daten_2$`mu[8]`
mu_9 <- Daten_2$`mu[9]`
mu_10 <- Daten_2$`mu[10]`
mu_11 <- Daten_2$`mu[11]`
mu_12 <- Daten_2$`mu[12]`
mu_13 <- Daten_2$`mu[13]`
mu_14 <- Daten_2$`mu[14]`
mu_15 <- Daten_2$`mu[15]`
mu_16 <- Daten_2$`mu[16]`
mu_17 <- Daten_2$`mu[17]`
mu_18 <- Daten_2$`mu[18]`
mu_19 <- Daten_2$`mu[19]`
mu_20 <- Daten_2$`mu[20]`
mu_21 <- Daten_2$`mu[21]`
mu_22 <- Daten_2$`mu[22]`
d_2 <- Daten_2$`d[2]`
sd <- Daten_2$sd
A <- Daten_2$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Effective Sample Size")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g


## => gespeichert unter "Bild_MCMC-Vergleich _ effectiveSize_coda.jpeg"



# Part 3: Rhat mit mcmcr    ----------------------------------------------------------------



Daten_3 <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Rhat_mcmcr")
#Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

#Parameter <- Daten$Verteilung 
#Parameter

head(Daten_3)


Sampler <- Daten_3$Sampler
Sampler
mu_1 <- Daten_3$`mu[1]`
#mu_1
mu_2 <- Daten_3$`mu[2]`
mu_3 <- Daten_3$`mu[3]`
mu_4 <- Daten_3$`mu[4]`
mu_5 <- Daten_3$`mu[5]`
mu_6 <- Daten_3$`mu[6]`
mu_7 <- Daten_3$`mu[7]`
mu_8 <- Daten_3$`mu[8]`
mu_9 <- Daten_3$`mu[9]`
mu_10 <- Daten_3$`mu[10]`
mu_11 <- Daten_3$`mu[11]`
mu_12 <- Daten_3$`mu[12]`
mu_13 <- Daten_3$`mu[13]`
mu_14 <- Daten_3$`mu[14]`
mu_15 <- Daten_3$`mu[15]`
mu_16 <- Daten_3$`mu[16]`
mu_17 <- Daten_3$`mu[17]`
mu_18 <- Daten_3$`mu[18]`
mu_19 <- Daten_3$`mu[19]`
mu_20 <- Daten_3$`mu[20]`
mu_21 <- Daten_3$`mu[21]`
mu_22 <- Daten_3$`mu[22]`
d_2 <- Daten_3$`d[2]`
sd <- Daten_3$sd
A <- Daten_3$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Rhat")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g


## => gespeichert unter "Bild_MCMC-Vergleich _ Rhat-mcmcr.jpeg"













# Part 4: ESS mit mcmcr    ----------------------------------------------------------------



Daten_4 <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "ess_mcmcr")
#Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

#Parameter <- Daten$Verteilung 
#Parameter

head(Daten_4)


Sampler <- Daten_4$Sampler
Sampler
mu_1 <- Daten_4$`mu[1]`
#mu_1
mu_2 <- Daten_4$`mu[2]`
mu_3 <- Daten_4$`mu[3]`
mu_4 <- Daten_4$`mu[4]`
mu_5 <- Daten_4$`mu[5]`
mu_6 <- Daten_4$`mu[6]`
mu_7 <- Daten_4$`mu[7]`
mu_8 <- Daten_4$`mu[8]`
mu_9 <- Daten_4$`mu[9]`
mu_10 <- Daten_4$`mu[10]`
mu_11 <- Daten_4$`mu[11]`
mu_12 <- Daten_4$`mu[12]`
mu_13 <- Daten_4$`mu[13]`
mu_14 <- Daten_4$`mu[14]`
mu_15 <- Daten_4$`mu[15]`
mu_16 <- Daten_4$`mu[16]`
mu_17 <- Daten_4$`mu[17]`
mu_18 <- Daten_4$`mu[18]`
mu_19 <- Daten_4$`mu[19]`
mu_20 <- Daten_4$`mu[20]`
mu_21 <- Daten_4$`mu[21]`
mu_22 <- Daten_4$`mu[22]`
d_2 <- Daten_4$`d[2]`
sd <- Daten_4$sd
A <- Daten_4$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Effective Sample Size")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g


## => gespeichert unter "Bild_MCMC-Vergleich _ ESS-mcmcr.jpeg"







# Part 5: MCMC Efficiency    ----------------------------------------------------------------



Daten_5 <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "MCMC efficiency")
#Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

#Parameter <- Daten$Verteilung 
#Parameter

head(Daten_5)


Sampler <- Daten_5$Sampler
Sampler
mu_1 <- Daten_5$`mu[1]`
#mu_1
mu_2 <- Daten_5$`mu[2]`
mu_3 <- Daten_5$`mu[3]`
mu_4 <- Daten_5$`mu[4]`
mu_5 <- Daten_5$`mu[5]`
mu_6 <- Daten_5$`mu[6]`
mu_7 <- Daten_5$`mu[7]`
mu_8 <- Daten_5$`mu[8]`
mu_9 <- Daten_5$`mu[9]`
mu_10 <- Daten_5$`mu[10]`
mu_11 <- Daten_5$`mu[11]`
mu_12 <- Daten_5$`mu[12]`
mu_13 <- Daten_5$`mu[13]`
mu_14 <- Daten_5$`mu[14]`
mu_15 <- Daten_5$`mu[15]`
mu_16 <- Daten_5$`mu[16]`
mu_17 <- Daten_5$`mu[17]`
mu_18 <- Daten_5$`mu[18]`
mu_19 <- Daten_5$`mu[19]`
mu_20 <- Daten_5$`mu[20]`
mu_21 <- Daten_5$`mu[21]`
mu_22 <- Daten_5$`mu[22]`
d_2 <- Daten_5$`d[2]`
sd <- Daten_5$sd
A <- Daten_5$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Effective Sample Size per Second")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g


## => gespeichert unter "Bild_MCMC-Vergleich _ MCMC Efficiency.jpeg"







# Part 6: MCMC Pace    ----------------------------------------------------------------



Daten_6 <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "MCMC Pace")
#Daten <- read_excel("Ergebnisse MCMC Vergleich.xlsx", sheet = "Verteilungen", skip=1)

#Parameter <- Daten$Verteilung 
#Parameter

head(Daten_6)


Sampler <- Daten_6$Sampler
Sampler
mu_1 <- Daten_6$`mu[1]`
#mu_1
mu_2 <- Daten_6$`mu[2]`
mu_3 <- Daten_6$`mu[3]`
mu_4 <- Daten_6$`mu[4]`
mu_5 <- Daten_6$`mu[5]`
mu_6 <- Daten_6$`mu[6]`
mu_7 <- Daten_6$`mu[7]`
mu_8 <- Daten_6$`mu[8]`
mu_9 <- Daten_6$`mu[9]`
mu_10 <- Daten_6$`mu[10]`
mu_11 <- Daten_6$`mu[11]`
mu_12 <- Daten_6$`mu[12]`
mu_13 <- Daten_6$`mu[13]`
mu_14 <- Daten_6$`mu[14]`
mu_15 <- Daten_6$`mu[15]`
mu_16 <- Daten_6$`mu[16]`
mu_17 <- Daten_6$`mu[17]`
mu_18 <- Daten_6$`mu[18]`
mu_19 <- Daten_6$`mu[19]`
mu_20 <- Daten_6$`mu[20]`
mu_21 <- Daten_6$`mu[21]`
mu_22 <- Daten_6$`mu[22]`
d_2 <- Daten_6$`d[2]`
sd <- Daten_6$sd
A <- Daten_6$A



### make Sampler an ordered factor
Sampler <- factor(Sampler, levels =Sampler)


### Data generation
x  <- Sampler
#x
y1 <- mu_1
y2 <- mu_2 
y3 <-  mu_3
y4 <-  mu_4
y5 <-  mu_5
y6 <-  mu_6
y7 <- mu_7
y8 <- mu_8
y9 <- mu_9
y10 <- mu_10

y11 <- mu_11
y12 <- mu_12
y13 <- mu_13
y14 <- mu_14
y15 <- mu_15
y16 <- mu_16
y17 <- mu_17
y18 <- mu_18
y19 <- mu_19
y20 <- mu_20

y21 <- mu_21 
y22 <- mu_22 

y23 <- d_2 
y24 <- sd 
y25 <- A 

### make x  an ordered factor
x <- factor(x, levels =x)
#y <- factor(y, levels =y)


df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25)
#df <- data.frame(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25, Parameter)

#df_1b <- data.frame(Parameter)
#mod <- c(x,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10, y11,y12,y13,y14,y15,y16,y17,y18,y19,y20, y21,y22,y23,y24,y25,Parameter )
#Parameter <- Daten$Verteilung
#Parameter
#example(colour)
#brewer.pal(11,"Spectral")


cols <- c("mu[1]"="green",   "mu[2]"="red",   "mu[3]"="hotpink",   "mu[4]"="yellow",   "mu[5]"="palegreen2", 
          "mu[6]"="blue",   "mu[7]"="grey",   "mu[8]"="tan1",   "mu[9]"="deeppink",   "mu[10]"="lawngreen",
          "mu[11]"="deepskyblue2",   "mu[12]"="chartreuse4",  "mu[13]"="turquoise1",  "mu[14]"="darkgreen",   "mu[15]"="mediumseagreen",
          "mu[16]"="azure4",   "mu[17]"="brown",   "mu[18]"="darkred",   "mu[19]"="darkviolet",   "mu[20]"="darkmagenta", 
          "mu[21]"="darkorchid",   "mu[22]"="gold4",   "d[2]"="sienna2",  "sd"= "olivedrab3",   "A"= "black")
#cols <- factor(cols, levels ="alphabetic") # egal wie: fuzt nicht



g <- ggplot(df,  aes(x=Sampler,group=7)  )  # geom_tile(fill = Parameter, aes(y=Parameter)) # 1, fill=y1
g <- g + geom_line(aes(y=y1, colour="mu[1]")) +geom_point(aes(y=y1, colour="mu[1]"))  #+geom_label(aes(y=y1, label="mu[1]"), colour="green") # , show.legend = T
g <- g + geom_line(aes(y=y2, colour="mu[2]"))  +geom_point(aes(y=y2, colour="mu[2]"))
g <- g + geom_line(aes(y=y3, colour="mu[3]")) +geom_point(aes(y=y3, colour="mu[3]"))
g <- g + geom_line(aes(y=y4, colour="mu[4]")) +geom_point(aes(y=y4, colour="mu[4]"))
g <- g + geom_line(aes(y=y5, colour="mu[5]")) +geom_point(aes(y=y5, colour="mu[5]"))
g <- g + geom_line(aes(y=y6, colour="mu[6]")) +geom_point(aes(y=y6, colour="mu[6]"))
g <- g + geom_line(aes(y=y7, colour="mu[7]")) +geom_point(aes(y=y7, colour="mu[7]"))
g <- g + geom_line(aes(y=y8, colour="mu[8]")) +geom_point(aes(y=y8, colour="mu[8]"))
g <- g + geom_line(aes(y=y9, colour="mu[9]")) +geom_point(aes(y=y9, colour="mu[9]"))
g <- g + geom_line(aes(y=y10, colour="mu[10]")) +geom_point(aes(y=y10, colour="mu[10]"))

g <- g + geom_line(aes(y=y11, colour="mu[11]")) +geom_point(aes(y=y11, colour="mu[11]"))
g <- g + geom_line(aes(y=y12, colour="mu[12]"))  +geom_point(aes(y=y12, colour="mu[12]"))
g <- g + geom_line(aes(y=y13, colour="mu[13]"))  +geom_point(aes(y=y13, colour="mu[13]"))
g <- g + geom_line(aes(y=y14, colour="mu[14]"))  +geom_point(aes(y=y14, colour="mu[14]"))
g <- g + geom_line(aes(y=y15, colour="mu[15]"))  +geom_point(aes(y=y15, colour="mu[15]"))
g <- g + geom_line(aes(y=y16, colour="mu[16]"))  +geom_point(aes(y=y16, colour="mu[16]"))
g <- g + geom_line(aes(y=y17, colour="mu[17]"))  +geom_point(aes(y=y17, colour="mu[17]"))
g <- g + geom_line(aes(y=y18, colour="mu[18]"))  +geom_point(aes(y=y18, colour="mu[18]"))
g <- g + geom_line(aes(y=y19, colour="mu[19]"))  +geom_point(aes(y=y19, colour="mu[19]"))
g <- g + geom_line(aes(y=y20, colour="mu[20]")) +geom_point(aes(y=y20, colour="mu[20]"))

g <- g + geom_line(aes(y=y21, colour="mu[21]")) +geom_point(aes(y=y21, colour="mu[21]"))
g <- g + geom_line(aes(y=y22, colour="mu[22]")) +geom_point(aes(y=y22, colour="mu[22]"))
g <- g + geom_line(aes(y=y23, colour="d[2]")) +geom_point(aes(y=y23, colour="d[2]"))
g <- g + geom_line(aes(y=y24, colour="sd")) +geom_point(aes(y=y24, colour="sd"))
g <- g + geom_line(aes(y=y25, colour="A")) +geom_point(aes(y=y25, colour="A"))


g <- g + scale_colour_manual(name="Parameter",values=cols) +  ylab("Seconds per Effective Size")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_x_discrete(limits=x)

g


## => gespeichert unter "Bild_MCMC-Vergleich _ MCMC Pace.jpeg"



