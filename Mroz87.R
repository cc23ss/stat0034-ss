rm(list=ls())

# Required R packages
library(ssmrob)
library(ssmodels)
library(sampleSelection)

# Required R functions
source("MLETobit2.R",local=FALSE)

#Data preparation: Female labour supply data set
data(Mroz87)
attach(Mroz87)
head(Mroz87)

Mroz87$kids<-(Mroz87$kids5 + Mroz87$kids618 > 0)
Mroz87$age2<-Mroz87$age^2
Mroz87$exper2<-Mroz87$exper^2
head(Mroz87)
attach(Mroz87)

########################################################################
#                          MLE (direct)                                #
########################################################################
# outcome
deso<-cbind(exper, exper2, educ, city)
# selection
dess<-cbind(age, age2, faminc, kids, educ)
# selection indicator
sel<-as.logical(lfp)
# outcome 
out<-wage
# Optimisation step
#init0<-rep(0,13)
init0<-c(3.10350,-0.13328,
         -1.9537242,0.0284295,-0.0001151,0.4562471,0.4451424,
         -4.120,0.1840,-0.002409,0.000005676,-0.4507,0.09533)
OPTlocal<-MLETobit2(init = round(init0,2), out = out, sel = sel, deso = deso, dess = dess, method = "nlminb", maxit = 10000)
OPT<-MLETobit2(init = round(init0), out = out, sel = sel, deso = deso, dess = dess, method = "nlminb", maxit = 10000)

#MLE
MLElocal<-OPTlocal$MLE
MLE<-OPT$MLE
MLElocal
MLE

########################################################################
#                         MLE (R packages)                             #
########################################################################
# Selection equation
selectEq<-lfp ~ age + I(age^2) + faminc + kids + educ
# Outcome equation
outcomeEq<-wage ~ exper + I(exper^2) + educ + city

# ssmodels R package
#HMLE<-HeckmanCL(selectEq, outcomeEq, data = Mroz87)
#summary(HMLE)

# sampleSelection
#ssMLElocal<-selection(selectEq, outcomeEq, data = Mroz87, maxMethod = "BHHH", iterlim = 500)
#summary(ssMLElocal) 

#ssMLE0<-selection(selectEq, outcomeEq, data = Mroz87,
#                  start = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.9),
#                  maxMethod = "SANN", parscale = 0.001)
#summary(ssMLE0)
#ssMLE<-selection(selectEq, outcomeEq, data = Mroz87,
#                 start = coef(ssMLE0))
#summary(ssMLE)

#coef(ssMLElocal)
#coef(ssMLE)
#MLE






