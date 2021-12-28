rm(list=ls())
library(tidyverse)
library(bbmle)
library(pracma)
library(fitdistrplus)

setwd('C:/Users/wen_x/Documents/main_documents/All Toulouse Files/TSE_Courses/TSE Deeqa/Policy Evaluation/Policy Evaluation Project 2')

data <- read.csv('casen_chile.txt', header = FALSE)

names(data)<- c('tu', 'wage', 'employ')

reserve_wage <- min(filter(data, employ!=1)$wage)
N <- dim(data)[1]
Nu <- dim(filter(data, employ==1))[1]
Ni <- dim(filter(data, employ==2))[1]
Nf <- dim(filter(data, employ==3))[1]
  
sum_tu <- sum(data$tu)

rho <- 0.05

data <- as.matrix(data)


## Inital parameter guesses

lambda <- 0.2
eta <- 0.1
mu_i <- 0.5
mu_f <- 0.6
sigma_i <- 1.0
sigma_f <- 1.2

initialparam <- c(lambda, eta, mu_i, mu_f, sigma_i, sigma_f)



negativeloglike <- function(lambda, eta, mu_i, mu_f, sigma_i, sigma_f,data){
  param <- c(lambda, eta, mu_i, mu_f, sigma_i, sigma_f)
  param <- exp(param)
  Gi <- gi(param)
  Gf <- gf(param)
  denom <- -N*log((param[1]*param[2])*(Gi+Gf)+param[2]^2)
  num_i_f <- (Ni+Nf)*log(param[1]*param[2])
  num_u <- Nu*log(param[2]^2*param[1]*(Gi+Gf))
  wage_mle <- sum(apply(data,FUN=singlelike, MARGIN = 1,param = param))
  time_mle <- -(param[1])*(Gf+Gi)*sum_tu
  return(-(denom+num_i_f+num_u+wage_mle+time_mle))
}
singlelike <- function(obs, param){
  if (obs[3] == 2){
    value <- log(dlnorm(obs[2], meanlog = param[3], sdlog = param[5]))
  }else if(obs[3]==3){
    value <- log(dlnorm(obs[2], meanlog = param[4], sdlog = param[6]))
  }else{
    value <- 0
  }
  return(value)
}
gi <- function(param){
  value <- 1-plnorm(reserve_wage, meanlog = param[3], sdlog = param[5])
  return(value)
}
gf <- function(param){
  value <- 1-plnorm(reserve_wage, meanlog = param[4], sdlog = param[6])
  return(value)
}

fitdist(data[data[,3]==3,2], "lnorm")
fitdist(data[data[,3]==2,2], "lnorm")
## Non-parametric fit to data
## Informal: mu_i = 1.02, sigma_i = 0.407
## Formal: mu_f = 0.81, sigma_f = 0.47
## These gives us inital guesses for distribution parameters


estimator <- mle2(negativeloglike, data = list(data = data), start = list(lambda=-0.82, eta=-0.62, mu_i=0.02, mu_f=-0.2, sigma_i=-0.9, sigma_f=-0.75))
summary(estimator)
try <- exp(coef(estimator))
se <- sqrt((exp(diag(vcov(estimator)))-1)*exp(2*coef(estimator)+diag(vcov(estimator))))
set.seed(1234)
gi_draws <- rlnorm(100000, meanlog = try[3], sdlog = try[5])
gf_draws <- rlnorm(100000, meanlog = try[4], sdlog = try[6])


int_i <- sum(gi_draws[gi_draws>reserve_wage]-reserve_wage)/length(gi_draws[gi_draws>reserve_wage])
int_f <- sum(gf_draws[gf_draws>reserve_wage]-reserve_wage)/length(gf_draws[gf_draws>reserve_wage])

b <- reserve_wage- try[1]/(rho+try[2])*(int_i+int_f)


guesses <- exp(-seq(0.1,2,0.25))
# 16 guesses for rates

guesslist <- vector('list', 16)
index <- seq(1,4,1)
k <- 1
for (i in index){
  for (j in index){
    guesslist[[k]] <- mle2(negativeloglike, data = list(data = data), start = list(lambda=log(guesses[i]), eta=log(guesses[j]), mu_i=0.02, mu_f=-0.2, sigma_i=-0.9, sigma_f=-0.75))
    k<-k+1
  }
}

