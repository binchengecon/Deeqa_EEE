#############Environment Setting################
################################################
# Load packages
rm(list=ls())
library(tidyverse)
library(stats4)
library(maxLik)
library(alr4)
library(bbmle)

setwd("C:/Users/33678/Desktop/Deeqa_EEE/Empirical methods for policy evaluation/Part 2/Takehome 2-20211227")
source("FunctionLib.R")
################################################
################################################
################################################

#Import Data

chile <- read.delim(file = "chile.txt", sep = ",", header = FALSE, col.names = c("duration","wage","status"))
argentina <- read.delim(file = "argentina.txt", sep = ",", header = FALSE, col.names = c("duration","wage","status"))
colombia  <-read.delim(file = "colombia.txt", sep = ",", header = FALSE, col.names = c("duration","wage","status"))
mexico    <-read.delim(file = "mexico.txt", sep = ",", header = FALSE, col.names = c("duration","wage","status"))

#Combine Data

chile$country =1
argentina$country =2
colombia$country =3
mexico$country =4

data <- bind_rows(chile,argentina,colombia,mexico)

#Setting Parameters

Rho <- 0.05

N <- rep(0,max(data$country))
Nu <- rep(0,max(data$country))
Nef <- rep(0,max(data$country))
Nei <- rep(0,max(data$country))
w_star <- rep(0,max(data$country))

# for (count in 1:4){

count <- 1

dataset <- data %>% filter(country == count)

#Market Stock
N[count] <- nrow(dataset)
Nu[count] <- nrow(dataset[dataset$status == 1,])
Nef[count] <- nrow(dataset[dataset$status == 2,])
Nei[count] <- nrow(dataset[dataset$status == 3,])






################################################

# Wage Estimate

w_star[count] = min(dataset$wage[dataset$status != 1])

# Parameter Estimate
param_init <- rep(0.5,6)
param_init2 <- c(-3.6969302,2.7128403,  0.9146251,  0.3604504,0.34,0.05)
param_optimal <- optim(
par = param_init2,
fn = LogLikelihood2,method = "BFGS", 
control = list(maxit = 10000, reltol = 1e-12),dataset=dataset,count = count ,Nu = Nu,Nei=Nei,Nef=Nef)

param_optimal$par
mle_results_alt <- maxLik(logLik = function(parameters){-LogLikelihood2(parameters)}, start = optimal_parameters$par)



## FOC that is lying around is fulfilled!
(G_tilde(mu = param_optimal$par[1], sigma = param_optimal$par[2], dataset = dataset)/Nei[1])/(G_tilde(mu = param_optimal$par[3], sigma = param_optimal$par[4], dataset = dataset)/Nef[1])


# Other initial values?
initial_parameters_alt <- c(1,0.4,-4.9,3.1,1,1)
optimal_parameters_alt <- optim(par = initial_parameters_alt,
                                fn = LogLikelihood2,method = "BFGS", 
                                control = list(maxit = 10000, reltol = 1e-12),dataset=dataset,count = count ,Nu = Nu,Nei=Nei,Nef=Nef)

## Estimate MaxLikelihood with package

mle_results <- maxLik(logLik = function(parameters){-LogLikelihood2(parameters)}, start = initial_parameters_alt)


mle_results_alt <- maxLik(logLik = function(parameters){-LogLikelihood(parameters)}, start = optimal_parameters$par)

# Formal sector G-distribution is well-identified, but algorithm has problems identifying informal wage distribution, because so much of it is unobserved!
# High variance, low mean is similar on the right tail then higher mean, lower variance. 

# But this does not satisfy the FOC!! So choose the other one!!
(G_tilde(mu = mle_results$estimate[1], sigma = mle_results$estimate[2], wage = w_star)/N_i)/(G_tilde(mu = mle_results$estimate[3], sigma = mle_results$estimate[4], wage = w_star)/N_f)

# So use mle_results_alt
summary(mle_results_alt)


lambda_estimate <- lambda(mu_i = mle_results_alt$estimate[1], sigma_i = mle_results_alt$estimate[2], mu_f = mle_results_alt$estimate[3], sigma_f = mle_results_alt$estimate[4])
eta_estimate <- eta(lambda = lambda_estimate, mu_i = mle_results_alt$estimate[1], sigma_i = mle_results_alt$estimate[2], mu_f = mle_results_alt$estimate[3], sigma_f = mle_results_alt$estimate[4])

# # Attempt to get delta method to work // but didn't work...
# deltaMethod(c(b0 = coef(mle_results_alt)[1],b1 = coef(mle_results_alt)[2],b2 = coef(mle_results_alt)[3], b3 = coef(mle_results_alt)[4]),
#             g = "lambda(mu_i = b0, sigma_i = b1, mu_f = b2, sigma_f = b3)", vcov. = vcov(mle_results_alt))


## Get the estimate of b

sample_G_i <- rlnorm(n = 100000, meanlog = mle_results_alt$estimate[1], sdlog = mle_results_alt$estimate[2])
conditional_mean_G_i <- mean(sample_G_i[sample_G_i > w_star])

sample_G_f <- rlnorm(n = 100000, meanlog = mle_results_alt$estimate[3], sdlog = mle_results_alt$estimate[4])
conditional_mean_G_f <- mean(sample_G_f[sample_G_f > w_star])

# Now compute
b <- w_star - (lambda_estimate/(0.05 + eta_estimate))*( 
  conditional_mean_G_i - w_star*G_tilde(mu = optimal_parameters$par[1], sigma = optimal_parameters$par[2], wage = w_star) + 
    conditional_mean_G_f - w_star*G_tilde(mu = optimal_parameters$par[3], sigma = optimal_parameters$par[4], wage = w_star))

# }