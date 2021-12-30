##### Max likelihood estimation EMPA course Matteo Bobba #####

## Load packages
rm(list=ls())

library(tidyverse)
library(stats4)
library(maxLik)
library(alr4)

setwd("C:/Users/alienware/Desktop/Deeqa_EEE/Empirical methods for policy evaluation/Part 2/Takehome 2-20211227")

## Load data
# chile <- read.delim(file = "casen_chile.txt", sep = ",", header = FALSE, col.names = c("unemp_duration","wage","status"))
chile <- read.delim(file = "chile.txt", sep = ",", header = FALSE, col.names = c("unemp_duration","wage","status"))

## define parameters
p <- 0.05
N <- nrow(chile)
N_u <- nrow(chile[chile$status == 1,])
N_f <- nrow(chile[chile$status == 2,])
N_i <- nrow(chile[chile$status == 3,])

# Take the estimate of the reservation wage
w_star <- min(chile$wage[chile$status != 1])

## Define likelihood function & functions to concentrate out some of the parameters

G_tilde <- function(mu, sigma, wage){
  return( 1 - plnorm(wage, meanlog = mu, sdlog = sigma, lower.tail = TRUE, log.p = FALSE) )
}

lambda <- function(mu_i,sigma_i,mu_f,sigma_f){
  return( ( N_u/(sum(chile$unemp_duration)) )*( 1/(G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star) + G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star)) ) )
}

eta <- function(lambda,mu_i,sigma_i,mu_f,sigma_f){
  return( lambda*( (N/N_i)*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star) - G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) - G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star) ) )
}

eta2 <- function(lambda,mu_i,sigma_i,mu_f,sigma_f){
  return( lambda*( (N/N_f)*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) - G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) - G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star) ) )
}

LogLikelihood <- function(parameters){
  
  print(parameters)
  
  # Assign parameters
  mu_i <- parameters[1]
  sigma_i <- parameters[2]
  
  mu_f <- parameters[3]
  sigma_f <- parameters[4]
  
  # get values for lambda and eta (concentrated out)
  lambda_value <- lambda(mu_i = mu_i, sigma_i = sigma_i, mu_f = mu_f, sigma_f  = sigma_f)
  eta_value <- eta(lambda = lambda_value, mu_i = mu_i, sigma_i = sigma_i, mu_f = mu_f, sigma_f = sigma_f)
  
  # Specify LL
  LL <- (
    -N*log(1 + lambda_value*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star)*(1/eta_value) + 
             lambda_value*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star)*(1/eta_value)
    ) + 
      # contribution of formally employed 
      sum(log(dlnorm(x = chile$wage[chile$status == 2], meanlog = mu_f, sdlog = sigma_f))) + 
      # contribution of informally employed
      sum(log(dlnorm(x = chile$wage[chile$status == 3], meanlog = mu_i, sdlog = sigma_i))) + 
      # 
      N_f*log(lambda_value/eta_value) + 
      N_i*log(lambda_value/eta_value) + 
      N_u*log(lambda_value*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) + lambda_value*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star)) -
      # 
      (lambda_value*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) + lambda_value*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star))*sum(chile$unemp_duration)
  )
  
  if(sigma_i > 0 & sigma_f > 0 & lambda_value > 0 & eta_value > 0){
    return(ifelse(LL == -Inf, 10^100, -LL))
  } else{
    return(10^100)
  }
  
}



LogLikelihood2 <- function(parameters){
  
  print(parameters)
  
  
  # Assign parameters
  mu_i <- parameters[1]
  sigma_i <- parameters[2]
  
  mu_f <- parameters[3]
  sigma_f <- parameters[4]
  
  lambda <- parameters[5]
  eta    <- parameters[6]
  # 
  # get values for lambda and eta 
  N <- N_u+N_i+N_f
  
  # Specify LL
  LL <- (
    -N*log(1 + lambda*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star)*(1/eta) + 
             lambda*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star)*(1/eta)
    )   
    # 
    +  N_f*log(lambda/eta)  
    +  N_i*log(lambda/eta)  
    +  N_u*log(lambda*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star)) 
    # 
    - (lambda*G_tilde(mu = mu_f, sigma = sigma_f, wage = w_star) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, wage = w_star))*sum(chile$unemp_duration)
    # formally employed 
    + sum(log(dlnorm(x = chile$wage[chile$status == 2], meanlog = mu_f, sdlog = sigma_f))) 
    # informally employed
    + sum(log(dlnorm(x = chile$wage[chile$status == 3], meanlog = mu_i, sdlog = sigma_i)))
  )
  
  if(sigma_i > 0 & sigma_f > 0 & lambda > 0 & eta > 0){
    return(ifelse(LL == -Inf, 10^100, -LL))
  } else{
    return(10^100)
  }
  
}
### Optimize over LL

initial_parameters2 <- rep(1.3,6)
optimal_parameters2 <- optim(par = initial_parameters2,
                             fn = LogLikelihood2,method = "BFGS", 
                             control = list(maxit = 10000, reltol = 1e-12))
 

Param <- optimal_parameters2$par
mu_i <- Param[1]
sigma_i<-Param[2]
mu_f<-Param[3]
sigma_f<-Param[4]

LogLikelihood2(parameters = Param)
G_tilde(mu_i,sigma_i,w_star)

# different algorithm converges to the same values
optimal_parameters3 <- optim(par = optimal_parameters$par,fn = LogLikelihood,method = "L-BFGS-B", 
                             lower = c(-10,0.75,0.0,0.1), upper = c(10,10,10,10), control = list(maxit = 10000, pgtol = 1e-12)
)

optimal_parameters4 <- optim(par = optimal_parameters2$par,fn = LogLikelihood2,method = "L-BFGS-B", 
                             lower = c(-10,0.75,0.0,0.1,0,0), upper = c(10,10,10,10,1,1), control = list(maxit = 10000, pgtol = 1e-12)
)

# FOC that is lying around is fulfilled!
(G_tilde(mu = optimal_parameters$par[1], sigma = optimal_parameters$par[2], wage = w_star)/N_i)/(G_tilde(mu = optimal_parameters$par[3], sigma = optimal_parameters$par[4], wage = w_star)/N_f)
(G_tilde(mu = optimal_parameters2$par[1], sigma = optimal_parameters2$par[2], wage = w_star)/N_i)/(G_tilde(mu = optimal_parameters2$par[3], sigma = optimal_parameters2$par[4], wage = w_star)/N_f)
(G_tilde(mu = optimal_parameters4$par[1], sigma = optimal_parameters4$par[2], wage = w_star)/N_i)/(G_tilde(mu = optimal_parameters4$par[3], sigma = optimal_parameters4$par[4], wage = w_star)/N_f)




# Other initial values?
initial_parameters_alt <- c(1,0.4,-4.9,3.1)
optimal_parameters_alt <- optim(par = initial_parameters_alt,
                                fn = LogLikelihood,method = "BFGS", 
                                control = list(maxit = 10000, reltol = 1e-12))



## Estimate MaxLikelihood with package

mle_results <- maxLik(logLik = function(parameters){-LogLikelihood(parameters)}, start = initial_parameters)
mle_results_alt <- maxLik(logLik = function(parameters){-LogLikelihood(parameters)}, start = optimal_parameters$par)

mle_results2 <- maxLik(logLik = function(parameters){-LogLikelihood2(parameters)}, start = initial_parameters2)
mle_results_alt2 <- maxLik(logLik = function(parameters){-LogLikelihood2(parameters)}, start = optimal_parameters2$par)

# Formal sector G-distribution is well-identified, but algorithm has problems identifying informal wage distribution, because so much of it is unobserved!
# High variance, low mean is similar on the right tail then higher mean, lower variance. 

# But this does not satisfy the FOC!! So choose the other one!!
(G_tilde(mu = mle_results$estimate[1], sigma = mle_results$estimate[2], wage = w_star)/N_i)/(G_tilde(mu = mle_results$estimate[3], sigma = mle_results$estimate[4], wage = w_star)/N_f)
(G_tilde(mu = mle_results2$estimate[1], sigma = mle_results2$estimate[2], wage = w_star)/N_i)/(G_tilde(mu = mle_results2$estimate[3], sigma = mle_results2$estimate[4], wage = w_star)/N_f)
(G_tilde(mu = mle_results_alt2$estimate[1], sigma = mle_results_alt2$estimate[2], wage = w_star)/N_i)/(G_tilde(mu = mle_results_alt2$estimate[3], sigma = mle_results_alt2$estimate[4], wage = w_star)/N_f)

# So use mle_results_alt
summary(mle_results_alt2)

lambda_estimate <- lambda(mu_i = mle_results_alt2$estimate[1], sigma_i = mle_results_alt2$estimate[2], mu_f = mle_results_alt2$estimate[3], sigma_f = mle_results_alt2$estimate[4])
eta_estimate <- eta(lambda = lambda_estimate, mu_i = mle_results_alt$estimate[1], sigma_i = mle_results_alt$estimate[2], mu_f = mle_results_alt$estimate[3], sigma_f = mle_results_alt$estimate[4])

eta_estimate2 <- eta2(lambda = lambda_estimate, mu_i = mle_results_alt$estimate[1], sigma_i = mle_results_alt$estimate[2], mu_f = mle_results_alt$estimate[3], sigma_f = mle_results_alt$estimate[4])

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





