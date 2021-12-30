###########Function Table of Content###############

##############Part 1: CDF and PDF########################
#####G_tilde: give lognormal complementary probability of wage
#####g: give lognormal density of wage vector

############### part 2: FOC checking ################
#####lambda: given input of estimated parameters (mu and sigma), 
# compute lambda according to equation 15 in solution document.
#####eta: given input of estimated parameters (mu and sigma), 
#         compute eta according to equation 16 in solution document.
#####ratioGandN: given input of estimated parameters (mu and sigma), 
#                compute ratio of probability and population size according to equation 17 in solution document.

############# part 3: Loglikelihood#################
#####LogLikelihood: given input of estimated parameters (mu and sigma), compute -loglikelihood of sample 
#########           with accelerated algo that use function lambda and eat to smooth search of lambda and eat.
#####LogLikelihood2: given input of estimated parameters (mu and sigma), compute -loglikelihood of sample 
#########           without accelerated algo and allow for std estimation in mle optimization.

########  Function Defintion#########################

#####G_tilde: give lognormal complementary probability of wage
#### Output: 1*1 scalar

### Input : 
### Mu: mean
### Sigma: var
### dataset: the country dataset

G_tilde <- function(mu, sigma, dataset){
  temp <- min(dataset$wage[dataset$status != 1])
  return( 1 - plnorm(temp, meanlog = mu, sdlog = sigma, lower.tail = TRUE, log.p = FALSE) )
}

#####g: give lognormal density of wage vector

#### Output: 1*N vector

### Input: 
### Mu: mean
### Sigma: var
### dataset: the country dataset

g       <- function(mu,sigma,dataset,employ){
  return(log(dlnorm(x = dataset$wage[dataset$status == employ], meanlog = mu, sdlog = sigma)))
}

#####lambda: given input of estimated parameters (mu and sigma), compute lambda according to equation 15 in solution document.

#### Output: 1*1 scalar

### Input: 
### Mu: mean of each job type
### Sigma: var of each job type
### dataset: the country dataset
### count: country number
### Nu: vector of number unemployed
### Nei: vector of number informal
### Nef: vector of number formal

lambda <- function(mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  return( ( Nu[count]/(sum(dataset$duration)) )*( 1/(G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset) + G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)) ) )
}

#####eta: given input of estimated parameters (mu and sigma), compute eta according to equation 16 in solution document.

#### Output: 1*1 scalar

### Input: 
### Lambda: lambda
### Mu: mean of each job type
### Sigma: var of each job type
### dataset: the country dataset
### count: country number
### Nu: vector of number unemployed
### Nei: vector of number informal
### Nef: vector of number formal

eta <- function(lambda,mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  N_temp <-Nu[count]+Nei[count]+Nef[count]
  Part1 <- lambda
  Part2 <- N_temp/Nei[count] *G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)
  Part3 <-  G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)
  return(  Part1*(Part2-Part3)) 
}

eta2 <- function(lambda,mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  N_temp <-Nu[count]+Nei[count]+Nef[count]
  Part1 <- lambda
  Part2 <- N_temp/Nef[count] *G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)
  Part3 <-  G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)
  return(  Part1*(Part2-Part3)) 
}
#####ratioGandN: given input of estimated parameters (mu and sigma), compute ratio of probability and population size according to equation 17 in solution document.

#### Output: 1*1 scalar

### Input: 
### Mu: mean of each job type
### Sigma: var of each job type
### dataset: the country dataset
### count: country number
### Nu: vector of number unemployed
### Nei: vector of number informal
### Nef: vector of number formal

ratioGandN <- function(mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  Temp1<-G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)/Nei[count]
  Temp2<-G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)/Nef[count]
  return(Temp1/Temp2)
}

#####LogLikelihood: given input of estimated parameters (mu and sigma), compute loglikelihood of sample 
#########           with accelerated algo that use function lambda and eat to smooth search of lambda and eat.

#### Output: 1*1 scalar

### Input: 
### Mu: mean of each job type
### Sigma: var of each job type
### dataset: the country dataset
### count: country number
### Nu: vector of number unemployed
### Nei: vector of number informal
### Nef: vector of number formal

LogLikelihood <- function(parameters,dataset,count,Nu,Nei,Nef){
  
  print(parameters)
  
  
  # Assign parameters
  mu_i <- parameters[1]
  sigma_i <- parameters[2]
  
  mu_f <- parameters[3]
  sigma_f <- parameters[4]
  
  # lambda <- parameters[5]
  # eta <- parameters[6]
  # 
  # get values for lambda and eta 
  lambda_current <- lambda(mu_i = mu_i, sigma_i = sigma_i, mu_f = mu_f, sigma_f  = sigma_f, dataset = dataset, count = count ,Nu = Nu,Nei=Nei,Nef=Nef)
  eta_current <- eta(lambda = lambda_current, mu_i = mu_i, sigma_i = sigma_i, mu_f = mu_f, sigma_f = sigma_f,dataset=dataset,count = count ,Nu = Nu,Nei=Nei,Nef=Nef)
  N <- Nu[count]+Nei[count]+Nef[count]
  
  # Specify LL
  LL <- (
    -N*log(
      1 + 
        lambda_current*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)*(1/eta_current)+ 
        lambda_current*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)*(1/eta_current)
    )+
      #
      Nef[count]*log(lambda_current/eta_current) + 
      Nei[count]*log(lambda_current/eta_current) + 
      #
      Nu[count]*log(lambda_current*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda_current*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)) -      
      (lambda_current*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda_current*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset))*sum(dataset$duration)+
      sum(g(mu_f,sigma_f,dataset,2)) + 
      sum(g(mu_i,sigma_i,dataset,3)) 
  )
  
  if(sigma_i > 0 & sigma_f > 0 & lambda_current > 0 & eta_current > 0){
    return(ifelse(LL == -Inf, 10^100, -LL))
  } else{
    return(10^100)
  }
  
}


#####LogLikelihood2: given input of estimated parameters (mu and sigma), compute loglikelihood of sample 
#########           without accelerated algo and allow for std estimation in mle optimization.

#### Output: 1*1 scalar

### Input: 
### Mu: mean of each job type
### Sigma: var of each job type
### dataset: the country dataset
### count: country number
### Nu: vector of number unemployed
### Nei: vector of number informal
### Nef: vector of number formal

LogLikelihood2 <- function(parameters,dataset,count,Nu,Nei,Nef){
  
  # print(parameters)
  
  
  # Assign parameters
  mu_i <- parameters[1]
  sigma_i <- parameters[2]
  
  mu_f <- parameters[3]
  sigma_f <- parameters[4]
  
  lambda <- parameters[5]
  eta    <- parameters[6]
  # 
  # get values for lambda and eta 
  N_temp <- Nu[count]+Nei[count]+Nef[count]
  
  # Specify LL
  LL <- (
    -N_temp*log(
      1
      +    lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)*(1/eta)
      +   lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)*(1/eta)
    )
    #
    + Nef[count]*log(lambda/eta)  
    +  Nei[count]*log(lambda/eta)  
    #
    +  Nu[count]*log(lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset))       
    #
    -  (lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset))*sum(dataset$duration)
    +  sum(g(mu_f,sigma_f,dataset,2))  
    +  sum(g(mu_i,sigma_i,dataset,3)) 
  )
  
  if(sigma_i > 0 & sigma_f > 0 & lambda > 0 & eta > 0){
    return(ifelse(LL == -Inf, 10^100, -LL))
  } else{
    return(10^100)
  }
  
}



FOCdiff <- function(LLK_maxLik){
    temp_index <- which.min(LLK_maxLik[,ncol(LLK_maxLik)])
    temp_param <- LLK_maxLik[temp_index,1:num_intervaldim]
    temp_std   <- LLK_maxLik[temp_index,(1+num_intervaldim):(2*num_intervaldim)]
    temp_prb   <- LLK_maxLik[temp_index,(1+3*num_intervaldim):(4*num_intervaldim)]
    FOC_Lambda[count] <- lambda(mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    FOC_Eta[count]  <- eta(lambda = temp_param[5],mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    FOC_Eta2[count]   <- eta2(lambda = temp_param[5],mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    FOC_Ratio[count] <- ratioGandN(mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    
}





