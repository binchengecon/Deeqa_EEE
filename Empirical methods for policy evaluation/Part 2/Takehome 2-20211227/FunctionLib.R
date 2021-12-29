
## Define likelihood function & functions to concentrate out some of the parameters

G_tilde <- function(mu, sigma, dataset){
  return( 1 - plnorm(min(dataset$wage), meanlog = mu, sdlog = sigma, lower.tail = TRUE, log.p = FALSE) )
}

g       <- function(mu,sigma,dataset,employ){
  return(log(dlnorm(x = dataset$wage[dataset$status == employ], meanlog = mu, sdlog = sigma)))
}

lambda <- function(mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  return( ( Nu[count]/(sum(dataset$duration)) )*( 1/(G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset) + G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)) ) )
}

eta <- function(lambda,mu_i,sigma_i,mu_f,sigma_f,dataset,count,Nu,Nei,Nef){
  return( 1/lambda * ( ((Nu[count]+Nei[count]+Nef[count]))/Nei[count])*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset) - G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) - G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset) ) 
}

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




LogLikelihood2 <- function(parameters,dataset,count,Nu,Nei,Nef){
  
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
  N <- Nu[count]+Nei[count]+Nef[count]
  
  # Specify LL
  LL <- (
    -N*log(
      1 + 
        lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset)*(1/eta)+ 
        lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)*(1/eta)
    )+
      #
      Nef[count]*log(lambda/eta) + 
      Nei[count]*log(lambda/eta) + 
      #
      Nu[count]*log(lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset)) -      
      (lambda*G_tilde(mu = mu_f, sigma = sigma_f, dataset = dataset) + lambda*G_tilde(mu = mu_i, sigma = sigma_i, dataset = dataset))*sum(dataset$duration)+
      sum(g(mu_f,sigma_f,dataset,2)) + 
      sum(g(mu_i,sigma_i,dataset,3)) 
  )
  
  if(sigma_i > 0 & sigma_f > 0 & lambda > 0 & eta > 0){
    return(ifelse(LL == -Inf, 10^100, -LL))
  } else{
    return(10^100)
  }
  
}









