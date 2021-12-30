#############Environment Setting################
################################################
# Load packages
rm(list=ls())
library(tidyverse)
library(stats4)
library(maxLik)
library(alr4)
library(bbmle)
library(xtable)

setwd("C:/Users/\alienware/Desktop/Deeqa_EEE/Empirical methods for policy evaluation/Part 2/Takehome 2-20211227")
# setwd("C:/Users/33678/Desktop/Deeqa_EEE/Empirical methods for policy evaluation/Part 2/Takehome 2-20211227")

source("FunctionLib2.R")
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

rho <- 0.05

N <- rep(0,max(data$country))
Nu <- rep(0,max(data$country))
Nef <- rep(0,max(data$country))
Nei <- rep(0,max(data$country))
w_star <- rep(0,max(data$country))
FOC_Lambda <- rep(0,max(data$country))
FOC_Eta <- rep(0,max(data$country))
FOC_Eta2 <- rep(0,max(data$country))
FOC_Ratio <- rep(0,max(data$country))

b <-rep(0,max(data$country))

name_country <- c("chile","argentina","colombia","mexico")

Table_Coefficient <- matrix(0,14,4)
colnames(Table_Coefficient) <- name_country
rownames(Table_Coefficient) <- c("mu_i"," ","sigma_i"," ","mu_f"," ","sigma_f"," ","lambda"," ","eta"," ","b","w^*")
Table_Coefficient <- as.table(Table_Coefficient)

Table_Population <- matrix(0,7,4)
colnames(Table_Population) <- name_country
rownames(Table_Population) <- c("N_u","Percentage","N_{e_i}","Percentage","N_{e_f}","Percentage","N")
Table_Population <- as.table(Table_Population)

for (count in 1:4){
  
  # count <- 1
  
  dataset <- data %>% filter(country == count)
  num_interval <- 100
  Lb <- c(-5,0,-5,0,0,0) 
  Ub <- c(10,10,10,10,1,1)
  
  num_intervaldim <- length(Lb)
  matrix_interval <- matrix(0,num_interval,num_intervaldim)
  
  for (i in 1:num_intervaldim){
    matrix_interval[1:num_interval,i] <- runif(num_interval,Lb[i],Ub[i])
  }
  
  #Market Stock
  N[count] <- nrow(dataset)
  Nu[count] <- nrow(dataset[dataset$status == 1,])
  Nef[count] <- nrow(dataset[dataset$status == 2,])
  Nei[count] <- nrow(dataset[dataset$status == 3,])
  
  
  
  
  
  
  
  ################################################
  
  # Wage Estimate
  
  w_star[count] = min(dataset$wage[dataset$status != 1])
  
  # Parameter Estimate
  assign(paste0("LLK_BFGS_",name_country[count]),matrix(0,num_interval,num_intervaldim+1))
  assign(paste0("LLK_LBFGSB_",name_country[count]),matrix(0,num_interval,num_intervaldim+1))
  assign(paste0("LLK_maxLik_",name_country[count]),matrix(0,num_interval,num_intervaldim*4+1))
  
  LLK_BFGS <- matrix(0,num_interval,num_intervaldim+1)
  LLK_LBFGSB <- matrix(0,num_interval,num_intervaldim+1)
  LLK_maxLik <- matrix(0,num_interval,num_intervaldim*4+1)
  
  
  for (i in 1:num_interval){
    # i<-1
    param_init <- matrix_interval[i,]
    
    ##minimization of -LogLikelihood
    
    param_optimal_BFGS <- optim(
      par = param_init,
      fn = LogLikelihood2,method = "BFGS", 
      control = list(maxit = 10000, reltol = 1e-12),
      dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    
    param_optimal_LBFGSB <- optim(
      par = param_init,
      fn = LogLikelihood2,method = "L-BFGS-B",
      lower = c(-7,0,-7, 0,0,0), upper = c(10,10,10,10,1,1),control = list(maxit = 10000, pgtol = 1e-12),
      dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    
    
    param_optimal_maxLik <- maxLik(
      logLik = function(parameters)
      {-LogLikelihood2(parameters,dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)},
      start = param_init)
    
    param_optimal_BFGS_par <- param_optimal_BFGS$par
    param_optimal_LBFGSB_par <- param_optimal_LBFGSB$par
    
    param_optimal_maxLik_class <- summary(param_optimal_maxLik)
    
    LLK_BFGS[i,1:num_intervaldim] <- param_optimal_BFGS_par
    LLK_LBFGSB[i,1:num_intervaldim] <- param_optimal_LBFGSB_par
    # Estimate
    LLK_maxLik[i,1:num_intervaldim] <- t(param_optimal_maxLik_class[["estimate"]][,1])
    # STD
    LLK_maxLik[i,(1+num_intervaldim):(2*num_intervaldim)] <- t(param_optimal_maxLik_class[["estimate"]][,2])
    # t-score
    LLK_maxLik[i,(1+2*num_intervaldim):(3*num_intervaldim)] <- t(param_optimal_maxLik_class[["estimate"]][,3])
    # probability
    LLK_maxLik[i,(1+3*num_intervaldim):(4*num_intervaldim)] <- t(param_optimal_maxLik_class[["estimate"]][,4])
    
    LLK_BFGS[i,num_intervaldim+1] <- LogLikelihood2(param_optimal_BFGS_par,dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    LLK_LBFGSB[i,num_intervaldim+1] <- LogLikelihood2(param_optimal_LBFGSB_par,dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    LLK_maxLik[i,4*num_intervaldim+1] <- LogLikelihood2(t(param_optimal_maxLik_class[["estimate"]][,1]),dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
    
    
    progress <- i/num_interval*100
    print(paste0("progress is currently ", progress , "% and country is ", name_country[count]," ranking ",count))
  }
  
  assign(paste0("LLK_BFGS_",name_country[count]),LLK_BFGS)
  assign(paste0("LLK_LBFGSB_",name_country[count]),LLK_LBFGSB)
  assign(paste0("LLK_maxLik_",name_country[count]),LLK_maxLik)
  
  
  temp_index <- which.min(LLK_maxLik[,ncol(LLK_maxLik)])
  temp_param <- LLK_maxLik[temp_index,1:num_intervaldim]
  temp_std   <- LLK_maxLik[temp_index,(1+num_intervaldim):(2*num_intervaldim)]
  temp_prb   <- LLK_maxLik[temp_index,(1+3*num_intervaldim):(4*num_intervaldim)]
  FOC_Lambda[count] <- lambda(mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
  FOC_Eta[count]  <- eta(lambda = temp_param[5],mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
  FOC_Eta2[count]   <- eta2(lambda = temp_param[5],mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
  FOC_Ratio[count] <- ratioGandN(mu_i = temp_param[1],sigma_i = temp_param[2],mu_f = temp_param[3],sigma_f = temp_param[4],dataset=dataset,count=count,Nu = Nu,Nei=Nei,Nef=Nef)
  
  
  ## Estimate b
  
  
  sample_wage_Gi <- rlnorm(n = 100000, meanlog = temp_param[1], sdlog = temp_param[2])
  wageHwstar_Gi <- mean(sample_wage_Gi[sample_wage_Gi > w_star[count]])
  
  sample_wage_Gf <- rlnorm(n = 100000, meanlog = temp_param[3], sdlog = temp_param[4])
  wageHwstar_Gf <- mean(sample_wage_Gf[sample_wage_Gf > w_star[count]])
  
  b[count] <- w_star[count] - (temp_param[5]/(rho + temp_param[6]))*(
    wageHwstar_Gi - w_star[count]*G_tilde(mu = temp_param[1], sigma = temp_param[2], dataset = dataset) +
      wageHwstar_Gf - w_star[count]*G_tilde(mu = temp_param[3], sigma = temp_param[4], dataset = dataset))
  
  # table part
  
  for (t in 1:num_intervaldim){
    Table_Coefficient[2*t-1,count] <- sprintf("%.2f",temp_param[t])
    if (temp_prb[t]<0.01){
      Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.2f",temp_std[t]),")^{***}")
      
    }else if (temp_prb[t]<0.05){
      Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.2f",temp_std[t]),")^{**}")
    }else if (temp_prb[t]<0.1){
      Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.2f",temp_std[t]),")^{*}")
      
    }else{
      Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.2f",temp_std[t]),")")
      
    }
    
    
  }
  
  Table_Coefficient[13,count] <- b[count]
  Table_Coefficient[14,count] <- w_star[count]
  
  Table_Population[1,count] <- paste0(Nu[count])
  Table_Population[3,count] <- Nei[count]
  Table_Population[5,count] <- Nef[count]
  Table_Population[7,count] <- N[count]
  Table_Population[2,count] <- paste0(sprintf("%.2f",Nu[count]/N[count]*100))
  Table_Population[4,count] <- paste0(sprintf("%.2f",Nei[count]/N[count]*100))
  Table_Population[6,count] <- paste0(sprintf("%.2f",Nef[count]/N[count]*100))
  
}

print(xtable(Table_Coefficient, type = "latex"), file = "Table_Coefficient.tex")
print(xtable(Table_Population, type = "latex"), file = "Table_Population.tex")






