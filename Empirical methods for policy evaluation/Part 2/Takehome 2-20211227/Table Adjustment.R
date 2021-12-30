
# load("test.RData")
# load("Result Collection.RData")

# 
# LLK_maxLik <- LLK_maxLik_chile
# count <- 1

# LLK_maxLik <- LLK_maxLik_argentina
# count <- 2
# 
# LLK_maxLik <- LLK_maxLik_colombia
# count <- 3
# 
LLK_maxLik <- LLK_maxLik_mexico
count <- 4

dataset <- data %>% filter(country == count)

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
  Table_Coefficient[2*t-1,count] <- sprintf("%.4f",temp_param[t])
  if (temp_prb[t]<0.01){
    Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.4f",temp_std[t]),")^{***}")
    
  }else if (temp_prb[t]<0.05){
    Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.4f",temp_std[t]),")^{**}")
  }else if (temp_prb[t]<0.1){
    Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.4f",temp_std[t]),")^{*}")
    
  }else{
    Table_Coefficient[2*t,count] <- paste0("(",sprintf("%.4f",temp_std[t]),")")
    
  }
  
  
}

Table_Coefficient[13,count] <-sprintf("%.4f", b[count])
Table_Coefficient[14,count] <- sprintf("%.4f", w_star[count])
# Ei
Table_Coefficient[15,count] <- sprintf("%.4f", exp(temp_param[1]+temp_param[2]^2/2))
# SDi  
Table_Coefficient[16,count] <- sprintf("%.4f", sqrt((exp(temp_param[2]^2)-1)*exp(2*temp_param[1]+temp_param[2]^2)))
#Ef
Table_Coefficient[17,count] <-sprintf("%.4f", exp(temp_param[3]+temp_param[4]^2/2))
# SDf
Table_Coefficient[18,count] <-sprintf("%.4f", sqrt((exp(temp_param[4]^2)-1)*exp(2*temp_param[3]+temp_param[4]^2)))

Table_Population[1,count] <- paste0(Nu[count])
Table_Population[3,count] <- Nei[count]
Table_Population[5,count] <- Nef[count]
Table_Population[7,count] <- N[count]
Table_Population[2,count] <- paste0(sprintf("%.2f",Nu[count]/N[count]*100))
Table_Population[4,count] <- paste0(sprintf("%.2f",Nei[count]/N[count]*100))
Table_Population[6,count] <- paste0(sprintf("%.2f",Nef[count]/N[count]*100))


# 
print(xtable(Table_Coefficient, type = "latex"), file = "Table_Coefficient.tex")
print(xtable(Table_Population, type = "latex"), file = "Table_Population.tex")
