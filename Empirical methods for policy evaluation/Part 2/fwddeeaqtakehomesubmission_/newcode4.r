rm(list=ls())
source("newfuns4.R")
setwd("D:/Dropbox/0 Development/project 2/code")


data <- read.csv("casen_chile.txt", header=F)



# time : unemployment duration.
# status : 1 = unemployed, 2 = formal, 3 = informal
data <- data %>% rename(time = V1, wage = V2, status = V3)
datau <- data %>% filter(status == 1)
dataf <- data %>% filter(status == 2)
datai <- data %>% filter(status == 3)
datae <- data %>% filter(status != 1) #employed




################## Estimation #####################

### Reservation wages for formal and informal
rwagef <- min(dataf$wage)
rwagei <- min(datai$wage)
rwage <- min(rwagef, rwagei)



### Find mu and sigma for formal and informal : LOG-NORMAL MLE

# formal
mlef <- mle2(likef, start = list(a = 5, b = 5))
summary(mlef)
paraf <- coef(mlef)

# informal
mlei <- mle2(likei, start = list(a = 5, b = 5)) 
summary(mlei)
parai <- coef(mlei)






### estimate Q4

# number of observations for each sector
n <- data %>% group_by(status) %>% count()
nu <- as.numeric(n[1,2])
nf <- as.numeric(n[2,2])
ni <- as.numeric(n[3,2])
n <- sum(n)


gf <- dlnorm(dataf[,2], meanlog = paraf[1], sdlog = paraf[2], log = T) #vector of log-pdf ln(g) : formal
gi <- dlnorm(datai[,2], meanlog = parai[1], sdlog = parai[2], log = T)
Gf <- 1-plnorm(rwage, meanlog = paraf[1], sdlog = paraf[2]) #G tilde : formal
Gi <- 1-plnorm(rwage, meanlog = parai[1], sdlog = parai[2])
tu <- datau[,1] # vector of time for unemployed only




mle <- mle2(like, start = list(a = 5, b = 5)) 
summary(mle)
para <- coef(mle)



## Find b
rho <- 0.05
lambda <- para[1]
eta <- para[2]
#draw random values from distribution
set.seed(9)
i_draw <-rlnorm(100000,mean=parai[1],sd=parai[2])
f_draw <-rlnorm(100000,mean=paraf[1],sd=paraf[2])
#truncated integral (expectation)
expe_i <- sum(i_draw[i_draw>rwage]-rwage) / length(i_draw[i_draw>rwage])
expe_f <- sum(f_draw[f_draw>rwage]-rwage) / length(f_draw[f_draw>rwage])

b <- rwage - lambda*inv(rho+eta)*(expe_i+expe_f)
b
