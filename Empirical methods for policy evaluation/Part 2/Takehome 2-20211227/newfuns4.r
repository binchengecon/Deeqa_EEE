library(tidyverse)
library(bbmle)





# like for log normal, formal and informal
likef <- function(a, b) {
  mu <- a
  sigma <- b
  -sum(dlnorm(dataf[,2], meanlog = mu, sdlog = sigma, log = T)) #log=T for log-sum
}



likei <- function(a, b) {
  mu <- a
  sigma <- b
  -sum(dlnorm(datai[,2], meanlog = mu, sdlog = sigma, log = T)) #log=T for log-sum
}





# inverse
inv <- function(x) {
  y <- x^-1
  return(y)
}




# estimate Q4
like <- function(a, b) {
  lambda <- a
  eta <- b
  h <- lambda*(Gf+Gi)
  -(-(n*log(lambda*Gi*eta+eta*eta+eta*lambda*Gf))
    +nf*log(eta*lambda*Gf)
    +ni*log(eta*lambda*Gi)
    +nu*log(eta*eta)+nu*log(h)
    +sum(gf)
    +sum(gi)
    -nf*log(Gf)
    -ni*log(Gi)
    -h*sum(tu))
}



