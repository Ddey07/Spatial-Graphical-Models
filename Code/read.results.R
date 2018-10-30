### Reading results file
rm(list=ls())
library(tidyverse)


setwd("Code/Results/Results")
df <- list.files() %>%
  purrr::map(readRDS)


f <- function(x,i){
  if(is.null(x[[i]])){
    NA
  } else {
    x[[i]]
  }
}


#Log likelihood plot
lik <- sapply(df,function(x){f(x,2)})
plot(seq(0,1,length=length(list.files())),lik)

SIGMA <- df[[1]]$truesigma

lik.max=which(lik==max(lik,na.rm = TRUE))
modelhat <- df[[lik.max]]
SIGMAhat <- modelhat$sigma

cov2cor(SIGMA)[c(1,201,401),c(1,201,401)]
cov2cor(SIGMAhat)[c(1,201,401),c(1,201,401)]

##load train test data
load("../../all.RData")

pred.matern <- function(loc=c(0,1),var=2, model=modelhat,train=Y[10,], n.var=3){
  d <- apply(coords,1,function(x){dist(rbind(loc,x))})
  cov.sp <- modelhat$sigmahat[var,var]*matern(d,phi= 1/modelhat$phihat[var,var], kappa= 0.5)
  
  
  n <- ncol(model$sigma)/n.var
  p <- n.var
  sigma.var <- modelhat$sigma[(n*(var-1)+1):(n*var),(n*(var-1)+1):(n*var)]
  
  Y.var <- train[(n*(var-1)+1):(n*var)]
  
  
  cond.sigma <- cov.sp %*% solve(sigma.var)
  pred <- mean(Y.var) + cond.sigma %*% (Y.var-rep(mean(Y.var),n))
  return(as.numeric(pred))
}


Y.pred <- matrix(ncol=length(Y.test)/n.test,nrow=n.test)

for(i in 1:nrow(coords.test)){
  
  Y.pred[i,] = unlist(lapply(c(1:3),function(x){pred.matern(loc=coords.test[i,],var=x)}))
}



