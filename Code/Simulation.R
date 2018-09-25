rm(list=ls())

library(gRim)
library(gRbase)
library(RandomFields)

#You might need bioconductor to install 'RBGL'
source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")

library(RBGL)
library(mvtnorm)
library(fields)
library(matrixcalc)
library(corrplot)
library(geoR)

if (!require("devtools")) install.packages("devtools")
#devtools::install_github("ArkajyotiSaha/BRISC")
library(BRISC)


##### Simulating from multivariate Matern model
cor.mat <- matrix(nc=3, c(1, 0.5, 0.2, 0.5, 1, 0.6, 0.2, 0.6, 1))
sigma.sq <- 2
phi <- 0.3
tau.sq <- 0

model <- RMparswmX(nudiag=c(0.5, 0.5, 0.5), rho=cor.mat, scale=1/phi, var=sigma.sq) 

set.seed(1)
#n <- 1000
coords <- cbind(seq(0,1,length.out = 10),seq(length.out = 10))

tau.sq <- 0
w <- RFsimulate(model = model, x=coords[,1], y=coords[,2]) 
e <- matrix(rnorm(ncol(w)*nrow(w),1,sd=sqrt(tau.sq)),ncol=ncol(w))

z <- w@data + e

#estimating marginal 
M <- list()
for(i in 1:ncol(z)){
M[[i]] <- likfit(coords=coordinates(w), data=z[,i], cov.model = "matern", fix.kappa=FALSE,ini = c(1, 0.3), fix.nugget = TRUE, nugget=tau.sq)
}


### Estimating marginal materns
#res1 <- BRISC_estimation(coordinates(w), x=as.matrix(rep(1,length(z[,1]))), y=z[,1])

d = as.matrix(dist(coordinates(w),method = "euclidean",diag = T,upper = T))

Diag.S <- list()

for(i in 1:ncol(z)){
Diag.S[[i]]<- matrix(Matern(as.numeric(d), alpha=1/M[[i]]$phi, smoothness=M[[i]]$kappa, phi=M[[i]]$sigmasq), ncol=nrow(d))
}


