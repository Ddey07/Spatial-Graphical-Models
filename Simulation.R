rm(list=ls())

library(gRim)
library(gRbase)
library(RandomFields)

#You might need bioconductor to install 'RBGL'
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")

library(RBGL)
library(mvtnorm)
library(fields)
library(matrixcalc)
library(corrplot)


if (!require("devtools")) install.packages("devtools")
devtools::install_github("ArkajyotiSaha/BRISC")
library(BRISC)


##### Simulating from multivariate Matern model
rho <- matrix(nc=3, c(1, 0.5, 0.2, 0.5, 1, 0.6, 0.2, 0.6, 1))

model <- RMparswmX(nudiag=c(1.3, 0.7, 2), rho=rho)
set.seed(1)
n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))

z <- RFsimulate(model = model, x=coords[,1], y=coords[,2])

### Estimating marginal materns
res1 <- BRISC_estimation(coords, x=as.matrix(rep(1,length(z@data[,1]))), y=z@data[,1], cov.model="matern")

d = as.matrix(dist(coords,method = "euclidean",diag = T,upper = T))
out1<- matrix(Matern(as.numeric(d), alpha=res1$Theta[3], smoothness=res1$Theta[4], phi=res1$Theta[1]), ncol=n)

