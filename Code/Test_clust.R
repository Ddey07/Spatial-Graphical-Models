
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

p=3
n=200

set.seed(10)
coords <- cbind(sort(runif(n,0,1)),sort(runif(n,0,1)))

nu.mat = matrix(0.5, ncol=p, nrow= p)
phi.diag= runif(p,0.1,5)
phi.mat = diag(phi.diag)
###Creating phi matrix
for(i in 1:(p-1)){
  for (j in (i+1): p){
    phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2)
    phi.mat[j,i] = phi.mat[i,j]
  }
}

##Creating smoothness matrix
sigma.diag=runif(p,0.2,2)
sigma.mat=diag(sigma.diag)

R_V=matrix(nc=p, c(1, 0.5, 0.2, 0.5, 1, 0.6, 0.2, 0.6, 1))

for(i in 1:(p-1)){
  for (j in (i+1): p){
    sigma.mat[i,j]= sqrt(sigma.mat[i,i] * sigma.mat[j,j]) *( (phi.mat[i,i]^nu.mat[i,i]) * (phi.mat[j,j]^nu.mat[j,j])) * ((1/ phi.mat[i,j]) ^(2*nu.mat[i,j])) * R_V[i,j] 
    sigma.mat[j,i] = sigma.mat[i,j]
  }
}

#Calculating distance
D <-  as.matrix(dist(coords))

#Using corollary 3 to get a valid correlation structure
SIGMA <- matrix(ncol=n*p, nrow=n*p)

for (i in 1:p){
  for(j in i:p){
    idx=c(((i-1)*n+1):(i*n))
    jdx=c(((j-1)*n+1):(j*n))
    
    SIGMA[idx,jdx] = sigma.mat[i,j]*geoR::matern(D,phi= 1/phi.mat[i,j], kappa= nu.mat[i,j])
    SIGMA[jdx,idx] = SIGMA[idx,jdx]
  }
}

save(SIGMA,file="sigma.RData")
quit('no')