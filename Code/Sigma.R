rm(list=ls())

#library(gRim)
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

ellK <- function (K, S, n) 
{
  value <- (n/2) * (sum(log(eigen(K)$values)) - sum(rowSums(K * S)))
  return(value)
}



ggmfitr.edit <- function (S, n.obs, glist, start = NULL, eps = 1e-12, iter = 1000, 
                          details = 0, ...) 
{
  if (is.null(start)) {
    K <- diag(1/diag(S))
  }
  else {
    K <- start
  }
  dimnames(K) <- dimnames(S)
  vn <- colnames(S)
  x <- lapply(glist, match, vn)
  varIndex = 1:nrow(K)
  itcount = 0
  if (length(x)) {
    my.complement <- function(C) return(setdiff(varIndex, 
                                                C))
    x.complements <- lapply(x, my.complement)
    if (length(x.complements[[1]]) == 0) {
      return(list(K = solve(S)))
    }
    logLvec <- NULL
    repeat {
      for (j in 1:length(x)) {
        C <- x[[j]]
        notC <- x.complements[[j]]
        K[C, C] <- solve(S[C, C, drop = FALSE]) + K[C, 
                                                    notC, drop = FALSE] %*% solve(K[notC, notC, 
                                                                                    drop = FALSE]) %*% K[notC, C, drop = FALSE]
      }
      logL <- ellK(K, S, n.obs)
      logLvec <- c(logLvec, logL)
      itcount <- itcount + 1
      if (itcount > 1) {
        if (logL - prevlogL < eps) {
          converged = TRUE
          break
        }
      }
      else {
        if (itcount == iter) {
          converged = FALSE
          break
        }
      }
      prevlogL <- logL
    }
  }
  df <- sum(K[upper.tri(K)] == 0)
  ans <- list(dev = -2 * logL, df = df, logL = logL, K = K, 
              S = S, n.obs = n.obs, itcount = itcount, converged = converged, 
              logLvec = logLvec)
  return(ans)
}


T0 <- Sys.time()
# ##### Simulating from multivariate Matern model
# cor.mat <- matrix(nc=3, c(1, 0.5, 0.2, 0.5, 1, 0.6, 0.2, 0.6, 1))
# sigma.sq <- 1
# phi <- 3
# tau.sq <- 0
# mu <- 1
# 
# model <- RMparswmX(nudiag=c(0.5, 0.5, 0.5), rho=cor.mat, scale=1/phi, var=sigma.sq) 
# 
# #set.seed(10)
# #n <- 200
# #coords <- cbind(runif(n,0,1),runif(n,0,1))
# #Setting up grid
# coords <- cbind(seq(0,1,length.out = 20), seq(0,1,length.out = 20))
# 
# #Setting nugget variance to be 0
# tau.sq <- 0
# w <- RFsimulate(model = model, x=coords[,1], y=coords[,2]) 
# ##Setting mean to be 1
# e <- matrix(rnorm(ncol(w)*nrow(w),mu,sd=sqrt(tau.sq)),ncol=ncol(w))
# 
# z <- w@data + e
# # 
# # #estimating marginal Materns with likfit
# # M <- list()
# # for(i in 1:ncol(z)){
# # M[[i]] <- likfit(coords=coordinates(w), data=z[,i], cov.model = "matern", fix.kappa=FALSE,ini = c(1.3, 3), fix.nugget = TRUE, nugget=tau.sq)
# # }
# 
# ### Estimating marginal materns wth
# M <- list()
# for(i in 1:ncol(z)){
# M[[i]] <- BRISC_estimation(coordinates(w), x=as.matrix(rep(1,length(z[,i]))), y=z[,i])
# }
# 
# d = as.matrix(dist(coordinates(w),method = "euclidean",diag = T,upper = T))
# 
# Diag.S <- list()
# 
# # for(i in 1:ncol(z)){
# # Diag.S[[i]]<- matrix(Matern(as.numeric(d), alpha=1/M[[i]]$phi, smoothness=M[[i]]$kappa, phi=M[[i]]$sigmasq), ncol=nrow(d))
# # }
# 
# ###From BRISC estimation
# for(i in 1:ncol(z)){
#   Diag.S[[i]]<- matrix(Matern(as.numeric(d), alpha=1/M[[i]]$Theta[3], smoothness=0.5, phi=M[[i]]$Theta[1]), ncol=nrow(d))
# }
# 
# 
# ######################################################################################################
# ######Trying one dimension to check the error#######################################################3#
# 
# sigma.sq <- 1
# phi <- 1
# 
# #model <- RMparswmX(nudiag=0.5, rho=1) 
# model <- RMmatern(nu=0.5)
# 
# n <- 1000
# 
# #Setting up grid
# coords <- cbind(seq(0,10,length.out = n), seq(0,10,length.out = n))
# 
# ###simulating using multivariate normal function
# R <- exp(-phi*as.matrix(dist(coords)))
# w <- rmvnorm(1, mean=rep(0,n), sigma=sigma.sq*R)
# 
# 
# #estimating marginal Materns with likfit
# fit.lik <- likfit(coords=coords, data=w, cov.model = "exponential", fix.kappa = FALSE, ini.cov.pars=c(1,1))
# fit.BRISC <- BRISC_estimation(coords, x=as.matrix(rep(1,length(w))), y=w)
# 
# #Simulating using RandomFields
# w <- RFsimulate(model=model, distance=dist(coords))
# 
# fit.lik <- likfit(coords=coords, data=w, cov.model = "matern", fix.kappa = FALSE, ini.cov.pars=c(1,1))
# 
# fit.BRISC <- BRISC_estimation(coords, x=as.matrix(rep(1,length(w))), y=w)
# 
# 

########## Simulation from Apanasovich paper corollary 3 #########
p=3
n=200

set.seed(30)
coords <- cbind(sort(runif(n,0,1)),sort(runif(n,0,1)))

nu.mat = matrix(0.5, ncol=p, nrow= p)
phi.diag= runif(p,2,6)
phi.mat = diag(phi.diag)
###Creating phi matrix
for(i in 1:(p-1)){
  for (j in (i+1): p){
    phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2)
    phi.mat[j,i] = phi.mat[i,j]
  }
}

##Creating smoothness matrix
sigma.diag=runif(p,1,10)
sigma.mat=diag(sigma.diag)

R_V=matrix(nc=p, c(1, -0.6, 0.4, -0.6, 1, 0.1, 0.4, 0.1, 1))

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

#####Simulating multivariate normal
Y <- rmvnorm(100, mean=rep(0,n*p), sigma=SIGMA) 

Y <- Y +1
###Creating n*p data matrix
Y.data <- matrix(Y[10,], ncol=p)


### Estimating marginal materns with BRISC
M <- list()
for(i in 1:ncol(Y.data)){
  M[[i]] <- BRISC_estimation(coords, x=as.matrix(rep(1,length(Y.data[,i]))), y=Y.data[,i], n.neighbors = 30, cov.model="exponential")
}

nuhat.mat = matrix(0.5, ncol=p, nrow= p)
phihat.diag= c(M[[1]]$Theta[3],M[[2]]$Theta[3], M[[3]]$Theta[3])
phihat.mat = diag(phihat.diag)

for(i in 1:(p-1)){
  for (j in (i+1): p){
    phihat.mat[i,j]= sqrt((phihat.mat[i,i]^2 + phihat.mat[j,j]^2)/2)
    phihat.mat[j,i] = phihat.mat[i,j]
  }
}

##Using estimates to initialize IPS
sigmahat.diag=c(M[[1]]$Theta[1],M[[2]]$Theta[1], M[[3]]$Theta[1])
sigmahat.mat=diag(sigmahat.diag)

colnames(SIGMA) <- c(1:(n*p))
rownames(SIGMA) <- c(1:(n*p))

n.test <- 50
coords.test <- cbind(sort(runif(n.test,0.1,0.9)),sort(runif(n.test,0.1,0.9)))
D.test <-  as.matrix(dist(coords.test))

SIGMA.test <- matrix(ncol=n.test*p, nrow=n.test*p)

for (i in 1:p){
  for(j in i:p){
    idx=c(((i-1)*n.test+1):(i*n.test))
    jdx=c(((j-1)*n.test+1):(j*n.test))
    
    SIGMA.test[idx,jdx] = sigma.mat[i,j]*geoR::matern(D.test,phi= 1/phi.mat[i,j], kappa= nu.mat[i,j])
    SIGMA.test[jdx,idx] = SIGMA.test[idx,jdx]
  }
}


Y.test <- rmvnorm(1, mean=rep(0,n.test*p), sigma=SIGMA.test) 

Y.test <- Y.test + 1

Y.test.data <- matrix(Y.test,ncol=p)



save.image("all.RData")
