library(gRim)
library(gRbase)

#You might need bioconductior to install 'RBGL'
library(RBGL)

library(mvtnorm)

#Set your own data
setwd("~/OneDrive - Johns Hopkins University/Spatial Graphical Models")

#Loading adjacency matrix
load("./DemoData/Graph.RData")
#Loading a covariance matrix
load("./DemoData/R2.RData")

#Simulating data
D <- rmvnorm(100000,sigma=R2)
colnames(D) <- c(1:32)
D <- as.data.frame(D)
R <- cov2cor(cov.wt(D, method="ML")$cov)
#names(G1) <- c(1:32)

#finding maximal cliques
cgens <- maxClique(as(G1,"graphNEL"))$maxCliques

T1 <- Sys.time()
carcfit2 <- ggmfit(R, n=100000, cgens)
Sys.time() - T1

#Estimated covariance matrix
Rhat <- solve(carcfit2$K)

#R=Rhat on a clique 
Rhat[cgens[[24]],cgens[[24]]]
R[cgens[[24]],cgens[[24]]]

