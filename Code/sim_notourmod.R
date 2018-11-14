library(gRim)
library(gRbase)
library(RandomFields)

#You might need bioconductor to install 'RBGL'
source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")
####

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



ggmfitr.edit <- function (S, n.obs, Y, glist, start = NULL,mu, eps = 1e-12, iter = 1000,
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
  Y <- as.matrix(Y)
  S_Y = (Y-mu) %*% t(Y-mu)
  ans <- list(dev = -2 * logL, df = df, logL = logL, K = K,
              S = S, n.obs = n.obs, itcount = itcount, converged = converged,
              logLvec = logLvec, logLtest=ellK(K=K,S=S_Y,n=ncol(Y)))
  return(ans)
}

ggmfit.edit <- function (S, n.obs, glist, mu, Y, start = NULL, eps = 1e-12, iter = 1000, 
                         details = 0, ...) 
{
  glist.save <- glist.num <- glist
  data.vn <- colnames(S)
  usevar <- unique.default(unlist(glist))
  zzz <- match(usevar, data.vn)
  if (any(is.na(zzz))) 
    stop("Variables ", usevar[is.na(zzz)], " not in data\n")
  usevar <- data.vn[sort(zzz)]
  S <- S[usevar, usevar, drop = FALSE]
  data.vn <- colnames(S)
  vn <- seq_along(data.vn)
  nvar <- length(vn)
  glist.num <- lapply(glist, match, data.vn)
  glen <- sapply(glist.num, length)
  ng <- length(glist.num)
  clist.num <- lapply(glist.num, function(x) vn[-x])
  clen <- sapply(clist.num, length)
  gg <- as.integer(unlist(glist.num) - 1)
  cc <- as.integer(unlist(clist.num) - 1)
  if (is.null(start)) {
    start <- diag(1/diag(S))
  }
  xxx <- .C("Cggmfit", S = S, n = as.integer(n.obs), K = start, 
            nvar = nvar, ngen = ng, glen = glen, glist = gg, clen = clen, 
            clist = cc, logL = numeric(1), eps = as.numeric(eps), 
            iter = as.integer(iter), converged = as.integer(1), details = as.integer(details), 
            PACKAGE = "gRim")
  xxx <- xxx[c("logL", "K", "iter")]
  dimnames(xxx$K) <- dimnames(S)
  detK <- det(xxx$K)
  dev <- -n.obs * log(det(S %*% xxx$K))
  df <- sum(xxx$K == 0)/2
  Y <- as.matrix(Y)
  S_Y = (Y-mu) %*% t(Y-mu)
  ans <- list(dev = dev, df = df, detK = detK, nvar = nvar, 
              S = S, n.obs = n.obs, logLtest=ellK(K=xxx$K,S=S_Y,n=ncol(Y)))
  ans <- c(ans, xxx)
  return(ans)
}


T0 <- Sys.time()

id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

seed <- seq(2,100,len=50)
sigma.phi <- c(1,2)
n.opt <- c(250,100)

setting  <- expand.grid(seed=seed,sigma.phi=sigma.phi,n.opt=n.opt)

########## Simulation from Apanasovich paper corollary 3 #########
p=3
n=setting$n.opt[id]
test.prop = 0.2

set.seed(setting$seed[id])
coords <- cbind(sort(runif(n,0,1)),sort(runif(n,0,1)))

nu.mat = matrix(0.5, ncol=p, nrow= p)
if(setting$sigma.phi[id]==1){
  phi.diag= c(3,3,3)
  sigma.diag=c(1,1,1)
} else {
  phi.diag= c(1,3,5)
  sigma.diag=c(1,3,5)
}

phi.mat = diag(phi.diag)
###Creating phi matrix
for(i in 1:(p-1)){
  for (j in (i+1): p){
    phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2)
    phi.mat[j,i] = phi.mat[i,j]
  }
}

##Creating smoothness matrix
sigma.mat=diag(sigma.diag)

R_V=matrix(nc=p, c(1, -0.5, 0.2, -0.5, 1, 0.6, 0.2, 0.6, 1))

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
Y <- rmvnorm(10, mean=rep(0,n*p), sigma=SIGMA) 

Y <- Y +1
###Creating n*p data matrix
Y.data <- matrix(Y[10,], ncol=p)

####resetting seed so that for a particular rho, we take same test locations
set.seed(12345+setting$sigma.phi)
test.sample <- sample(1:n,size=n*test.prop)

### Splitting into train test
Y.train <- Y.data[-test.sample,]
Y.test <- Y.data[test.sample,]
coords.test <- coords[test.sample,]
coords.train <- coords[-test.sample,]
exc <- rep(test.sample,each=p) + rep(seq(0,n*p,len=p+1)[-(p+1)],length(test.sample))
SIGMA.train <- SIGMA[-exc,-exc]


### Estimating marginal materns with BRISC
M <- list()
for(i in 1:ncol(Y.train)){
  M[[i]] <- BRISC_estimation(coords.train, x=as.matrix(rep(1,length(Y.train[,i]))), y=Y.train[,i], n.neighbors = 15, cov.model="exponential")
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


SIGMA.create <- function(theta,phihat, sigmahat,Distance,Y){
  
  p <- length(phihat)
  n <- length(Y)/p
  d <- 2
  D_B= theta[1]
  r_B= diag(p)
  r_B[r_B==0]=theta[2]
  r_V= diag(p)
  r_V[upper.tri(r_V)] <- theta[3:(2+(p*(p-1))/2)]
  r_V[lower.tri(r_V)] = t(r_V)[lower.tri(r_V)]
  phi.mat = diag(phihat)
  ###Creating phi matrix
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2) + D_B*(1-r_B[i,j])
      phi.mat[j,i] = phi.mat[i,j]
    }
  }
  nu.mat=matrix(nc=p,nr=p,0.5)
  W = sqrt((sigmahat * phihat^(2*diag(nu.mat)))/(gamma(diag(nu.mat))))
  ##Creating smoothness matrix
  sigma.mat=diag(sigmahat)
  
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      sigma.mat[i,j]= W[i]*W[j] * r_V[i,j] * phi.mat[i,j] ^(-(nu.mat[i,i] + nu.mat[j,j])) * gamma((nu.mat[i,i] + nu.mat[j,j])/2 + d/2) * gamma(nu.mat[i,j]) * 1/(gamma(nu.mat[i,j] + d/2)) 
      sigma.mat[j,i] = sigma.mat[i,j]
    }
  }
  
  SIGMA <- matrix(ncol=n*p, nrow=n*p)
  
  for (i in 1:p){
    for(j in i:p){
      idx=c(((i-1)*n+1):(i*n))
      jdx=c(((j-1)*n+1):(j*n))
      
      SIGMA[idx,jdx] = sigma.mat[i,j]*geoR::matern(Distance,phi= 1/phi.mat[i,j], kappa= nu.mat[i,j])
      SIGMA[jdx,idx] = SIGMA[idx,jdx]
    }
  }
  
  return(list(phihat=phi.mat,sigmahat=sigma.mat,sigma=SIGMA))
  
}


###### Apanasovich algorithm for data with 3 variables #########
lik <- function(theta, phihat, sigmahat, Y, Distance, mu){
  p <- length(phihat)
  n <- length(Y)/p
  d <- 2
  D_B= theta[1]
  r_B= diag(p)
  r_B[r_B==0]=theta[2]
  r_V= diag(p)
  r_V[upper.tri(r_V)] <- theta[3:(2+(p*(p-1))/2)]
  r_V[lower.tri(r_V)] = t(r_V)[lower.tri(r_V)]
  phi.mat = diag(phihat)
  ###Creating phi matrix
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2) + D_B*(1-r_B[i,j])
      phi.mat[j,i] = phi.mat[i,j]
    }
  }
  nu.mat=matrix(nc=p,nr=p,0.5)
  W = sqrt((sigmahat * phihat^(2*diag(nu.mat)))/(gamma(diag(nu.mat))))
  ##Creating smoothness matrix
  sigma.mat=diag(sigmahat)
  
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      sigma.mat[i,j]= W[i]*W[j] * r_V[i,j] * phi.mat[i,j] ^(-(nu.mat[i,i] + nu.mat[j,j])) * gamma((nu.mat[i,i] + nu.mat[j,j])/2 + d/2) * gamma(nu.mat[i,j]) * 1/(gamma(nu.mat[i,j] + d/2)) 
      sigma.mat[j,i] = sigma.mat[i,j]
    }
  }
  
  SIGMA <- matrix(ncol=n*p, nrow=n*p)
  
  for (i in 1:p){
    for(j in i:p){
      idx=c(((i-1)*n+1):(i*n))
      jdx=c(((j-1)*n+1):(j*n))
      
      SIGMA[idx,jdx] = sigma.mat[i,j]*geoR::matern(Distance,phi= 1/phi.mat[i,j], kappa= nu.mat[i,j])
      SIGMA[jdx,idx] = SIGMA[idx,jdx]
    }
  }
  #SIGMA <- SIGMA.create(theta,phihat, sigmahat,Distance,Y)
  Y <- as.matrix(Y)
  mu0 <- rep(mu,each=n)
  S_Y = (Y-mu0) %*% t(Y-mu0)
  K <- solve(SIGMA)
  #loglik <- -ellK(S=S_Y,K=K,n=1)
  loglik <- ifelse(is.positive.semi.definite(SIGMA),-log(dmvnorm(Y[,1],mean=mu0,sigma=SIGMA)),1000)
  return(loglik)
}

### model with r_V equicorrelated

lik.parsi=function(theta, phihat, sigmahat, Y, Distance, mu){
  p <- length(phihat)
  n <- length(Y)/p
  d <- 2
  D_B= theta[1]
  r_B= diag(p)
  r_B[r_B==0]=theta[2]
  #r_V= diag(p)
  #r_V[r_V==0]=theta[3]
  r_V= (1-theta[2])*diag(p) + theta[2]*cor(Y.train)
  phi.mat = diag(phihat)
  ###Creating phi matrix
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2) + D_B*(1-r_B[i,j])
      phi.mat[j,i] = phi.mat[i,j]
    }
  }
  nu.mat=matrix(nc=p,nr=p,0.5)
  W = sqrt((sigmahat * phihat^(2*diag(nu.mat)))/(gamma(diag(nu.mat))))
  ##Creating smoothness matrix
  sigma.mat=diag(sigmahat)
  
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      sigma.mat[i,j]= W[i]*W[j] * r_V[i,j] * phi.mat[i,j] ^(-(nu.mat[i,i] + nu.mat[j,j])) * gamma((nu.mat[i,i] + nu.mat[j,j])/2 + d/2) * gamma(nu.mat[i,j]) * 1/(gamma(nu.mat[i,j] + d/2)) 
      sigma.mat[j,i] = sigma.mat[i,j]
    }
  }
  
  SIGMA <- matrix(ncol=n*p, nrow=n*p)
  
  for (i in 1:p){
    for(j in i:p){
      idx=c(((i-1)*n+1):(i*n))
      jdx=c(((j-1)*n+1):(j*n))
      
      SIGMA[idx,jdx] = sigma.mat[i,j]*geoR::matern(Distance,phi= 1/phi.mat[i,j], kappa= nu.mat[i,j])
      SIGMA[jdx,idx] = SIGMA[idx,jdx]
    }
  }
  #SIGMA <- SIGMA.create(theta,phihat, sigmahat,Distance,Y)
  Y <- as.matrix(Y)
  mu0 <- rep(mu,each=n)
  S_Y = (Y-mu0) %*% t(Y-mu0)
  K <- solve(SIGMA)
  #loglik <- -ellK(S=S_Y,K=K,n=1)
  loglik <- ifelse(is.positive.semi.definite(SIGMA),-log(dmvnorm(Y[,1],mean=mu0,sigma=SIGMA)),1000)
  return(loglik)
}


mu.est <- c(M[[1]]$Beta,M[[2]]$Beta,M[[3]]$Beta)

S <- Sys.time()
opt <- optim(par=c(0,0.3,0.1,0.2,0.3),lik,phihat=diag(phihat.mat),sigmahat=diag(sigmahat.mat),Y=as.numeric(Y.train),Dist=as.matrix(dist(coords.train)),mu=mu.est, method="L-BFGS-B",control = list(trace=1),lower = c(0,rep(-0.98,4))) 
opt.parsi <- optim(par=c(0,0.3,0.1),lik.parsi,phihat=diag(phihat.mat),sigmahat=diag(sigmahat.mat),Y=as.numeric(Y.train),Dist=as.matrix(dist(coords.train)),mu=mu.est,method="L-BFGS-B",control = list(trace=1),lower = c(0,rep(-0.98,4))) 

t <- Sys.time()-S


model.apag= SIGMA.create(theta=opt$par,phihat =diag(phihat.mat),sigmahat=diag(sigmahat.mat),Y=as.numeric(Y.train),Dist=as.matrix(dist(coords.train)))
model.apag.parsi=SIGMA.create(theta=c(opt.parsi$par[1:2],rep(opt.parsi$par[3],3)),phihat =diag(phihat.mat),sigmahat=diag(sigmahat.mat),Y=as.numeric(Y.train),Dist=as.matrix(dist(coords.train)))

colnames(SIGMA) <- c(1:(n*p))
rownames(SIGMA) <- c(1:(n*p))

##taking an initial rho to get initial Sigma

results <- list()

model.apag$mu=mu.est
model.apag.parsi$mu=mu.est


  
  pred.apag <- function(loc=c(0,1),var=2, model=model1,train=list(coords=coords.train,Y=as.numeric(Y.train)), n.var=p){
    d <- apply(train$coords,1,function(x){dist(rbind(loc,x))})
    n <- ncol(model$sigma)/n.var
    p <- n.var
    
    cov.sp <- numeric(length(train$Y))
    for(i in 1:length(train$Y)){
      cov.sp[i] = model$sigmahat[var,ceiling(i/n)] * geoR::matern(d[ifelse(i %% n==0,n,i %% n)],phi= 1/model$phihat[var,ceiling(i/n)], kappa= 0.5)
    }
    
    sigma.var <- model$sigma
    Y.var <- train$Y
    
    cond.sigma <- t(as.matrix(cov.sp)) %*% solve(sigma.var)
    pred <- model$mu[var] + as.matrix(cond.sigma) %*% (Y.var-rep(model$mu,each=n))
    return(as.numeric(pred))
  }
  
  
  Y.pred <- matrix(ncol=ncol(Y.test),nrow=nrow(Y.test))
  Y.pred.apag <- Y.pred
  Y.pred.apag.parsi <- Y.pred
  
  for(j in 1:nrow(coords.test)){
   # Y.pred[j,] = unlist(lapply(c(1:3),function(x){pred.matern(loc=coords.test[j,],var=x)}))
    Y.pred.apag[j,]=unlist(lapply(c(1:3),function(x){pred.apag(loc=coords.test[j,],var=x,model=model.apag)}))
    Y.pred.apag.parsi[j,]=unlist(lapply(c(1:3),function(x){pred.apag(loc=coords.test[j,],var=x,model=model.apag.parsi)}))
  }
  
  results$pred.mse.apag <- mean((Y.test-Y.pred.apag)^2, na.rm=TRUE)
  results$pred.mse.apag.parsi <- mean((Y.test-Y.pred.apag.parsi)^2, na.rm=TRUE)




results$test <- Y.test
elapsed <- Sys.time() - T0
units(elapsed) <- "secs"
results$elapsed <- elapsed

n <- n*(1-test.prop)

####### Comparing methods data frame
compare.results=list(true=list(),apag=list(),apag.parsi=list())
compare.results$true = list(corr=cov2cor(SIGMA)[-exc,-exc][c(1,(n+1),(2*n+1)),c(1,(n+1),(2*n+1))],phi=phi.mat,sigma=sigma.mat,mse=0,out=Y.test, sigma=SIGMA, mu=rep(1,3))
#compare.results$ourmod = list(corr=cov2cor(model1$sigma)[c(1,(n+1),(2*n+1)),c(1,(n+1),(2*n+1))],phi=phihat.mat,sigma=sigmahat.mat,mse=results$pred.mse,out=Y.pred, sigma=model1$sigma, mu=mu.est, lik=results$lik)
compare.results$apag = list(corr=cov2cor(model.apag$sigma)[c(1,(n+1),(2*n+1)),c(1,(n+1),(2*n+1))],phi=model.apag$phihat,sigma=model.apag$sigmahat,mse=results$pred.mse.apag,out=Y.pred.apag, sigma=model.apag$sigma,mu=mu.est)
compare.results$apag.parsi = list(corr=cov2cor(model.apag.parsi$sigma)[c(1,(n+1),(2*n+1)),c(1,(n+1),(2*n+1))],phi=model.apag.parsi$phihat,sigma=model.apag.parsi$sigmahat,mse=results$pred.mse.apag.parsi,out=Y.pred.apag.parsi, sigma=model.apag.parsi$sigma,mu=mu.est)


## phi directory create
dir.path <- file.path("Results_notourmod", paste0("run-set_(",setting$sigma.phi[id],"_",setting$n.opt[id],")"))
dir.create(dir.path,showWarnings = FALSE)

results.file <- file.path(dir.path,paste0("seed_",setting$seed[id],".rds"))
saveRDS(compare.results, results.file)
# results.file <- file.path(dir.path,paste0("n_",n),paste0("seed_",setting$seed[id],".rds"))
# saveRDS(compare.results, results.file)
# 


