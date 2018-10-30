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

ellK <- function (K, S, n) 
{
  value <- (n/2) * (sum(log(eigen(K)$values)) - sum(rowSums(K * S)))
  return(value)
}



ggmfitr.edit <- function (S, n.obs, Y, glist, start = NULL, eps = 1e-12, iter = 1000, 
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
  S_Y = (Y-mean(Y)) %*% t(Y-mean(Y))
  ans <- list(dev = -2 * logL, df = df, logL = logL, K = K, 
              S = S, n.obs = n.obs, itcount = itcount, converged = converged, 
              logLvec = logLvec, logLtest=ellK(K=K,S=S_Y,n=ncol(Y)))
  return(ans)
}


T0 <- Sys.time()

id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

seed <- seq(2,100,len=50)
phi.lo <- c(0.1,0.7)
phi.hi <- c(5,7)
rho.seq <- seq(0,1,len=50)

setting  <- expand.grid(rho=rho.seq,seed=seed,phi.lo=phi.lo,phi.hi=phi.hi)

########## Simulation from Apanasovich paper corollary 3 #########
p=3
n=250
test.prop = 0.2

set.seed(setting$seed[id])
coords <- cbind(sort(runif(n,0,1)),sort(runif(n,0,1)))

nu.mat = matrix(0.5, ncol=p, nrow= p)
phi.diag= seq(setting$phi.lo[id],setting$phi.hi[id],len=p)
phi.mat = diag(phi.diag)
###Creating phi matrix
for(i in 1:(p-1)){
  for (j in (i+1): p){
    phi.mat[i,j]= sqrt((phi.mat[i,i]^2 + phi.mat[j,j]^2)/2)
    phi.mat[j,i] = phi.mat[i,j]
  }
}

##Creating smoothness matrix
sigma.diag=c(3*(1:p))
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
Y <- rmvnorm(100, mean=rep(0,n*p), sigma=SIGMA) 

Y <- Y +1
###Creating n*p data matrix
Y.data <- matrix(Y[10,], ncol=p)

####resetting seed so that for a particular rho, we take same test locations
set.seed(12345+setting$phi.lo[id] + setting$phi.hi[id])
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

colnames(SIGMA) <- c(1:(n*p))
rownames(SIGMA) <- c(1:(n*p))


##taking an initial rho to get initial Sigma

pseudo.est <- function(rho,...){
  R_V= (1-rho)*diag(p) + rho*cor(Y.train)
  for(i in 1:(p-1)){
    for (j in (i+1): p){
      sigmahat.mat[i,j]= sqrt(sigmahat.mat[i,i] * sigmahat.mat[j,j]) * R_V[i,j] 
      sigmahat.mat[j,i] = sigmahat.mat[i,j]
    }
  }
  
  SIGMAhat <- matrix(ncol=n*p, nrow=n*p)
  colnames(SIGMAhat) <- c(1:(n*p))
  rownames(SIGMAhat) <- c(1:(n*p))
  
  for (i in 1:p){
    for(j in i:p){
      idx=c(((i-1)*n+1):(i*n))
      jdx=c(((j-1)*n+1):(j*n))
      SIGMAhat[idx,jdx] = sigmahat.mat[i,j]*geoR::matern(D,phi= 1/phihat.mat[i,j], kappa= nuhat.mat[i,j])
      SIGMAhat[jdx,idx] = SIGMAhat[idx,jdx]
    }
  }
  
  if(is.positive.semi.definite(SIGMAhat)){
    #### Creating the adjacency matrix for conditional independency structure (graph) #######
    A <- matrix(0,nrow=n*p,ncol=n*p)
    for(i in (1:(n*p-1))){
      for(j in ((i+1): (n*p))){
        j.q <- (j-1) %/% n
        i.q <- (i-1) %/% n
        if((j-i) %% n ==0 | ((j-i) < n & (j.q-i.q)==0)){
          A[i,j]=1
          A[j,i]=A[i,j]
        }
      }
    }
    
    
    #finding maximal cliques
    cgens <- maxClique(as(A,"graphNEL"))$maxCliques
    
    
    #save(SIGMAhat,file="Sigmahat.RData")
    
    #T1 <- Sys.time()
    
    #T2 <- T1-T0
    carcfit2 <- ggmfitr.edit(SIGMAhat, n=1, glist=cgens,...)
    #elapsed <- Sys.time() - T1
    
    # #Estimated covariance matrix
    Rhat <- solve(carcfit2$K)
    return(list(fit=carcfit2,lik=carcfit2$logLtest,sigma=Rhat))
  } else{
    return(list(fit=NULL,lik=NULL,sigma=NULL))
  }
}

results <- pseudo.est(setting$rho[id],eps=10^(-3), Y=as.numeric(Y.data))
results$truephi <- phi.mat
results$truesig <- sigma.mat

if(!is.null(results$fit)){
model1 <- list(phihat=phihat.mat,sigmahat=sigmahat.mat,sigma=results$sigma[-exc,-exc])

pred.matern <- function(loc=c(0,1),var=2, model=model1,train=list(coords=coords.train,Y=as.numeric(Y.train)), n.var=p){
  d <- apply(train$coords,1,function(x){dist(rbind(loc,x))})
  cov.sp <- model$sigmahat[var,var]*geoR::matern(d,phi= 1/model$phihat[var,var], kappa= 0.5)
  
  
  n <- ncol(model$sigma)/n.var
  p <- n.var
  sigma.var <- model$sigma[(n*(var-1)+1):(n*var),(n*(var-1)+1):(n*var)]
  
  Y.var <- train$Y[(n*(var-1)+1):(n*var)]
  
  
  cond.sigma <- t(as.matrix(cov.sp)) %*% solve(sigma.var)
  pred <- mean(Y.var) + as.matrix(cond.sigma) %*% (Y.var-rep(mean(Y.var),n))
  return(as.numeric(pred))
}


Y.pred <- matrix(ncol=ncol(Y.test),nrow=nrow(Y.test))

for(j in 1:nrow(coords.test)){
  
  Y.pred[j,] = unlist(lapply(c(1:3),function(x){pred.matern(loc=coords.test[j,],var=x)}))
}

results$pred.mse <- mean((Y.test-Y.pred)^2, na.rm=TRUE)
} else { 
  results$pred.mse <- NULL
  }

results$rho <- setting$rho[id]
results$philim <- c(setting$phi.lo[id],setting$phi.hi[id])
results$test <- Y.test
results$pred <- Y.pred

## phi directory create
dir.path <- file.path("Results", paste0("run-phi_(",setting$phi.lo[id],",",setting$phi.hi[id],")"))
dir.create(dir.path,showWarnings = FALSE)

##rho directory create
rho.ind=which(rho.seq == setting$rho[id])
dir.create(file.path(dir.path,paste0("rho_",rho.ind)),showWarnings = FALSE)

results.file <- file.path(dir.path,paste0("rho_",rho.ind),paste0("seed_",setting$seed[id],".rds"))
saveRDS(results, results.file)



