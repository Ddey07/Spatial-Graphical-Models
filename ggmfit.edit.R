S= SIGMAhat
n.obs=1
glist=cgens

start = NULL
eps = 1e-12
iter = 1000
details = 0

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
ans <- list(dev = dev, df = df, detK = detK, nvar = nvar, 
            S = S, n.obs = n.obs)
ans <- c(ans, xxx)