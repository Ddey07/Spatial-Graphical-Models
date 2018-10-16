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