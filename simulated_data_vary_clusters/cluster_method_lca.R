#
#
### Multimorbidity project
#
#
#
#
# LCA function
#
#
#




### LCA


run_lca <- function(df, Kval, ...)  {
  on.exit(sink(file=NULL))

  # recode presence/absence of condition as 1 = absence, 2 = presence as poLCA package needs positive integer for category labels
  varlist <- names(df)[2:27]
  df[varlist] <- df[varlist] + 1

  f <- with(df, cbind(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z) ~ 1) 

  sink(tmpfile <- tempfile())
  ## fit LCA with K + 1 class
  lca_mm <- poLCA(f, df, nclass=(Kval), ...)
                  #maxiter=1000, nrep=25
                  #tol=1e-5, na.rm=FALSE, 
                  #verbose=FALSE, calc.se=FALSE) 

  unlink(tmpfile)
  
  # add most appropriate class based on highest posterior probability of membership
  #outdf <- cbind(df, "predclust" = lca_mm$predclass)
  #adjri <- adjustedRandIndex(outdf$group, outdf$predclust)
  #return(adjri)
  
  allocate <- data.frame("lca" = lca_mm$predclass)
  return(allocate)
  
}


