#
#
### Multimorbidity project
#
#
#
#
# kmodes function
#
#
#


#rv1 <- make_fake_data( K = 6, allow_overlap = TRUE, balanced = TRUE, N_PATIENTS = 3000, N_NULL = 100, NOISE_PARAMETER = -Inf, MEAN_VAL = 2.0, COR_VAL = 0.9)


#dat2 <- rv1$data


### kmodes

run_kmode <- function(df, Kval)  {
  
  kmd <- kmodes( df[, 2:27], modes=Kval + 1, iter.max=10, fast=TRUE)
  outdf <- cbind(df, "predclust" = kmd$cluster)
  
  adjri <- adjustedRandIndex(outdf$group, outdf$predclust)
 
  #print(table(outdf$group, outdf$predclust))
  #print(adjri)
  
  return(adjri)
  
}

