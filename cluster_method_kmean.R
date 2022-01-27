#
#
### Multimorbidity project
#
#
#
#
# kmeans
#
#
#


#rv1 <- make_fake_data( K = 3, allow_overlap = TRUE, balanced = TRUE, N_PATIENTS = 3000, N_NULL = 1000, NOISE_PARAMETER = -4.0, MEAN_VAL = -1.0, COR_VAL = 0.9)

#dat2 <- rv1$data


### kmeans


run_kmean <- function(df, Kval)  {

  #set.seed(123)

  kmn <- kmeans( df[, 2:27], centers=Kval + 1, nstart=10)
  outdf <- cbind(df, "predclust" = kmn$cluster)

  adjri <- adjustedRandIndex(outdf$group, outdf$predclust)

  #print(table(outdf$group, outdf$predclust))
  #print(adjri)

  return(adjri)

}

