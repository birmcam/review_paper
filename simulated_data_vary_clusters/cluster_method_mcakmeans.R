#
#
### Multimorbidity project
#
#
#
#
# MCA-kmeans function
#
#
#


#rv1 <- make_fake_data( K = 3, allow_overlap = TRUE, balanced = TRUE, N_PATIENTS = 3000, N_NULL = 100, NOISE_PARAMETER = -Inf, MEAN_VAL = 2.0, COR_VAL = 0.9)


#dat2 <- rv1$data


### MCA-kmeans


run_mcakmeans <- function(df, Kval, ...) {
  on.exit(sink(file=NULL))

  sink(tmpfile <- tempfile())
  mcak <- clusmca(df[ , 2:27], nclus=Kval, ndim=2, method="MCAk", ...)
  unlink(tmpfile)

  #outdf <- cbind(df, "predclust" = mcak$cluster)
  #adjri <- adjustedRandIndex(outdf$group, outdf$predclust)
  #return(adjri)

  allocate <- data.frame("mcak" = mcak$cluster)
  return(allocate)

}


