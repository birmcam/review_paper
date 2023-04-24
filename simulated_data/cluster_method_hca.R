#
#
### Multimorbidity project
#
#
#
#
#  Hierarchical cluster analysis
#
#
#


#rv1 <- make_fake_data( K = 3, allow_overlap = TRUE, balanced = TRUE, N_PATIENTS = 3000, N_NULL = 100, NOISE_PARAMETER = -Inf, MEAN_VAL = 2.0, COR_VAL = 0.9)


#dat2 <- rv1$data


### HCA.  Distance measure: binary

run_hca <- function(df, Kval)  {
  
  dis_ind_bin <- dist(df[ , 2:27], method="binary", by_rows=TRUE)

#  hca.indb.c <- hclust(dis_ind_bin, method = "complete")
#  hca.indb.a <- hclust(dis_ind_bin, method = "average")
  hca.indb.w <- hclust(dis_ind_bin, method = "ward.D")

  hca.indb.w.cut <- cutree(hca.indb.w, k = Kval + 1) # cut tree at Kval + 1 cluster
  #outdf <- cbind(df, "predclust" = hca.indb.w.cut)
  #adjri <- adjustedRandIndex(outdf$group, outdf$predclust)
  #return(adjri)
  
  allocate <- data.frame("hca" = hca.indb.w.cut)
  return(allocate)
  
}


