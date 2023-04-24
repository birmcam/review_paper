#
### 01/06/2021
### runs clustering algorithms on full simulated dataset, then runs on dataset of observations with 2 or more conditions
### for each simulated dataset, number of clusters found by algorithm is varied
#
#
#
#


#setwd("U:/Multimorbid/R/Feb2020")
setwd("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters")

library("poLCA")
library("mclust")
library("mvtnorm")
library("klaR")
library("proxy")
library("factoextra")
library("dplyr")
library("clustrd")
library("cba")
library("apcluster")
library("matrixStats")
library("stringr")


source("C:/Users/linda/Documents/UoB/R/linda_sim_data_round3_tidied.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters/cluster_method_lca.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters/cluster_method_kmode.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters/cluster_method_hca.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters/cluster_method_mcakmeans.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/vary_clusters/cluster_method_kmean.R", echo=TRUE, max.deparse.length=5000)


runalgs <- function(inputdf, ncl) {
  
  intmat <- rep(NA, 10) 
  
  # run clustering algorithms on full dataset 
  intmat[[ 1 ]] <- run_lca(inputdf, Kval = ncl, verbose=FALSE, tol=1e-5, calc.se=TRUE)
  
  kmodes_output <- NA
  # 20 restarts
  for(i in 1:20){
    kmodes_output <- tryCatch({
      run_kmode(inputdf, Kval = ncl)
    }, error=function(e) NA)
    if(!is.na(kmodes_output)) break;
  }
  intmat[[ 2 ]] <- kmodes_output
  
  intmat[[ 3 ]] <- run_hca(inputdf, Kval = ncl)
  intmat[[ 4 ]] <- run_mcakmeans(inputdf, Kval = ncl, nstart=1)
  intmat[[ 5 ]] <- run_kmean(inputdf, Kval = ncl)
  
  names(intmat[[ 1 ]]) <- paste0("lca", ncl)
  names(intmat[[ 2 ]]) <- paste0("kmode", ncl)
  names(intmat[[ 3 ]]) <- paste0("hca", ncl)
  names(intmat[[ 4 ]]) <- paste0("mcak", ncl)
  names(intmat[[ 5 ]]) <- paste0("kmean", ncl)
  
  
  # run clustering algorithms on dataset of observations with two or more conditions 
  infile_2plus <- inputdf[ rowSums(inputdf[ , 2:27])>1, ]
  
  intmat[[ 6]] <- run_lca(infile_2plus, Kval = ncl, verbose=FALSE, tol=1e-5, calc.se=TRUE)
  
  kmodes_output <- NA
  # 20 restarts
  for(i in 1:20){
   kmodes_output <- tryCatch({
     run_kmode(infile_2plus, Kval = ncl)
   }, error=function(e) NA)
  if(!is.na(kmodes_output)) break;
  }
  intmat[[ 7]] <- kmodes_output
   
  intmat[[ 8]] <- run_hca(infile_2plus, Kval = ncl)
  intmat[[ 9]] <- run_mcakmeans(infile_2plus, Kval = ncl, nstart=1)
  intmat[[10]] <- run_kmean(infile_2plus, Kval = ncl)

  names(intmat[[ 6 ]]) <- paste0("lca", ncl)
  names(intmat[[ 7 ]]) <- paste0("kmode", ncl)
  names(intmat[[ 8 ]]) <- paste0("hca", ncl)
  names(intmat[[ 9 ]]) <- paste0("mcak", ncl)
  names(intmat[[10 ]]) <- paste0("kmean", ncl)
  
  return(intmat)
}


#rv1 <- make_fake_data(K = 4, allow_overlap = TRUE, balanced = TRUE, N_PATIENTS = 2000, N_NULL = 1, NOISE_PARAMETER = -4.0, MEAN_VAL = -1.5, COR_VAL = 0.9)
#infile <- rv1$data

#zz <- runalgs(infile, 2)


runmeth <- function(Kx = 3,
                     allow_overlapx = TRUE, 
                     balancedx = TRUE, 
                     N_PATIENTSx = 3000, 
                     N_NULLx = 100,
                     NOISE_PARAMETERx = -Inf, 
                     MEAN_VALx = 2.0,
                     COR_VALx = 0.9) {
  
  #intmat <- rep(NA, 10) 
  
      # generate simulated data
      rv1 <- make_fake_data( K = Kx, allow_overlap = allow_overlapx, balanced = balancedx, N_PATIENTS = N_PATIENTSx, N_NULL = N_NULLx,
                           NOISE_PARAMETER = NOISE_PARAMETERx, MEAN_VAL = MEAN_VALx, COR_VAL = COR_VALx )
      

      # run clustering algorithms 
      infile <- rv1$data
      infile_2plus <- infile[ rowSums(infile[ , 2:27])>1, ]
      
      intmat <- runalgs(infile, ncl = Kx - 2)
      outdf_all <- cbind(infile, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(infile_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      intmat <- runalgs(infile, ncl = Kx - 1)
      outdf_all <- cbind(outdf_all, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(outdf_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      intmat <- runalgs(infile, ncl = Kx)
      outdf_all <- cbind(outdf_all, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(outdf_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      intmat <- runalgs(infile, ncl = Kx + 1)
      outdf_all <- cbind(outdf_all, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(outdf_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      intmat <- runalgs(infile, ncl = Kx + 2)
      outdf_all <- cbind(outdf_all, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(outdf_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      intmat <- runalgs(infile, ncl = Kx + 3)
      outdf_all <- cbind(outdf_all, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outdf_2plus <- cbind(outdf_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      
      
      outri_all   <- apply(outdf_all[ , 29:58], 2, function(x) adjustedRandIndex(outdf_all$group, x))
      outri_2plus <- apply(outdf_2plus[ , 29:58], 2, function(x) adjustedRandIndex(outdf_2plus$group, x))
  
      endlist <- list(outdf_all, outri_all, outdf_2plus, outri_2plus)

  return(endlist)
}

#zz <- runmeth(Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 1, NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.9) 





## START CLUSTER CODE

require(parallel)
runmeth_wrap <- function(i, ...){
	runmeth(...)
	}



#intdf <- data.frame(LCA=0, kmodes=0, HCA=0, "MCAkmeans"=0) 
detectCores() # this should be 8
cl <- makeCluster(mc <- getOption("cl.cores", 8))

# export the workspace to each cluster
clusterExport(cl, ls())

set.seed(0)
clusterSetRNGStream(cl, 123)

#system.time(a <- parLapply(cl, seq_len(8), runmeth_wrap))

# We need these invocations for the cores to know what we're talking about
clusterEvalQ(cl,{
library("poLCA")
library("mclust")
library("mvtnorm")
library("klaR")
library("proxy")
library("factoextra")
library("dplyr")
library("clustrd")
library("cba")
library("matrixStats")
}
)

clusterExport(cl, ls())

run_all <- function(NSIMS=NSIMS,
                    Kx = Kx,
                    allow_overlapx = allow_overlapx, 
                    balancedx = balancedx, 
                    N_PATIENTSx = N_PATIENTSx, 
                    N_NULLx = N_NULLx,
                    NOISE_PARAMETERx = NOISE_PARAMETERx, 
                    MEAN_VALx = MEAN_VALx,
                    COR_VALx = COR_VALx){
  
  # Here's the main code
  system.time(a <- parLapply(cl, seq(NSIMS), runmeth_wrap,
  	Kx = Kx,
  	allow_overlapx = allow_overlapx, 
  	balancedx = balancedx, 
  	N_PATIENTSx = N_PATIENTSx, 
  	N_NULLx = N_NULLx,
  	NOISE_PARAMETERx = NOISE_PARAMETERx, 
  	MEAN_VALx = MEAN_VALx,
  	COR_VALx = COR_VALx))
  	
  outmat <- do.call(rbind, a)
  
  #return(outmat)
  
  df_all    <- outmat[ 1 : NSIMS]
  df_2plus  <- outmat[ ((2 * NSIMS) + 1) : (3 * NSIMS)]
  
  ri_all    <- data.frame(matrix(ncol = 30, nrow = 0))
  ri_2plus  <- data.frame(matrix(ncol = 30, nrow = 0))
  for (i in 1:NSIMS) {
    ri_all   <- rbind(ri_all, outmat[[NSIMS + i]])
    ri_2plus <- rbind(ri_2plus, outmat[[(3 * NSIMS) + i]])
  }
  cnames <- c("lca2", "kmode2", "hca2", "mcak2", "kmean2", "lca3", "kmode3", "hca3", "mcak3", "kmean3", "lca4", "kmode4", "hca4", "mcak4", "kmean4",
              "lca5", "kmode5", "hca5", "mcak5", "kmean5", "lca6", "kmode6", "hca6", "mcak6", "kmean6", "lca7", "kmode7", "hca7", "mcak7", "kmean7")
  colnames(ri_all)   <- cnames
  colnames(ri_2plus) <- cnames
  
  ri_all$NumClusters      <- Kx
  ri_all$Overlap          <- allow_overlapx
  ri_all$Balanced         <- balancedx
  ri_all$Noise            <- NOISE_PARAMETERx
  ri_all$Prevalence       <- MEAN_VALx
  ri_all$Corr             <- COR_VALx
  ri_all$TotalPatients    <- N_PATIENTSx + N_NULLx
  ri_all$NumInClusters    <- N_PATIENTSx
  ri_all$NumNull          <- N_NULLx
  
  ri_2plus$NumClusters    <- Kx
  ri_2plus$Overlap        <- allow_overlapx
  ri_2plus$Balanced       <- balancedx
  ri_2plus$Noise          <- NOISE_PARAMETERx
  ri_2plus$Prevalence     <- MEAN_VALx
  ri_2plus$Corr           <- COR_VALx
  ri_2plus$TotalPatients  <- N_PATIENTSx + N_NULLx
  ri_2plus$NumInClusters  <- N_PATIENTSx
  ri_2plus$NumNull        <- N_NULLx
  
  outlist <- list(df_all, df_2plus, ri_all, ri_2plus)
  return(outlist)

}

# the cluster invocation lines above need to be run first


#set.seed(0)


# 6 - 24% prev, noise 0.5%
system.time(out30 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 20, NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.0, 
                              COR_VALx = 0.9) )
save(out30,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out30.RData")


# 7 - 24% prev, noise 8%
system.time(out31 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 20, NOISE_PARAMETERx = -2.0, MEAN_VALx = -1.0, 
                              COR_VALx = 0.9) )
save(out31,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out31.RData")
rm(out31)


# 8 - 14% prev, noise 1.5%
system.time(out32 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 20, NOISE_PARAMETERx = -3.0, MEAN_VALx = -1.5, 
                              COR_VALx = 0.9) )
save(out32,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out32.RData")
rm(out32)


# 9 - 14% prev, noise 1.5%, 5 true clusters
system.time(out33 <-  run_all(NSIMS=1000, Kx = 5, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 3000, N_NULLx = 20, NOISE_PARAMETERx = -3.0, MEAN_VALx = -1.5, 
                              COR_VALx = 0.9) )
save(out33,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out33.RData")
rm(out33)


# 10 - 50% prev, noise 1.5%, 4 true clusters
system.time(out34 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 20, NOISE_PARAMETERx = -3.0, MEAN_VALx = 0.0, 
                              COR_VALx = 0.9) )
save(out34,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out34.RData")
rm(out34)


# 11 - 75% prev, noise 1.5%, 4 true clusters
system.time(out35 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 2000, N_NULLx = 20, NOISE_PARAMETERx = -3.0, MEAN_VALx = 1.0, 
                              COR_VALx = 0.9) )
save(out35,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out35.RData")
rm(out35)

# run this afterwards

stopCluster(cl)

## END CLUSTER CODE



# find how many observations in each simulated dataset have > 1 conditions

# load(file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out35.RData")
# 
# n2plus <- vector(mode="numeric", length=1000)
# 
# for (i in 1:1000) {
#   temp <- out35[[2]][[i]]
#   n2plus[i] <- nrow( temp)
# }
# out35[[4]] <- cbind(out35[[4]], n2plus)
#  
# save(out35, file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out35.RData")








# to get summary stats on adj RIs
NSIMS <- 1000
j <- ceiling( NSIMS*0.5 - 1.96*sqrt(NSIMS*0.5*0.5) )
k <- ceiling( NSIMS*0.5 + 1.96*sqrt(NSIMS*0.5*0.5) )

# loop over a block of files (files are too big to load all of them at once)

# set up empty data frame to copy summary results into
simout_0506_vary_clusters <- setNames( data.frame( matrix(, , ncol=22, ) ), c("Algorithm", "AlgorithmR", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", 
                                                                         "TotalPatients", "NumInClusters", "NumNull", "Num2plus", "Alg_NumClus", "MeanAdjRI", "SEAdjRI",
                                                                         "CIMeanLower", "CIMeanUpper", "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile",
                                                                         "UpperQuartile") )

#simout_full_0506_vary_clusters <- setNames( data.frame( matrix(, , ncol=21, ) ), c("Algorithm", "AlgorithmR", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", 
#                                                                                   "TotalPatients", "NumInClusters", "NumNull", "Alg_NumClus", "MeanAdjRI", "SEAdjRI", 
#                                                                                   "CIMeanLower", "CIMeanUpper", "MedianAdjRI", "CIMedianLower", "CIMedianUpper", 
#                                                                                   "LowerQuartile", "UpperQuartile") )



file.names <- c("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out30.RData", 
                "C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out31.RData", 
                "C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out32.Rdata",
                "C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out33.Rdata",
                "C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out34.RData",
                "C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out35.RData")
lapply(file.names, load, .GlobalEnv)
big.list <- list(out30[[4]], out31[[4]], out32[[4]], out33[[4]], out34[[4]], out35[[4]])  # use this for dataset >=2 conditions
#big.list <- list(out30[[3]], out31[[3]], out32[[3]], out33[[3]], out34[[3]], out35[[3]]) # use this for full dataset

  
  
for (i in 1:6) {

  temp <- big.list[[i]]
  
  outdf <- setNames( data.frame( matrix(, nrow=30 , ncol=22, ) ), c("Algorithm", "AlgorithmR", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", "TotalPatients", 
                                                                   "NumInClusters", "NumNull", "Num2plus", "Alg_NumClus", "MeanAdjRI", "SEAdjRI", "CIMeanLower", 
                                                                   "CIMeanUpper", "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )

  #outdf <- setNames( data.frame( matrix(, nrow=30 , ncol=21, ) ), c("Algorithm", "AlgorithmR", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", "TotalPatients", 
  #                                                                  "NumInClusters", "NumNull", "Alg_NumClus", "MeanAdjRI", "SEAdjRI", "CIMeanLower", "CIMeanUpper", 
  #                                                                  "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )
 
  outdf$MeanAdjRI     <- colMeans(temp[ , 1:30], na.rm=TRUE) 
  outdf$SEAdjRI       <- colSds(as.matrix(temp[ , 1:30])) / sqrt(nrow(temp))
  outdf$CIMeanLower   <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(0.025), type=4)
  outdf$CIMeanUpper   <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(0.975), type=4)
  outdf$MedianAdjRI   <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(0.5), type=4)
  outdf$CIMedianLower <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(j/NSIMS), type=4)
  outdf$CIMedianUpper <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(k/NSIMS), type=4)
  outdf$LowerQuartile <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(0.25), type=4)
  outdf$UpperQuartile <- colQuantiles(as.matrix(temp[ , 1:30]), probs=c(0.75), type=4)
  outdf$Algorithm     <- colnames(temp[ , 1:30])
  outdf$AlgorithmR    <- substr(outdf$Algorithm, 1, nchar(outdf$Algorithm) - 1)
  outdf$NumClusters   <- temp$NumClusters[1]
  outdf$Overlap       <- temp$Overlap[1]
  outdf$Balanced      <- temp$Balanced[1]
  outdf$Noise         <- temp$Noise[1]
  outdf$Prevalence    <- temp$Prevalence[1]
  outdf$Corr          <- temp$Corr[1]
  outdf$TotalPatients <- temp$TotalPatients[1]
  outdf$NumInClusters <- temp$NumInClusters[1]
  outdf$NumNull       <- temp$NumNull[1]
  outdf$Num2plus      <- round(mean(temp$n2plus), 0) # comment this out when running on full simulated files
  outdf$Alg_NumClus   <- str_sub(outdf$Algorithm, -1)
  
  simout_0506_vary_clusters    <- rbind(simout_0506_vary_clusters, outdf)
  #simout_full_0506_vary_clusters    <- rbind(simout_full_0506_vary_clusters, outdf)
  
}

save(simout_0506_vary_clusters, file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_0506_vary_clusters.RData")
#save(simout_full_0506_vary_clusters, file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_full_0506_vary_clusters.RData")

                 




### plots

library(ggplot2)

load("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_0506_vary_clusters.RData")
simout <- simout_0506_vary_clusters
simout <- simout[ substr(simout$Algorithm, 1, 5)!="kmode", ]

load("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_full_0506_vary_clusters.RData")
simout_full <- simout_full_0506_vary_clusters
simout_full <- simout_full[ substr(simout_full$Algorithm, 1, 5)!="kmode", ]


pd <- position_dodge(0.05) # to move error bars .05 to the left and right




# 6: 4 true clusters, noise=-4.0 (0.5%), prev=-1.0 (24%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==-1.0 & simout_full$Noise==-4.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==2000 
#                      & simout_full$NumClusters==4, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==-1.0 & simout$Noise==-4.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==2000 & simout$NumClusters==4, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 4a") +
  labs(x="Number of clusters found by algorithm", y="Median \naRI", color="Algorithm") +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.5), "cm"), 
        axis.title.y=element_text(angle=0, hjust=0))
#dev.copy(jpeg, 'plot6_full_prev1_noise4.jpg')
#dev.copy(jpeg, 'plot6_prev1_noise4.jpg')
#dev.off()



# 7: 4 true clusters, noise=-2.0 (8%), prev=-1.0 (24%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==-1.0 & simout_full$Noise==-2.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==2000 
#                       & simout_full$NumClusters==4, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==-1.0 & simout$Noise==-2.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==2000 & simout$NumClusters==4, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 7a: 2000 patients in 4 clusters \n prevalence ~ 24%, noise ~ 8% \n all patients") +
  labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 7b: 2000 patients in 4 clusters \n prevalence ~ 24%, noise ~ 8% \n >=2 conditions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.copy(jpeg, 'plot7_full_prev1_noise2.jpg')
dev.copy(jpeg, 'plot7_prev1_noise2.jpg')
dev.off()



# 8: 4 true clusters, noise=-3.0 (1.5%), prev=-1.5 (14%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==-1.5 & simout_full$Noise==-3.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==2000 
#                       & simout_full$NumClusters==4, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==-1.5 & simout$Noise==-3.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==2000 & simout$NumClusters==4, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 8a: 2000 patients in 4 clusters \n prevalence ~ 14%, noise ~ 1.5% \n all patients") +
  labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 8b: 2000 patients in 4 clusters \n prevalence ~ 14%, noise ~ 1.5% \n >=2 conditions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.copy(jpeg, 'plot8_full_4_prev15_noise3.jpg')
dev.copy(jpeg, 'plot8_4_prev15_noise3.jpg')
dev.off()



# 9: 5 true clusters, noise=-3.0 (1.5%), prev=-1.5 (14%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==-1.5 & simout_full$Noise==-3.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==3000 
#                       & simout_full$NumClusters==5, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==-1.5 & simout$Noise==-3.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==3000 & simout$NumClusters==5, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 9a: 3000 patients in 5 clusters \n prevalence ~ 14%, noise ~ 1.5% \n all patients") +
  labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 9b: 3000 patients in 5 clusters \n prevalence ~ 14%, noise ~ 1.5% \n >=2 conditions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.copy(jpeg, 'plot9_full_5_prev15_noise3.jpg')
dev.copy(jpeg, 'plot9_5_prev15_noise3.jpg')
dev.off()



# 10: 4 true clusters, noise=-3.0 (1.5%), prev=0.0 (50%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==0 & simout_full$Noise==-3.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==2000 
#                       & simout_full$NumClusters==4, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==0 & simout$Noise==-3.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==2000 & simout$NumClusters==4, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 10a: 2000 patients in 4 clusters \n prevalence ~ 50%, noise ~ 1.5% \n all patients") +
  labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 10b: 2000 patients in 4 clusters \n prevalence ~ 50%, noise ~ 1.5% \n >=2 conditions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.copy(jpeg, 'plot10_full_prev0_noise3.jpg')
dev.copy(jpeg, 'plot10_prev0_noise3.jpg')
dev.off()



# 11: 4 true clusters, noise=-3.0 (1.5%), prev=1.0 (75%)

# complete dataset
#plotdf <- simout_full[ simout_full$Prevalence==1.0 & simout_full$Noise==-3.0 & simout_full$Balanced==TRUE & simout_full$Overlap==TRUE & simout_full$NumInClusters==2000 
#                       & simout_full$NumClusters==4, ]

# >1 condition
plotdf <- simout[ simout$Prevalence==1.0 & simout$Noise==-3.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==2000 & simout$NumClusters==4, ]

# plot median RI
ggplot(plotdf, aes(x=Alg_NumClus, y = MedianAdjRI, group=AlgorithmR)) +
  geom_line(aes(color = AlgorithmR), position=pd) +
  geom_point(aes(color = AlgorithmR), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = AlgorithmR), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 11a: 2000 patients in 4 clusters \n prevalence ~ 75%, noise ~ 1.5% \n all patients") +
  labs(x="Number of clusters found by algorithm", y="Median adjusted Rand Index", color="Algorithm", title="Figure 11b: 2000 patients in 4 clusters \n prevalence ~ 75%, noise ~ 1.5% \n >=2 conditions") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.copy(jpeg, 'plot11_full_prev1_noise3.jpg')
dev.copy(jpeg, 'plot11_prev1_noise3.jpg')
dev.off()
