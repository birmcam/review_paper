#
### 11/3/2020 - added in additional clustering algorithm (kmeans)
#
#
#


#setwd("U:/Multimorbid/R/Feb2020")
setwd("C:/Users/linda/Documents/UoB/R")

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


source("linda_sim_data_round3_tidied.R", echo=TRUE, max.deparse.length=5000)
source("cluster_method_lca.R", echo=TRUE, max.deparse.length=5000)
source("cluster_method_kmode.R", echo=TRUE, max.deparse.length=5000)
source("cluster_method_hca.R", echo=TRUE, max.deparse.length=5000)
source("cluster_method_mcakmeans.R", echo=TRUE, max.deparse.length=5000)
source("cluster_method_kmean.R", echo=TRUE, max.deparse.length=5000)



runmeth <- function(Kx = 3,
                     allow_overlapx = TRUE, 
                     balancedx = TRUE, 
                     N_PATIENTSx = 3000, 
                     N_NULLx = 100,
                     NOISE_PARAMETERx = -Inf, 
                     MEAN_VALx = 2.0,
                     COR_VALx = 0.9) {
  intmat <- rep(NA, 5) # intmat is matrix of RIs - one row per simulated dataset
  
    #tryCatch({
    
      rv1 <- make_fake_data( K = Kx, allow_overlap = allow_overlapx, balanced = balancedx, N_PATIENTS = N_PATIENTSx, N_NULL = N_NULLx,
                           NOISE_PARAMETER = NOISE_PARAMETERx, MEAN_VAL = MEAN_VALx, COR_VAL = COR_VALx )
    
      infile <- rv1$data
    
      #intmat[ j, 1] <- run_lca(infile, Kval = Kx)
      #intmat[ j, 2] <- run_kmode(infile, Kval = Kx)
      #intmat[ j, 3] <- run_hca(infile, Kval = Kx)
      #intmat[ j, 4] <- run_mcakmeans(infile, Kval = Kx)
      
      system.time(intmat[ 1] <- run_lca(infile, Kval = Kx, verbose=FALSE, tol=1e-5, calc.se=TRUE))
      #system.time(intmat[ 2] <- run_kmode(infile, Kval = Kx))
      
      kmodes_output <- NA
      # 20 restarts
      for(i in 1:20){
        kmodes_output <- tryCatch({
          run_kmode(infile, Kval = Kx)
        }, error=function(e) NA)
        if(!is.na(kmodes_output)) break;
      }
      intmat[ 2] <- kmodes_output
    
      
      system.time(intmat[ 3] <- run_hca(infile, Kval = Kx))
      system.time(intmat[ 4] <- run_mcakmeans(infile, Kval = Kx, nstart=1))
      system.time(intmat[ 5] <- run_kmean(infile, Kval = Kx))
       
    #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  intmat
}






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
                    COR_VALx = COR_VALx
  ){
  
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
  
  	
  intmat <- do.call(rbind, a)
  colnames(intmat) <- c("LCA", "kmode", "HCA", "MCAkm", "kmean")
  
  #return(intmat) # for initial checking that summary figs looked ok and RIs were sensible
  
  j <- ceiling( NSIMS*0.5 - 1.96*sqrt(NSIMS*0.5*0.5) )
  k <- ceiling( NSIMS*0.5 + 1.96*sqrt(NSIMS*0.5*0.5) )
  
   outdf <- setNames( data.frame( matrix(, nrow=5, ncol=19, ) ), c("Algorithm", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr",
                                                                   "TotalPatients", "NumInClusters", "NumNull", "MeanAdjRI", "SEAdjRI", "CIMeanLower", "CIMeanUpper",
                                                                   "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )
   
   outdf$MeanAdjRI     <- colMeans(intmat, na.rm=TRUE)
   outdf$SEAdjRI       <- colSds(intmat) / sqrt(nrow(intmat))
   outdf$CIMeanLower   <- colQuantiles(intmat, probs=c(0.025), type=4)
   outdf$CIMeanUpper   <- colQuantiles(intmat, probs=c(0.975), type=4)
   outdf$MedianAdjRI   <- colQuantiles(intmat, probs=c(0.5), type=4)
   outdf$CIMedianLower <- colQuantiles(intmat, probs=c(j/NSIMS), type=4)
   outdf$CIMedianUpper <- colQuantiles(intmat, probs=c(k/NSIMS), type=4)
   outdf$LowerQuartile <- colQuantiles(intmat, probs=c(0.25), type=4)
   outdf$UpperQuartile <- colQuantiles(intmat, probs=c(0.75), type=4)
   outdf$Algorithm     <- colnames(intmat)
   outdf$NumClusters   <- Kx + 1
   outdf$Overlap       <- allow_overlapx
   outdf$Balanced      <- balancedx
   outdf$Noise         <- NOISE_PARAMETERx
   outdf$Prevalence    <- MEAN_VALx
   outdf$Corr          <- COR_VALx
   outdf$TotalPatients <- N_PATIENTSx + N_NULLx
   outdf$NumInClusters <- N_PATIENTSx
   outdf$NumNull       <- N_NULLx
   
   outdf
}

# the cluster invocation lines above need to be run first


set.seed(0)


# 1b - vary correlation with prev=-1.5 and fixed noise (=-4.0)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.9) ) 
system.time(out2 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.7) ) 
system.time(out3 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out4 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.3) ) 

simout_v2 <- rbind(simout_v2, out1, out2, out3, out4)
save(simout_v2, file="C:/Users/linda/Documents/UoB/R/simout_v2.RData")



# 2b - vary noise with prev=-1.5 and fixed corr (=0.5)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -Inf, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out2 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -5.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out3 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -3.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out4 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out5 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -1.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 

simout_v2 <- rbind(simout_v2, out1, out2, out3, out4, out5)
save(simout_v2, file="C:/Users/linda/Documents/UoB/R/simout_v2.RData")



# 3 - vary proportion of patients not in a cluster, with fixed prev (=-1.5), noise (=-4.0) and corr (=0.5)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7920, N_NULLx = 80, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 1%
system.time(out2 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7840, N_NULLx = 160, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 2%
system.time(out3 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7600, N_NULLx = 400, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 5%
system.time(out4 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7200, N_NULLx = 800, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 10%
system.time(out5 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6400, N_NULLx = 1600, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 20%
system.time(out6 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 5600, N_NULLx = 2400, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 30%
system.time(out7 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 4800, N_NULLx = 3200, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 40%
system.time(out8 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 4000, N_NULLx = 4000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 50%

simout_v2 <- rbind(simout_v2, out1, out2, out3, out4, out5, out6, out7, out8)
save(simout_v2, file="C:/Users/linda/Documents/UoB/R/simout_v2.RData")



# 4 - vary prevalence with fixed noise (=-2.5) and fixed corr (=0.5)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -3.0, COR_VALx = 0.5) )
system.time(out2 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -2.0, COR_VALx = 0.5) ) 
system.time(out3 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -1.0, COR_VALx = 0.5) ) 
system.time(out4 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx =  0.0, COR_VALx = 0.5) ) 
system.time(out5 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx =  1.0, COR_VALx = 0.5) ) 
system.time(out6 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx =  2.0, COR_VALx = 0.5) )

simout_v2 <- rbind(simout_v2, out1, out2, out3, out4, out5, out6)
save(simout_v2, file="C:/Users/linda/Documents/UoB/R/simout_v2.RData")



# 5 - vary number of clusters, with fixed prev (=-1.5), noise (=-4.0) and corr (=0.5)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 2, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out2 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )
system.time(out3 <-  run_all(NSIMS=1000, Kx = 5, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )
system.time(out4 <-  run_all(NSIMS=1000, Kx = 6, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )
system.time(out5 <-  run_all(NSIMS=1000, Kx = 7, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out6 <-  run_all(NSIMS=1000, Kx = 8, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 

simout_v2 <- rbind(simout_v2, out5, out6)
save(simout_v2, file="C:/Users/linda/Documents/UoB/R/simout_v2.RData")




# run this afterwards

stopCluster(cl)

## END CLUSTER CODE



