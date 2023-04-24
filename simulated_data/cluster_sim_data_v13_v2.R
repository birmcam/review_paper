#
### 05/05/2021
### runs clustering algorithms on full simulated dataset, then runs on dataset of observations with 2 or more conditions
#
#
#
#


#setwd("U:/Multimorbid/R/Feb2020")
setwd("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy")

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


source("C:/Users/linda/Documents/UoB/R/linda_sim_data_round3_tidied.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/cluster_method_lca.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/cluster_method_kmode.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/cluster_method_hca.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/cluster_method_mcakmeans.R", echo=TRUE, max.deparse.length=5000)
source("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/cluster_method_kmean.R", echo=TRUE, max.deparse.length=5000)


runmeth <- function(Kx = 3,
                     allow_overlapx = TRUE, 
                     balancedx = TRUE, 
                     N_PATIENTSx = 3000, 
                     N_NULLx = 100,
                     NOISE_PARAMETERx = -Inf, 
                     MEAN_VALx = 2.0,
                     COR_VALx = 0.9) {
  intmat <- rep(NA, 10) 
  
    #tryCatch({
  
      # generate simulated data
      rv1 <- make_fake_data( K = Kx, allow_overlap = allow_overlapx, balanced = balancedx, N_PATIENTS = N_PATIENTSx, N_NULL = N_NULLx,
                           NOISE_PARAMETER = NOISE_PARAMETERx, MEAN_VAL = MEAN_VALx, COR_VAL = COR_VALx )
      
      # run clustering algorithms on full dataset  
      infile <- rv1$data
      
      intmat[[ 1]] <- run_lca(infile, Kval = Kx, verbose=FALSE, tol=1e-5, calc.se=TRUE)
      
      kmodes_output <- NA
      # 20 restarts
      for(i in 1:20){
        kmodes_output <- tryCatch({
          run_kmode(infile, Kval = Kx)
        }, error=function(e) NA)
        if(!is.na(kmodes_output)) break;
      }
      intmat[[ 2]] <- kmodes_output
    
      intmat[[ 3]] <- run_hca(infile, Kval = Kx)
      intmat[[ 4]] <- run_mcakmeans(infile, Kval = Kx, nstart=1)
      intmat[[ 5]] <- run_kmean(infile, Kval = Kx)
       
    #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      outdf_all <- cbind(infile, intmat[[1]], intmat[[2]], intmat[[3]], intmat[[4]], intmat[[5]])
      outri_all <- apply(outdf_all[ , 29:33], 2, function(x) adjustedRandIndex(outdf_all$group, x))
  
  
  
      # run clustering algorithms on dataset of observations with two or more conditions 
      infile_2plus <- rv1$data
      infile_2plus <- infile_2plus[ rowSums(infile_2plus[ , 2:27])>1, ]
  
      intmat[[ 6]] <- run_lca(infile_2plus, Kval = Kx, verbose=FALSE, tol=1e-5, calc.se=TRUE)
  
      kmodes_output <- NA
      # 20 restarts
      for(i in 1:20){
        kmodes_output <- tryCatch({
          run_kmode(infile_2plus, Kval = Kx)
        }, error=function(e) NA)
      if(!is.na(kmodes_output)) break;
      }
      intmat[[ 7]] <- kmodes_output
  
      intmat[[ 8]] <- run_hca(infile_2plus, Kval = Kx)
      intmat[[ 9]] <- run_mcakmeans(infile_2plus, Kval = Kx, nstart=1)
      intmat[[10]] <- run_kmean(infile_2plus, Kval = Kx)
  
      outdf_2plus <- cbind(infile_2plus, intmat[[6]], intmat[[7]], intmat[[8]], intmat[[9]], intmat[[10]])
      outri_2plus <- apply(outdf_2plus[ , 29:33], 2, function(x) adjustedRandIndex(outdf_2plus$group, x))
  
  
  endlist <- list(outdf_all, outri_all, outdf_2plus, outri_2plus)

  return(endlist)
}

#zz <- runmeth(Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.9) 
  

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
  
  #return(intmat)
  
  df_all    <- intmat[ 1 : NSIMS]
  df_2plus  <- intmat[ ((2 * NSIMS) + 1) : (3 * NSIMS)]
  
  ri_all    <- data.frame(matrix(ncol = 6, nrow = 0))
  ri_2plus  <- data.frame(matrix(ncol = 6, nrow = 0))
  for (i in 1:NSIMS) {
   ri_all   <- rbind(ri_all, intmat[[NSIMS + i]])
   ri_2plus <- rbind(ri_2plus, intmat[[(3 * NSIMS) + i]])
  }
  cnames <- c("LCA", "kmode", "HCA", "MCAkm", "kmean")
  colnames(ri_all)   <- cnames
  colnames(ri_2plus) <- cnames
  
  ri_all$NumClusters      <- Kx + 1
  ri_all$Overlap          <- allow_overlapx
  ri_all$Balanced         <- balancedx
  ri_all$Noise            <- NOISE_PARAMETERx
  ri_all$Prevalence       <- MEAN_VALx
  ri_all$Corr             <- COR_VALx
  ri_all$TotalPatients    <- N_PATIENTSx + N_NULLx
  ri_all$NumInClusters    <- N_PATIENTSx
  ri_all$NumNull          <- N_NULLx
  
  ri_2plus$NumClusters    <- Kx + 1
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


# 1b - vary correlation with prev=-1.5 and fixed noise (=-4.0)
system.time(out1 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.9) )
system.time(out2 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.7) ) 
system.time(out3 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out4 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.3) ) 

save(out1,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out1.RData")
save(out2,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out2.RData")
save(out3,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out3.RData")
save(out4,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out4.RData")

rm(out1, out2, out3, out4)



# 2b - vary noise with prev=-1.5 and fixed corr (=0.5)
#system.time(out5 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
#                             NOISE_PARAMETERx = -Inf, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out6 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -5.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out7 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -3.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out8 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out9 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -1.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 

save(out5,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out5.RData")
save(out6,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out6.RData")
save(out7,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out7.RData")
save(out8,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out8.RData")
save(out9,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out9.RData")

rm(out5, out6, out7, out8, out9)



# 3 - vary proportion of patients not in a cluster, with fixed prev (=-1.5), noise (=-4.0) and corr (=0.5)
system.time(out10 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7920, N_NULLx = 80, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 1%
system.time(out11 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7840, N_NULLx = 160, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 2%
system.time(out12 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7600, N_NULLx = 400, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 5%
system.time(out13 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 7200, N_NULLx = 800, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 10%
system.time(out14 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6400, N_NULLx = 1600, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 20%
system.time(out15 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 5600, N_NULLx = 2400, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 30%
system.time(out16 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 4800, N_NULLx = 3200, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 40%
system.time(out17 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 4000, N_NULLx = 4000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) # 50%

save(out10,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out10.RData")
save(out11,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out11.RData")
save(out12,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out12.RData")
save(out13,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out13.RData")
save(out14,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out14.RData")
save(out15,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out15.RData")
save(out16,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out16.RData")
save(out17,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out17.RData")

rm(out10, out11, out12, out13, out14, out15, out16, out17)



# 4 - vary prevalence with fixed noise (=-2.5) and fixed corr (=0.5)
system.time(out18 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -3.0, COR_VALx = 0.5) )
system.time(out19 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -2.0, COR_VALx = 0.5) ) 
system.time(out20 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx = -1.0, COR_VALx = 0.5) ) 
system.time(out21 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -2.5, MEAN_VALx =  0.0, COR_VALx = 0.5) )
system.time(out22 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                              NOISE_PARAMETERx = -2.5, MEAN_VALx =  1.0, COR_VALx = 0.5) ) 
system.time(out23 <-  run_all(NSIMS=1000, Kx = 3, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                              NOISE_PARAMETERx = -2.5, MEAN_VALx =  2.0, COR_VALx = 0.5) )

save(out18,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out18.RData")
save(out19,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out19.RData")
save(out20,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out20.RData")
save(out21,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out21.RData")
save(out22,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out22.RData")
save(out23,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out23.RData")

rm(out18, out19, out20, out21, out22, out23)



# 5 - vary number of clusters, with fixed prev (=-1.5), noise (=-4.0) and corr (=0.5)
system.time(out24 <-  run_all(NSIMS=1000, Kx = 2, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out25 <-  run_all(NSIMS=1000, Kx = 4, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )

save(out24,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out24.RData")
save(out25,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out25.RData")

rm(out24, out25)


system.time(out26 <-  run_all(NSIMS=1000, Kx = 5, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                              NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )
system.time(out27 <-  run_all(NSIMS=1000, Kx = 6, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                              NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )
system.time(out28 <-  run_all(NSIMS=1000, Kx = 7, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) ) 
system.time(out29 <-  run_all(NSIMS=1000, Kx = 8, allow_overlapx = TRUE, balancedx = TRUE, N_PATIENTSx = 6000, N_NULLx = 2000, 
                             NOISE_PARAMETERx = -4.0, MEAN_VALx = -1.5, COR_VALx = 0.5) )

save(out26,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out26.RData")
save(out27,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out27.RData")
save(out28,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out28.RData")
save(out29,  file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/sim_20210506_out29.RData")

rm(out26, out27, out28, out29)


# run this afterwards

stopCluster(cl)

## END CLUSTER CODE



# find how many observations in each simulated dataset have > 1 conditions

# load(file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/n2plus_sim_20210429_out23.RData")
# 
# n2plus <- vector(mode="numeric", length=1000)
# 
# for (i in 1:1000) {
#   temp <- out23[[1]][[i]]
#   n2plus[i] <- nrow( temp)
# }
# out23[[2]] <- cbind(out23[[2]], n2plus)
# 
# save(out23, file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/n2plus_sim_20210429_out23.RData")



# to get summary stats on adj RIs
NSIMS <- 1000
j <- ceiling( NSIMS*0.5 - 1.96*sqrt(NSIMS*0.5*0.5) )
k <- ceiling( NSIMS*0.5 + 1.96*sqrt(NSIMS*0.5*0.5) )

# loop over a block of files (files are too big to load all of them at once)

# # set up empty data frame to copy summary results into
#simout_0506 <- setNames( data.frame( matrix(, , ncol=20, ) ), c("Algorithm", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", "TotalPatients", 
#                                                                  "NumInClusters", "NumNull", "Num2plus", "MeanAdjRI", "SEAdjRI", "CIMeanLower", "CIMeanUpper",
#                                                                  "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )

simout_0506 <- setNames( data.frame( matrix(, , ncol=20, ) ), c("Algorithm", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", "TotalPatients", 
                                                                "NumInClusters", "NumNull", "MeanAdjRI", "SEAdjRI", "CIMeanLower", "CIMeanUpper",
                                                                "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )



file.names <- c("outputs/sim_20210506_out25.RData", "outputs/sim_20210506_out26.RData", "outputs/sim_20210506_out27.Rdata")

lapply(file.names, load, .GlobalEnv)
big.list <- list(out25[[4]], out26[[4]], out27[[4]])

for (i in 1:3) {
  temp <- big.list[[i]]
  outdf <- setNames( data.frame( matrix(, nrow=5 , ncol=20, ) ), c("Algorithm", "NumClusters", "Overlap", "Balanced", "Noise", "Prevalence", "Corr", "TotalPatients", 
                                                                   "NumInClusters", "NumNull", "MeanAdjRI", "SEAdjRI", "CIMeanLower", "CIMeanUpper",
                                                                   "MedianAdjRI", "CIMedianLower", "CIMedianUpper", "LowerQuartile", "UpperQuartile") )
  
  outdf$MeanAdjRI     <- colMeans(temp[ , 1:5], na.rm=TRUE) 
  outdf$SEAdjRI       <- colSds(as.matrix(temp[ , 1:5])) / sqrt(nrow(temp))
  outdf$CIMeanLower   <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(0.025), type=4)
  outdf$CIMeanUpper   <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(0.975), type=4)
  outdf$MedianAdjRI   <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(0.5), type=4)
  outdf$CIMedianLower <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(j/NSIMS), type=4)
  outdf$CIMedianUpper <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(k/NSIMS), type=4)
  outdf$LowerQuartile <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(0.25), type=4)
  outdf$UpperQuartile <- colQuantiles(as.matrix(temp[ , 1:5]), probs=c(0.75), type=4)
  outdf$Algorithm     <- c("LCA", "kmode", "HCA", "MCAk", "kmeans")
  outdf$NumClusters   <- temp$NumClusters[1]
  outdf$Overlap       <- temp$Overlap[1]
  outdf$Balanced      <- temp$Balanced[1]
  outdf$Noise         <- temp$Noise[1]
  outdf$Prevalence    <- temp$Prevalence[1]
  outdf$Corr          <- temp$Corr[1]
  outdf$TotalPatients <- temp$TotalPatients[1]
  outdf$NumInClusters <- temp$NumInClusters[1]
  outdf$NumNull       <- temp$NumNull[1]
  #outdf$Num2plus     <- round(mean(temp$n2plus), 0)
  
  simout_0506         <- rbind(simout_0506, outdf)

}

save(simout_0506, file="C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_0506.RData")

                 




### plots

library(ggplot2)

load("C:/Users/linda/Documents/UoB/R/ZipFileCode - Copy/outputs/simout_0506.RData")
simout <- simout_0506
simout <- simout[ simout$Algorithm!="kmode", ]


pd <- position_dodge(0.05) # to move error bars .05 to the left and right

# vary correlation at fixed prevalence and noise
plotdf <- simout[ simout$Prevalence==-1.5 & simout$Noise==-4.0 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==6000 & simout$NumClusters==4, ]
plotdf$lcl <- plotdf$MeanAdjRI - (1.96 * plotdf$SEAdjRI)
plotdf$ucl <- plotdf$MeanAdjRI + (1.96 * plotdf$SEAdjRI)
# # plot mean RI
# ggplot(plotdf, aes(x=Corr, y = MeanAdjRI)) +
#      geom_line(aes(color = Algorithm), position=pd) +
#      geom_point(aes(color = Algorithm), position=pd) +
#      geom_errorbar(aes(ymin = lcl, ymax = ucl, color = Algorithm), width=0.1, position=pd) +
#      #geom_ribbon(data=plotdf, aes(ymin = CIMeanLower, ymax = CIMeanUpper, color = Algorithm), alpha=0.15) +
#      xlim(0, 1) +
#      ylim(0, 1) + 
#      labs(x="Correlation", y="Mean Adjusted Rand Index", title="8000 patients, 4 clusters, overlap,\n balanced, noise ~ 0.5%, prevalence ~ 15%") +
#      theme(plot.title = element_text(hjust = 0.5))
# plot median RI
ggplot(plotdf, aes(x=Corr, y = MedianAdjRI)) +
  geom_line(aes(color = Algorithm), position=pd) +
  geom_point(aes(color = Algorithm), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = Algorithm), width=0.1, position=pd) +
  xlim(0, 1) +
  ylim(0, 1) + 
  #labs(x="Correlation", y="Median adjusted Rand Index", title="8000 patients, 4 clusters, overlap,\n balanced, noise ~ 0.5%, prevalence ~ 15%") +
  labs(x="Correlation", y="Median \naRI") +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.5), "cm"), 
        axis.title.y=element_text(angle=0, hjust=0))
#dev.copy(jpeg, 'plot1_corr.jpg')
#dev.off()



# vary noise with fixed prevalence (=-1.5) and correlation (=0.5)
plotdf <- simout[ simout$Prevalence==-1.5 & (simout$Noise==-Inf | (simout$Noise>=-5.0 & simout$Noise<=-1.0)) & simout$NumInClusters==6000 & simout$NumClusters==4 & simout$Corr==0.5, ]
plotdf[plotdf$Noise==-Inf, "Noise"] <- "none"
plotdf[plotdf$Noise==-5.0, "Noise"] <- "0.1%"
plotdf[plotdf$Noise==-4.0, "Noise"] <- "0.5%"
plotdf[plotdf$Noise==-3.0, "Noise"] <- "1.5%"
plotdf[plotdf$Noise==-2.0, "Noise"] <- "8%"
plotdf[plotdf$Noise==-1.0, "Noise"] <- "25%"
plotdf$Noise <- as.factor(plotdf$Noise)
plotdf$Noise <- ordered(plotdf$Noise, levels = c("none", "0.1%", "0.5%", "1.5%", "8%", "25%"))
# plot median RI
ggplot(plotdf, aes(x=Noise, y = MedianAdjRI, group=Algorithm)) +
  geom_point(aes(color = Algorithm), position=pd) +
  geom_line(aes(color = Algorithm), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = Algorithm), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Noise", y="Median adjusted Rand Index", title="8000 patients, 4 clusters, overlap,\n balanced, correlation=0.5, prevalence ~ 15%") +
  labs(x="Noise", y="Median \naRI") +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.5), "cm"), 
        axis.title.y=element_text(angle=0, hjust=0))
#dev.copy(jpeg, 'plot2_noise.jpg')
#dev.off()



# vary prevalence with fixed noise (=-2.5) and correlation (=0.5)
plotdf <- simout[ simout$Noise==-2.5 & (simout$Prevalence>=-3.0 & simout$Prevalence<=2.0) & simout$NumInClusters==6000 & simout$NumClusters==4 & simout$Corr==0.5, ]
plotdf[plotdf$Prevalence==-3.0, "Prevalence"] <- "1.5%"
plotdf[plotdf$Prevalence==-2.0, "Prevalence"] <- "8%"
plotdf[plotdf$Prevalence==-1.0, "Prevalence"] <- "25%"
plotdf[plotdf$Prevalence== 0.0, "Prevalence"] <- "50%"
plotdf[plotdf$Prevalence== 1.0, "Prevalence"] <- "75%"
plotdf[plotdf$Prevalence== 2.0, "Prevalence"] <- "90%"
plotdf$Prevalence <- as.factor(plotdf$Prevalence)
plotdf$Prevalence <- ordered(plotdf$Prevalence, levels = c("1.5%", "8%", "25%", "50%", "75%", "90%"))
# plot median RI
ggplot(plotdf, aes(x=Prevalence, y = MedianAdjRI, group=Algorithm)) +
  geom_line(aes(color = Algorithm), position=pd) +
  geom_point(aes(color = Algorithm), position=pd) +
  geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = Algorithm), width=0.1, position=pd) +
  ylim(0, 1) + 
  #labs(x="Prevalence", y="Median adjusted Rand Index", title="8000 patients, 4 clusters, overlap,\n balanced, correlation=0.5, noise ~ 4%") +
  labs(x="Within cluster disease prevalence", y="Median \naRI") +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.5), "cm"), 
        axis.title.y=element_text(angle=0, hjust=0))
#dev.copy(jpeg, 'plot3_prev.jpg')
#dev.off()



# vary proportion of patients not in a cluster
#plotdf <- simout[ simout$Prevalence==-1.5 & simout$Noise==-4.0 & simout$Corr==0.5 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$TotalPatients==8000 & simout$NumClusters==4, ]
#plotdf$PropNull <- round( ( plotdf$NumNull / (plotdf$TotalPatients) ) * 100, 2)
# plot median RI
#ggplot(plotdf, aes(x=PropNull, y = MedianAdjRI)) +
#   geom_line(aes(color = Algorithm), position=pd) +
#   geom_point(aes(color = Algorithm), position=pd) +
#   geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = Algorithm), width=0.1, position=pd) +
#   xlim(0, 51) +
#   ylim(0, 1) + 
#   labs(x="Proportion of patients not in a cluster (%)", y="Median adjusted Rand Index", title="8000 patients, 4 clusters, overlap,\n balanced, noise ~ 0.5%, prevalence ~ 15%") +
#   theme(plot.title = element_text(hjust = 0.5))
# dev.copy(jpeg, 'plot4_propcluster.jpg')
# dev.off()



# vary number of clusters with fixed prevalence (=-1.5), noise (=-4.0) and corr (=0.5)
# plotdf <- simout[ simout$Prevalence==-1.5 & simout$Noise==-4.0 & simout$Corr==0.5 & simout$Balanced==TRUE & simout$Overlap==TRUE & simout$NumInClusters==6000 & simout$NumNull==2000, ]
# # plot median RI
# ggplot(plotdf, aes(x=NumClusters, y = MedianAdjRI)) +
#   geom_line(aes(color = Algorithm), position=pd) +
#   geom_point(aes(color = Algorithm), position=pd) +
#   geom_errorbar(aes(ymin = LowerQuartile, ymax = UpperQuartile, color = Algorithm), width=0.1, position=pd) +
#   ylim(0, 1) + 
#   #labs(x="Number of clusters", y="Median adjusted Rand Index", title="8000 patients, overlap,\n balanced, noise ~ 0.5%, prevalence ~ 15%") +
#   labs(x="Number of clusters found by algorithm", y="Median \naRI") +
#   theme_bw() +
#   theme(legend.position="bottom", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
#         plot.margin = unit(c(0.5, 0.2, 0.2, 0.5), "cm"), 
#         axis.title.y=element_text(angle=0, hjust=0))
#dev.copy(jpeg, 'plot5_numcluster.jpg')
#dev.off()




