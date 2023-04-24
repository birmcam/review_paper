#setwd("U:/Multimorbid/R/Feb2020/thin")
setwd("C:/Users/linda/Documents/UoB/R/THIN/finaldata")


library("poLCA")
library("mclust")
library("mvtnorm")
library("klaR")
library("proxy")
library("factoextra")
library("plyr")
library("dplyr")
library("clustrd")
library("cba")
library("NbClust")
library("entropy")
library("fpc")
library("boot")
library("reshape")
library("matrixStats")
library("caret")

options(scipen=999)



# load sample of thin data
load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/thinsample.RData")

# keep patients with two or more conditions (multimorbid)
thin.sample <- thin.sample[ rowSums(thin.sample[ , 8:66])>1, ]


# generate stratified data 
# sex: 0 = male, 1 = female
condnames <- c("pracid", "patid", "age_at_index", "townsend", "ethnicity", "ca.all", "str.tia", "arr.af", 
               "hypertension", "hf", "ihd", "valve.dis", "pvd", "aortic.an", "eczema", "psoriasis", "rhin.conjunc", 
               "sinusitis", "deaf", "blind", "cataract", "glaucoma", "amd", "peptic", "ibs", "divertic", "pancreatitis",
               "endometriosis", "pcos", "vte", "pulm.emb", "coagulopathy", "depression", "anxiety", "smi", "eating.dis", 
               "parkinsons", "migraine", "epilepsy", "hemiplegia", "oa", "osteoporosis", "gout", "ckd", "asthma", "copd", 
               "osa", "bronch", "pulm.fib", "hyperthyroid", "hypothyroid", "hiv", "bph", "erec.dys", "neurodev", "somat",
               "chron.liver", "alz.dem", "other.heart", "addict", "auto.bowel", "auto.tissue", "auto.other", "diabetes")

#tb_1844m    <- subset(thin.sample, sex==0 & ageg=="18-44", select=condnames)
#tb_1844f    <- subset(thin.sample, sex==1 & ageg=="18-44", select=condnames)
#tb_4564m    <- subset(thin.sample, sex==0 & ageg=="45-64", select=condnames)
#tb_4564f    <- subset(thin.sample, sex==1 & ageg=="45-64", select=condnames)
tb_6584m    <- subset(thin.sample, sex==0 & ageg=="65-84", select=condnames)
#tb_6584f    <- subset(thin.sample, sex==1 & ageg=="65-84", select=condnames)
#tb_85plusm  <- subset(thin.sample, sex==0 & ageg==">=85",  select=condnames)
#tb_85plusf  <- subset(thin.sample, sex==1 & ageg==">=85",  select=condnames)





### Males, aged 65-84y subset

allconds <- colnames(tb_6584m[ , 6:64])

# check if any morbidity is NOT present in the data
colSums(tb_6584m[ , allconds]==1) 


# prevalence of each disease in whole dataset
colcount <- data.frame(colSums(tb_6584m[ , allconds]==1))
colcount[2] <- round( colcount[1]/nrow(tb_6584m), 5)
names(colcount) <- c("Freq.","TotalProp")
prev_6584m <- tibble::rownames_to_column(colcount, "Morbidity")
prev_6584m <- prev_6584m[ order(-prev_6584m$TotalProp), ]


# exclude female only conditions: endometriosis, pcos
# exclude conditions with prevalence <1%
excludeconds <- c("pcos", "endometriosis", "eating.dis", "hiv", "coagulopathy", "pancreatitis", "hemiplegia", 
                  "pulm.fib", "somat", "neurodev")
colmorbidlist <- allconds[ !(allconds %in% excludeconds) ]


# create vector of short name for conditions to merge in with prevalence, profiles
Morbidity <- c("ca.all", "str.tia", "arr.af", "hypertension", "hf", "ihd", "valve.dis", "pvd", "aortic.an", 
               "eczema", "psoriasis", "rhin.conjunc", "sinusitis", "deaf", "blind", "cataract", "glaucoma", 
               "amd", "peptic", "ibs", "divertic", "vte", "pulm.emb", "depression", "anxiety", "smi", "parkinsons",
               "migraine", "epilepsy", "oa", "osteoporosis", "gout", "ckd", "asthma", "copd", "osa", "bronch", 
               "hyperthyroid", "hypothyroid", "bph", "erec.dys", "chron.liver", "alz.dem", "other.heart", "addict", 
               "auto.bowel", "auto.tissue", "auto.other", "diabetes") 

Morb <- c("Ca", "stroke", "arr.af", "hyperten", "HF", "IHD", "valve", "PVD", "aneu", "ecz", "psor", "rhin", "sinus",
          "deaf", "blind", "cat", "glau", "AMD", "peptic", "IBS", "div", "VTE", "PE", "dep", "anx", "SMI", "park", 
          "mig", "epi", "OA", "osteo", "gout", "CKD", "asth", "COPD", "OSA", "bronch", "hyperthy", "hypothy", "BPH",
          "ED", "liver", "dem", "oth.heart", "addict", "bowel", "tissue", "auto.oth", "diab") 

shortcode <- cbind(Morbidity, Morb)

prev_6584m <- merge(prev_6584m, shortcode, by="Morbidity")
prev_6584m <- prev_6584m[ , c(1, 4, 3)]
prev_6584m <- prev_6584m[order(-prev_6584m$TotalProp), ]




set.seed(20200908)



# function to run LCA algorithm
runlca <- function(inputfile) {
  
  infile_lca <- inputfile
  infile_lca[colmorbidlist] <- infile_lca[colmorbidlist] + 1
  f <- with(infile_lca, cbind(ca.all, str.tia, arr.af, hypertension, hf, ihd, valve.dis, pvd, aortic.an, eczema, 
                              psoriasis, rhin.conjunc, sinusitis, deaf, blind, cataract, glaucoma, amd, peptic, ibs,
                              divertic, vte, pulm.emb, depression, anxiety, smi, parkinsons, migraine, epilepsy, oa, 
                              osteoporosis, gout, ckd, asthma, copd, osa, bronch, hyperthyroid, hypothyroid, bph, 
                              erec.dys, chron.liver, alz.dem, other.heart, addict, auto.bowel, auto.tissue, auto.other,
                              diabetes) ~ 1)
  lca_mm <- poLCA(f, infile_lca[ , colmorbidlist], nclass=Kx, verbose=FALSE, maxiter=3000, tol=1e-5, nrep=5)
  lca_out  <- cbind(inputfile, "cluster"=lca_mm$predclass)
  return(lca_out)
}


# function to run MCA-kmeans algorithm
runmcak <- function(inputfile) {
  mcak_mm  <- clusmca(inputfile[ , colmorbidlist], nclus=Kx, ndim=2, method="MCAk", nstart=5)
  mcak_out <- cbind(inputfile, "cluster"=mcak_mm$cluster)
  return(mcak_out)
}


# function to run kmeans algorithm
runkmean <- function(inputfile) {
  kmn_mm <- kmeans(inputfile[ , colmorbidlist], centers=Kx, nstart=150, iter.max=50)
  kmn_out  <- cbind(inputfile, "cluster"=kmn_mm$cluster)
  return(kmn_out)
}


# function to run kmeans-HCA algorithm
runkmhca <- function(inputfile) {
  
  # k-means
  step1 <- kmeans(inputfile[ , colmorbidlist], centers=50, nstart=150, iter.max=50)
  step1clus <- data.frame(step1$cluster)
  colnames(step1clus) <- c("kmcluster")
  kmhca_out <- cbind(inputfile, step1clus)
  
  # then HCA
  dis <- dist(step1$centers, method="euclidean", by_rows=TRUE)
  hcares <- hclust(dis, method = "ward.D")
  step2 <- data.frame( "cluster"=cutree(hcares, k=Kx) )
  step2$kmcluster <- as.numeric( rownames(step2) )
  
  kmhca_out <- merge(kmhca_out, step2, by="kmcluster", sort=FALSE, all.x=TRUE)
  kmhca_out <- kmhca_out[ order(kmhca_out$pracid, kmhca_out$patid), ]
  kmhca_out$kmcluster <- NULL
  
  return(kmhca_out)
}



# function to calculate prevalence of disease in each cluster. Needs to have disease coded as 0/1 where 1=presence
prevcalc <- function(inputfile) {
  prevout <- data.frame( apply(inputfile[ , colmorbidlist], 2, FUN=function(z) round( tapply(z , inputfile$cluster, FUN=mean), 5) ) )
  prevout$cluster <- as.numeric( rownames(prevout) )
  prevout <- prevout[ , c(50, 1:49)]
  return(prevout)
}



# function to label clusters and calculate prevalence of disease in each cluster
clusterprev <- function(infile) {
  
  df_orig_cluster_prev <- prevcalc(infile)
  
  orig_cluster_prev <- melt(df_orig_cluster_prev, id.vars="cluster")
  orig_cluster_prev <- orig_cluster_prev[ , c(2, 1, 3)]
  colnames(orig_cluster_prev) <- c("Morbidity", "Cluster", "ClustProp")
  
  orig_cluster_total_prev <- merge(orig_cluster_prev, prev_6584m, by='Morbidity',all.x=TRUE)
  
  temp1 <- orig_cluster_total_prev
  #temp1$value1 <- (temp1$ClustProp > temp1$TotalProp) * temp1$ClustProp # only posterior prob. bigger than overall prevalence is counted
  temp1$value2 <- (temp1$ClustProp > temp1$TotalProp) * (temp1$ClustProp - temp1$TotalProp) # difference in relative prevalence
  
  temp1 <- temp1 %>%
    group_by(Cluster) %>%
    arrange(Cluster, desc(value2)) %>% # within each cluster sort by prevalence of each disease (descending order). Use difference in prevalence between cluster and overall dataset!
    filter(value2>0) %>% # select top 3 rows where cluster prevalence is greater than overall prevalence
    top_n(3, value2) %>% 
    ungroup
  
  # Report top 3 morbidities for each cluster with greatest difference between prevalencein cluster and prevalence in whole dataset
  temp2 <- temp1 %>%
    group_by(Cluster) %>%
    mutate(top3morb_3 = paste0(Morbidity, collapse =", ")) %>%
    top_n(1, value2) %>%
    select(Cluster, top3morb_3)
  
  orig_cluster_labels <- merge(orig_cluster_prev, temp2, by='Cluster')
  orig_cluster_labels$ClusterLabel <- paste(orig_cluster_labels$Cluster, ": " , orig_cluster_labels$top3morb_3)
  orig_cluster_labels <- orig_cluster_labels[ , c(2, 3, 5)] 
  
  df_orig_cluster_prev <- t(df_orig_cluster_prev)
  df_orig_cluster_prev <- df_orig_cluster_prev[ -1, ]
  
  outlist <- list(df_orig_cluster_prev, temp2)
  return(outlist)
}


# function to calculate PCC matrix between disease profiles in each cluster from two clustering algorithms.  
# finds most similar (highest PCC>=0.5) cluster between the two methods. If PCC<0.5  implies there is no similar cluster  
calcpcc <- function(alg1, alg2) {
  
  pccmat <- round( cor(alg1, alg2), 8)
  
  maxpos <- matrix(, nrow=5, ncol=2)
  rowInd <- vector(mode="numeric", length=Kx)
  colInd <- vector(mode="numeric", length=Kx)
  pccvec <- vector(mode="numeric", length=Kx)
  
  a <- pccmat
  print(a)
  
  for (i in 1:Kx) {
    if (i==1) {
      maxpos[i, ] <- which(a==max(a), arr.ind=TRUE) 
      rowInd[i]   <- maxpos[i, 1]
      colInd[i]   <- maxpos[i, 2]
      pccvec[i]   <- a[ rowInd[i], colInd[i] ]
      #print(paste("Row: ", rownames(a)[ rowInd[i]], ", Col: ",  colnames(a)[colInd[i]], ", PCC=", a[ rowInd[i], colInd[i] ] ) )
    }   else {
      maxpos[i, ] <- which( a==max( a[ -rowInd,  -colInd ] ), arr.ind=TRUE ) 
      rowInd[i]   <- maxpos[i, 1]
      colInd[i]   <- maxpos[i, 2]
      pccvec[i]   <- a[ rowInd[i], colInd[i] ]
      #print(paste("Row: ", rownames(a)[rowInd[i]], ", Col: ",  colnames(a)[colInd[i]], ", PCC=", a[ rowInd[i], colInd[i] ] ) )
    }  
  }
  matchout <- data.frame(cbind(rowInd, colInd, pccvec))
  colnames(matchout) <- c("Cluster_orig", "Cluster_bstrap", "PCC")
  #matchout <- merge(matchout, labels1, by.x="Cluster_Alg1", by.y="Cluster")
  #matchout <- merge(matchout, labels2, by.x="Cluster_Alg2", by.y="Cluster")
  #colnames(matchout)[ colnames(matchout) %in% c("top3morb_3.x", "top3morb_3.y") ] <- c("Labels_Alg1", "Labels_Alg2")
  #matchout <- matchout[ , c(2, 4, 1, 5, 3)]
  #matchout$PCC <- ifelse(matchout$PCC<0.5, NA, matchout$PCC) 
  matchout <- matchout[order(matchout$Cluster_orig), ]
  return(matchout)
}





Kx = 4 # 4 clusters looks most appropriate from LCA entropy results


lca_orig <- runlca(tb_6584m)
round(prop.table(table(lca_orig$cluster)), 4)
lca_orig_profile <- clusterprev(lca_orig)[[1]]
lca_orig_labels  <- clusterprev(lca_orig)[[2]]

mcak_orig <- runmcak(tb_6584m)
round(prop.table(table(mcak_orig$cluster)), 4)
mcak_orig_profile <- clusterprev(mcak_orig)[[1]]
mcak_orig_labels  <- clusterprev(mcak_orig)[[2]]

kmean_orig <- runkmean(tb_6584m)
round(prop.table(table(kmean_orig$cluster)), 4)
kmean_orig_profile <- clusterprev(kmean_orig)[[1]]
kmean_orig_labels  <- clusterprev(kmean_orig)[[2]]

kmhca_orig <- runkmhca(tb_6584m)
round(prop.table(table(kmhca_orig$cluster)), 4)
kmhca_orig_profile <- clusterprev(kmhca_orig)[[1]]
kmhca_orig_labels  <- clusterprev(kmhca_orig)[[2]]




# matrix of PCC comparing two clustering methods
#zz <- calcpcc(lca_orig_profile, lca_orig_labels, mcak_orig_profile, mcak_orig_labels)
#zz




# function to draw b bootstrap samples, run algorithms on each bootstrap sample dataset and calculate
# Pearson correlation (PCC) between each cluster in bootstrap sample and clusters found in original dataset 
# (using same method)
# uses greedy matching algorithm to find cluster in bootstrap sample which has most similar disease profile 
# (highest PCC) with cluster in original dataset

run_bstrap <- function(infile, b) {
  
  lca_outlist   <- list()
  mcak_outlist  <- list()
  kmean_outlist <- list()
  kmhca_outlist <- list()
  
  for(i in 1:b) {
    
    
    # draw bootstrap sample
    boot_infile <- infile[ sample(x=seq_len(nrow(infile)), size=nrow(infile), replace=TRUE), ]
    print(i)
    
    
    # run LCA (method=1) on bootstrap sample
    lca_bs <- runlca(boot_infile)
    lca_bs_profile <- clusterprev(lca_bs)[[1]]
    lca_pcc <- calcpcc(lca_orig_profile, lca_bs_profile)
    lca_pcc$bstrap <- i
    lca_outlist[[i]] <- lca_pcc
    
    
    # run MCA-kmeans (method=2)
    mcak_bs <- runmcak(boot_infile)
    mcak_bs_profile <- clusterprev(mcak_bs)[[1]]
    mcak_pcc <- calcpcc(mcak_orig_profile, mcak_bs_profile)
    mcak_pcc$bstrap <- i
    mcak_outlist[[i]] <- mcak_pcc
    
    
    # run kmeans (method=3)
    kmean_bs <- runkmean(boot_infile)
    kmean_bs_profile <- clusterprev(kmean_bs)[[1]]
    kmean_pcc <- calcpcc(kmean_orig_profile, kmean_bs_profile)
    kmean_pcc$bstrap <- i
    kmean_outlist[[i]] <- kmean_pcc
    
     
    # run kmeans-HCA (method=4)
    kmhca_bs <- runkmhca(boot_infile)
    kmhca_bs_profile <- clusterprev(kmhca_bs)[[1]]
    kmhca_pcc <- calcpcc(kmhca_orig_profile, kmhca_bs_profile)
    kmhca_pcc$bstrap <- i
    kmhca_outlist[[i]] <- kmhca_pcc
     
  }
  
  lca_out   <- as.data.frame(do.call(rbind, lca_outlist))
  mcak_out  <- as.data.frame(do.call(rbind, mcak_outlist))
  kmean_out <- as.data.frame(do.call(rbind, kmean_outlist))
  kmhca_out <- as.data.frame(do.call(rbind, kmhca_outlist))
  
  outlist <- list(lca_out, mcak_out, kmean_out, kmhca_out)
  return(outlist)
}

bootfile <- run_bstrap(tb_6584m,                   400) 
bootfile



save(bootfile, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/bootfile.RData")





load("C:/Users/linda/Documents/UoB/R/THIN/bootfile.RData")

ave_out <- list() 

for (i in 1:4) {
  bootfile[[i]]$PCC <- ifelse(bootfile[[i]]$PCC<0.5, NA, bootfile[[i]]$PCC)
  ave_out[[i]] <- bootfile[[i]] %>%
    group_by(Cluster_orig) %>%
    summarise(avepcc=round(mean(PCC, na.rm=TRUE), 4), 
              sdpcc =round(sd(PCC, na.rm=TRUE), 4),
              pccct =sum(!is.na(PCC)))
}

lca_ave <- ave_out[[1]]
lca_ave_labels <- merge(lca_ave, lca_orig_labels, by.x="Cluster_orig", by.y="Cluster")

mcak_ave <- ave_out[[2]]
mcak_ave_labels <- merge(mcak_ave, mcak_orig_labels, by.x="Cluster_orig", by.y="Cluster")

kmean_ave <- ave_out[[3]]
kmean_ave_labels <- merge(kmean_ave, kmean_orig_labels, by.x="Cluster_orig", by.y="Cluster")

kmhca_ave <- ave_out[[4]]
kmhca_ave_labels <- merge(kmhca_ave, kmhca_orig_labels, by.x="Cluster_orig", by.y="Cluster")



