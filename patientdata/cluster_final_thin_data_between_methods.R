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
library("ggpubr")

options(scipen=999)

# run this because sex column is all zeros in Multimorbid_mm_agesexethn.csv!!
# thinagesex <- read.csv("C:/Users/linda/Documents/UoB/R/THIN/finaldata/Multimorbid_mm_agesex.csv", stringsAsFactors = FALSE)
# thin <- read.csv("C:/Users/linda/Documents/UoB/R/THIN/finaldata/Multimorbid_mm_agesexethn.csv", stringsAsFactors = FALSE)
# 
# thin$SEX <- NULL
# thin <- cbind(thin, SEX=thinagesex$SEX)
# thin <- thin[ , c(100, 1:99) ]
# 
# save(thin, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/thin.RData")



# load thin data
load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/thin.RData")

thin$AGE_CATG <- NULL

# change all column names to lower case
colnames(thin) <- tolower(colnames(thin))

names(thin)[names(thin)=="practice_id"] <- "pracid"
# add patient ID variable and create age group 
thin$patid <- as.numeric( rownames(thin) )
# check minimum age
#min(thin$age_at_index)
thin$ageg <- cut(thin$age_at_index, breaks=c(-Inf, 44, 64, 84, Inf), labels=c("18-44", "45-64", "65-84", ">=85"))
thin <- thin[ , c(2, 100, 1, 5, 101, 3, 4, 6:99)]


# tidy column names
names(thin)[names(thin)=="cprd_pri_lung"]   <- "ca.lung"
names(thin)[names(thin)=="cprd_pri_breast"] <- "ca.breast"
names(thin)[names(thin)=="cprd_pri_bowel"]  <- "ca.bowel"
names(thin)[names(thin)=="cprd_pri_prost"]  <- "ca.prostate"
names(thin)[names(thin)=="tas_lymphoma"]    <- "lymphoma"
names(thin)[names(thin)=="cprd_leukaemia"]  <- "leukaemia"
names(thin)[names(thin)=="cprd_pri_skin"]   <- "ca.skin"
names(thin)[names(thin)=="cprd_pri_melanoma"] <- "melanoma"
names(thin)[names(thin)=="charlsonmetastaticcancer"] <- "metastatic"
names(thin)[names(thin)=="hfincidenceprevkoo"]  <- "hf"
names(thin)[names(thin)=="incidenceprevihd"]    <- "ihd"
names(thin)[names(thin)=="heartvalvedisease"]   <- "valve.dis"
names(thin)[names(thin)=="congenitalheartdseases"] <- "cong.heart.dis"
names(thin)[names(thin)=="pvdprevincidkoo"] <- "pvd"
names(thin)[names(thin)=="aortic_aneurysm"] <- "aortic.an"
names(thin)[names(thin)=="cprd_tia"] <- "tia"
names(thin)[names(thin)=="cprd_isch_stroke"]    <- "isch.stroke"
names(thin)[names(thin)=="haemorrhagicstroke"]  <- "haem.stroke"
names(thin)[names(thin)=="cprd_stroke_nos"]     <- "stroke.nos"
names(thin)[names(thin)=="stroketiaincidprevkoo"] <- "stroke.tia"
names(thin)[names(thin)=="eczemalong2018"]  <- "eczema"
names(thin)[names(thin)=="cprd_psoriasis"]  <- "psoriasis"
names(thin)[names(thin)=="cprd_vitiligo"]   <- "vitiligo"
names(thin)[names(thin)=="cprd_alopecia_areata"] <- "alopecia"
names(thin)[names(thin)=="allergicrhinconjunc"] <- "rhin.conjunc"
names(thin)[names(thin)=="cprd_sinusitis"]  <- "sinusitis"
names(thin)[names(thin)=="cprd_deaf"]       <- "deaf"
names(thin)[names(thin)=="cprd_blind"]      <- "blind"
names(thin)[names(thin)=="cprd_cataract"]   <- "cataract"
names(thin)[names(thin)=="cprd_glaucoma"]   <- "glaucoma"
names(thin)[names(thin)=="amd_alastair"]    <- "amd"
names(thin)[names(thin)=="stretinopathy"]   <- "retinop"
names(thin)[names(thin)=="tas_alluveitis"]  <- "uveitis"
names(thin)[names(thin)=="tas_scleritis_jan2020"]   <- "scleritis"
names(thin)[names(thin)=="cprd_ulcer_peptic"] <- "peptic"
names(thin)[names(thin)=="tas_ibd_jan2020"] <- "ibd"
names(thin)[names(thin)=="irritablebowelsyndrome2019"] <- "ibs"
names(thin)[names(thin)=="chronicliverdisease"]   <- "liver.dis"
names(thin)[names(thin)=="nashandnafl"]           <- "nash.nafl"
names(thin)[names(thin)=="cprd_diverticuli"]      <- "divertic"
names(thin)[names(thin)=="cprd_coeliac"]          <- "coeliac"
names(thin)[names(thin)=="chronicpancreatitis"]   <- "pancreatitis"
names(thin)[names(thin)=="endometriosiskoo"]      <- "endometriosis"
names(thin)[names(thin)=="pcoskoo"] <- "pcos"
names(thin)[names(thin)=="perniciousanaemia"] <- "pern.anaem"
names(thin)[names(thin)=="cprd_vte_ex_pe"] <- "vte"
names(thin)[names(thin)=="pulmonaryembolismjscv1"] <- "pulm.emb"
names(thin)[names(thin)=="depressionrr"]  <- "depression"
names(thin)[names(thin)=="anxietyrr"]     <- "anxiety"
names(thin)[names(thin)=="smirr"]         <- "smi"
names(thin)[names(thin)=="cprd_substance_misuse"] <- "sub.misuse"
names(thin)[names(thin)=="cam_alcoholproblem"]   <- "alcohol"
names(thin)[names(thin)=="cprd_adhd"] <- "adhd"
names(thin)[names(thin)=="cprd_eating_dz"] <- "eating.dis"
names(thin)[names(thin)=="cam_learningdisability"] <- "ld"
names(thin)[names(thin)=="alzheimersdisease"] <- "alzheimers"
names(thin)[names(thin)=="vasculardementiav"] <- "vasc.dementia"
names(thin)[names(thin)=="cprd_parkinsons"]   <- "parkinsons"
names(thin)[names(thin)=="cprd_migraine"]     <- "migraine"
names(thin)[names(thin)=="multiplesclerosis"] <- "ms"
names(thin)[names(thin)=="cprd_epilepsy"] <- "epilepsy"
names(thin)[names(thin)=="charlsonhemiplegia"]          <- "hemiplegia"
names(thin)[names(thin)=="chronicfatiguesyndrome2019"]  <- "cfs"
names(thin)[names(thin)=="fibromyalgia2019"]  <- "fibromyal"
names(thin)[names(thin)=="cprd_oa"] <- "oa"
names(thin)[names(thin)=="cprd_osteoporosis"] <- "osteoporosis"
names(thin)[names(thin)=="polymyalgiarheumaticadraftjscv1"] <- "polymyal"
names(thin)[names(thin)=="cprd_rha"] <- "ra"
names(thin)[names(thin)=="cprd_sjogren"]  <- "sjogren"
names(thin)[names(thin)=="cprd_sle"] <- "sle"
names(thin)[names(thin)=="cprd_sys_sclerosis"] <- "sys.sclerosis"
names(thin)[names(thin)=="psoriaticarthritis"] <- "psor.arthritis"
names(thin)[names(thin)=="cprd_ank_spond"]  <- "ank.spond"
names(thin)[names(thin)=="cprd_gout"]       <- "gout"
names(thin)[names(thin)=="ckdstage3to5"]    <- "ckd"
names(thin)[names(thin)=="asthmalonglist2018"] <- "asthma"
names(thin)[names(thin)=="cprd_copd"] <- "copd"
names(thin)[names(thin)=="osafinal"] <- "osa"
names(thin)[names(thin)=="cprd_bronchiectasis"] <- "bronch"
names(thin)[names(thin)=="cprd_pulm_fibrosis"]  <- "pulm.fib"
names(thin)[names(thin)=="hyperthyroidism"]   <- "hyperthyroid"
names(thin)[names(thin)=="hypothyroidism"]    <- "hypothyroid"
names(thin)[names(thin)=="type1diabetes2018"] <- "t1dm"
names(thin)[names(thin)=="type2diabetes"]     <- "t2dm"
names(thin)[names(thin)=="cprd_hiv"] <- "hiv"
names(thin)[names(thin)=="cprd_bph"] <- "bph"
names(thin)[names(thin)=="erectiledysfunction"] <- "erec.dys"



# go through columns and flag any pairs where one is nested inside the other
#thintest <- thin[ , 8:101]
#for(i in 1:(ncol(thintest)-1)) {
#  for (j in (i+1):ncol(thintest)){
#    di <- which(thintest[ ,i]==1)
#    dj <- which(thintest[, j]==1)
#    if(all(dj %in% di)) cat("All", colnames(thintest)[j], "in", colnames(thintest)[i], "\n")
#    if(all(di %in% dj)) cat("All", colnames(thintest)[i], "in", colnames(thintest)[j], "\n")
#  }   
#}

# drop conditions which are nested within others 
# All haem.stroke in stroke.tia 
# All psor.arthritis in psoriasis 
# All vasc.dementia in dementia
dropcols <- c("psor.arthritis", "vasc.dementia",  "haem.stroke")
thin.work <- thin[ , !(names(thin) %in% dropcols) ]


# look at specific pairs of conditions, how frequently they occur together to see if there is almost overlap
#checkvars <- thin.work[ , c("blind", "cataract", "glaucoma", "amd", "retinop", "alluveitis", "scleritis") ]
#checkvars <- thin.work[, c("arrhythmia", "af", "hypertension", "hf", "ihd", "valve.dis", "cardiomyopathy", "cong.heart.dis") ]
#checkvars <- thin.work[ , c("nash.nafl", "liver.dis") ]
#checkvars <- thin.work[ , c("alzheimers", "dementia", "parkinsons") ]
#checkvars <- thin.work[ , c("polymyal", "ra") ]
#checkvars <- thin.work[ , c("tia", "isch.stroke", "stroke.nos", "stroke.tia") ]

#for(i in 1:(ncol(checkvars)-1)) {
#  for (j in (i+1):ncol(checkvars)){
#    print(i)
#    print(j)
#    print(table(checkvars[ , i] , checkvars[ , j]))
#  }   
#}
table(thin.work$arrhythmia, thin.work$af)


# make combined arrhythmia/AF variable
thin.work$arr.af <- ifelse(rowSums(thin.work[ , c("arrhythmia", "af")])>0, 1 , 0)
dropcols <- c("arrhythmia", "af")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# make stroke/TIA variable
thin.work$str.tia <- ifelse(rowSums(thin.work[ , c("tia", "isch.stroke", "stroke.nos", "stroke.tia")])>0, 1 , 0)
dropcols <- c("tia", "isch.stroke", "stroke.nos", "stroke.tia")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# make all cancer variable
thin.work$ca.all <- ifelse(rowSums(thin.work[ , c("cancerall", "ca.lung", "ca.breast", "ca.bowel", "ca.prostate",
                                                  "lymphoma", "leukaemia", "ca.skin", "melanoma", "metastatic")])>0, 1 , 0)
dropcols <- c("cancerall", "ca.lung", "ca.breast", "ca.bowel", "ca.prostate", "lymphoma", "leukaemia", "ca.skin", 
              "melanoma", "metastatic")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]# keep only ca.all variable


# create combined categories for those conditions with prevalence < 1% (in males 65-84y). Using categories specificed by Ing
# neurodevelopmental disorder
thin.work$neurodev <- ifelse(rowSums(thin.work[ , c("adhd", "ld")])>0, 1 , 0)
dropcols <- c("adhd", "ld")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# somatoform disorder
thin.work$somat <- ifelse(rowSums(thin.work[ , c("cfs", "fibromyal")])>0, 1 , 0)
dropcols <- c("cfs", "fibromyal")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

#chronic liver disease
thin.work$chron.liver <- ifelse(rowSums(thin.work[ , c("liver.dis", "nash.nafl")])>0, 1 , 0)
dropcols <- c("liver.dis", "nash.nafl")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# dementia
thin.work$alz.dem <- ifelse(rowSums(thin.work[ , c("alzheimers", "dementia")])>0, 1 , 0)
dropcols <- c("alzheimers", "dementia")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# other heart disease
thin.work$other.heart <- ifelse(rowSums(thin.work[ , c("cardiomyopathy", "cong.heart.dis")])>0, 1 , 0)
dropcols <- c("cardiomyopathy", "cong.heart.dis")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# addiction
thin.work$addict <- ifelse(rowSums(thin.work[ , c("alcohol", "sub.misuse")])>0, 1 , 0)
dropcols <- c("alcohol", "sub.misuse")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# autoimmune bowel disease
thin.work$auto.bowel <- ifelse(rowSums(thin.work[ , c("ibd", "coeliac", "pern.anaem")])>0, 1 , 0)
dropcols <- c("ibd", "coeliac", "pern.anaem")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# autoimmune connective tissue disorders and related disease
thin.work$auto.tissue <- ifelse(rowSums(thin.work[ , c("ra", "sle", "sys.sclerosis", "sjogren", "ank.spond", "polymyal")])>0, 1 , 0)
dropcols <- c("ra", "sle", "sys.sclerosis", "sjogren", "ank.spond", "polymyal")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# other autoimmune
thin.work$auto.other <- ifelse(rowSums(thin.work[ , c("alopecia", "vitiligo", "uveitis", "scleritis", "ms")])>0, 1 , 0)
dropcols <- c("alopecia", "vitiligo", "uveitis", "scleritis", "ms")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# diabetes
thin.work$diabetes <- ifelse(rowSums(thin.work[ , c("t1dm", "t2dm")])>0, 1 , 0)
dropcols <- c("t1dm", "t2dm")
thin.work <- thin.work[ , !(names(thin.work) %in% dropcols) ]

# drop retinopathy variable as lookup needs reviewing
thin.work$retinop <- NULL


# choose random sample of 50 practices for analysis

set.seed(20200828)

pracs <- as.vector(unique(thin.work$pracid))
pracs.sample <- sample( pracs, 50, replace=FALSE)

thin.sample <- thin.work[ thin.work$pracid %in% pracs.sample, ]
thin.sample <- thin.sample[ order(thin.sample$pracid, thin.sample$patid), ]

# save sample of thin data
save(thin.sample, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/thinsample.RData")

### end of data prep







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


# create variable for count of num of morbidities
tb_6584m$ctmorb <- rowSums(tb_6584m[ , allconds] )


# check if any morbidity is NOT present in the data
colSums(tb_6584m[ , allconds]==1) 


# prevalence of each disease in whole dataset
colcount <- data.frame(colSums(tb_6584m[ , allconds]==1))
colcount[2] <- round( colcount[1]/nrow(tb_6584m), 5)
names(colcount) <- c("TotalFreq","TotalProp")
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
prev_6584m <- prev_6584m[ , c(1, 4, 3, 2)]
prev_6584m <- prev_6584m[order(-prev_6584m$TotalProp), ]




set.seed(20200908)



if(FALSE){
  ### run LCA finding AIC, BIC to get number of clusters
  infile <- tb_6584m
  
  infile[colmorbidlist] <- infile[colmorbidlist] + 1
  f <- with(infile, cbind(ca.all, str.tia, arr.af, hypertension, hf, ihd, valve.dis, pvd, aortic.an, eczema, 
                          psoriasis, rhin.conjunc, sinusitis, deaf, blind, cataract, glaucoma, amd, peptic, ibs,
                          divertic, vte, pulm.emb, depression, anxiety, smi, parkinsons, migraine, epilepsy, oa, 
                          osteoporosis, gout, ckd, asthma, copd, osa, bronch, hyperthyroid, hypothyroid, bph, 
                          erec.dys, chron.liver, alz.dem, other.heart, addict, auto.bowel, auto.tissue, auto.other,
                          diabetes) ~ 1)

  for(i in 1:12) { # note class=1 -> log-linear independence model
    print(i)
    temp <- poLCA(f, infile[ , colmorbidlist], nclass=i, maxiter=5000, tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE) 
    assign(paste0("lc_mm", i), temp)
  }
  
  #   Step 2: compute a set of model-selection statistics and choose the best model
  results_mm <- data.frame(Modell=1,
                           log_likelihood=lc_mm1$llik,
                           df = lc_mm1$resid.df,
                           BIC=lc_mm1$bic,
                           SABIC=  (-2*lc_mm1$llik) + (log((lc_mm1$N + 2)/24)) * (lc_mm1$bic-(-2*lc_mm1$llik))/log(lc_mm1$N),
                           likelihood_ratio=lc_mm1$Gsq)
  
  #   results_mm$R2_entropy
  results_mm$Modell<-as.integer(results_mm$Modell)
  results_mm[1,7]<-c("-")
  
  #   build lists of latent class objects
  my.lcalist <- lapply(paste('lc_mm', seq(1, 12, 1), sep=''), get)
  
  # 	for all results, i in 1:12
  for(i in 2:12) {
    results_mm[i, 1] <- i
    results_mm[i, 2] <- my.lcalist[[i]]$llik
    results_mm[i, 3] <- my.lcalist[[i]]$resid.df
    results_mm[i, 4] <- my.lcalist[[i]]$bic
    results_mm[i, 5] <- (-2*my.lcalist[[i]]$llik) + (log((my.lcalist[[i]]$N + 2)/24)) * (my.lcalist[[i]]$bic-(-2*my.lcalist[[i]]$llik))/log(my.lcalist[[i]]$N) #SAbic
    results_mm[i, 6] <- my.lcalist[[i]]$Gsq  
    error_prior <- entropy(my.lcalist[[i]]$P) # class proportions model i
    error_post <- mean(apply(my.lcalist[[i]]$posterior, 1, entropy), na.rm = TRUE)
    results_mm[i, 7] <- round(((error_prior - error_post) / error_prior), 3)
  }
  
  # 	combining results to a dataframe
  colnames(results_mm) <- c("Model", "log-likelihood", "resid. df", "BIC", "SABIC", "likelihood-ratio", "Entropy")
  lca_results_mm <- results_mm
  
  # plot BIC, SABIC, entropy for each number of classes tried
  
  # Order categories of results$model in order of appearance
  results_mm$Model <- as.factor(results_mm$Model) 
  results_mm$Entropy <- as.numeric(results_mm$Entropy)
  
  #convert to long format
  results2 <- tidyr::gather(results_mm, Meas, Value, 4:7)
  results2
  
  #plot
  
  fit.plot <- ggplot(results2) + 
    geom_point(aes(x=Model, y=Value), size=3) +
    geom_line(aes(Model, Value, group = 1)) +
    theme_bw() +
    labs(x = "", y="", title = "") + 
    facet_grid(Meas ~. , scales = "free") +
    theme_bw(base_size = 16, base_family = "") +   
    theme(panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(colour="grey", size=0.5),
          legend.title = element_text(size = 16, face = 'bold'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text=  element_text(size=16),
          axis.line = element_line(colour = "black"))
  jpeg(filename="C:/Users/linda/Documents/UoB/R/THIN/finaldata/LCA_6584m.jpg")
  plot(fit.plot)
  dev.off()
  
  fit.plot
}

#  most appropriate number of classes = 4 (entropy is increased)




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



# function to calculate frequency of disease in each cluster. Needs to have disease coded as 0/1 where 1=presence
freqcalc <- function(inputfile) {
  freqout <- data.frame( apply(inputfile[ , colmorbidlist], 2, FUN=function(z) round( tapply(z , inputfile$cluster, FUN=sum), 5) ) )
  freqout$cluster <- as.numeric( rownames(freqout) )
  freqout <- freqout[ , c(50, 1:49)]
  return(freqout)
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



# function to calculate frequency and prevalence of disease in each cluster, exclusivity and O/E ratio
clusterfreqprev <- function(infile) {

  # frequency of disease within each cluster
  df_orig_cluster_freq_prev <- freqcalc(infile)
  
  # find total number of patients in each cluster
  clustn <- infile %>% 
    group_by(Cluster=cluster) %>%
    summarise(clustn = length(cluster))
  
  orig_cluster_freq_prev <- melt(df_orig_cluster_freq_prev, id.vars="cluster")
  colnames(orig_cluster_freq_prev) <- c("Cluster", "Morbidity", "ClustFreq")
  orig_cluster_total_freq_prev <- merge(orig_cluster_freq_prev, prev_6584m, by="Morbidity", all.x=TRUE)
  orig_cluster_total_freq_prev <- merge(orig_cluster_total_freq_prev, clustn, by="Cluster", all.x=TRUE)
  orig_cluster_total_freq_prev <- orig_cluster_total_freq_prev[ , c(1, 4, 2, 3, 7, 6, 5)]

  # prevalence of disease within each cluster
  orig_cluster_total_freq_prev$ClustProp <- round(orig_cluster_total_freq_prev$ClustFreq/orig_cluster_total_freq_prev$clustn, 4)
  
  # calculate exclusivity and O/E ratio
  orig_cluster_total_freq_prev$exclusivity <- round((orig_cluster_total_freq_prev$ClustFreq/orig_cluster_total_freq_prev$TotalFreq) * 100, 2) 
  orig_cluster_total_freq_prev$oeratio     <- round(orig_cluster_total_freq_prev$ClustProp/orig_cluster_total_freq_prev$TotalProp, 2)
  
  return(orig_cluster_total_freq_prev)

}



# function to calculate PCC matrix between disease profiles in each cluster from two clustering algorithms.  
# finds most similar (highest PCC>=0.5) cluster between the two methods. If PCC<0.5  implies there is no similar cluster  
calcpcc <- function(alg1, labels1, alg2, labels2) {
  
  pccmat <- round( cor(alg1, alg2), 4)
  
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
  matchout <- cbind(rowInd, colInd, pccvec)
  colnames(matchout) <- c("Cluster_Alg1", "Cluster_Alg2", "PCC")
  matchout <- merge(matchout, labels1, by.x="Cluster_Alg1", by.y="Cluster")
  matchout <- merge(matchout, labels2, by.x="Cluster_Alg2", by.y="Cluster")
  colnames(matchout)[ colnames(matchout) %in% c("top3morb_3.x", "top3morb_3.y") ] <- c("Labels_Alg1", "Labels_Alg2")
  matchout <- matchout[ , c(2, 4, 1, 5, 3)]
  matchout$PCC <- ifelse(matchout$PCC<0.5, NA, matchout$PCC) 
  matchout <- matchout[order(matchout$Cluster_Alg1), ]
  return(matchout)
}


# function to calculate summary statistics for each cluster: median age, median number of conditions, distribution of smoking, ethnicity
cluststat <- function(inputfile) {
  med.num.morb <- tapply(inputfile$ctmorb, inputfile$cluster, FUN=median)
  med.age  <- tapply(inputfile$age_at_index, inputfile$cluster, FUN=function(x) round(median(x),2))
  eth.dist <- tapply(inputfile$ethnicity, inputfile$cluster, FUN=function(x) round( prop.table(table(x)) *100, 2 ) )
  townsend.dist <- tapply(inputfile$townsend, inputfile$cluster, FUN=function(x) round( prop.table(table(x)) *100, 2 ) )
  outlist <- list( med.num.morb, med.age, eth.dist, townsend.dist)
  return(outlist)
}


Kx = 4 # 4 clusters looks most appropriate from LCA entropy results



### Latent class analysis
lca_orig <- runlca(tb_6584m)
round(prop.table(table(lca_orig$cluster)), 4)
cluststat(lca_orig)
lca_orig_profile <- clusterprev(lca_orig)[[1]]
lca_orig_labels  <- clusterprev(lca_orig)[[2]]


# bubble plots
lca_exc_oe <- clusterfreqprev(lca_orig)
lca_exc_oe_melt <- lca_exc_oe[ , c(1:2, 9:10)]
lca_exc_oe_melt <- melt(lca_exc_oe_melt, id = c("Morb", "Cluster"))

# exclusivity                        
lca_exc <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="exclusivity", ]                      
lca_exc <- lca_exc[ , c(1, 2, 4)]
ggplot(lca_exc, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) + 
  labs(x="Cluster", y="Morbidity", title="Latent class analysis: exclusivity") +
  theme_bw() +
  scale_size_area(max_size=5, name="Exclusivity %") 
dev.copy(jpeg, 'bubble_lca_exc.jpg')
dev.off()

# O/E ratio
lca_oeratio <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="oeratio", ]                      
lca_oeratio <- lca_oeratio[ , c(1, 2, 4)]
ggplot(lca_oeratio, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) +
  labs(x="Cluster", y="Morbidity", title="Latent class analysis: O/E ratio") +
  theme_bw() +
  scale_size_area(max_size=10, name="O/E ratio")
dev.copy(jpeg, 'bubble_lca_oeratio.jpg')
dev.off()



# try running LCA with varying number of clusters
if(FALSE){
  Kx = 5
  lca_orig5 <- runlca(tb_6584m)
  table(lca_orig5$cluster)
  round(prop.table(table(lca_orig5$cluster)), 4)
  cluststat(lca_orig5)
  lca_orig_profile5 <- clusterprev(lca_orig5)[[1]]
  lca_orig_labels5  <- clusterprev(lca_orig5)[[2]]

  # bubble plots
  lca_exc_oe <- clusterfreqprev(lca_orig5)
  lca_exc_oe_melt <- lca_exc_oe[ , c(1:2, 9:10)]
  lca_exc_oe_melt <- melt(lca_exc_oe_melt, id = c("Morb", "Cluster"))
  
  # exclusivity                        
  lca_exc <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="exclusivity", ]                      
  lca_exc <- lca_exc[ , c(1, 2, 4)]
  ggplot(lca_exc, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) + 
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: exclusivity") +
    theme_bw() +
    scale_size_area(max_size=5, name="Exclusivity %") 
  dev.copy(jpeg, 'bubble_lca5_exc.jpg')
  dev.off()
  
  # O/E ratio
  lca_oeratio <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="oeratio", ]                      
  lca_oeratio <- lca_oeratio[ , c(1, 2, 4)]
  ggplot(lca_oeratio, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) +
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: O/E ratio") +
    theme_bw() +
    scale_size_area(max_size=10, name="O/E ratio")
  dev.copy(jpeg, 'bubble_lca5_oeratio.jpg')
  dev.off()
  
  
  
  Kx = 6
  lca_orig6 <- runlca(tb_6584m)
  table(lca_orig6$cluster)
  round(prop.table(table(lca_orig6$cluster)), 4)
  cluststat(lca_orig6)
  lca_orig_profile6 <- clusterprev(lca_orig6)[[1]]
  lca_orig_labels6  <- clusterprev(lca_orig6)[[2]]

  # bubble plots
  lca_exc_oe <- clusterfreqprev(lca_orig6)
  lca_exc_oe_melt <- lca_exc_oe[ , c(1:2, 9:10)]
  lca_exc_oe_melt <- melt(lca_exc_oe_melt, id = c("Morb", "Cluster"))
  
  # exclusivity                        
  lca_exc <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="exclusivity", ]                      
  lca_exc <- lca_exc[ , c(1, 2, 4)]
  ggplot(lca_exc, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) + 
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: exclusivity") +
    theme_bw() +
    scale_size_area(max_size=5, name="Exclusivity %") 
  dev.copy(jpeg, 'bubble_lca6_exc.jpg')
  dev.off()
  
  # O/E ratio
  lca_oeratio <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="oeratio", ]                      
  lca_oeratio <- lca_oeratio[ , c(1, 2, 4)]
  ggplot(lca_oeratio, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) +
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: O/E ratio") +
    theme_bw() +
    scale_size_area(max_size=10, name="O/E ratio")
  dev.copy(jpeg, 'bubble_lca6_oeratio.jpg')
  dev.off()
  
  
  
  Kx = 7
  lca_orig7 <- runlca(tb_6584m)
  table(lca_orig7$cluster)
  round(prop.table(table(lca_orig7$cluster)), 4)
  cluststat(lca_orig7)
  lca_orig_profile7 <- clusterprev(lca_orig7)[[1]]
  lca_orig_labels7  <- clusterprev(lca_orig7)[[2]]
  
  # bubble plots
  lca_exc_oe <- clusterfreqprev(lca_orig7)
  lca_exc_oe_melt <- lca_exc_oe[ , c(1:2, 9:10)]
  lca_exc_oe_melt <- melt(lca_exc_oe_melt, id = c("Morb", "Cluster"))
  
  # exclusivity                        
  lca_exc <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="exclusivity", ]                      
  lca_exc <- lca_exc[ , c(1, 2, 4)]
  ggplot(lca_exc, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) + 
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: exclusivity") +
    theme_bw() +
    scale_size_area(max_size=5, name="Exclusivity %") 
  dev.copy(jpeg, 'bubble_lca7_exc.jpg')
  dev.off()
  
  # O/E ratio
  lca_oeratio <- lca_exc_oe_melt[ lca_exc_oe_melt$variable=="oeratio", ]                      
  lca_oeratio <- lca_oeratio[ , c(1, 2, 4)]
  ggplot(lca_oeratio, aes(x = Cluster, y = Morb, size=value)) +
    geom_point(alpha=0.7) +
    labs(x="Cluster", y="Morbidity", title="Latent class analysis: O/E ratio") +
    theme_bw() +
    scale_size_area(max_size=10, name="O/E ratio")
  dev.copy(jpeg, 'bubble_lca7_oeratio.jpg')
  dev.off()
  
  

  # add cluster allocation from LCA with differing number of classes into original data
  names(lca_orig)[names(lca_orig)=="cluster"] <- "lca4"
  lca_orig <- merge(lca_orig, lca_orig5[, c("pracid", "patid", "cluster")], by=c("pracid", "patid"))
  names(lca_orig)[names(lca_orig)=="cluster"] <- "lca5"
  lca_orig <- merge(lca_orig, lca_orig6[, c("pracid", "patid", "cluster")], by=c("pracid", "patid"))
  names(lca_orig)[names(lca_orig)=="cluster"] <- "lca6"
  lca_orig <- merge(lca_orig, lca_orig7[, c("pracid", "patid", "cluster")], by=c("pracid", "patid"))
  names(lca_orig)[names(lca_orig)=="cluster"] <- "lca7"
  
  # investigate which classes each patient goes into when differing number of classes used
  table(lca_orig$lca4, lca_orig$lca5)
  table(lca_orig$lca4, lca_orig$lca6)
  table(lca_orig$lca4, lca_orig$lca7)
}




### MCA-kmeans
mcak_orig <- runmcak(tb_6584m)
round(prop.table(table(mcak_orig$cluster)), 4)
cluststat(mcak_orig)
mcak_orig_profile <- clusterprev(mcak_orig)[[1]]
mcak_orig_labels  <- clusterprev(mcak_orig)[[2]]


# bubble plots
mcak_exc_oe <- clusterfreqprev(mcak_orig)
mcak_exc_oe_melt <- mcak_exc_oe[ , c(1:2, 9:10)]
mcak_exc_oe_melt <- melt(mcak_exc_oe_melt, id = c("Morb", "Cluster"))

# exclusivity                        
mcak_exc <- mcak_exc_oe_melt[ mcak_exc_oe_melt$variable=="exclusivity", ]                      
mcak_exc <- mcak_exc[ , c(1, 2, 4)]
ggplot(mcak_exc, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) + 
  labs(x="Cluster", y="Morbidity", title="MCA-kmeans: exclusivity") +
  theme_bw() +
  scale_size_area(max_size=5, name="Exclusivity %") 
dev.copy(jpeg, 'bubble_mcak_exc.jpg')
dev.off()

# O/E ratio
mcak_oeratio <- mcak_exc_oe_melt[ mcak_exc_oe_melt$variable=="oeratio", ]                      
mcak_oeratio <- mcak_oeratio[ , c(1, 2, 4)]
ggplot(mcak_oeratio, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) +
  labs(x="Cluster", y="Morbidity", title="MCA-kmeans: O/E ratio") +
  theme_bw() +
  scale_size_area(max_size=10, name="O/E ratio")
dev.copy(jpeg, 'bubble_mcak_oeratio.jpg')
dev.off()




### kmeans
kmean_orig <- runkmean(tb_6584m)
round(prop.table(table(kmean_orig$cluster)), 4)
cluststat(kmean_orig)
kmean_orig_profile <- clusterprev(kmean_orig)[[1]]
kmean_orig_labels  <- clusterprev(kmean_orig)[[2]]


# bubble plots
kmean_exc_oe <- clusterfreqprev(kmean_orig)
kmean_exc_oe_melt <- kmean_exc_oe[ , c(1:2, 9:10)]
kmean_exc_oe_melt <- melt(kmean_exc_oe_melt, id = c("Morb", "Cluster"))

# exclusivity                        
kmean_exc <- kmean_exc_oe_melt[ kmean_exc_oe_melt$variable=="exclusivity", ]                      
kmean_exc <- kmean_exc[ , c(1, 2, 4)]
ggplot(kmean_exc, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) + 
  labs(x="Cluster", y="Morbidity", title="Kmeans: exclusivity") +
  theme_bw() +
  scale_size_area(max_size=5, name="Exclusivity %") 
dev.copy(jpeg, 'bubble_kmean_exc.jpg')
dev.off()

# O/E ratio
kmean_oeratio <- kmean_exc_oe_melt[ kmean_exc_oe_melt$variable=="oeratio", ]                      
kmean_oeratio <- kmean_oeratio[ , c(1, 2, 4)]
ggplot(kmean_oeratio, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) +
  labs(x="Cluster", y="Morbidity", title="Kmeans: O/E ratio") +
  theme_bw() +
  scale_size_area(max_size=10, name="O/E ratio")
dev.copy(jpeg, 'bubble_kmean_oeratio.jpg')
dev.off()




### Kmeans-HCA
kmhca_orig <- runkmhca(tb_6584m)
round(prop.table(table(kmhca_orig$cluster)), 4)
cluststat(kmhca_orig)
kmhca_orig_profile <- clusterprev(kmhca_orig)[[1]]
kmhca_orig_labels  <- clusterprev(kmhca_orig)[[2]]

#kmn_mm <- kmeans(tb_6584m[ , colmorbidlist], centers=Kx, nstart=150, iter.max=50)
#kmn_mm$withinss


# bubble plots
kmhca_exc_oe <- clusterfreqprev(kmhca_orig)
kmhca_exc_oe_melt <- kmhca_exc_oe[ , c(1:2, 9:10)]
kmhca_exc_oe_melt <- melt(kmhca_exc_oe_melt, id = c("Morb", "Cluster"))

# exclusivity                        
kmhca_exc <- kmhca_exc_oe_melt[ kmhca_exc_oe_melt$variable=="exclusivity", ]                      
kmhca_exc <- kmhca_exc[ , c(1, 2, 4)]
ggplot(kmhca_exc, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) + 
  labs(x="Cluster", y="Morbidity", title="Kmeans-HCA: exclusivity") +
  theme_bw() +
  scale_size_area(max_size=5, name="Exclusivity %") 
dev.copy(jpeg, 'bubble_kmhca_exc.jpg')
dev.off()

# O/E ratio
kmhca_oeratio <- kmhca_exc_oe_melt[ kmhca_exc_oe_melt$variable=="oeratio", ]                      
kmhca_oeratio <- kmhca_oeratio[ , c(1, 2, 4)]
ggplot(kmhca_oeratio, aes(x = Cluster, y = Morb, size=value)) +
  geom_point(alpha=0.7) +
  labs(x="Cluster", y="Morbidity", title="Kmeans-HCA: O/E ratio") +
  theme_bw() +
  scale_size_area(max_size=10, name="O/E ratio")
dev.copy(jpeg, 'bubble_kmhca_oeratio.jpg')
dev.off()





# add cluster allocations by each method into one data file
all_clust <- lca_orig
names(all_clust)[names(all_clust)=="cluster"]   <- "LCA"
all_clust <- merge(all_clust, mcak_orig[ , c(1:2,66)], by=c("pracid", "patid"))
names(all_clust)[names(all_clust)=="cluster"]   <- "MCAk"
all_clust <- merge(all_clust, kmean_orig[ , c(1:2,66)], by=c("pracid", "patid"))
names(all_clust)[names(all_clust)=="cluster"]   <- "kmean"
all_clust <- merge(all_clust, kmhca_orig[ , c(1:2,66)], by=c("pracid", "patid"))
names(all_clust)[names(all_clust)=="cluster"]   <- "kmHCA"




clusterlabels <- unique(all_clust$kmHCA)
withinssvector <- vector(mode="numeric", length=length(clusterlabels) * 2)
count <- 1

for (k in clusterlabels) {
  currentClusterData <- all_clust[all_clust$kmHCA==k, colmorbidlist]
  withinssvector[count] <- sum(as.matrix(dist(currentClusterData)^2)) / (2 * nrow(currentClusterData))
  withinssvector[count + 4] <- withinssvector[count] / nrow(currentClusterData)
  count <- count + 1
}

clusterlabels
withinssvector



# adjustedRandIndex(all_clust$LCA,  all_clust$MCAk)
# table(all_clust$LCA,  all_clust$MCAk)
# table(all_clust$LCA,  all_clust$kmean)
# table(all_clust$LCA,  all_clust$kmHCA)
# table(all_clust$MCAk, all_clust$kmean)
# table(all_clust$MCAk, all_clust$kmHCA)
# table(all_clust$kmean,all_clust$kmHCA)


# Jaccard similarity function
jaccard <- function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

# table(all_clust$LCA)
# lca_c1 <- all_clust[ all_clust$LCA==1, colmorbidlist ]
# lca_c1 <- lca_c1[ 1:250, ]

# jaccmat <- matrix(, nrow=nrow(lca_c1), ncol=nrow(lca_c1))
# 
# system.time(
# for(i in 1:(nrow(lca_c1)-1)) {
#     for (j in (i+1):nrow(lca_c1)){
#       jaccmat[ i, j] <- jaccard( lca_c1[ i , ], lca_c1[ j, ] )
#   }   
# }
# )
# 
# mean(jaccmat, na.rm=TRUE)


#system.time(jaccmat3 <- apply(combn(nrow(lca_c1), 2), 2, function(x) jaccard(lca_c1[ x[1], ], lca_c1[ x[2], ] )) )


aim_ss <- 1000

jacc_clust <- function(alg, m) {
  
  jaccmat <- matrix( , nrow=aim_ss, ncol=aim_ss)
  hammmat <- matrix( , nrow=aim_ss, ncol=aim_ss)
  outlist <- vector(mode="list", length=m*6)
  
  for (i in 1:m) {
    
    clus <- all_clust[all_clust[[alg]]==i, ]
    
    jacc_ss <- ifelse(nrow(clus)<aim_ss, nrow(clus), aim_ss)
    
    ind  <- sample( 1:nrow(clus), jacc_ss, replace=FALSE) 
    tmp <- clus[ ind, ]

    for(j in 1:(nrow(tmp)-1)) {
      for (k in (j+1):nrow(tmp)){
        jaccmat[ j, k] <- jaccard( tmp[ j , colmorbidlist ], tmp[ k, colmorbidlist ] )
        hammmat[ j, k] <- sum( tmp[ j , colmorbidlist ] != tmp[ k , colmorbidlist ] ) 
      }   
    }
    
    outlist[[i]] <- jaccmat
    outlist[[i+m]] <- tmp
    outlist[[i+(2*m)]] <- mean(jaccmat, na.rm=TRUE)
    outlist[[i+(3*m)]] <- hammmat
    outlist[[i+(4*m)]] <- mean(hammmat, na.rm=TRUE)
    outlist[[i+(5*m)]] <- median(hammmat, na.rm=TRUE)
  }
 
  return(outlist)
  
}

system.time(lca_jacc   <- jacc_clust("LCA", 4))
system.time(mcak_jacc  <- jacc_clust("MCAk", 4))
system.time(kmean_jacc <- jacc_clust("kmean", 4))
system.time(kmhca_jacc <- jacc_clust("kmHCA", 4))


save(lca_jacc, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/lca_jacc.RData")
save(mcak_jacc, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/mcak_jacc.RData")
save(kmean_jacc, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/kmean_jacc.RData")
save(kmhca_jacc, file="C:/Users/linda/Documents/UoB/R/THIN/finaldata/kmhca_jacc.RData")


load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/lca_jacc.RData")
load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/mcak_jacc.RData")
load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/kmean_jacc.RData")
load("C:/Users/linda/Documents/UoB/R/THIN/finaldata/kmhca_jacc.RData")


# other summary stats for Jaccard and Hamming distance
x1 <- kmean_jacc[[4]]
x1v <- c(x1)
sd(x1v, na.rm=TRUE)/sqrt(length(na.omit(x1v)))

y1 <- kmean_jacc[[16]]
y1v <- c(y1)
sd(y1v, na.rm=TRUE)/sqrt(length(na.omit(y1v)))

quantile(y1v, na.rm=TRUE, probs=c(0.25, 0.5, 0.75))





# matrix of PCC comparing two clustering methods
zz <- calcpcc(kmean_orig_profile, kmean_orig_labels, kmhca_orig_profile, kmhca_orig_labels)
zz




# bar chart: prevalence of each condition in full dataset
#ggplot(prev_6584) +
#  geom_bar(aes(x = Morb, y = TotalProp), stat="identity") +
#  theme(axis.text.x = element_text(angle=90, size=10),
#        axis.ticks  = element_blank()) +
#  ggtitle("Prevalence of conditions: age 65-84y") +
#  xlab("Condition") + ylab("Proportion with disease")



# bar chart: prevalence of each condition in LCA cluster 1 (diabetes-ED-hyperten) compared with overall data
lca_orig_profilex <- tibble::rownames_to_column(data.frame(lca_orig_profile), "Morbidity")
lca_orig_profilex <- merge(lca_orig_profilex, prev_6584m, by="Morbidity")
lca_cluster_prev  <- melt(lca_orig_profilex, id.vars=c("Morbidity", "Morb"))
colnames(lca_cluster_prev) <- c("Morbidity", "Morb", "Cluster", "Prop")

lca_cluster_x1 <- lca_cluster_prev[ lca_cluster_prev$Cluster=="X4" | lca_cluster_prev$Cluster=="TotalProp" , ]
bar.colors <- c(X4 = "red", TotalProp = "black")
ggplot(lca_cluster_x1) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("LCA cluster 4 (PVD-aortic an-IHD)", "Full dataset")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("Prevalence of conditions: LCA cluster 4 (PVD-aortic an-IHD) compared to full dataset") +
  xlab("Condition") + ylab("Proportion with disease")



# bar chart: prevalence of each condition in MCA-kmenas cluster 3 (IHD-AF-COPD) compared with overall data
mcak_orig_profilex <- tibble::rownames_to_column(data.frame(mcak_orig_profile), "Morbidity")
mcak_orig_profilex <- merge(mcak_orig_profilex, prev_6584m, by="Morbidity")
mcak_cluster_prev  <- melt(mcak_orig_profilex, id.vars=c("Morbidity", "Morb"))
colnames(mcak_cluster_prev) <- c("Morbidity", "Morb", "Cluster", "Prop")

mcak_cluster_x3 <- mcak_cluster_prev[ mcak_cluster_prev$Cluster=="X4" | mcak_cluster_prev$Cluster=="TotalProp" , ]
bar.colors <- c(X4 = "red", TotalProp = "black")
ggplot(mcak_cluster_x3) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("MCA-kmeans cluster 4 (HF-IHD-AF)", "Full dataset")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("Prevalence of conditions: MCA-kmeans cluster 4 (HF-IHD-AF) compared to full dataset") +
  xlab("Condition") + ylab("Proportion with disease")




# bar chart: prevalence of each condition in LCA cluster 1 (Diabetes-ED-hyperten) compared with each cluster from kmeans
lca_orig_profilex  <- tibble::rownames_to_column(data.frame(lca_orig_profile), "Morbidity")
kmean_orig_profilex <- tibble::rownames_to_column(data.frame(kmean_orig_profile), "Morbidity")
lca_orig_profilex  <- lca_orig_profilex[ , c(1,2)] # select cluster 1 (coulmn 2) from LCA profile
colnames(lca_orig_profilex) <- c("Morbidity", "LCA_1")
lca_kmean_cluster_prev <- merge(lca_orig_profilex, kmean_orig_profilex, by="Morbidity")
lca_kmean_cluster_prev <- merge(lca_kmean_cluster_prev, shortcode, by="Morbidity")
lca_kmean_cluster_prev <- melt(lca_kmean_cluster_prev, id.vars=c("Morbidity", "Morb"))
colnames(lca_kmean_cluster_prev) <- c("Morbidity", "Morb", "Cluster", "Prop")

# kmeans cluster 1
lca_kmean_cluster_x1 <- lca_kmean_cluster_prev[ lca_kmean_cluster_prev$Cluster=="X1" | lca_kmean_cluster_prev$Cluster=="LCA_1" , ]
bar.colors <- c(X1 = "red", LCA_1 = "black")
plotx1 <- ggplot(lca_kmean_cluster_x1) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("LCA cluster 1", "kmeans cluster 1")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("kmeans cluster 1 (hypertension-CKD-gout), PCC=0.51") +
  ylab("Proportion with disease")
plotx1

# kmeans cluster 2
lca_kmean_cluster_x2 <- lca_kmean_cluster_prev[ lca_kmean_cluster_prev$Cluster=="X2" | lca_kmean_cluster_prev$Cluster=="LCA_1" , ]
bar.colors <- c(X2 = "red", LCA_1 = "black")
plotx2 <- ggplot(lca_kmean_cluster_x2) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("LCA cluster 1", "kmeans cluster 2")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("kmeans cluster 2 (Diabetes-ED-hypertension), PCC=1.00")
plotx2

# kmeans cluster 3
lca_kmean_cluster_x3 <- lca_kmean_cluster_prev[ lca_kmean_cluster_prev$Cluster=="X3" | lca_kmean_cluster_prev$Cluster=="LCA_1" , ]
bar.colors <- c(X3 = "red", LCA_1 = "black")
plotx3 <- ggplot(lca_kmean_cluster_x3) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("LCA cluster 1", "kmeans cluster 3")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("kmeans cluster 3 (OA-BPH-deaf), PCC=0.39")
plotx3

# kmeans cluster 4
lca_kmean_cluster_x4 <- lca_kmean_cluster_prev[ lca_kmean_cluster_prev$Cluster=="X4" | lca_kmean_cluster_prev$Cluster=="LCA_1" , ]
bar.colors <- c(X4 = "red", LCA_1 = "black")
plotx4 <- ggplot(lca_kmean_cluster_x4) +
  geom_bar(aes(x = Morb, y = Prop, fill=Cluster), position="dodge", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=10, vjust=0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=11, face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.85)) +
  scale_fill_manual(values=bar.colors, labels = c("LCA cluster 1", "kmeans cluster 4")) +
  scale_y_continuous(breaks=seq(0, 1,, by = 0.25), limits=c(0,1)) +
  ggtitle("kmeans cluster 4 (Cancer-depression-COPD), PCC=0.24") +
  ylab("Proportion with disease")
plotx4


plotall <- ggarrange(plotx1, plotx2, plotx3, plotx4)
plotall

