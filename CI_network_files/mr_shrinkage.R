# Script to try shrinkage of correlation values to make MR matrix positive definite

path<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

load("../Results/MR/MR_results.rds")
effect_estimate<-read.csv("../Results/MR/effect_estimate_matrix.csv")
rownames(effect_estimate)<-effect_estimate[,1]
effect_estimate<-effect_estimate[,-1]
effect_estimate<-as.matrix(effect_estimate)

mr_ci<-mr_table[,c("Exposure","Outcome","Estimate","SE")]
mr_ci$zscore<-mr_ci$Estimate/mr_ci$SE

# Note: In the matrix, the columns are the outcomes and rows are the exposures
mr_zscores<-read.csv("../Results/MR/mr_zscores.csv")
rownames(mr_zscores)<-mr_zscores[,1]
mr_zscores<-mr_zscores[,-1]
mr_zscores<-as.matrix(mr_zscores)

#install.packages("DescTools")
library(DescTools)
mr_cors_ci<-FisherZInv(mr_zscores)

for (i in 1:nrow(mr_cors_ci)){
  mr_cors_ci[i,i]=0
}

mr_cors_ci.no_ms<-mr_cors_ci[-5,-5] # Remove MS
mr_cors_ci.no_epi<-mr_cors_ci[-3,-3]# Remove epilepsy
mr_cors_limited<-mr_cors_ci[c(-3,-5),c(-3,-5)] #Remove both

lambda=seq(0.00001,0.1,by=0.00001)

for (i in 1:length(lambda)){
  corr_shrink = sign(mr_cors_ci.no_ms) * (min((abs(mr_cors_ci.no_ms) - lambda[i]), 0))
  verdict=is.positive.definite(corr_shrink)
  if (verdict==TRUE){
    print(paste(lambda[i]," = smallest lambda to get a positive definite matrix"))
  } else {
    x=i
  }
} # All 10,000 were false

for (i in 1:length(lambda)){
  corr_shrink = sign(mr_cors_ci.no_epi) * (min((abs(mr_cors_ci.no_epi) - lambda[i]), 0))
  verdict=is.positive.definite(corr_shrink)
  if (verdict==TRUE){
    print(paste(lambda[i]," = smallest lambda to get a positive definite matrix"))
  } else {
    x=i
  }
} # All 10,000 were false

for (i in 1:length(lambda)){
  corr_shrink = sign(mr_cors_limited) * (min((abs(mr_cors_limited) - lambda[i]), 0))
  verdict=is.positive.definite(corr_shrink)
  if (verdict==TRUE){
    print(paste(lambda[i]," = smallest lambda to get a positive definite matrix"))
  } else {
    x=i
  }
} # All 10,000 were false