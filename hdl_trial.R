rm(list=ls())
setwd("~/Summer_project/HDL")

.libPaths(c("/rds/general/user/sr4515/home/R/x86_64-redhat-linux-gnu-library/3.6",
            "/usr/lib64/R/library",
            "/usr/share/R/library" ))
start_time <- Sys.time()

library(devtools)
#install_github("zhenin/HDL/HDL")
#install.packages("doSNOW")
library(doSNOW)
library(HDL)
library(GenomicSEM)

#How the function works
#source("~/Summer_project/HDL/HDL.run.R")
gwas2.example <- readRDS("~/Summer_project/HDL/gwas2.array.example.rds")
gwas2.example=gwas2.example[1:10,]
gwas1.example <- readRDS("~/Summer_project/HDL/gwas1.array.example.rds")
gwas1.example=gwas1.example[1:10,]
#untar("uk_biobank/UKB_imputed_SVD_eigen99_extraction.tar.gz")

LD.path <- "~/Desktop/Summer_project/Data/UKB_imputed_SVD_eigen99_extraction"
#list.files(LD.path)
res.HDL <- HDL.rg.parallel(gwas1.example, gwas2.example, LD.path, numCores = 2)
#res.HDL

#Do the HDL function on the actual datasets

#Load and format the GWAS datasets to be compatible with the HDL function
source("/rds/general/project/hda_students_data/live/Group2/Shinil/Summer_project/Scripts/load_significant_datasets.R")

alz.data<-alz.sig
p.data<-p.sig

#alz.data
#als.data
#epilepsy.data
#lewy.data
#ms.data
#p.data
#stroke.data

#Put it in a loop so it goes over each combination automatically 
alz.data<-alz.data[1:10,]
p.data<-p.data[1:10,]

#stroke.data<-stroke.data[1:100,]

#The format of the dataset should be snp,a1, a2,n,beta,se,(Z-score)
#head(alz.data)
alz.data<-alz.data[,c(3:8)]
#head(alz.data)
colnames(alz.data)<-c("SNP","a1","a2","or","se","p")

#head(p.data)
for (i in 1:nrow(p.data)){
  p.data$N[i]<-sum(p.data$N_cases[i],p.data$N_controls[i])
}
colnames(p.data)<-c("SNP","A1","A2","maf","b","se","p","N_cases","N_controls","N")
p.data$SNP<- gsub(".*:", "rs", p.data$SNP)

LD.path <- "~/Summer_project/HDL/UKB_imputed_SVD_eigen99_extraction"
output="/rds/general/project/hda_students_data/live/Group2/Shinil/Summer_project/processed/hdl_trial.Rout"
result.HDL <- HDL.rg.parallel(gwas1, gwas2, LD.path, numCores = 2)
result.HDL

data_list<-c("alz.data","als.data","epilepsy.data","lewy.data","ms.data","p.data","stroke.data")
combn(data_list,2)

end_time <- Sys.time()
end_time - start_time