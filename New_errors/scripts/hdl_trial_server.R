rm(list=ls())


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

#LD.path <- "~/Desktop/Summer_project/Data/UKB_imputed_SVD_eigen99_extraction"
#list.files(LD.path)
#res.HDL <- HDL.rg.parallel(gwas1.example, gwas2.example, LD.path, numCores = 2)
#res.HDL

#Do the HDL function on the actual datasets

#Load and format the GWAS datasets to be compatible with the HDL function
#This script does both those things

source("~/Summer_project/Scripts/process_colnames.R")

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

#The format of the dataset should be snp,a1, a2,n,beta,se,(Z-score)

LD.path <- "~/Summer_project/HDL/UKB_imputed_SVD_eigen99_extraction"

#output="/rds/general/project/hda_students_data/live/Group2/Shinil/Summer_project/processed/hdl_trial.Rout"

#These were the three versions I tried running

#Version 1: Two different datasets, GWAS significant SNPs
#result.HDL1 <- HDL.rg.parallel(alz.data,p.data, LD.path,numCores = 2)
#result.HDL1

#Version 2: Same as above, the only change one dataset was subsetted to be the same
# size and the SNPs were forced to be the same
#"Simulated" SNPs for Alzheimer's dataset
#p.sim<-p.sig[1:nrow(alz.sig),]
#p.sim$SNP<-alz.sig$SNP
#result.HDL <- HDL.rg.parallel(alz.sig,p.sim, LD.path,numCores = 2)
#result.HDL

#Version 3: The same dataset against each other
result.HDL2 <- HDL.rg.parallel(p.sig,p.sig, LD.path,numCores = 2)
result.HDL2

data_list<-c("alz.data","als.data","epilepsy.data","lewy.data","ms.data","p.data","stroke.data")
combn(data_list,2)

end_time <- Sys.time()
end_time - start_time