# Script to load all the datasets ------------------
rm(list=ls())

library(data.table)
data.path = "~/Summer_project/processed"
setwd(data.path)

#Alzheimer's
alz.sig = readRDS("alz.sig.rds")

#Parkinson's
p.sig=readRDS("p.sig.rds")
p.sig$SNP<- gsub(".*:", "rs", p.sig$SNP)
#p.sig=p.sig[complete.cases(p.sig)]

#Stroke
stroke.sig=readRDS("stroke.sig.rds")
#stroke.sig=stroke.sig[complete.cases(stroke.sig)]

#ALS
als.sig=readRDS("als.sig.rds")
#als.sig=als.sig[complete.cases(als.sig)]

#Lewy body dementia
lewy.sig=readRDS("lewy.sig.rds")
#lewy.sig=lewy.sig[complete.cases(lewy.sig)]

#Epilepsy
epilepsy.sig=readRDS("epilepsy.sig.rds")
#epilepsy.sig=epilepsy.sig[complete.cases(epilepsy.sig)]

#MS
ms.sig=readRDS("ms.sig.rds")
#ms.sig=ms.sig[complete.cases(ms.sig)]
