#Harmonize column names

library(doSNOW)
library(HDL)
library(GenomicSEM)

source("~/Summer_project/Scripts/load_significant_datasets.R")

head(alz.sig)
alz.sig<-alz.sig[,c(3:8)]
#alz.sig$Beta<-exp(alz.sig$Beta)
colnames(alz.sig)<-c("SNP","A1","A2","b","se","p")
alz.sig$N<-as.numeric(rep("94437",times=nrow(alz.sig)))

head(als.sig)
als.sig<-als.sig[,c(2,4,5,7,8,9)]
head(als.sig)
colnames(als.sig)<-c("SNP","A1","A2","b","se","p")
als.sig$N<-as.numeric(rep("36052",times=nrow(als.sig)))
als.sig$A1<-toupper(als.sig$A1)
als.sig$A2<-toupper(als.sig$A2)

for (i in 1:nrow(p.sig)){
  p.sig$N[i]<-sum(p.sig$N_cases[i],p.sig$N_controls[i])
}
colnames(p.sig)<-c("SNP","A1","A2","maf","b","se","p","N_cases","N_controls","N")
p.sig$SNP<- gsub(".*:", "rs", p.sig$SNP)
p.sig<-p.sig[,c(1,2,3,5,6,7,10)]

head(epilepsy.sig)
epilepsy.sig<-epilepsy.sig[,c(3,4,5,9,10)]
colnames(epilepsy.sig)<-c("SNP","A1","A2","zscore","p")
epilepsy.sig$N<-as.numeric(rep("44889",times=nrow(epilepsy.sig)))
epilepsy.sig$A1<-toupper(epilepsy.sig$A1)
epilepsy.sig$A2<-toupper(epilepsy.sig$A2)

head(lewy.sig)
lewy.sig<-lewy.sig[,c(1,2,5,6,9,10)]
colnames(lewy.sig)<-c("SNP","p","A1","A2","b","se")
lewy.sig$N<-as.numeric(rep("7372",times=nrow(lewy.sig)))

head(ms.sig)
ms.sig<-ms.sig[,c("variant_id","other_allele","effect_allele","p_value","beta","standard_error")]
colnames(ms.sig)<-c("SNP","A1","A2","p","b","se")
ms.sig$N<-as.numeric(rep("24091",times=nrow(ms.sig)))

head(stroke.sig)
stroke.sig<-stroke.sig[,c(1,2,3,5,6,7)]
colnames(stroke.sig)<-c("SNP","A1","A2","b","se","p")
stroke.sig$N<-as.numeric(rep("521612",times=nrow(stroke.sig)))
stroke.sig$A1<-toupper(stroke.sig$A1)
stroke.sig$A2<-toupper(stroke.sig$A2)

#stroke.sig$b<-exp(stroke.sig$b)
#colnames(stroke.sig)[4]<-"or"


