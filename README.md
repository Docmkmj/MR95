
PATH<-"E:\\sum"
setwd(PATH)
library(ggplot2)
library(plyr)
library(data.table)
library(devtools)
library(MendelianRandomization)
library(TwoSampleMR)

EXP1<-data.table::fread("exposure.F.csv")
EXP1<-data.frame(EXP1)
head(EXP1)
colnames(EXP1)

colnames(EXP1)[c(1,5,6,7,12,10,11)]<-
  c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure",
    "pval.exposure","beta.exposure","se.exposure")

EXP1$id.exposure<-""
EXP1$exposure<-""

colnames(EXP1)[13] <- "samplesize.exposure" 
head(EXP1,6)
EXP1_IV<-EXP1

EXP1_IV<-subset(EXP1,pval.exposure<5e-08)

EXP1_IV<-clump_data(EXP1_IV,clump_kb = 10000,clump_r2 = 0.01,pop = "EUR")#如果是亚洲则EUR改为EAS

OUT<-data.table::fread("")
OUT<-data.frame(OUT)
head(OUT)
colnames(OUT)

colnames(OUT)[c(4,3,5,7,9,10,11)]<-c("effect_allele.outcome","other_allele.outcome", "SNP","pval.outcome","beta.outcome","se.outcome","eaf.outcome")

OUT$id.outcome<-""
OUT$outcome<-""
OUT$samplesize.outcome<-
head(OUT)#


total<-merge(OUT,EXP1_IV,by.x="SNP",by.y="SNP",all = F)

total<-subset(total,pval.outcome>5e-08)
total<-total[!duplicated(total$SNP),]


EXP3<-total[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure", "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure","samplesize.exposure")]
OUT3<-total[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome","samplesize.outcome")]

data_h<-harmonise_data(exposure_dat=EXP3,outcome_dat=OUT3,action=2)
View(data_h)
write.table(data_h, file="data.xls",sep="\t",quote=F)

R2<-NULL
attach(data_h)
R2<-data_h$beta.exposure*data_h$beta.exposure/(data_h$beta.exposure*data_h$beta.exposure+data_h$se.exposure*data_h$se.exposure*data_h$samplesize.exposure)
F_statistics<-R2*data_h$samplesize.exposure/(1-R2)
data_h<-cbind(data_h,R2,F_statistics)
write.table(data_h, file="F&R.xls",sep="\t",row.names = F)

data_h<-data_h[-c(),]

library(phenoscanner)
dim(data_h)[1]
PhenoScan=phenoscanner(snpquery=data_h$SNP[1:8],pvalue = 5e-08)
#导出数据
write.csv(PhenoScan$result,file="PhenoScan.csv")

data_h<-data_h[-c(3,5,10,11),]
write.table(data_h, file="混杂因素剔出以后的SNP.xls",sep="\t",quote=F)


mr<-mr(data_h)

mr_OR<-generate_odds_ratios(mr) 
mr_OR

write.csv(mr_OR,file="mr_OR11.csv")


p1 <- mr_scatter_plot(mr[1:3,], data_h)
p1[[1]]
print(p1)

mr_outcome_loo <- mr_leaveoneout(data_h)
p3 <- mr_leaveoneout_plot(mr_outcome_loo)
p3[[1]]

mr_outcome_single <- mr_singlesnp(data_h)
p2 <- mr_forest_plot(mr_outcome_single)
p2[[1]]

mr_outcome_single <- mr_singlesnp(data_h)
p4 <- mr_funnel_plot(mr_outcome_single)
p4[[1]]

H<-mr_heterogeneity(data_h)
H
write.csv(H,file="h.csv")

ple <- mr_pleiotropy_test(data_h)
ple
write.csv(ple,file="p.csv")


library(MRPRESSO)
presso<-run_mr_presso(data_h, NbDistribution = 1000, SignifThreshold = 0.05)
presso
or1<-exp(presso[[1]][["Main MR results"]][["Causal Estimate"]][1])
or1
or2<-exp(presso[[1]][["Main MR results"]][["Causal Estimate"]][1]-1.96*presso[[1]][["Main MR results"]][["Sd"]][1])
or2
or3<-exp(presso[[1]][["Main MR results"]][["Causal Estimate"]][1]+1.96*presso[[1]][["Main MR results"]][["Sd"]][1])
or3
