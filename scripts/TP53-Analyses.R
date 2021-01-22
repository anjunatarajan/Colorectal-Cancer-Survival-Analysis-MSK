library(dplyr)
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)

# Generating a Kaplan Meier curve with just TP53 Amino Acid Change status
# 0 = WT (Wild Type)
# 1 = TP53_Other (Non Hotspot mutations and also mutations < 15)
# 2 = G245
# 3 = R175
# 4 = R248
# 5 = R273
# 6 = R282

alldata <-  read_excel('X1CRC_FINAL_GENEDATA.xls')
KM_TP53_AA <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_TP53_AA
KM_TP53_AA_FIT <- survfit(KM_TP53_AA ~ alldata$HS_Status, data = alldata)
summary(KM_TP53_AA_FIT)
ggsurvplot(KM_TP53_AA_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)


#Creating 2 CoxPH Box plots

tp53_aa <- alldata
time_tp53 <- tp53_aa$`Overall Survival (Months)`
status_tp53 <- tp53_aa$`Overall Survival Status (Num)`

tp53_within <- within(tp53_aa, {
  Amino_Acid_Change <- factor(HS_Status, labels = c("0", "1", "2", "3", "4", "5", "6"))
})

#1: Creating an overall Plot
bigmodel_tp53 <- coxph(Surv(time_tp53, status_tp53) ~ WT_TP53 + TP53_Other + G245 + R175 + R248 + R273 + R282, data = tp53_aa)
ggforest(bigmodel_tp53)
#2: Creating a factored Box Plot using WT as the reference
factored_bigmodel_tp53 <- coxph(Surv(time_tp53, status_tp53) ~ Amino_Acid_Change, data = tp53_within)
ggforest(factored_bigmodel_tp53)

#This code analyses the correlation between TP53 mutations and Metastatic sites
# fisher test - TP53 and metastatic sites
#Fisher Test: Correlation between TP53 mutations and Mets site
TP53MetsSite <- matrix(c(8,3,25,367,77,87,14,67, 4,2,13,136,26,51,6,49, 1,0,2,20,3,7,1,3, 2,0,3,48,12,6,2,6, 0,0,1,42,8,6,2,6, 2,0,5,53,11,8,6,14, 0,0,0,17,3,7,1,4), 7, 8,
                       dimnames = list(income = c("TP53 Other", "WT_TP53", "G245", "R175","R248","R273","R282"),
                                       metsite = c("Bone", "Brain", "Gynecological", "Liver", "Lymph", "Lung","Pelvis","PeriOmiAbd")))
fisher.test(TP53MetsSite,simulate.p.value=TRUE)


#Pairwise analysis of TP53 hotspot sites
library(readxl)
genedata <- read_excel("MERGED_TP53_Hotspot.xls")
genedata <- genedata[c('HS_Status','Time','Event')]

data1<-genedata[!(is.na(genedata$Time) | genedata$Time==""), ]
# Pairwise survdiff
res <- pairwise_survdiff(Surv(Time, Event) ~ HS_Status,
     data = data1, p.adjust.method ='fdr')
res
