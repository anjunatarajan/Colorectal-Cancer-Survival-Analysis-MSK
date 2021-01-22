library(dplyr)
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(readxl)

#This script creates 3 different CoxPH Box Plots

# 1: Clinical and gene Cox model
df <- read_excel('X1CRC_FINAL_GENEDATA.xls')

time <- df$`Overall Survival (Months)`
status <- df$`Overall Survival Status (Num)`

crcwithin <- within(df, {
  Sex <- factor(Sex, labels = c("Female", "Male"))
  Primary_Tumor_Location <- factor(Primary_Tumor_Location, labels = c("Left", "Right"))
  Stage_at_Diagnosis <- factor(Stage_At_Diagnosis, labels = c("I", "II", "III", "IV"))
})

modelwithgene <- coxph(Surv(time, status) ~ Sex + Primary_Tumor_Location + Stage_at_Diagnosis + Metastasectomy + APC + KRAS + BRAF + BCL2L1 + SMAD4, data = crcwithin )
ggforest(modelwithgene)

# 2: Just gene Cox model
alldata <- read_excel('X1CRC_FINAL_GENEDATA.xls')

time_gene <- alldata$`Overall Survival (Months)`
status_gene <- alldata$`Overall Survival Status (Num)`

bigmodel_gene <- coxph(Surv(time_gene, status_gene) ~ APC + KRAS + TP53 + PIK3CA + SMAD4 + FBXW7 + BRAF + TCF7L2 + SOX9 + ATM + PTPRT + BCL2L1, data = alldata)
ggforest(bigmodel_gene)

#3: Just clinical
crcwithin <- within(df, {
  Sex <- factor(Sex, labels = c("Female", "Male"))
  Primary_Tumor_Location <- factor(Primary_Tumor_Location, labels = c("Left", "Right"))
  Stage_at_Diagnosis <- factor(Stage_At_Diagnosis, labels = c("I", "II", "III", "IV"))
})

bigmodel1 <- coxph(Surv(time, status) ~ Sex + Primary_Tumor_Location + Stage_at_Diagnosis + Metastasectomy, data = crcwithin )
ggforest(bigmodel1)
