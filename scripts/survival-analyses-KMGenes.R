#statistical analyses using KM curves for the genomic profile
# Making Kaplan-Meier Curves

library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(ranger)
library(ggplot2)
library(ggfortify)

alldata <- read_excel("X1CRC_FINAL_GENEDATA.xls")

# APC
KM_APC <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_APC

KM_APC_FIT <- survfit(KM_APC ~ alldata$APC, data = alldata)
summary(KM_APC_FIT)
ggsurvplot(KM_APC_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# TP53
KM_TP53 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_TP53

KM_TP53_FIT <- survfit(KM_TP53 ~ alldata$TP53, data = alldata)
summary(KM_TP53_FIT)
ggsurvplot(KM_TP53_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# KRAS
KM_KRAS <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_KRAS

KM_KRAS_FIT <- survfit(KM_KRAS ~ alldata$KRAS, data = alldata)
summary(KM_KRAS_FIT)
ggsurvplot(KM_KRAS_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# PIK3CA
KM_PIK3CA <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_PIK3CA

KM_PIK3CA_FIT <- survfit(KM_PIK3CA ~ alldata$PIK3CA, data = alldata)
summary(KM_PIK3CA_FIT)
ggsurvplot(KM_PIK3CA_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# SMAD4
KM_SMAD4 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_SMAD4

KM_SMAD4_FIT <- survfit(KM_SMAD4 ~ alldata$SMAD4, data = alldata)
summary(KM_SMAD4_FIT)
ggsurvplot(KM_SMAD4_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# FBXW7
KM_FBXW7 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_FBXW7

KM_FBXW7_FIT <- survfit(KM_FBXW7 ~ alldata$FBXW7, data = alldata)
summary(KM_FBXW7_FIT)
ggsurvplot(KM_FBXW7_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# BRAF
KM_BRAF <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_BRAF

KM_BRAF_FIT <- survfit(KM_BRAF ~ alldata$BRAF, data = alldata)
summary(KM_BRAF_FIT)
ggsurvplot(KM_BRAF_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# TCF7L2
KM_TCF7L2 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_TCF7L2

KM_TCF7L2_FIT <- survfit(KM_TCF7L2 ~ alldata$TCF7L2, data = alldata)
summary(KM_TCF7L2_FIT)
ggsurvplot(KM_TCF7L2_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# SOX9
KM_SOX9 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_SOX9

KM_SOX9_FIT <- survfit(KM_SOX9 ~ alldata$SOX9, data = alldata)
summary(KM_SOX9_FIT)
ggsurvplot(KM_SOX9_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# ATM
KM_ATM <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_ATM

KM_ATM_FIT <- survfit(KM_ATM ~ alldata$ATM, data = alldata)
summary(KM_ATM_FIT)
ggsurvplot(KM_ATM_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# PTPRT
KM_PTPRT <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_PTPRT

KM_PTPRT_FIT <- survfit(KM_PTPRT ~ alldata$PTPRT, data = alldata)
summary(KM_PTPRT_FIT)
ggsurvplot(KM_PTPRT_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# BCL2L1
KM_BCL2L1 <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_BCL2L1

KM_BCL2L1_FIT <- survfit(KM_BCL2L1 ~ alldata$BCL2L1, data = alldata)
summary(KM_BCL2L1_FIT)
ggsurvplot(KM_BCL2L1_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)
