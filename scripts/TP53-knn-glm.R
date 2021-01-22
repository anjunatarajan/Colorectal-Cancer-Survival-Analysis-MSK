#KNN model

#Load the data, normalize it and splice it into training and test dataset
library(readxl)

genedata <- read_excel("MERGED_TP53_Hotspot.xls")

genedata <- genedata[c('TP53','Mets_Site_First_Bone','Mets_Site_First_Brain','Mets_Site_First_Gyn','Mets_Site_First_Liver','Mets_Site_First_Ln','Mets_Site_First_Lung','Mets_Site_First_Pelvis','Mets_Site_First_Pt_oment_abd','Fraction_Genome_Altered','Overall Survival Status (Num)','Sex_Num','Sample_Type_Num')]


#function for normalizing data
normalize <- function(x) { return((x - min(x)) / (max(x) - min(x)))}

head(genedata)
print('Total datapoints ')
NROW(genedata)
genedata_norm <- as.data.frame(genedata, normalize)

#Splice the data into training and test dataset
set.seed(123)
dat.d <- sample(1:nrow(genedata_norm),size=nrow(genedata_norm)*0.7, replace = FALSE) #random 70% sample

train.genedata <- genedata_norm[dat.d,2:12] #70% training data
test.genedata <- genedata_norm[-dat.d,2:12] #Remaining 30% test data

#Let's create the expected results vector

train.TP53 <- genedata_norm[dat.d,1] #creating a frame with answer to compare
test.TP53 <- genedata_norm[-dat.d,1]

print('training sample info: row, col of feature matrix')
NROW(train.genedata)
NCOL(train.genedata)
print('training sample info: row, col of results vector')
NROW(train.TP53)
NCOL(train.TP53)

print('test sample info: row, col of feature matrix')
NROW(test.genedata)
NCOL(test.genedata)
print('test sample info: row, col of results vector')
NROW(test.TP53)
NCOL(test.TP53)

#Builing a model to predict whether there is TP53 mutation or not given Mets site
#We think this is possbile since the Fisher Test showed a close association
#between gene mutation and Metastatic sites
library(class)
library(dplyr)
library("RSNNS")
#square root of 692 is 26.3. Let's pick 27 for K value

#knn.27 <- knn(train = train.genedata[complete.cases(train.genedata), ], test = test.genedata[complete.cases(test.genedata), ], cl = train.TP53, k=27, prob= TRUE)
knn.27 <- knn(train=train.genedata, test=test.genedata, cl=train.TP53, k=10)
#check the prediction using confusion table against the test data
table(knn.27, test.TP53)
summary(knn.27)
ACC.27 <- 100 * sum(test.TP53 == knn.27)/NROW(test.TP53)
print(ACC.27)

#ROC AUC Curve
library(ROCR)
knn_isolet <- knn(train.genedata, test.genedata, train.TP53, k=10, prob=TRUE)
prob <- attr(knn_isolet, "prob")
prob <- 2*ifelse(knn_isolet == "-1", 1-prob, prob) - 1
pred_knn <- prediction(prob, test.TP53)

pred_knn <- performance(pred_knn, "tpr", "fpr")

plot(pred_knn, avg= "threshold", colorize=T, lwd=3, main="ROC Curve for TP53 Prediction")

#Same dataset applied to Logistic Regression model

lin_model = glm(I(train.TP53=="0")~., data=train.genedata, family="binomial")
formula<- (train.TP53=="0")~.
logit <- glm(formula, data = train.genedata, family = 'binomial')
summary(logit)
predict <- predict(logit, test.genedata, type = 'response')
# confusion matrix
table_mat <- table(test.TP53, predict < 0.5)
table_mat
accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test
library(ROCR)
ROCRpred <- prediction(predict, test.TP53)
ROCRperf <- performance(ROCRpred, 'tpr', 'fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7))
