
# Colorectal Cancer Survival Analysis

# Description

This R Repository is a comprehensive statistical analysis of the clinical and genomic attributes based on samples collected from 1134 colorectal cancer patients treated at the Memorial Sloan Kettering Cancer Institute. This data is available publicly through cBioPortal for cancer genomics. In addition to statistical analyses, further machine learning models were built to predict the presence or absence of TP53 gene mutations as a novel study since TP53 mutations were present in 77% of the cancer patients within the cohort, and no similar studies have been conducted on TP53. The intent behind the statistical analyses and machine learning models is to enable clinicians who are providing care for cancer patients to provide a more targeted pattern based treatments and build a more useful genomic and clinical profile for each patient.

## R packages used
The following R packages were instrumental for this project:
readxl
class
dplyr
RSMNS
ROCR
tidyverse
ggplot2
survival
survminer
lubridate
lranger
ggfortify

## Sample Kaplan Meier Survival Curves for univariate analysis of Clinical Attributes
![Median Age](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%201AGE%20AT%20DIAG%20KM%20CURVE.jpeg)
![Masectomy](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%20METASTASECTOMY%20KM%20CURVE.jpeg)
![Sex](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%20SEX%20KM%20CURVE.jpeg)


## Sample Kaplan Meier Survival Curves for multivariate analysis of genomic attributes
![Gene Alterations](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%20TP53PROTEIN%20CHANGE%20KM%20CURVE.jpeg)

## CoxPH Survival Analysis Box Plots to perform comparative reference based multivariate analysis on TP53 Gene Mutations
![TP53 Mutations](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%20AACHANGE%20COX%20MODEL.jpeg)

## CoxPH Survival Analysis Box Plots to perform comparative reference based multivariate analysis on all attributes
![CoxPH Clinical and Genomic ](https://github.com/anjunatarajan/Colorectal-Cancer-Survival-Analysis-MSK/blob/main/images/Rplot%20NEW%20gene%2Bclindata%20Cox%20Model.jpeg)


## Repository Contents
Non curated public data in Excel format containing clinical and genomic attributes for the patient cohort studied. Curated Text file and Excel file for performing pairwise analysis of certain insightful attributes. R code that executes statistical and machine learning methods on the Excel data


## Usage

Install R Studio or feel free to use an online tool like Kaggle to execute the R code. Download the Excel and Text files, alter the source code to load the data from your local path. Execute the R code and have fun! It is recommended to gain some basis understanding of Survival Analysis using Kaplan Meier, CoxPH models and the usage of KNN and GLM machine learning models prior to trying out this script.

## Roadmap

Future releases to this script will be made based on the following vision of the author:

 - Gain access to more data in order to develop the existing statistical analyses and to improve the Accuracy and AUC for the machine learning models
 - Obtain data about samples of primary tumor along with the information of which ones later on metastasized. Such a model will be vital for cancer care clinicians since they can use this information for preventive and more targeted treatment measures by monitoring specific patients and also by observing patients with similar prognosis based on how the Machine Learning models clustered them.
 - Submit the models for experimental use for clinicians providing cancer care based on reviews and approvals from bioinformaticians

# Support

Any one running into issues or have questions can open a Github issue and the author will respond to the issues within 2 weeks of opening the issue.
