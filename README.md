# Persistence-Inducing-Score
This is the R code to calculate persister inducing score for a drug in high throughput drug screening experiments
PIS is a computational pipeline for quantification of drug capacity to growth inhibit the cells without killing them. The standard drug metrics, such as DSS, AUC and IC50 for both CTG and CTX assay can be computed in the Breeze application (https://github.com/potdarswapnil/Breeze). The persister inducing score (or PIS) is calculated as difference between the CTG-based DSS/AUC and CTX-based DSS/AUS. Given enough samples, one can estimate the mean, median and standard deviation of compound-specific PIS distribution over the samples. When combining drug response profiles from different sources, the ComBat batch effect correction method should be used when there appear visible batch effects.

**Citation**: If you use this work, please cite the paper Tripathi et al. Robust scoring of selective drug responses for patient-tailored therapy selection Nature Protocols (2024).

**Instructions**

R version 3.5.1 or newer is required.

The way to install packages 'sva' and 'pcaMethods' in Bioconductor differs from other packages in R:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace(c("sva", "pcaMethods"), quietly = TRUE))
    BiocManager::install(c("sva", "pcaMethods"))

**Under Linux/Unix**

# download the example data and R scripts
git clone https://github.com/anilkumargiri/Persistence-Inducing-Score.git

# change the directory
cd ./PIS

# start the R program
R

Start by loading libraries and functions

lapply(c("matrixStats","dplyr","reshape","reshape2", "scales", "drc", "caTools","ggplot2", "data.table", "stringr","MESS", "BiocManager","svMisc", "egg", "pheatmap", "sva", "pcaMethods"), library, character.only = T)

source('./PIS.R');
source('./HelperFunctions.R')

 PIS calculation example
Load the example data and compute PIS with CALC_METRICS function updated from Breeze


# load the ex vivo dose-response profiles (cell viability at five drug concentrations)
df_dose.responses <- read.csv('./exampleData_procedure1.csv', header = T,sep = ',',check.names = F)

head(df_dose.responses)

##    DRUG_NAME CONCENTRATION_nM SCREEN_NAME    CELL_VIABILITY
## 1 Nelarabine         10000  AML_013_01            0.3012
## 2 Nelarabine          1000  AML_013_01            0.5565
## 3 Nelarabine           100  AML_013_01            0.7578
## 4 Nelarabine            10  AML_013_01            0.9133
## 5 Nelarabine             1  AML_013_01            0.8610
## 6 Decitabine         10000  AML_013_01            0.4829

# calculate the percentage of growth inhibition 
df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = TRUE)

# calculate DSS metrics (DSS1, DSS2, DSS3), AUC and relative IC50
df.metrics <- CALC_METRICS(df_dose.responses.list[[1]], df_dose.responses.list[[2]], graph = FALSE)
# calculate the percentage of growth inhibition 
df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = TRUE)

Import the control sample DSS profiles
controls.dss <- read.csv('./controls/File_1_Drugname_response_DSS_10Healthy.txt', header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)

Calculate the selective drug response scores(sDSS, zDSS, and rDSS)

# compute the descriptive statistics of DSSs for each drug over 10 controls
controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))

# define names of statistics
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')

# let's set DSS2 as the drug response metrics of patient samples (as an example)
patients.dss <- as.data.frame(acast(df.metrics,df.metrics$Patient.num ~ df.metrics$drug , value.var  = 'DSS2'))

patients.dss[, 1:3]

##                    1-methyl-D-tryptophan 4-hydroxytamoxifen 8-amino-adenosine
## AML_003_01                     0                9.2              25.4
## AML_004_01                     0                2.3              22.6
## AML_013_01                     0                4.3              23.9

# normalize and scale patient-specific responses to drugs with control DSS profiles
# patient sDSS
patients.sdss <- patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss)))

# patient zDSS
patients.zdss <- (patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['sd', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)

# patient rDSS
patients.rdss <- (patients.dss - slice(controls.summary['median', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['mad', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)

