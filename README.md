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
df_dose.responses_CTG <- read.csv('./exampleData_CTG.csv', header = T,sep = ',',check.names = F)

# load the ex vivo dose-response profiles (cell toxicity at five drug concentrations)
df_dose.responses_CTX <- read.csv('./exampleData_CTX.csv', header = T,sep = ',',check.names = F)

head(df_dose.responses_CTG)

##    DRUG_NAME CONCENTRATION_nM SCREEN_NAME    CELL_VIABILITY
## 1 Nelarabine         10000  AML_013_01            0.3012
## 2 Nelarabine          1000  AML_013_01            0.5565
## 3 Nelarabine           100  AML_013_01            0.7578
## 4 Nelarabine            10  AML_013_01            0.9133
## 5 Nelarabine             1  AML_013_01            0.8610
## 6 Decitabine         10000  AML_013_01            0.4829

# calculate the percentage of growth inhibition 
df_dose.responses.list_CTG <- DOSE_RESPONSE_PROCESS(df_dose.responses_CTG, viability = TRUE)
# calculate the percentage of cell death
df_dose.responses.list_CTX <- DOSE_RESPONSE_PROCESS(df_dose.responses_CTX, viability = TRUE)

# calculate DSS metrics (DSS1, DSS2, DSS3), AUC and relative IC50 for growth inhibition (CTG assay)
df.metrics_CTG <- CALC_METRICS(df_dose.responses.list_CTG[[1]], df_dose.responses.list_CTG[[2]], graph = FALSE)
# calculate DSS metrics (DSS1, DSS2, DSS3), AUC and relative IC50 for cell death (CTX assay)
df.metrics_CTX <- CALC_METRICS(df_dose.responses.list_CTX [[1]], df_dose.responses.list_CTX[[2]], graph = FALSE)



# compute the descriptive statistics of DSSs for each drug over CTG
df.metrics_CTG_summary <- as.data.frame(rbind(colMeans(as.matrix(df.metrics_CTG)),colSds(as.matrix(df.metrics_CTG)),colMedians(as.matrix(df.metrics_CTG)), colMads(as.matrix(df.metrics_CTG))))

# define names of statistics
rownames(df.metrics_CTG_summary  ) <- c('mean', 'sd', 'median', 'mad')

df.metrics_CTX_summary <- as.data.frame(rbind(colMeans(as.matrix(df.metrics_CTX)),colSds(as.matrix(df.metrics_CTX)),colMedians(as.matrix(df.metrics_CTX)), colMads(as.matrix(df.metrics_CTX))))

# define names of statistics
rownames(df.metrics_CTX_summary  ) <- c('mean', 'sd', 'median', 'mad')


# patient PIS
patients.PIS <- slice(df.metrics_CTG_summary['mean', rep(1:n())]) - slice(df.metrics_CTX_summary['mean', colnames(df.metrics_CTG_summary)],rep(1:n(), each = nrow(df.metrics_CTG_summary)))
