# another run of PCA, this time using the 'mixOmics' package

#BiocManager::install('mixOmics')
library(mixOmics)

# load data
#data(liver.toxicity)
df <- read.csv("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end/normalised_cpm.csv", row.names = 1)
df_clinic <- read.csv("/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary.csv")


# select imput data (data.frame, where rows = samples and cols = genes)
X <- liver.toxicity$gene


# 1 Run the PCA
MyResult.pca <- pca(df_ano_final)  

# PCA plot
plotIndiv(MyResult.pca)

# Plot the variables
plotVar(MyResult.pca)

# 2 Run sparse PCA
# (sparse PCA can be applied to select the top 5 variables contributing to each of the two components in PCA. 
# The user specifies the number of variables to selected on each component, 
# for example, here 5 variables are selected on each of the first two components (keepX=c(5,5)):)
MyResult.spca <- spca(X, keepX=c(5,5)) # 1 Run the method
plotIndiv(MyResult.spca)

# 3 Plot the variables
plotVar(MyResult.spca)








