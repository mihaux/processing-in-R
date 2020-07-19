# this script performs PCA analysis on read counts data (output from featureCounts)

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)


#args <- commandArgs(trailingOnly = TRUE)

args <- c("/cloud/project/normalised_cpm.csv",
          "/cloud/project/cic_clinical_data_v2_summary.csv",
          "/cloud/project/outputs")

cat("Example of usage: \n Rscript downstream.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .csv annotation file 
       \nand (3 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[3], sep="\n")

setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)

# if running for mode_II then, the colnames need to be changed
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
IDs_final <- sub("X*", "", IDs)
colnames(df) <- IDs_final

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

#Create working folder and get data
dir.create(paste(workFolder, "std_DeSeq2", sep="/"), recursive = TRUE)
setwd(paste(workFolder, "std_DeSeq2", sep="/"))
counts <- read.delim(rSubReadFile)

#subset data and filter it for empty rows
counts <- as.matrix(counts[, columnsFromDataFile])
storage.mode(counts)<-"integer"
counts <-counts[rowSums(counts>3)>2,]

#Create experiment design model
sampleConditions=conditionGroupsForeachBamFile
sampleNames=SampleNamesForEachOfTheBamFiles
sampleTable<-data.frame(sample=sampleNames, condition=sampleConditions)
#sampleTable


