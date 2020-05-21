# this script performs PCA analysis on read counts data (output from featureCounts)

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)


#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/R_local/processing-in-R/extraction/counts_merged.csv",
          "/Users/ummz/R_local/processing-in-R/extraction/cic_clinical_data_v2_summary.csv",
          "/Users/ummz/R_local/processing-in-R/extraction/outputs")

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
counts <- read.delim(args[1], sep = ",", row.names = 1)

# filter data for empty rows
counts <-counts[rowSums(counts>3)>2,]


# run for genders
condycje <- anno$gender..1.male..2.female.
condycje_1 <- sub(1, "male", condycje)
conditionGroupsForeachBamFile <- sub(2, "female", condycje_1)

conditionGroupsForDataFirstIsReference=c("male","female")

#Create experiment design model
sampleConditions=conditionGroupsForeachBamFile
sampleNames = rownames(anno)
sampleTable<-data.frame(sample=sampleNames, condition=sampleConditions)
#sampleTable

#Set column names to those in model and not from the BAMs
colnames(counts)<-sampleTable$sample
#head(counts)
#Se row names in model design to sample names
rownames(sampleTable)<-sampleTable$sample
sampleTable

#Check all is well with names
all (rownames(sampleTable) %in% colnames(counts))
all (rownames(sampleTable) == colnames(counts))

#Put data from counts file in to correct format and set factors and levels
dds <- DESeqDataSetFromMatrix(countData=counts, colData=sampleTable, design=~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=conditionGroupsForDataFirstIsReference)

#set up the dds object
#mcols(dds)
featureData<-data.frame(gene=rownames(counts))
#head(featureData)
mcols(dds)<-DataFrame(mcols(dds), featureData)
#mcols(dds)
#does the work of the differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

fileBaseName="michal"	

#write.csv(resultsNames(dds), file=paste(fileBaseName,"Analysis.txt", sep="_"))
#save(dds,file=paste(fileBaseName,"DESeqDataSet.Rda", sep="_"))

nt<- normTransform(dds)
plotPCA(nt, intgroup=c("condition", "sample"))
dev.copy(png, paste(fileBaseName,"PCA_plot_1_Deseq2.png", sep="_"))
dev.off()#

