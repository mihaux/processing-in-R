# this script plots distribution plots of counts
# output files:  

# interesting source: http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf

# install (if necessary) and load package
if (!requireNamespace("DESeq", quietly = TRUE)) BiocManager::install("DESeq"); library(DESeq)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
  stop("2 argument must be supplied: 
        \n(1 - input) path to .csv file with count data and
        \n(2 - output) path to folder for output files.", call.=FALSE)
}

# Example of usage: 
# Rscript run_Distribution_plot.R /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/run_1/distribution_plot

# INPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_nodups_run1_SE_mod.csv

# OUTPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/run_1/distribution_plot/[TO BE COMPLETED]

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)

# histogram of raw counts data
ggplot(rawCountTable, aes(x = control1)) + geom_histogram(fill = "#525252", binwidth = 2000)



# Possible ways of transforming RNA-seq data:
# logarithm transformation => it will get rid of some extreme values. 
# variance-stabilizing transformation (VST), implemented in the DESeq package (Anders and Huber, 2010)


#vst <- function(countdata){
#  condition <- factor(rep("Tumour", ncol(countdata)))
#  countdata <- newCountDataSet(countdata,condition )
#  countdata <- estimateSizeFactors( countdata )
#  cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
#  vstdata <- varianceStabilizingTransformation( cdsBlind )
#  return(exprs(vstdata))
#}
 
#load("rnaseq_lusc_example_SeqQC.Rda")
#data.log2 <- log2(data+1)
#data.vst <- vst(data)

# write plots
#write.csv(df, str_replace(args[1], ".csv", "_mod.csv"))

#cat("Finished!", "\nCreated:", str_replace(args[1], ".csv", "_mod.csv"))

