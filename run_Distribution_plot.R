# this script plots distribution plots of counts
# output files:  

# interesting source: http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf

# install (if necessary) and load package
if (!requireNamespace("DESeq", quietly = TRUE)) BiocManager::install("DESeq"); library(DESeq)
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse"); library(tidyverse)
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra"); library(gridExtra)
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid"); library(grid)
library(ggplot2)

# create a shortcut for the OneDrive directory where all files are stored
main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"      # on my mac
# main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
  stop("2 argument must be supplied: 
        \n(1 - input) path to .csv file with count data and
        \n(2 - output) path to folder for output files.", call.=FALSE)
}

# Example of usage: 
# Rscript run_Distribution_plot.R /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/run_1/distribution_plots

# INPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_nodups_run1_SE_mod.csv

# OUTPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/run_1/distribution_plot/[TO BE COMPLETED]

args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/distribution_plots/"))

# extract info about data from input filename
info_1 <- tail(unlist(str_split(args[1], "/")), n=1)
info_2 <- unlist(str_split(info_1, "_"))[3:5]
info_3 <- paste(info_2[1], info_2[2], info_2[3], sep = "_")

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)

# histogram of raw counts data per sample
histo_list <- list()
for (i in 1:ncol(df)) {
  histo_list[[i]] <- ggplot(df, aes_string(x = colnames(df)[i])) + geom_histogram(binwidth = 5000) 
}

# define "shortcut" for removing y axis title
rem_y <- theme(axis.title.y=element_blank())

# create multiple plots 4 x 3
pdf(file=paste0(args[2], "raw_counts_histo_", info_3, ".pdf"), width=12, height=12)
grid.arrange(histo_list[[1]], histo_list[[2]]+rem_y, histo_list[[3]]+rem_y, histo_list[[4]]+rem_y, 
             histo_list[[5]], histo_list[[6]]+rem_y, histo_list[[7]]+rem_y, histo_list[[8]]+rem_y, 
             histo_list[[9]], histo_list[[10]]+rem_y, histo_list[[11]]+rem_y, histo_list[[12]]+rem_y, nrow = 3,
             top = textGrob(paste("Histogram of raw counts data per sample:", info_3, ".Part 1 of 4", sep = " "), 
                            gp=gpar(fontsize=20,font=3)))

grid.arrange(histo_list[[13]], histo_list[[14]]+rem_y, histo_list[[15]]+rem_y, histo_list[[16]]+rem_y, 
             histo_list[[17]], histo_list[[18]]+rem_y, histo_list[[19]]+rem_y, histo_list[[20]]+rem_y, 
             histo_list[[21]], histo_list[[22]]+rem_y, histo_list[[23]]+rem_y, histo_list[[24]]+rem_y,nrow = 3,
             top = textGrob(paste("Histogram of raw counts data per sample:", info_3, ".Part 2 of 4", sep = " "), 
                            gp=gpar(fontsize=20,font=3)))

grid.arrange(histo_list[[25]], histo_list[[26]]+rem_y, histo_list[[27]]+rem_y, histo_list[[28]]+rem_y, 
             histo_list[[29]], histo_list[[30]]+rem_y, histo_list[[31]]+rem_y, histo_list[[32]]+rem_y, 
             histo_list[[33]], histo_list[[34]]+rem_y, histo_list[[35]]+rem_y, histo_list[[36]]+rem_y, nrow = 3,
             top = textGrob(paste("Histogram of raw counts data per sample:", info_3, ".Part 3 of 4", sep = " "), 
                            gp=gpar(fontsize=20,font=3)))

grid.arrange(histo_list[[37]], histo_list[[38]], histo_list[[39]], histo_list[[40]], histo_list[[41]], nrow = 2,
             top = textGrob(paste("Histogram of raw counts data per sample:", info_3, ".Part 4 of 4", sep = " "), 
                            gp=gpar(fontsize=20,font=3)))
dev.off()

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

