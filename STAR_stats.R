# STAR processing script that allows to merge all [sample_ID]_Log.final.out files and generate a summary table

# TODO:

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

args <- commandArgs(trailingOnly = TRUE)

#args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam", 
#          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement")

if (length(args)!=2) {
  stop("2 arguments must be supplied: \n(1 - input) path to directory with data and \n(2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[2])





