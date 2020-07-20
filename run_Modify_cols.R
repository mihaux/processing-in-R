# this script changes default colnames from "X11026.Aligned.sortedByCoord.out.bam" to "ID_11026"
# output files have "_mod" added in their names and are placed in the same directory as input files

# install (if necessary) and load package
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)

if (length(args)!=1) {
  stop("1 argument must be supplied: 
       \n(1 - input) path to .csv file with count data , call.=FALSE)
}

cat("Example of usage: \nRscript run_Modify_cols.R /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE.csv

  # INPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE.csv
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_nodups_run1_SE.csv

# OUTPUT examples
# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv

args <- commandArgs(trailingOnly = TRUE)

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)

# if running for mode_II then, the colnames need to be changed
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
IDs_final <- sub("X*", "ID_", IDs)
colnames(df) <- IDs_final
       
# write modified count tables as .csv with "_mod" added in the new file name
write.csv(df, str_replace(args[1], ".csv", "_mod.csv"))

