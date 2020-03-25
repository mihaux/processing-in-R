# featureCounts processing script that allows to merge all .csv files with counts and stats

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

args <- commandArgs(trailingOnly = TRUE)

# for running in R-studio
# args <- c("/Users/ummz/R_local/input", "/Users/ummz/R_local/output")

if (length(args)!=2) {
  stop("2 arguments must be supplied: \n(1 - input) path to directory with data and \n(2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[1])

temp_counts = list.files(pattern="counts*")
temp_stats = list.files(pattern="stats*")

all_files_counts = lapply(temp_counts, function(x) read.csv(x,sep = ",", header = TRUE))
all_files_stats = lapply(temp_stats, function(x) read.csv(x,sep = ",", header = TRUE, row.names = 1))

setwd(args[2])

# merge all all_files
merged_counts <- do.call("cbind", all_files_counts)
merged_stats <- do.call("cbind", all_files_stats)

# remove every second column starting from the 3rd one
to_delete <- seq(3, length(merged_counts), 2)
merged_counts_final <- merged_counts[,-to_delete]
merged_stats_final <- merged_stats[,-to_delete]

# remove unnecessary elements from colnames
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(merged_counts))
IDs_final <- sub("X*", "", IDs)
  
colnames(merged_counts_final) <- IDs_final
colnames(merged_stats_final) <- IDs_final

# write output table as .csv [rows = genes | cols = IDs]
write.csv(merged_counts_final, file="counts_merged.csv")
write.csv(merged_stats_final, file="counts_merged.csv")

