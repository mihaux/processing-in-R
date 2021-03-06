# STAR processing script that allows to merge all [sample_ID]_Log.final.out files and generate a summary table

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
library(lubridate)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - mode) running parameter (SE or PE),
       \n(2 - input) path to directory with data and 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

#args <- c("SE", 
#          paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/alignment_SE/bam/"),
#          paste0(main_dir, "/ANALYSES/RNA-seq_pipeline_QC/4_alignment/SE"))

args <- c("PE",
          paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/alignment_PE/bam/"),
          paste0(main_dir, "/ANALYSES/RNA-seq_pipeline_QC/4_alignment/PE"))

# Example of usage: 
# Rscript get_STARstats.R [arg1] [arg2] [args3]

cat("Directories with data (IN): "); cat(args[2], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")

setwd(args[2])

temp = list.files(pattern="*Log.final.out")
all_files = lapply(temp, function(x) read.csv(x,sep = "\t", header = FALSE))

setwd(args[3])

# extract IDs from file names 
IDs <- sub("_.*", "", temp)

# add colnames to each table
for (i in 1:length(all_files)) {
  names(all_files[[i]]) <- c("V1", IDs[i])
}

# merge all all_files
merged <- do.call("cbind", all_files)

# remove every second column starting from the 3rd one
to_delete <- seq(3, length(merged), 2)
merged_final <- merged[,-to_delete]

# extract 2nd & 3rd rows to calculate running time
df_time <- merged_final[2:3, -1]

RunningTime <- function(x){     # x -> data frame with 2nd & 3rd rows extracted from the main table
  
  # remove date, to keep time only (format HH:MM:SS)
  dropped_date <- sapply(df_time, function(x) substr(x, 8,16))
  
  # split into starting & end time
  df_start <- dropped_date[1,]
  df_stop <- dropped_date[2,]
  
  # convert to seconds
  df_start_converted <- sapply(df_start, function(x) period_to_seconds(hms(x)))
  df_stop_converted <- sapply(df_stop, function(x) period_to_seconds(hms(x)))

  # calculate running time (and add one second)
  run_time <- (df_stop_converted - df_start_converted) + 1

  # convert running time to seconds
  run_time_converted <- seconds_to_period(run_time)
  return(run_time_converted)
}

running_time <- RunningTime(df_time)

# replace the first row with calculated running times
df <- as.matrix(merged_final)
df[1,2:42] <- sprintf('%02d:%02d:%02d', hour(running_time), minute(running_time), second(running_time))

df_final <- data.frame(df[,-1], row.names = df[,1], stringsAsFactors = FALSE)

# remove 2nd, 3rd, 5th and 20th row
df_fin <- df_final[-c(2,3,7,22,27,34),]

# remove "X" from each colname
colnames(df_fin) <- colnames(merged_final)[-1]
  
# remove last 2 characters from each rowname and change first row and empty spaces
rownames(df_fin) <- gsub('.{2}$', '', rownames(df_fin))
rownames(df_fin)[1] <- "Running time of mapping"
rownames(df_fin) <- trimws(rownames(df_fin))

# generate a table (all together & split in 5)
write.csv(df_fin, file=paste0("STAR_stats_all_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_all_",  args[1], ".csv"))

# general info
df_sub1 <- df_fin[1:4,]
write.csv(df_sub1, file=paste0("STAR_stats_sub1_general_info_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_sub1_general_info_",  args[1], ".csv"))

#UNIQUE READS section
df_sub2 <- df_fin[5:18,]
write.csv(df_sub2, file=paste0("STAR_stats_sub2_unique_reads_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_sub2_unique_reads_",  args[1], ".csv"))

#MULTI-MAPPING READS section
df_sub3 <- df_fin[19:22,]
write.csv(df_sub3, file=paste0("STAR_stats_sub3_multi-mapping_reads_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_sub3_multi-mapping_reads_",  args[1], ".csv"))

#UNMAPPED READS section
df_sub4 <- df_fin[23:28,]
write.csv(df_sub4, file=paste0("STAR_stats_sub4_unmapped_reads_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_sub4_unmapped_reads_",  args[1], ".csv"))

#CHIMERIC READS section
df_sub5 <- df_fin[29:30,]
write.csv(df_sub5, file=paste0("STAR_stats_sub5_chimeric_reads_",  args[1], ".csv"))
cat("Created: ", paste0("STAR_stats_sub5_chimeric_reads_",  args[1], ".csv"))


