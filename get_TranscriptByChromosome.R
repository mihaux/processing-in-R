# the scripts check where all mapped transcripts (colnames of count matrix) are located (i.e. on which chromosome)

if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)

# create a shortcut for the OneDrive directory where all files are stored
main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"      # on my mac
# main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .gtf annotation file and 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

args <- c(paste0(main_dir, "/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_dups_rr5_x.csv"),
          paste0(main_dir, "/gtf_files/hg38_Ian.gtf"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/"))

df <- read.csv(args[1], row.names = 1, header = TRUE) # data.frame with counts only



gtf <- read.delim(args[2], header = FALSE)

# gtf[1] => chromosome info
# gtf[9] => transcript info (to be split)

transcript_info <- gtf[9]

transcript_ID <- matrix(apply(transcript_info, 1, function(x) str_split_fixed(x,pattern = " ", n=4)),ncol = 4, byrow = TRUE)

# transcript_ID[,1] <- "gene_id"
# transcript_ID[,2] <- "NM_001276352;"     => in case of duplicates, the ID is just repeated
# transcript_ID[,3] <- "transcript_id"
# transcript_ID[,4] <- "NM_001276352; "    => in case of duplicates, "_dupX" is added (where X is a number)
# NOTE: transcript_ID[,4] needs to be modified to remove the space at the end => trimws(kot[,4])

# check which transcript_ID[,2] does not correspond totranscript_ID[,4]
length(which(transcript_ID[,2] != trimws(transcript_ID[,4]))) 

# remove ";" from the end of each element
transcript_ID[,2] <- gsub(";", "", transcript_ID[,2])

# => there are 10 707 rows, because they have "_dupX" included at the end but it always matches
# get these rows and remove "_dup1" from them
#gsub("_dup1", transcript_ID[which(transcript_ID[,2] != trimws(transcript_ID[,4])),4]) 
#unlist(strsplit(transcript_ID[which(transcript_ID[,2] != trimws(transcript_ID[,4])),4], "_dup1"))

# create a table with chr & transcript_ID[,2] 
tab_1 <- cbind(as.vector(gtf[1]), transcript_ID[,2])
tab_1_uniq <- unique(tab_1)

names(tab_1_uniq) <- c("chromosomes", "IDs")

# compare tab_1_uniq with rownames(df)

# sort tab_1_uniq by IDs  (dim = 58875     2)
tab_1_uniq_sorted <- tab_1_uniq[order(tab_1_uniq$IDs),] 

# sort rownames(df)       (length = 58608)
rows_df_sorted <- sort(rownames(df)) 

# tab_1_uniq_sorted is slightly longer, because of non-uniqness
# get duplicated elements and skip them
which(duplicated(x = tab_1_uniq_sorted$IDs))

tab_2_uniq_sorted <- tab_1_uniq_sorted[-which(duplicated(x = tab_1_uniq_sorted$IDs)),]

# sort by chromosome
tab_2_uniq_sorted_byChr <- tab_2_uniq_sorted[order(tab_2_uniq_sorted$chromosomes),] 

# save results as .csv
write.csv2(x=tab_2_uniq_sorted_byChr, file=paste0(args[3], "table_TranscriptsByChromosome.csv"))

