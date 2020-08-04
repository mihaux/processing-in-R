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
transcript_pos <- gtf[4]

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
tab_1 <- cbind(as.vector(gtf[1]), transcript_ID[,2], transcript_pos$V4)
names(tab_1) <- c("chromosomes", "IDs", "positions")

# subset for each chromosome and then sort it
tab_1_chr1 <- tab_1[which(tab_1$chromosomes == "chr1"),]
tab_1_chr2 <- tab_1[which(tab_1$chromosomes == "chr2"),]
tab_1_chr3 <- tab_1[which(tab_1$chromosomes == "chr3"),]
tab_1_chr4 <- tab_1[which(tab_1$chromosomes == "chr4"),]
tab_1_chr5 <- tab_1[which(tab_1$chromosomes == "chr5"),]
tab_1_chr6 <- tab_1[which(tab_1$chromosomes == "chr6"),]
tab_1_chr7 <- tab_1[which(tab_1$chromosomes == "chr7"),]
tab_1_chr8 <- tab_1[which(tab_1$chromosomes == "chr8"),]
tab_1_chr9 <- tab_1[which(tab_1$chromosomes == "chr9"),]
tab_1_chr10 <- tab_1[which(tab_1$chromosomes == "chr10"),]
tab_1_chr11 <- tab_1[which(tab_1$chromosomes == "chr11"),]
tab_1_chr12 <- tab_1[which(tab_1$chromosomes == "chr12"),]
tab_1_chr13 <- tab_1[which(tab_1$chromosomes == "chr13"),]
tab_1_chr14 <- tab_1[which(tab_1$chromosomes == "chr14"),]
tab_1_chr15 <- tab_1[which(tab_1$chromosomes == "chr15"),]
tab_1_chr16 <- tab_1[which(tab_1$chromosomes == "chr16"),]
tab_1_chr17 <- tab_1[which(tab_1$chromosomes == "chr17"),]
tab_1_chr18 <- tab_1[which(tab_1$chromosomes == "chr18"),]
tab_1_chr19 <- tab_1[which(tab_1$chromosomes == "chr19"),]
tab_1_chr20 <- tab_1[which(tab_1$chromosomes == "chr20"),]
tab_1_chr21 <- tab_1[which(tab_1$chromosomes == "chr21"),]
tab_1_chr22 <- tab_1[which(tab_1$chromosomes == "chr22"),]
tab_1_chrM <- tab_1[which(tab_1$chromosomes == "chrM"),]
tab_1_chrX <- tab_1[which(tab_1$chromosomes == "chrX"),]
tab_1_chrY <- tab_1[which(tab_1$chromosomes == "chrY"),]

# sort within chr only
tab_1_chr1_sorted_pos <- tab_1_chr1[order(tab_1_chr1$positions),]
tab_1_chr2_sorted_pos <- tab_1_chr2[order(tab_1_chr2$positions),]
tab_1_chr3_sorted_pos <- tab_1_chr3[order(tab_1_chr3$positions),]
tab_1_chr4_sorted_pos <- tab_1_chr4[order(tab_1_chr4$positions),]
tab_1_chr5_sorted_pos <- tab_1_chr5[order(tab_1_chr5$positions),]
tab_1_chr6_sorted_pos <- tab_1_chr6[order(tab_1_chr6$positions),]
tab_1_chr7_sorted_pos <- tab_1_chr7[order(tab_1_chr7$positions),]
tab_1_chr8_sorted_pos <- tab_1_chr8[order(tab_1_chr8$positions),]
tab_1_chr9_sorted_pos <- tab_1_chr9[order(tab_1_chr9$positions),]
tab_1_chr10_sorted_pos <- tab_1_chr10[order(tab_1_chr10$positions),]
tab_1_chr11_sorted_pos <- tab_1_chr11[order(tab_1_chr11$positions),]
tab_1_chr12_sorted_pos <- tab_1_chr12[order(tab_1_chr12$positions),]
tab_1_chr13_sorted_pos <- tab_1_chr13[order(tab_1_chr13$positions),]
tab_1_chr14_sorted_pos <- tab_1_chr14[order(tab_1_chr14$positions),]
tab_1_chr15_sorted_pos <- tab_1_chr15[order(tab_1_chr15$positions),]
tab_1_chr16_sorted_pos <- tab_1_chr16[order(tab_1_chr16$positions),]
tab_1_chr17_sorted_pos <- tab_1_chr17[order(tab_1_chr17$positions),]
tab_1_chr18_sorted_pos <- tab_1_chr18[order(tab_1_chr18$positions),]
tab_1_chr19_sorted_pos <- tab_1_chr19[order(tab_1_chr19$positions),]
tab_1_chr20_sorted_pos <- tab_1_chr20[order(tab_1_chr20$positions),]
tab_1_chr21_sorted_pos <- tab_1_chr21[order(tab_1_chr21$positions),]
tab_1_chr22_sorted_pos <- tab_1_chr22[order(tab_1_chr22$positions),]
tab_1_chrM_sorted_pos <- tab_1_chrM[order(tab_1_chrM$positions),]
tab_1_chrX_sorted_pos <- tab_1_chrX[order(tab_1_chrX$positions),]
tab_1_chrY_sorted_pos <- tab_1_chrY[order(tab_1_chrY$positions),]

# get uniq IDs
chr1_positions <- unique(tab_1_chr1_sorted_pos$IDs)
chr2_positions <- unique(tab_1_chr2_sorted_pos$IDs)
chr3_positions <- unique(tab_1_chr3_sorted_pos$IDs)
chr4_positions <- unique(tab_1_chr4_sorted_pos$IDs)
chr5_positions <- unique(tab_1_chr5_sorted_pos$IDs)
chr6_positions <- unique(tab_1_chr6_sorted_pos$IDs)
chr7_positions <- unique(tab_1_chr7_sorted_pos$IDs)
chr8_positions <- unique(tab_1_chr8_sorted_pos$IDs)
chr9_positions <- unique(tab_1_chr9_sorted_pos$IDs)
chr10_positions <- unique(tab_1_chr10_sorted_pos$IDs)
chr11_positions <- unique(tab_1_chr11_sorted_pos$IDs)
chr12_positions <- unique(tab_1_chr12_sorted_pos$IDs)
chr13_positions <- unique(tab_1_chr13_sorted_pos$IDs)
chr14_positions <- unique(tab_1_chr14_sorted_pos$IDs)
chr15_positions <- unique(tab_1_chr15_sorted_pos$IDs)
chr16_positions <- unique(tab_1_chr16_sorted_pos$IDs)
chr17_positions <- unique(tab_1_chr17_sorted_pos$IDs)
chr18_positions <- unique(tab_1_chr18_sorted_pos$IDs)
chr19_positions <- unique(tab_1_chr19_sorted_pos$IDs)
chr20_positions <- unique(tab_1_chr20_sorted_pos$IDs)
chr21_positions <- unique(tab_1_chr21_sorted_pos$IDs)
chr22_positions <- unique(tab_1_chr22_sorted_pos$IDs)
chrM_positions <- unique(tab_1_chrM_sorted_pos$IDs)
chrX_positions <- unique(tab_1_chrX_sorted_pos$IDs)
chrY_positions <- unique(tab_1_chrY_sorted_pos$IDs)

# save outputs
write.csv(as.vector(chr1_positions), file = "chr1_positions.csv")
write.csv(as.vector(chr2_positions), file = "chr2_positions.csv")
write.csv(as.vector(chr3_positions), file = "chr3_positions.csv")
write.csv(as.vector(chr4_positions), file = "chr4_positions.csv")
write.csv(as.vector(chr5_positions), file = "chr5_positions.csv")
write.csv(as.vector(chr6_positions), file = "chr6_positions.csv")
write.csv(as.vector(chr7_positions), file = "chr7_positions.csv")
write.csv(as.vector(chr8_positions), file = "chr8_positions.csv")
write.csv(as.vector(chr9_positions), file = "chr9_positions.csv")
write.csv(as.vector(chr10_positions), file = "chr10_positions.csv")
write.csv(as.vector(chr11_positions), file = "chr11_positions.csv")
write.csv(as.vector(chr12_positions), file = "chr12_positions.csv")
write.csv(as.vector(chr13_positions), file = "chr13_positions.csv")
write.csv(as.vector(chr14_positions), file = "chr14_positions.csv")
write.csv(as.vector(chr15_positions), file = "chr15_positions.csv")
write.csv(as.vector(chr16_positions), file = "chr16_positions.csv")
write.csv(as.vector(chr17_positions), file = "chr17_positions.csv")
write.csv(as.vector(chr18_positions), file = "chr18_positions.csv")
write.csv(as.vector(chr19_positions), file = "chr19_positions.csv")
write.csv(as.vector(chr20_positions), file = "chr20_positions.csv")
write.csv(as.vector(chr21_positions), file = "chr21_positions.csv")
write.csv(as.vector(chr22_positions), file = "chr22_positions.csv")
write.csv(as.vector(chrM_positions), file = "chrM_positions.csv")
write.csv(as.vector(chrX_positions), file = "chrX_positions.csv")
write.csv(as.vector(chrY_positions), file = "chrY_positions.csv")

# past all files
all_positions <- c(as.vector(chr1_positions), as.vector(chr2_positions), as.vector(chr3_positions),
                   as.vector(chr4_positions), as.vector(chr5_positions), as.vector(chr6_positions),
                   as.vector(chr7_positions), as.vector(chr8_positions), as.vector(chr9_positions),
                   as.vector(chr10_positions), as.vector(chr11_positions), as.vector(chr12_positions),
                   as.vector(chr13_positions), as.vector(chr14_positions), as.vector(chr15_positions),
                   as.vector(chr16_positions), as.vector(chr17_positions), as.vector(chr18_positions),
                   as.vector(chr19_positions), as.vector(chr20_positions), as.vector(chr21_positions),
                   as.vector(chr22_positions), as.vector(chrM_positions), as.vector(chrX_positions),
                   as.vector(chrY_positions))

write.csv(as.vector(unique(all_positions)), file = "all_positions_unique.csv")

#sum(length(chr1_positions), length(chr2_positions), length(chr3_positions), length(chr4_positions), 
#    length(chr5_positions), length(chr6_positions), length(chr7_positions), length(chr8_positions), 
#    length(chr9_positions), length(chr10_positions), length(chr11_positions), length(chr12_positions),
#    length(chr13_positions), length(chr14_positions), length(chr15_positions), length(chr16_positions), 
#    length(chr17_positions), length(chr18_positions), length(chr19_positions), length(chr20_positions), 
#    length(chr21_positions), length(chr22_positions), length(chrM_positions), length(chrX_positions),
#    length(chrY_positions))


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

