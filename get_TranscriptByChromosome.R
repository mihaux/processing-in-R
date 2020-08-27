# the scripts check where all mapped transcripts (colnames of count matrix) are located (i.e. on which chromosome)

if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .gtf annotation file and 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# USE ONLY raw - PE - all_chr DATA AS THEY ARE ALL THE SAME
args <- c(paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts"),
          paste0(main_dir, "/RNA-Sequencing/gtf_files/annotation/processed_hg38.ncbiRefSeq.gtf"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/"))

load(paste0(args[1], "/Raw_DESeq_dataset_all.Rda"))

df <- assay(dds_all)        # data.frame with counts only

gtf <- read.delim(args[2], header = FALSE)

# gtf[1] => chromosome info
# gtf[9] => transcript info (to be split)

#transcript_info <- gtf[9]
#transcript_pos <- gtf[4]

#transcript_ID <- matrix(apply(gtf, 1, function(x) str_split_fixed(x,pattern = "\t", n=4)),ncol = 4, byrow = TRUE)

# TODO: these lines need to be modified to make them up to date for the current .gtf file
# transcript_ID[,1] <- "gene_id"
# transcript_ID[,2] <- "NM_001276352;"     => in case of duplicates, the ID is just repeated
# transcript_ID[,3] <- "transcript_id"
# transcript_ID[,4] <- "NM_001276352; "    => in case of duplicates, "_dupX" is added (where X is a number)
# NOTE: transcript_ID[,4] needs to be modified to remove the space at the end => trimws(kot[,4])

# TODO: need to make sure that the column with 'gene_id' always equals the column with 'gene_name'

# check which transcript_ID[,2] does not correspond totranscript_ID[,4]
#length(which(transcript_ID[,2] != trimws(transcript_ID[,4]))) 

# remove ";" from the end of each element
#transcript_ID[,2] <- gsub(";", "", transcript_ID[,2])

# => there are 10 707 rows, because they have "_dupX" included at the end but it always matches
# get these rows and remove "_dup1" from them
#gsub("_dup1", transcript_ID[which(transcript_ID[,2] != trimws(transcript_ID[,4])),4]) 
#unlist(strsplit(transcript_ID[which(transcript_ID[,2] != trimws(transcript_ID[,4])),4], "_dup1"))

# create a table with chr & transcript_ID[,2] 
#tab_1 <- cbind(as.vector(gtf[1]), transcript_ID[,2], transcript_pos$V4)
#names(tab_1) <- c("chromosomes", "IDs", "positions")

names(gtf) <- c("chromosomes", "IDs")   # => 4 131 447

uniq_chr <- as.vector(unique(gtf$chromosomes))

chr_min_coords <- list()
chr_max_coords <- list()

for(i in 1:length(uniq_chr)){
  print(uniq_chr[i])
  chr_min_coords[i] <- min(which(gtf$chromosomes == uniq_chr[i]))
  chr_max_coords[i] <- max(which(gtf$chromosomes == uniq_chr[i]))
}

names(chr_min_coords) <- uniq_chr
names(chr_max_coords) <- uniq_chr

# create a new gtf matrix (with columns of 'chromosomes' and 'IDs' ordered correctly)

# keep only chr1 - chr22 and chrX, chrY and chrM

# subset for each chromosome and then sort it
gtf_chr1 <- gtf[which(gtf$chromosomes == "chr1"),]
gtf_chr2 <- gtf[which(gtf$chromosomes == "chr2"),]
gtf_chr3 <- gtf[which(gtf$chromosomes == "chr3"),]
gtf_chr4 <- gtf[which(gtf$chromosomes == "chr4"),]
gtf_chr5 <- gtf[which(gtf$chromosomes == "chr5"),]
gtf_chr6 <- gtf[which(gtf$chromosomes == "chr6"),]
gtf_chr7 <- gtf[which(gtf$chromosomes == "chr7"),]
gtf_chr8 <- gtf[which(gtf$chromosomes == "chr8"),]
gtf_chr9 <- gtf[which(gtf$chromosomes == "chr9"),]
gtf_chr10 <- gtf[which(gtf$chromosomes == "chr10"),]
gtf_chr11 <- gtf[which(gtf$chromosomes == "chr11"),]
gtf_chr12 <- gtf[which(gtf$chromosomes == "chr12"),]
gtf_chr13 <- gtf[which(gtf$chromosomes == "chr13"),]
gtf_chr14 <- gtf[which(gtf$chromosomes == "chr14"),]
gtf_chr15 <- gtf[which(gtf$chromosomes == "chr15"),]
gtf_chr16 <- gtf[which(gtf$chromosomes == "chr16"),]
gtf_chr17 <- gtf[which(gtf$chromosomes == "chr17"),]
gtf_chr18 <- gtf[which(gtf$chromosomes == "chr18"),]
gtf_chr19 <- gtf[which(gtf$chromosomes == "chr19"),]
gtf_chr20 <- gtf[which(gtf$chromosomes == "chr20"),]
gtf_chr21 <- gtf[which(gtf$chromosomes == "chr21"),]
gtf_chr22 <- gtf[which(gtf$chromosomes == "chr22"),]
gtf_chrM <- gtf[which(gtf$chromosomes == "chrM"),]
gtf_chrX <- gtf[which(gtf$chromosomes == "chrX"),]
gtf_chrY <- gtf[which(gtf$chromosomes == "chrY"),]

gtf_new <- rbind(gtf_chr1, gtf_chr2, gtf_chr3, gtf_chr4, gtf_chr5, gtf_chr6, gtf_chr7, gtf_chr8,
                 gtf_chr9, gtf_chr10, gtf_chr11, gtf_chr12, gtf_chr13, gtf_chr14, gtf_chr15,
                 gtf_chr16, gtf_chr17, gtf_chr18, gtf_chr19, gtf_chr20, gtf_chr21, gtf_chr22,
                 gtf_chrM, gtf_chrX, gtf_chrY)

dim(gtf_new)    # [1] 3 968 352       2

write.csv(x=gtf_new, file=paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/", "table_TranscriptsByChromosome_modified.csv"))

# difference in length when only main chromosomes kept
dim(gtf)[1] - dim(gtf_new)[1]   # [1] 163 095

# percentage:
((dim(gtf)[1] - dim(gtf_new)[1])/dim(gtf)[1]) * 100   # [1] 3.95% of lines excluded

# sort within chr only
#tab_1_chr1_sorted_pos <- tab_1_chr1[order(tab_1_chr1$positions),]
#tab_1_chr2_sorted_pos <- tab_1_chr2[order(tab_1_chr2$positions),]
#tab_1_chr3_sorted_pos <- tab_1_chr3[order(tab_1_chr3$positions),]
#tab_1_chr4_sorted_pos <- tab_1_chr4[order(tab_1_chr4$positions),]
#tab_1_chr5_sorted_pos <- tab_1_chr5[order(tab_1_chr5$positions),]
#tab_1_chr6_sorted_pos <- tab_1_chr6[order(tab_1_chr6$positions),]
#tab_1_chr7_sorted_pos <- tab_1_chr7[order(tab_1_chr7$positions),]
#tab_1_chr8_sorted_pos <- tab_1_chr8[order(tab_1_chr8$positions),]
#tab_1_chr9_sorted_pos <- tab_1_chr9[order(tab_1_chr9$positions),]
#tab_1_chr10_sorted_pos <- tab_1_chr10[order(tab_1_chr10$positions),]
#tab_1_chr11_sorted_pos <- tab_1_chr11[order(tab_1_chr11$positions),]
#tab_1_chr12_sorted_pos <- tab_1_chr12[order(tab_1_chr12$positions),]
#tab_1_chr13_sorted_pos <- tab_1_chr13[order(tab_1_chr13$positions),]
#tab_1_chr14_sorted_pos <- tab_1_chr14[order(tab_1_chr14$positions),]
#tab_1_chr15_sorted_pos <- tab_1_chr15[order(tab_1_chr15$positions),]
#tab_1_chr16_sorted_pos <- tab_1_chr16[order(tab_1_chr16$positions),]
#tab_1_chr17_sorted_pos <- tab_1_chr17[order(tab_1_chr17$positions),]
#tab_1_chr18_sorted_pos <- tab_1_chr18[order(tab_1_chr18$positions),]
#tab_1_chr19_sorted_pos <- tab_1_chr19[order(tab_1_chr19$positions),]
#tab_1_chr20_sorted_pos <- tab_1_chr20[order(tab_1_chr20$positions),]
#tab_1_chr21_sorted_pos <- tab_1_chr21[order(tab_1_chr21$positions),]
#tab_1_chr22_sorted_pos <- tab_1_chr22[order(tab_1_chr22$positions),]
#tab_1_chrM_sorted_pos <- tab_1_chrM[order(tab_1_chrM$positions),]
#tab_1_chrX_sorted_pos <- tab_1_chrX[order(tab_1_chrX$positions),]
#tab_1_chrY_sorted_pos <- tab_1_chrY[order(tab_1_chrY$positions),]

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

