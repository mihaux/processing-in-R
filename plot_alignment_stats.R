# script to plot alignment stats (from STAR)
# calculate average for 8 selected metrics
# create a barplot for 8 selected metrics

# install (if necessary) and load package
library(ggplot2)
library(stringr)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"      # on uni mac    
} else {
  print("Unrecognised machine.")
}

# load table with trimming stats
df_align_se <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/4_alignment/SE/STAR_stats_all_SE.csv"), row.names = 1)
df_align_pe <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/4_alignment/PE/STAR_stats_all_PE.csv"), row.names = 1)

# check if rownames match between df_X_se and df_X_pe
any(rownames(df_align_se) == rownames(df_align_pe))

setwd(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/4_alignment/plots"))

# ALL METRICS (only important were selected):
# "Uniquely mapped reads number"                                    | "Uniquely mapped reads %"
#  NA                                                               | "Mismatch rate per base, %"                    
# "Number of reads mapped to multiple loci"                         | "% of reads mapped to multiple loci"           
# "Number of reads mapped to too many loci"                         | "% of reads mapped to too many loci"           
# [0 - all samples] "Number of reads unmapped: too many mismatches" | "% of reads unmapped: too many mismatches"
# "Number of reads unmapped: too short"                             | "% of reads unmapped: too short"               
# "Number of reads unmapped: other"                                 | "% of reads unmapped: other"                   
# [0 - all samples] "Number of chimeric reads"                      | "% of chimeric reads"  

# create a df with selected metrics only (numbers and percentages separately)
df_selected_nb_se <- data.frame(a=as.numeric(levels(droplevels(unlist(matrix(df_align_se[5,]))))),
                                b=as.numeric(levels(droplevels(unlist(matrix(df_align_se[19,]))))),
                                c=as.numeric(levels(droplevels(unlist(matrix(df_align_se[21,]))))),
                                d=as.numeric(levels(droplevels(unlist(matrix(df_align_se[23,]))))),
                                e=as.numeric(levels(droplevels(unlist(matrix(df_align_se[25,]))))),
                                f=as.numeric(levels(droplevels(unlist(matrix(df_align_se[27,]))))),
                                g=as.numeric(levels(droplevels(unlist(matrix(df_align_se[29,]))))))

df_selected_nb_pe <- data.frame(a=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[5,]))))),
                                b=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[19,]))))),
                                c=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[21,]))))),
                                d=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[23,]))))),
                                e=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[25,]))))),
                                f=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[27,]))))),
                                g=as.numeric(levels(droplevels(unlist(matrix(df_align_pe[29,]))))))


names_numbers <- c("Uniquely mapped reads number",
                   "Number of reads mapped to multiple loci",      
                   "Number of reads mapped to too many loci",
                   "Number of reads unmapped: too many mismatches",
                   "Number of reads unmapped: too short",      
                   "Number of reads unmapped: other",           
                   "Number of chimeric reads")

names(df_selected_nb_se) <- names_numbers 
names(df_selected_nb_pe) <- names_numbers 
  
df_selected_percent_se <- data.frame(a=as.character(unlist(matrix(df_align_se[6,]))),
                                     b=as.character(unlist(matrix(df_align_se[14,]))),
                                     c=as.character(unlist(matrix(df_align_se[20,]))),
                                     d=as.character(unlist(matrix(df_align_se[22,]))),
                                     e=as.character(unlist(matrix(df_align_se[24,]))),
                                     f=as.character(unlist(matrix(df_align_se[26,]))),
                                     g=as.character(unlist(matrix(df_align_se[28,]))),
                                     h=as.character(unlist(matrix(df_align_se[30,]))))

df_selected_percent_pe <- data.frame(a=as.character(unlist(matrix(df_align_pe[6,]))),
                                     b=as.character(unlist(matrix(df_align_pe[14,]))),
                                     c=as.character(unlist(matrix(df_align_pe[20,]))),
                                     d=as.character(unlist(matrix(df_align_pe[22,]))),
                                     e=as.character(unlist(matrix(df_align_pe[24,]))),
                                     f=as.character(unlist(matrix(df_align_pe[26,]))),
                                     g=as.character(unlist(matrix(df_align_pe[28,]))),
                                     h=as.character(unlist(matrix(df_align_pe[30,]))))

names_percent <- c("Uniquely mapped reads %", 
                   "Mismatch rate per base, %",
                   "% of reads mapped to multiple loci",         
                   "% of reads mapped to too many loci",           
                   "% of reads unmapped: too many mismatches",     
                   "% of reads unmapped: too short",               
                   "% of reads unmapped: other",                   
                   "% of chimeric reads")

names(df_selected_percent_se) <- names_percent
names(df_selected_percent_pe) <- names_percent

# calculate averages for SE and PE and save as table
tab_summary_se_pe <- data.frame(single_end=c(round(mean(as.numeric(str_replace_all(df_selected_percent_se[,1], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,2], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,3], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,4], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,5], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,6], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,7], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_se[,8], "[%]", ""))), digits = 2)),
                                paired_end=c(round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,1], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,2], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,3], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,4], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,5], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,6], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,7], "[%]", ""))), digits = 2),
                                             round(mean(as.numeric(str_replace_all(df_selected_percent_pe[,8], "[%]", ""))), digits = 2)))

rownames(tab_summary_se_pe) <- c("Uniquely mapped reads %", 
                                 "Mismatch rate per base, %",
                                 "% of reads mapped to multiple loci",         
                                 "% of reads mapped to too many loci",           
                                 "% of reads unmapped: too many mismatches",     
                                 "% of reads unmapped: too short",               
                                 "% of reads unmapped: other",                   
                                 "% of chimeric reads")

write.csv(tab_summary_se_pe, "alignment_percentages_se_pe.csv")

# create boxplots to compare percentages
df_boxplot_uniquely <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                  value=c(as.numeric(str_replace_all(df_selected_percent_se[,1], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,1], "[%]", ""))))

df_boxplot_mismatch <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                  value=c(as.numeric(str_replace_all(df_selected_percent_se[,2], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,2], "[%]", ""))))

df_boxplot_multiple <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                  value=c(as.numeric(str_replace_all(df_selected_percent_se[,3], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,3], "[%]", ""))))

df_boxplot_toomanyloci <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                     value=c(as.numeric(str_replace_all(df_selected_percent_se[,4], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,4], "[%]", ""))))

df_boxplot_toomanymismatches <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                           value=c(as.numeric(str_replace_all(df_selected_percent_se[,5], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,5], "[%]", ""))))

df_boxplot_tooshort <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                  value=c(as.numeric(str_replace_all(df_selected_percent_se[,6], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,6], "[%]", ""))))

df_boxplot_other <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                               value=c(as.numeric(str_replace_all(df_selected_percent_se[,7], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,7], "[%]", ""))))

df_boxplot_chimeric <- data.frame(name=c(rep("SE", 41), rep("PE", 41)),
                                  value=c(as.numeric(str_replace_all(df_selected_percent_se[,8], "[%]", "")), as.numeric(str_replace_all(df_selected_percent_pe[,8], "[%]", ""))))

png("uniquely_mapped_reads.png")
df_boxplot_uniquely %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("Uniquely mapped reads %") + xlab("")
dev.off()

png("mismatch_rate_per_base.png")
df_boxplot_mismatch %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("Mismatch rate per base, %") + xlab("")
dev.off()

png("reads_mapped_to_multiple_loci.png")
df_boxplot_multiple %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of reads mapped to multiple loci") + xlab("")
dev.off()

png("reads_mapped_to_too_many_loci.png")
df_boxplot_toomanyloci %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of reads mapped to too many loci") + xlab("")
dev.off()

png("reads_unmapped_too_many_mismatches.png")
df_boxplot_toomanymismatches %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of reads unmapped: too many mismatches") + xlab("")
dev.off()

png("reads_unmapped_too_short.png")
df_boxplot_tooshort %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of reads unmapped: too short") + xlab("")
dev.off()

png("reads_unmapped_other.png")
df_boxplot_other %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of reads unmapped: other" ) + xlab("")
dev.off()

png("chimeric_reads.png")
df_boxplot_chimeric %>%
  ggplot(aes(x=name, y=value, fill=name)) + geom_boxplot() + ggtitle("% of chimeric reads") + xlab("")
dev.off()

