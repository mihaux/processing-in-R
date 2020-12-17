# script to plot counting stats (from featureCounts)
# calculate average for ...
# create a barplot for ...

# install (if necessary) and load package
library(ggplot2)
library(stringr)
library(hrbrthemes)

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
df_count_se <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/6_featCounts/SE/all_stats_dups_run1_SE.csv"), row.names = 2)
df_count_pe <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/6_featCounts/PE/all_stats_dups_run1_PE.csv"), row.names = 2)

# check if rownames match between df_X_se and df_X_pe
any(rownames(df_count_se) == rownames(df_count_pe))

setwd(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/6_featCounts/plots"))

# ALL METRICS (only these that have results):
# Assigned
# Unassigned_NoFeatures

# remove first column from df_count_se and df_count_pe
df_count_se$X <- NULL
df_count_pe$X <- NULL

# correct colnames
colnames(df_count_se) <- str_replace(colnames(df_count_se), "X", "ID_")
colnames(df_count_pe) <- str_replace(colnames(df_count_pe), "X", "ID_")

colnames(df_count_se) <- str_replace(colnames(df_count_se), ".Aligned.sortedByCoord.out.bam", "")
colnames(df_count_pe) <- str_replace(colnames(df_count_pe), ".Aligned.sortedByCoord.out.bam", "")

# create a df with selected metrics only (numbers and percentages separately)
df_selected_se <- data.frame(assigned=unlist(c(df_count_se[1,])), 
                             unassigned=unlist(c(df_count_se[12,])),
                             sum_both=unlist(c(df_count_se[1,])) + unlist(c(df_count_se[12,])))

df_selected_pe <- data.frame(assigned=unlist(c(df_count_pe[1,])),
                             unassigned=unlist(c(df_count_pe[12,])),
                             sum_both=unlist(c(df_count_pe[1,])) + unlist(c(df_count_pe[12,])))

# calculate percentages
df_percent_se <- data.frame(percent_assigned=(df_selected_se$assigned/df_selected_se$sum_both),
                            percent_unassigned=(df_selected_se$unassigned/df_selected_se$sum_both))

df_percent_pe <- data.frame(percent_assigned=(df_selected_pe$assigned/df_selected_pe$sum_both),
                            percent_unassigned=(df_selected_pe$unassigned/df_selected_pe$sum_both))

# create a dataset
ids_se <- c(colnames(df_count_se), colnames(df_count_se))
condition_se <- c(rep("assigned", 41), rep("unassigned", 41))
percentage_se <- c(df_percent_se$percent_assigned, df_percent_se$percent_unassigned)
data_se <- data.frame(ids_se, condition_se, percentage_se)

ids_pe <- c(colnames(df_count_pe), colnames(df_count_pe))
condition_pe <- c(rep("assigned", 41), rep("unassigned", 41))
percentage_pe <- c(df_percent_pe$percent_assigned, df_percent_pe$percent_unassigned)
data_pe <- data.frame(ids_pe, condition_pe, percentage_pe)

# get average percentages of assigned and unassigned for SE and PE
ave_assigned_se <- mean(data_se$percentage_se[1:41])      # [1] 0.8293393
ave_unassigned_se <- mean(data_se$percentage_se[42:82])   # [1] 0.1706607

ave_assigned_pe <- mean(data_pe$percentage_pe[1:41])      # [1] 0.876267
ave_unassigned_pe <- mean(data_pe$percentage_pe[42:82])   # [1] 0.123733

png("se_assigne-unassigned.png")
ggplot(data_se, aes(fill=condition_se, y=percentage_se, x=ids_se)) + 
  geom_bar(position="fill", stat="identity") + coord_flip() + ggtitle("Single-end data") + xlab("sample IDs")
dev.off()

png("pe_assigne-unassigned.png")
ggplot(data_pe, aes(fill=condition_pe, y=percentage_pe, x=ids_pe)) + 
  geom_bar(position="fill", stat="identity") + coord_flip() + ggtitle("Paired-end data") + xlab("sample IDs")
dev.off()


