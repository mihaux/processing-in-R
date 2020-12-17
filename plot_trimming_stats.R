# script to plot trimming stats (from Trimmomatic)
# calculate average for SE: "Input.Reads." | "Surviving." | "Dropped."    
# calculate average for PE: "Input.Read.Pairs." | "Both.Surviving." | "Forward.Only.Surviving." | "Reverse.Only.Surviving." "Dropped."   

# create a barplot for "Input.Reads." and "Input.Read.Pairs." (only one, as they are both the same)
# create a barplot for "Surviving." (SE) and "Both.Surviving." (PE)
# create a barplot for "Dropped." (SE) and "Dropped." (PE)

# install (if necessary) and load package
library(ggplot2)
library(stringr)
library(tidyr)

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
df_trim_se <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/2_trimming/trimming_stats_SE.csv"), row.names = 1)
df_trim_pe <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/2_trimming/trimming_stats_PE.csv"), row.names = 1)

# check if rownames match between df_X_se and df_X_pe
any(rownames(df_trim_se) == rownames(df_trim_pe))

setwd(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/2_trimming/plots"))

# split columns with percentages
split_se_Surviving <- str_split_fixed(as.character(df_trim_se$Surviving.),' ', n=2)
split_se_Surviving[,2] <- str_replace_all(split_se_Surviving[,2], "[(%)]", "")

split_se_Dropped <- str_split_fixed(as.character(df_trim_se$Dropped.),' ', n=2)
split_se_Dropped[,2] <- str_replace_all(split_se_Dropped[,2], "[(%)]", "")

split_pe_Both.Surviving <- str_split_fixed(as.character(df_trim_pe$Both.Surviving.),' ', n=2)
split_pe_Both.Surviving[,2] <- str_replace_all(split_pe_Both.Surviving[,2], "[(%)]", "")

split_pe_Forward <- str_split_fixed(as.character(df_trim_pe$Forward.Only.Surviving.),' ', n=2)
split_pe_Forward[,2] <- str_replace_all(split_pe_Forward[,2], "[(%)]", "")

split_pe_Reverse <- str_split_fixed(as.character(df_trim_pe$Reverse.Only.Surviving.),' ', n=2)
split_pe_Reverse[,2] <- str_replace_all(split_pe_Reverse[,2], "[(%)]", "")

split_pe_Dropped <- str_split_fixed(as.character(df_trim_pe$Dropped.),' ', n=2)
split_pe_Dropped[,2] <- str_replace_all(split_pe_Dropped[,2], "[(%)]", "")

# calculate mean values and save them as a table
mean(df_trim_se$Input.Reads.)
mean(as.numeric(split_se_Surviving[,2]))
mean(as.numeric(split_se_Dropped[,2]))

mean(df_trim_pe$Input.Read.Pairs.)
mean(as.numeric(split_pe_Both.Surviving[,2]))
mean(as.numeric(split_pe_Forward[,2]))
mean(as.numeric(split_pe_Reverse[,2]))
mean(as.numeric(split_pe_Dropped[,2]))

tab_summary <- data.frame(SE=c(round(mean(df_trim_se$Input.Reads.)), round(mean(as.numeric(split_se_Surviving[,2])), 2), "NA", "NA", round(mean(as.numeric(split_se_Dropped[,2])), 2)),
                          PE=c(round(mean(df_trim_pe$Input.Read.Pairs.)), round(mean(as.numeric(split_pe_Both.Surviving[,2])),2), round(mean(as.numeric(split_pe_Forward[,2])),2), round(mean(as.numeric(split_pe_Reverse[,2])),2), round(mean(as.numeric(split_pe_Dropped[,2])),2)))

rownames(tab_summary) <- c("Total number of reads",
                           "Read_Surviving_%",
                           "Forward_reads_surviving_%",
                           "Reverse_reads_surviving_%",
                           "Dropped_%")

write.csv(tab_summary, "trimming_summary_mean.csv")

# create a barplot for total number of obtained reads (same for SE and PE)
df_nb_of_reads <- data.frame(id=unlist(lapply(str_replace_all(rownames(df_trim_se), "_R1", ""), function(x) paste0("ID_", x))),
                             nb=df_trim_se$Input.Reads.)

png("total_number_of_reads.png")
ggplot(df_nb_of_reads, aes(x=id, y=nb)) + geom_bar(stat="identity") + coord_flip() + xlab("sample ID") + ylab("Number of reads")
dev.off()


# create a barplot for "Surviving." (SE) and "Both.Surviving." (PE)
dat_surviving <- data.frame(surviving_se=c(split_se_Surviving[,2]),
                            surviving_pe=c(split_pe_Both.Surviving[,2]),
                            ID = as.factor(c(unlist(lapply(str_replace_all(rownames(df_trim_se), "_R1", ""), function(x) paste0("ID_", x))))))

dat_surviving_long <- dat_surviving %>%
  gather("Stat", "Value", -ID)

# switch values to numeric
dat_surviving_long$Value <- as.numeric(dat_surviving_long$Value)

png("barplot_surviving_se_pe.png")
ggplot(dat_surviving_long, aes(x = ID, y = Value, fill = Stat)) + geom_col(position = "dodge") + coord_flip()
dev.off()

# create a barplot for "Dropped." (SE) and "Dropped." (PE)
dat_dropped <- data.frame(dropped_se = c(split_se_Dropped[,2]),
                          dropped_pe = c(split_pe_Dropped[,2]),
                          ID = as.factor(c(unlist(lapply(str_replace_all(rownames(df_trim_se), "_R1", ""), function(x) paste0("ID_", x))))))

dat_dropped_long <- dat_dropped %>%
  gather("Stat", "Value", -ID)

# switch values to numeric
dat_dropped_long$Value <- as.numeric(dat_dropped_long$Value)

png("barplot_dropped_se_pe.png")
ggplot(dat_dropped_long, aes(x = ID, y = Value, fill = Stat)) + geom_col(position = "dodge") + coord_flip()
dev.off()

# plot SE against PE
df_se_vs_pe <- data.frame(se=as.numeric(split_se_Dropped[,2]),
                          pe=as.numeric(split_pe_Dropped[,2]))

png("dropped_reads_SE-vs-PE.png")
ggplot(df_se_vs_pe, aes(se, pe)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentages of dropped reads - SE vs. PE") +
  xlim(c(0,3)) + ylim(c(0,3))
dev.off()







