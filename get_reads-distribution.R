# script to plot reads distribution and other plots for trimming stats (from Trimmomatic)

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

if (length(args)!=2) {
  stop("2 arguments must be supplied: 
       \n(1 - input) path to _.csv file with trimming stats, 
       \n(2 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/RNA-seq_pipeline_QC/2_trimming/trimming_stats_merged.csv"),
          paste0(main_dir, "/ANALYSES/RNA-seq_pipeline_QC/2_trimming/"))

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[2], sep="\n")
setwd(args[2])

# load normalised counts from DESeq
df <- read.csv(args[1], row.names = 1)

library(ggplot2)
# Basic barplot
p <- ggplot(data=df, aes(x= rownames(df), y=Input_Reads_SE)) + geom_bar(stat="identity", width=0.5)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Samples", y = "Number of reads obtained") +
  geom_hline(yintercept=median(mean(df$Input_Reads_SE) + 2*sd(df$Input_Reads_S)))
# geom_hline(yintercept=median(mean(df$Input_Reads_SE) + 2*sd(df$Input_Reads_S)))

# samples with smaller nb of reads than mean - 2 standard deviations => ID_8546
# samples greater number of reads than mean + 2 standard deviations => ID_16187 and ID_16142

pdf(file="plot_nb_of_reads_per_sample.pdf", width=9, height=5)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Samples", y = "Number of reads obtained")
# + geom_hline(yintercept=median(df$Input_Reads_SE))
dev.off()

