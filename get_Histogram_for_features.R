# this script creates histograms for continuous features to see their distribution

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

# args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
  stop("2 arguments must be supplied: 
       \n(1 - input) path to .csv with clinical features and, 
       \n(2 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3
args <- c(paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/histograms/"))

# Example of usage: 
# Rscript run_DESeq2.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/            

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[2], sep="\n")
setwd(args[2])

# load clinical data 
df_meta <- read.csv(args[1], row.names = 1, header = TRUE)

# add "ID_" to all rownames
rownames(df_meta) <- paste0("ID_", rownames(df_meta))

png(file = "histogram_nb_days_on_steroids.png")
hist(df_meta$number.of.days.on.steroids.at.TAB,
     xlab = "Number of days on steroids",
     col = "blue",
     xlim = c(0,20), 
     ylim = c(0,14),
     main = "Histogram of 'Number of days on steroids'")
dev.off()

png(file = "histogram_age.png")
hist(df_meta$age.at.BL,
     xlab = "Age",
     col = "gray",
     xlim = c(50,100), 
     ylim = c(0,13),
     main = "Histogram of 'Age'")
dev.off()

png(file = "histogram_year_TAB_sample_was_collected.png")
hist(df_meta$year.TAB.sample.was.collected,
     xlab = "year.TAB.sample.was.collected",
     col = "red",
     xlim = c(2004,2016), 
     ylim = c(0,20),
     main = "Histogram of 'year.TAB.sample.was.collected'")
dev.off()

png(file = "histogram_number_of_days_between_TAB._and_BL_blood_sample.png")
hist(as.numeric(df_meta$number.of.days.between.TAB.and.BL.blood.sample),
     xlab = "number.of.days.between.TAB.and.BL.blood.sample",
     col = "green",
     xlim = c(0,14), 
     ylim = c(0,10),
     main = "Histogram of 'number.of.days.between.TAB.and.BL.blood.sample'")
dev.off()

