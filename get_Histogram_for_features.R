# this script creates histograms and barplots for continuous features to see their distribution

### plots included: ### 
# number.of.days.on.steroids.at.TAB (all / male only / female only)
# age (all / male only / female only)
# year.TAB.sample.was.collected (all / male only / female only)
# number.of.days.between.TAB.and.BL.blood.sample (all / male only / female only)

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

# define wheather output files should be saved or not [TRUE / FALSE]
output_save <- TRUE

# args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
  stop("2 arguments must be supplied: 
       \n(1 - input) path to .csv with clinical features and, 
       \n(2 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3
args <- c(paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/data/histograms/"))

# Example of usage: 
# Rscript run_DESeq2.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/            

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[2], sep="\n")
setwd(args[2])

# load clinical data 
df_meta <- read.csv(args[1], row.names = 1, header = TRUE)

# add "ID_" to all rownames
rownames(df_meta) <- paste0("ID_", rownames(df_meta))

# split dataset by gender
df_male <- df_meta[which(df_meta$gender == 1),]
df_female <- df_meta[which(df_meta$gender == 2),]

################################################
### print how many patients in each variable ###
#for (i in sort(unique(df_meta$age.at.BL))) {
#  cat("Age",i)
#  cat(" Male:", length(which(df_male$age.at.BL == i)))
#  cat(" Female:", length(which(df_female$age.at.BL == i)),"\n")
#}
################################################

### number.of.days.on.steroids.at.TAB ###
if(output_save==TRUE){ png(file = "histogram_nb_days_on_steroids_all.png") }
hist(df_meta$number.of.days.on.steroids.at.TAB,
     col = "blue", xlim = c(0,20), ylim = c(0,14), 
     xlab = "Number of days on steroids",
     main = "Histogram of 'Number of days on steroids' - all data")
if(output_save==TRUE){ dev.off() }

# male
if(output_save==TRUE){ png(file = "histogram_nb_days_on_steroids_male.png") }
hist(df_male$number.of.days.on.steroids.at.TAB,
     col = "blue", xlim = c(0,20), ylim = c(0,10), 
     xlab = "Number of days on steroids",
     main = "Histogram of 'Number of days on steroids' - male only")
if(output_save==TRUE){ dev.off() }

# female
if(output_save==TRUE){ png(file = "histogram_nb_days_on_steroids_female.png") }
hist(df_female$number.of.days.on.steroids.at.TAB,
     col = "blue", xlim = c(0,20), ylim = c(0,10), 
     xlab = "Number of days on steroids",
     main = "Histogram of 'Number of days on steroids' - female only")
if(output_save==TRUE){ dev.off() }

### age ###
if(output_save==TRUE){ png(file = "histogram_age_all.png") }
hist(df_meta$age.at.BL,
     col = "gray", xlim = c(50,100), ylim = c(0,13), 
     xlab = "Age",
     main = "Histogram of 'Age' - all data")
if(output_save==TRUE){ dev.off() }

# male
if(output_save==TRUE){ png(file = "histogram_age_male.png") }
hist(df_male$age.at.BL,
     col = "gray", xlim = c(50,100), ylim = c(0,13), 
     xlab = "Age",
     main = "Histogram of 'Age' - male only")
if(output_save==TRUE){ dev.off() }

# female
if(output_save==TRUE){ png(file = "histogram_age_female.png") }
hist(df_female$age.at.BL,
     col = "gray", xlim = c(50,100), ylim = c(0,13), 
     xlab = "Age",
     main = "Histogram of 'Age' - female only")
if(output_save==TRUE){ dev.off() }


### year.TAB.sample.was.collected ###
if(output_save==TRUE){ png(file = "histogram_year_TAB_sample_was_collected_all.png") }
hist(df_meta$year.TAB.sample.was.collected,
     col = "red", xlim = c(2004,2016), ylim = c(0,20),
     xlab = "year.TAB.sample.was.collected",
     main = "Histogram of 'year.TAB.sample.was.collected' - all data")
if(output_save==TRUE){ dev.off() }

# male
if(output_save==TRUE){ png(file = "histogram_year_TAB_sample_was_collected_male.png") }
hist(df_male$year.TAB.sample.was.collected,
     col = "red", xlim = c(2004,2016), ylim = c(0,20),
     xlab = "year.TAB.sample.was.collected",
     main = "Histogram of 'year.TAB.sample.was.collected' - male only")
if(output_save==TRUE){ dev.off() }

# female
if(output_save==TRUE){ png(file = "histogram_year_TAB_sample_was_collected_female.png") }
hist(df_female$year.TAB.sample.was.collected,
     col = "red", xlim = c(2004,2016), ylim = c(0,20),
     xlab = "year.TAB.sample.was.collected",
     main = "Histogram of 'year.TAB.sample.was.collected' - female only")
if(output_save==TRUE){ dev.off() }


### number.of.days.between.TAB.and.BL.blood.sample ###
if(output_save==TRUE){ png(file = "histogram_number_of_days_between_TAB._and_BL_blood_sample_all.png") }
hist(as.numeric(df_meta$number.of.days.between.TAB.and.BL.blood.sample),
     col = "green", xlim = c(0,14), ylim = c(0,10),
     xlab = "number.of.days.between.TAB.and.BL.blood.sample",
     main = "Histogram of 'number.of.days.between.TAB.and.BL.blood.sample' \n - all data")
if(output_save==TRUE){ dev.off() }

# male
if(output_save==TRUE){ png(file = "histogram_number_of_days_between_TAB._and_BL_blood_sample_male.png") }
hist(as.numeric(df_male$number.of.days.between.TAB.and.BL.blood.sample),
     col = "green", xlim = c(0,14), ylim = c(0,10),
     xlab = "number.of.days.between.TAB.and.BL.blood.sample",
     main = "Histogram of 'number.of.days.between.TAB.and.BL.blood.sample' \n - male only")
if(output_save==TRUE){ dev.off() }

# female
if(output_save==TRUE){ png(file = "histogram_number_of_days_between_TAB._and_BL_blood_sample_female.png") }
hist(as.numeric(df_female$number.of.days.between.TAB.and.BL.blood.sample),
     col = "green", xlim = c(0,14), ylim = c(0,10),
     xlab = "number.of.days.between.TAB.and.BL.blood.sample",
     main = "Histogram of 'number.of.days.between.TAB.and.BL.blood.sample' \n - female only")
if(output_save==TRUE){ dev.off() }

