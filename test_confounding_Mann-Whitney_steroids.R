# this script checks the confounding influence of 'number_of_days_on_steroids'
# => using Mann-Whitney test (split the data in 2 approximately equal groups in terms of the number of patients to get comparison groups for Mann-Whitney test using gene expression data as input)
# => ADDED LATER: split the data in 2 but use the quartiles Q1 vs Q4

# GENERAL NOTE ABOUT CONFOUNDING:
# Confounding assessment was performed for the following features: age, gender and duration of steroid treatment to make sure that they do not cause any distortion in the association between the exposure and the outcome (e.g. visual loss). 
# To do so, a paired samples t-test was run for each of these features to see if they were associated with visual loss, which could suggest possible confounding.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(ggplot2)

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

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to the file with counts matrix
       \n(2 - metadata) path to .csv with clinical features and, 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3
args <- c(paste0(main_dir,"/data/count_matrices/outliers_excluded/counts_rlog_no_outliers.csv"),
          paste0(main_dir, "/data/metadata/outliers_excluded/cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv"),
          paste0(main_dir, "/ANALYSES/Dec20_steroids_age_rerun/steroids/rlog/"))

# Example of usage: 
# Rscript test_confounding_Mann-Whitney_steroids.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load counts matrix
dds <- read.csv(args[1], row.names = 1, header = TRUE)

# define running ID (either "raw" / "vst" / "rlog")
run_id <- "vst"

# load clinical data 
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)

# add "ID_" to all rownames
#rownames(df_meta) <- paste0("ID_", rownames(df_meta))
  
# change data format to matrix and integer (PROBABLY NOT NEEDED)
dat <- as.matrix(dds)

# add 0.0001 - 0.0041 to columns (not sure this works as expected)
#to_be_added <- seq(0.0001,0.0041, 0.0001)
#for(i in 1:40) {
#    dat[,i] <- dat[,i]+to_be_added[i]
#}

#dat <- as.matrix(assay(dds))
# storage.mode(dat) <- "integer"          # class: matrix | type: integer

###########################################################################################
# => (1) using Mann-Whitney test (split the data in 2 approximately equal groups in terms of the number of patients to get comparison groups for Mann-Whitney test using gene expression data as input)
###########################################################################################
# df_meta$number.of.days.on.steroids.at.TAB
# NOTE: there are 2 samples with 16 days, can consider them as too extreme and exclude
# split as:
# group_1 => 0 - 5 days
# group_2 => 6 - 16 days
# group_2_bis => 6 - 11 days (when both samples with 16 days are excluded)

gr_1 <- df_meta[which(df_meta$number.of.days.on.steroids.at.TAB < 6),]
gr_2 <- df_meta[which(df_meta$number.of.days.on.steroids.at.TAB >= 6),]

gr_2_bis <- gr_2[-which(gr_2$number.of.days.on.steroids.at.TAB == 16),]

# add new column to define comparison groups 1 and 2
steroids <- c(1:40)
steroids[which(df_meta$number.of.days.on.steroids.at.TAB < 6)] <- 1
steroids[which(df_meta$number.of.days.on.steroids.at.TAB >= 6)] <- 2

df_meta$number.of.days.on.steroids.at.TAB <- steroids

# define groups for comparison
group = df_meta$number.of.days.on.steroids.at.TAB

# initiate list for results
res <- list()
for(j in 1:nrow(dat)){
  
  # get data frame for one gene at a time
  temp = dat[j, ]
  
  # perform statistical testing
  res[[j]] = wilcox.test(temp[group==1], temp[group==2],
                         alternative = "two.sided", exact = TRUE)  
  
  # exact = FALSE => to suppress the warning of “cannot compute exact p-value with tie”
  # which is caused by the assumption of a Wilcoxon test that the responses are continuous. 
}

# retrieve p-value for each gene (transcript)
all_pval = unlist(lapply(res, function(x) x$p.value))

# Multiplicity adjustment
# After obtaining the p-values, you can make multiplicity correction to the pvalue. 
all_adjusted <- p.adjust(all_pval, "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

# create a table as R data frame
res_table = data.frame(ID=rownames(dat), 
                       p.value=all_pval, 
                       fdr.pvalue=all_adjusted)

if(FALSE){
# check other methods
all_adjusted_holm       <- p.adjust(all_pval, "holm")       # "holm"
all_adjusted_hochberg   <- p.adjust(all_pval, "hochberg")   # "hochberg"
all_adjusted_hommel     <- p.adjust(all_pval, "hommel")     # "hommel"
all_adjusted_bonferroni <- p.adjust(all_pval, "bonferroni") # "bonferroni"
all_adjusted_by         <- p.adjust(all_pval, "BY")         # "BY"
}

### histogram of p-values ###
# create a df with IDs and p-values
df_pval <- data.frame(ID=rownames(dat), pval=all_pval)

# make a histogram of pvalues
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_", run_id, ".png")) }
ggplot(df_pval, aes(x=pval)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_20bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_50bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# NOTE: all the p-adjusted values equal 1

# create a table as R data frame
res_table = data.frame(ID=rownames(dat), 
                       pvalue=all_pval, 
                       fdr.pvalue=all_adjusted)

# sort the table by p-values
res_table_sorted <- res_table[order(all_pval),]

# sort the table by p-adjusted
#res_table_sorted_bis <- res_table[order(all_adjusted),]

# save tables
write.csv2(res_table_sorted, file=paste0("table_sorted_pvalues_", run_id, ".csv"))

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval <- length(which(res_table_sorted$pvalue < 0.05))
significant_padjusted <- length(which(res_table_sorted$fdr.pvalue < 0.05))

# create a histogram for significant 
# p-values
if(output_save==TRUE){ png(file = paste0("histogram_p-value_significant_only_", run_id, ".png")) }
ggplot(res_table_sorted[1:752,], aes(x=pvalue)) + 
  geom_histogram(binwidth=0.0001) + 
  labs(title=paste0("Histogram of significant p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# p-adjusted
if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_significant_only_", run_id, ".png")) }
ggplot(res_table_sorted[1:752,], aes(x=fdr.pvalue)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of significant p-adjusted (fdr): ", run_id), x="p-ajdusted")
if(output_save==TRUE){ dev.off() }

# check correlation between p-val and p-adjusted
if(output_save==TRUE){ png(file = paste0("correlation_pval-padj_", run_id, ".png")) }
ggplot(res_table, aes(pvalue, fdr.pvalue)) + geom_point()
if(output_save==TRUE){ dev.off() }

# NOTES:
# total length => 26 486
# nb of NA values => 5 121 (same genes in p-value and p-adj)

########################################################
### (2) same run, but with 2 samples of 16 days excluded ###
########################################################

# need to reload meta data
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)
rownames(df_meta) <- paste0("ID_", rownames(df_meta))

# exclude these 2 samples from count datasets and clinical spreadsheet
dat_excluded <- dat[,-which(gr_2$number.of.days.on.steroids.at.TAB == 16)]
df_meta_excluded <- df_meta[-which(gr_2$number.of.days.on.steroids.at.TAB == 16),]

# add new column to define comparison groups 1 and 2
steroids_excluded <- c(1:39)
steroids_excluded[which(df_meta_excluded$number.of.days.on.steroids.at.TAB < 6)] <- 1
steroids_excluded[which(df_meta_excluded$number.of.days.on.steroids.at.TAB >= 6)] <- 2

df_meta_excluded$number.of.days.on.steroids.at.TAB <- steroids_excluded

# define groups for comparison
group_excluded = df_meta_excluded$number.of.days.on.steroids.at.TAB

# initiate list for results
res_excluded <- list()
for(j in 1:nrow(dat_excluded)){
  
  # get data frame for one gene at a time
  temp_excluded = dat_excluded[j, ]
  
  # perform statistical testing
  res_excluded[[j]] = wilcox.test(temp_excluded[group_excluded==1], temp_excluded[group_excluded==2], 
                             alternative = "two.sided", exact = FALSE)               
  
  # exact = FALSE => to suppress the warning of “cannot compute exact p-value with tie”
  # which is caused by the assumption of a Wilcoxon test that the responses are continuous. 
}

# retrieve p-value for each gene (transcript)
all_pval_excluded = unlist(lapply(res_excluded, function(x) x$p.value))

# get corrected p-values
all_adjusted_excluded <- p.adjust(all_pval_excluded, "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

### histogram of p-values ###
# create a df with IDs and p-values
df_pval_excluded <- data.frame(ID=rownames(dat_excluded), pval=all_pval_excluded)

# make a histogram of pvalues
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_excluded_", run_id, ".png")) }
ggplot(df_pval_excluded, aes(x=pval)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# NOTE: all the p-adjusted values equal 1

# create a table as R data frame
res_table_excluded = data.frame(ID=rownames(dat_excluded), 
                                pvalue=all_pval_excluded, 
                                fdr.pvalue=all_adjusted_excluded)

# sort the table by p-values
res_table_sorted_excluded <- res_table[order(all_pval_excluded),]

# sort the table by p-adjusted
#res_table_sorted_bis <- res_table[order(all_adjusted),]

# save tables
if(output_save==TRUE){  
write.csv2(res_table_sorted_excluded, file=paste0("table_sorted_pvalues_excluded_", run_id, ".csv"))
}

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval_excluded <- length(which(res_table_sorted_excluded$pvalue < 0.05))
significant_padjusted_excluded <- length(which(res_table_sorted_excluded$fdr.pvalue < 0.05))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
########################################################
### (3) same run, but with ID_8546 sample excluded ###
########################################################

# need to reload meta data
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)
rownames(df_meta) <- paste0("ID_", rownames(df_meta))

# exclude these ID_8546 from count datasets and clinical spreadsheet
df_meta_no_outliers <- df_meta[-which(rownames(df_meta) == "ID_8546"),]
dat_no_outliers <- dat[,-which(colnames(dat) == "ID_8546")]

gr_1_no_outliers <- df_meta_no_outliers[which(df_meta_no_outliers$number.of.days.on.steroids.at.TAB < 6),]
gr_2_no_outliers <- df_meta_no_outliers[which(df_meta_no_outliers$number.of.days.on.steroids.at.TAB >= 6),]

# add new column to define comparison groups 1 and 2
steroids <- c(1:40)
steroids[which(df_meta_no_outliers$number.of.days.on.steroids.at.TAB < 6)] <- 1
steroids[which(df_meta_no_outliers$number.of.days.on.steroids.at.TAB >= 6)] <- 2

# define groups for comparison
group_no_outliers = steroids

# initiate list for results
res_no_outliers <- list()
for(j in 1:nrow(dat_no_outliers)){
  
  # get data frame for one gene at a time
  temp_no_outliers = dat_no_outliers[j, ]
  
  # perform statistical testing
  res_no_outliers[[j]] = wilcox.test(temp_no_outliers[group_no_outliers==1], temp_no_outliers[group_no_outliers==2], 
                                  alternative = "two.sided", exact = FALSE)               
  
  # exact = FALSE => to suppress the warning of “cannot compute exact p-value with tie”
  # which is caused by the assumption of a Wilcoxon test that the responses are continuous. 
}

# retrieve p-value for each gene (transcript)
all_pval_no_outliers = unlist(lapply(res_no_outliers, function(x) x$p.value))

# get corrected p-values
all_adjusted_no_outliers <- p.adjust(all_pval_no_outliers, "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

### histogram of p-values ###
# create a df with IDs and p-values
df_pval_no_outliers <- data.frame(ID=rownames(dat_no_outliers), pval=all_pval_no_outliers)

# make a histogram of pvalues
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_no_outliers_", run_id, ".png")) }
ggplot(df_pval_no_outliers, aes(x=pval)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# NOTE: all the p-adjusted values equal 1

# create a table as R data frame
res_table_no_outliers = data.frame(ID=rownames(dat_no_outliers), 
                                   pvalue=all_pval_no_outliers, 
                                   fdr.pvalue=all_adjusted_no_outliers)

# sort the table by p-values
res_table_sorted_no_outliers <- res_table_no_outliers[order(all_pval_no_outliers),]

# sort the table by p-adjusted
#res_table_sorted_bis <- res_table[order(all_adjusted),]

# save tables
if(output_save==TRUE){  
  write.csv2(res_table_sorted_no_outliers, file=paste0("table_sorted_pvalues_no_outliers_", run_id, ".csv"))
}

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval_no_outliers <- length(which(res_table_sorted_no_outliers$pvalue < 0.05))
significant_padjusted_no_outliers <- length(which(res_table_sorted_no_outliers$fdr.pvalue < 0.05))



###########################################################################################
# => (4) using Mann-Whitney test (Q1 vs Q3)
###########################################################################################

quantile(df_meta$number.of.days.on.steroids.at.TAB)
#  0%  25%  50%  75% 100% 
#  0    3    6    7   16 

# group_1 => 0 - 3 days
length(which(df_meta$number.of.days.on.steroids.at.TAB <= 3)) # 12 samples

# group_2 => 7 - 16 days
length(which(df_meta$number.of.days.on.steroids.at.TAB >= 7)) # 14 samples

# create a subset of df_meta that will include only samples of interest (12 + 15)
dat_new <- dat[,c(which(df_meta$number.of.days.on.steroids.at.TAB <= 3), which(df_meta$number.of.days.on.steroids.at.TAB >= 7))]
df_meta_new <- df_meta[c(which(df_meta$number.of.days.on.steroids.at.TAB <= 3), which(df_meta$number.of.days.on.steroids.at.TAB >= 7)),]

any(colnames(dat_new) == rownames(df_meta_new))

gr_1 <- df_meta_new[which(df_meta_new$number.of.days.on.steroids.at.TAB <= 3),]
gr_2 <- df_meta_new[which(df_meta_new$number.of.days.on.steroids.at.TAB >= 7),]

# add new column to define comparison groups 1 and 2
steroids <- c(1:26)
steroids[which(df_meta_new$number.of.days.on.steroids.at.TAB <= 3)] <- 1
steroids[which(df_meta_new$number.of.days.on.steroids.at.TAB >= 7)] <- 2

df_meta_new$number.of.days.on.steroids.at.TAB <- steroids

# define groups for comparison
group = df_meta_new$number.of.days.on.steroids.at.TAB

# initiate list for results
res <- list()
for(j in 1:nrow(dat_new)){
  
  # get data frame for one gene at a time
  temp = dat_new[j, ]
  
  # perform statistical testing
  res[[j]] = wilcox.test(temp[group==1], temp[group==2],
                         alternative = "two.sided", exact = FALSE)  
  
  # exact = FALSE => to suppress the warning of “cannot compute exact p-value with tie”
  # which is caused by the assumption of a Wilcoxon test that the responses are continuous. 
}

# retrieve p-value for each gene (transcript)
all_pval = unlist(lapply(res, function(x) x$p.value))

# Multiplicity adjustment
# After obtaining the p-values, you can make multiplicity correction to the pvalue. 
all_adjusted <- p.adjust(all_pval, "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

# create a table as R data frame
res_table = data.frame(ID=rownames(dat), 
                       p.value=all_pval, 
                       fdr.pvalue=all_adjusted)

if(FALSE){
  # check other methods
  all_adjusted_holm       <- p.adjust(all_pval, "holm")       # "holm"
  all_adjusted_hochberg   <- p.adjust(all_pval, "hochberg")   # "hochberg"
  all_adjusted_hommel     <- p.adjust(all_pval, "hommel")     # "hommel"
  all_adjusted_bonferroni <- p.adjust(all_pval, "bonferroni") # "bonferroni"
  all_adjusted_by         <- p.adjust(all_pval, "BY")         # "BY"
}

### histogram of p-values ###
# create a df with IDs and p-values
df_pval <- data.frame(ID=rownames(dat_new), pval=all_pval)

# make a histogram of pvalues
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_", run_id, ".png")) }
ggplot(df_pval, aes(x=pval)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_20bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_50bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# NOTE: all the p-adjusted values equal 1

# create a table as R data frame
res_table = data.frame(ID=rownames(dat_new), 
                       pvalue=all_pval, 
                       fdr.pvalue=all_adjusted)

# sort the table by p-values
res_table_sorted <- res_table[order(all_pval),]

# sort the table by p-adjusted
#res_table_sorted_bis <- res_table[order(all_adjusted),]

# save tables
write.csv2(res_table_sorted, file=paste0("table_sorted_pvalues_", run_id, ".csv"))

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval <- length(which(res_table_sorted$pvalue < 0.05))
significant_padjusted <- length(which(res_table_sorted$fdr.pvalue < 0.05))

# create a histogram for significant 
# p-values
if(output_save==TRUE){ png(file = paste0("histogram_p-value_significant_only_", run_id, ".png")) }
ggplot(res_table_sorted[1:752,], aes(x=pvalue)) + 
  geom_histogram(binwidth=0.0001) + 
  labs(title=paste0("Histogram of significant p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# p-adjusted
if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_significant_only_", run_id, ".png")) }
ggplot(res_table_sorted[1:752,], aes(x=fdr.pvalue)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of significant p-adjusted (fdr): ", run_id), x="p-ajdusted")
if(output_save==TRUE){ dev.off() }

# check correlation between p-val and p-adjusted
if(output_save==TRUE){ png(file = paste0("correlation_pval-padj_", run_id, ".png")) }
ggplot(res_table, aes(pvalue, fdr.pvalue)) + geom_point()
if(output_save==TRUE){ dev.off() }


