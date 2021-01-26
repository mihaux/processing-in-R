# this script checks the confounding influence of 'number_of_days_on_steroids'
# => using Mann-Whitney test (split the data in 2 approximately equal groups in terms of the number of patients to get comparison groups for Mann-Whitney test using gene expression data as input)

### List of outputs: [8 files] ###
# histogram of all p-values                               (.png)
# histogram of all p-adjusted                             (.png)
# table with IDs, p-val and p-adj, sorted_by_p-values     (.csv)
# table with IDs, p-val and p-adj, sorted_by_p-adjusted   (.csv)
# table with summary of nb of significant results         (.csv)
# histogram of p-values significant only                  (.png) [if there are some]
# histogram of p-adjusted significant only                (.png) [if there are some]
# correlation plot p-values agains p-adjusted             (.png)

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
       \n(1 - input) path to the directory with counts data
       \n(2 - metadata) path to .csv with clinical features and, 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3
args <- c(paste0(main_dir,"/data/count_matrices/outliers_excluded/counts_vst_no_outliers.csv"),
          paste0(main_dir, "/data/metadata/outliers_excluded/cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv"),
          paste0(main_dir, "/ANALYSES/Dec20_steroids_age_rerun/age/vst/"))

# Example of usage: 
 

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

#dat <- as.matrix(assay(dds))
#storage.mode(dat) <- "integer"          # class: matrix | type: integer

###########################################################################################
# => (1) using Mann-Whitney test (split the data in 2 approximately equal groups in terms of the number of patients to get comparison groups for Mann-Whitney test using gene expression data as input)
###########################################################################################

# df_meta$age.at.BL
# NOTE: tthe samples have a nice normal distribution
# split as:
# group_1 => less than 75 y/o
# group_2 => 75 y/o and more

gr_1 <- df_meta[which(df_meta$age.at.BL < 75),]
gr_2 <- df_meta[which(df_meta$age.at.BL >= 75),]

# add new column to define comparison groups 1 and 2
age_binary <- c(1:40)
age_binary[which(df_meta$age.at.BL < 75)] <- 1
age_binary[which(df_meta$age.at.BL >= 75)] <- 2

df_meta$age.at.BL <- age_binary

# define groups for comparison
group = df_meta$age.at.BL

# initiate list for results
res <- list()
for(j in 1:nrow(dat)){
  
  # get data frame for one gene at a time
  temp = dat[j, ]
  
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

if(FALSE){
# check other methods
all_adjusted_holm       <- p.adjust(all_pval, "holm")       # "holm"
all_adjusted_hochberg   <- p.adjust(all_pval, "hochberg")   # "hochberg"
all_adjusted_hommel     <- p.adjust(all_pval, "hommel")     # "hommel"
all_adjusted_bonferroni <- p.adjust(all_pval, "bonferroni") # "bonferroni"
all_adjusted_by         <- p.adjust(all_pval, "BY")         # "BY"
}

# create a table as R data frame
res_table = data.frame(ID=rownames(dat), 
                       p.value=all_pval, 
                       fdr.pvalue=all_adjusted)

# make a histogram of p-values
if(output_save==TRUE){ png(file = paste0("histogram_p-values_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_20bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_25bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=25) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_50bins_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-values: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

# make a histogram of p-adjusted
if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_", run_id, ".png")) }
ggplot(res_table, aes(x=fdr.pvalue)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-adjusted: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_25bins_", run_id, ".png")) }
ggplot(res_table, aes(x=fdr.pvalue)) + geom_histogram(bins=25) + 
  labs(title=paste0("Histogram of p-adjusted: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

# sort the table by p-values
res_table_sorted <- res_table[order(res_table$p.value),]

# sort the table by p-adjusted
res_table_sorted_bis <- res_table[order(res_table$fdr.pvalue),]

# save tables
if(output_save==TRUE){ 
  write.csv2(res_table_sorted, file=paste0("table_sorted_by_pvalues_", run_id, ".csv"))
  write.csv2(res_table_sorted_bis, file=paste0("table_sorted_by_padjusted_", run_id, ".csv")) 
}

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval <- length(which(res_table_sorted$p.value < 0.05))
significant_padjusted <- length(which(res_table_sorted$fdr.pvalue < 0.05))

# write a summary table with numbers of significant results
summary_significant = data.frame(V1=run_id, 
                                 V2=significant_pval,
                                 V3=significant_padjusted)

colnames(summary_significant) <- c("", "p-valiue < 0.05", "p-adjusted < 0.05")

if(output_save==TRUE){ write.csv2(summary_significant, file=paste0("table_summary_significant_", run_id, ".csv")) }

# create a histogram for significant 
# p-values
if(output_save==TRUE){ print("kot");png(file = paste0("histogram_p-value_significant_only_", run_id, ".png")) }
  if(significant_pval != 0) {  
    ggplot(res_table_sorted[1:significant_pval,], aes(x=p.value)) + 
    geom_histogram(binwidth=0.001) + 
    labs(title=paste0("Histogram of significant p-values: ", run_id), x="p-values")
  } else {
    cat("No significant results with p-value < 0.05")
  } 
if(output_save==TRUE){ print("kot");dev.off() }

# p-adjusted

if(output_save==TRUE){ print("kot");png(file = paste0("histogram_p-adjusted_significant_only_", run_id, ".png")) }
  if(significant_padjusted != 0) {
    ggplot(res_table_sorted[1:significant_padjusted,], aes(x=fdr.pvalue)) + geom_histogram(binwidth=0.01) + 
    labs(title=paste0("Histogram of significant p-adjusted (fdr): ", run_id), x="p-ajdusted")
  } else {
    cat("No significant results with p-adjusted < 0.05")
  }
if(output_save==TRUE){ print("kot");dev.off() }


# check correlation between p-val and p-adjusted
if(output_save==TRUE){ png(file = paste0("correlation_pval-padj_", run_id, ".png")) }
ggplot(res_table, aes(p.value, fdr.pvalue)) + geom_point()
if(output_save==TRUE){ dev.off() }

# NOTES:
# total length => 26 486
# nb of NA values => 0

# TODO: merge manually all three "table_summary_significant_" .csv files

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###########################################################################################
# => (4) using Mann-Whitney test (Q1 vs Q3)
###########################################################################################

quantile(df_meta$age.at.TAB)
# 0%      25%       50%       75%       100% 
#60.00    71.00     75.50     80.25     92.00

# group_1 => 60 - 71 years old
length(which(df_meta$age.at.BL <= 71)) # 11 samples

# group_2 => 80 - 92 years old
length(which(df_meta$age.at.BL >= 80)) # 11 samples

# create a subset of df_meta that will include only samples of interest (12 + 15)
dat_new <- dat[,c(which(df_meta$age.at.BL <= 71), which(df_meta$age.at.BL >= 80))]
df_meta_new <- df_meta[c(which(df_meta$age.at.BL <= 71), which(df_meta$age.at.BL >= 80)),]

any(colnames(dat_new) == rownames(df_meta_new))

gr_1 <- df_meta_new[which(df_meta_new$age.at.BL <= 71),]
gr_2 <- df_meta_new[which(df_meta_new$age.at.BL >= 80),]

# add new column to define comparison groups 1 and 2
steroids <- c(1:22)
steroids[which(df_meta_new$age.at.BL <= 71)] <- 1
steroids[which(df_meta_new$age.at.BL >= 80)] <- 2

df_meta_new$age.at.BL <- steroids

# define groups for comparison
group = df_meta_new$age.at.BL

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




