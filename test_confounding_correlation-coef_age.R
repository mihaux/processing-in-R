# this script checks the confounding influence of 'number_of_days_on_steroids'
# => using spearman correlation coefficients between the gene expressions and the values of the number of days on steroids; use the cor.test() function in R to obtain p-values as well; can also try using Pearson (second method)

### List of outputs: [15 files] ###
# table pearson / spearman
# summary of nb of significant results
# histogram p-values pearson 20 / 50 bins
# histogram p-values spearman 20 / 50 bins
# histogram p-adjusted pearson 50 / 20 bins
# histogram p-adjusted spearman 50 / 20 bins
# histogram pearson coefficients 50 / 20 bins
# histogram spearman coefficients 50 / 20 bins

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
args <- c(paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/oct20_confounding/age/method_2_correlation_coef/vst/"))

# Example of usage: 
# Rscript test_confounding_steroids.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load data RAW | VST | rlog (run one at the time)
#load(paste0(args[1], "/Raw_DESeq_dataset_all.Rda"), verbose = TRUE); dds <- dds_all
load(paste0(args[1], "/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
#load(paste0(args[1], "/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# define running ID (either "raw", "vst" pr "rlog")
run_id <- "vst"

# load clinical data 
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)

# add "ID_" to all rownames
rownames(df_meta) <- paste0("ID_", rownames(df_meta))
  
# change data format to matrix and integer
dat <- as.matrix(assay(dds))
storage.mode(dat) <- "integer"          # class: matrix | type: integer

###########################################################################################
# => (2) using spearman correlation coefficients between the gene expressions and the values of the number of days on steroids; use the cor.test() function in R to obtain p-values as well; can also try using Pearson (second method)
###########################################################################################

# initiate list for results
res_pearson <- list()
res_spearman <- list()

# iterate over genes 
for(j in 1:nrow(dat)){
  
  # get data frame for one gene at a time
  temp = dat[j, ]                                           # x => 41 values for gene n
  age_vector <- df_meta$age.at.BL   # y => 41 values of the number of days on steroids
  
  # compute correlation coefficients 
  res_pearson[[j]] <- cor.test(temp, age_vector, alternative = c("two.sided"), method = c("pearson"), conf.level = 0.95, exact=FALSE)
  res_spearman[[j]] <- cor.test(temp, age_vector, alternative = c("two.sided"), method = c("spearman"), conf.level = 0.95, exact=FALSE)
}
  
# Multiple testing correction
p_adj_pearson <- p.adjust(unlist(lapply(res_pearson, function(x) x$p.value)), "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")
p_adj_spearman <- p.adjust(unlist(lapply(res_spearman, function(x) x$p.value)), "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

# create table with results
res_table_pearson <- data.frame(ID=rownames(dat),
                                p.value=unlist(lapply(res_pearson, function(x) x$p.value)),
                                p.adj=p_adj_pearson,
                                pearson_coef=unlist(lapply(res_pearson, function(x) x$estimate)),
                                CI_1=unlist(lapply(res_pearson, function(x) x$conf.int[1])),
                                CI_2=unlist(lapply(res_pearson, function(x) x$conf.int[2])))


res_table_spearman <- data.frame(ID=rownames(dat),
                                 p.value=unlist(lapply(res_spearman, function(x) x$p.value)),
                                 p.adj=p_adj_spearman,
                                 spearman_coef=unlist(lapply(res_spearman, function(x) x$estimate)))

# write a summary table with numbers of significant results
summary_significant = data.frame(V1=run_id, 
                                 V2=length(which(res_table_pearson$p.value < 0.05)),
                                 V3=length(which(res_table_pearson$p.adj < 0.05)),
                                 V4=length(which(res_table_spearman$p.value < 0.05)),
                                 V5=length(which(res_table_spearman$p.adj < 0.05)))

colnames(summary_significant) <- c("", "p-value < 0.05 (Pearson)", "p-adjusted < 0.05 (Pearson)",
                                       "p-value < 0.05 (Spearman)", "p-adjusted < 0.05 (Spearman)")

if(output_save==TRUE){  
write.csv2(summary_significant, file=paste0("table_summary_significant_", run_id, ".csv"))
write.csv2(res_table_pearson, file=paste0("table_pearson_", run_id, ".csv"))
write.csv2(res_table_spearman, file=paste0("table_spearman_", run_id, ".csv"))
}

### make a histogram of p-values ###
if(output_save==TRUE){ png(file = paste0("histogram_p-values_pearson_20bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=p.value)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-values pearson: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_pearson_50bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=p.value)) + 
  geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-values pearson: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_spearman_20_bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=p.value)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-values spearman: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-values_spearman_50_bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=p.value)) + 
  geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-values spearman: ", run_id), x="p-values")
if(output_save==TRUE){ dev.off() }


### make a histogram of p-adj ###
if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_pearson_50bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=p.adj)) + 
  geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-adjusted pearson: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_pearson_20bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=p.adj)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-adjusted pearson: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_spearman_50bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=p.adj)) + 
  geom_histogram(bimns=50) + 
  labs(title=paste0("Histogram of p-adjusted spearman: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_spearman_20bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=p.adj)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-adjusted spearman: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

### make a histogram of coefficients ###
if(output_save==TRUE){ png(file = paste0("histogram_pearson_coefficients_50bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=pearson_coef)) + 
  geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of pearson coefficients: ", run_id), x="pearson coefficients")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_pearson_coefficients_20bins_", run_id, ".png")) }
ggplot(res_table_pearson, aes(x=pearson_coef)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of pearson coefficients: ", run_id), x="pearson coefficients")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_spearman_coefficients_50bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=spearman_coef)) + 
  geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of spearman coefficients: ", run_id), x="spearman coefficients")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_spearman_coefficients_20bins_", run_id, ".png")) }
ggplot(res_table_spearman, aes(x=spearman_coef)) + 
  geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of spearman coefficients: ", run_id), x="spearman coefficients")
if(output_save==TRUE){ dev.off() }
