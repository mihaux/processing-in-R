# this script runs Mann-Whitney testing between the outlier samples vs. non-outliers 

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
          paste0(main_dir, "/ANALYSES/nov20_outliers_testing/rlog/"))

# Example of usage: 
# Rscript test_confounding_Mann-Whitney_steroids.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load data RAW | VST | rlog (run one at the time)
#load(paste0(args[1], "/Raw_DESeq_dataset_all.Rda"), verbose = TRUE); dds <- dds_all
#load(paste0(args[1], "/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
load(paste0(args[1], "/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# define running ID (either "raw", "vst" pr "rlog")
run_id <- "rlog"

# load clinical data 
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)

# add "ID_" to all rownames
rownames(df_meta) <- paste0("ID_", rownames(df_meta))
  
# change data format to matrix and integer
dat <- as.matrix(assay(dds))
storage.mode(dat) <- "integer"          # class: matrix | type: integer

# define both groups (outliers and non-outliers)
gr_1 <- as.vector(c("ID_12330", "ID_12331", "ID_14930", "ID_16137", "ID_16142", "ID_16153", "ID_16187"))
to_be_excluded <- which(colnames(dat) %in% gr_1)
gr_2 <- colnames(dat)[-to_be_excluded]

# add new column to define comparison groups 1 - outliers and 2 - non-outliers
outlier_status <- c(rep(2, 41))
outlier_status[which(colnames(dat) %in% gr_1)] <- 1   # outliers

# add "outlier_status" to the metadata spreadsheet
df_meta <- cbind(df_meta, outlier_status)

# define groups for comparison
group = df_meta$outlier_status

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

# make a histogram of p-values
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_", run_id, ".png")) }
ggplot(res_table, aes(x=p.value)) + geom_histogram(binwidth=0.01) + 
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

# make a histogram of p-adjusted
if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_", run_id, ".png")) }
ggplot(res_table, aes(x=fdr.pvalue)) + geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of p-adjusted: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_20bins_", run_id, ".png")) }
ggplot(res_table, aes(x=fdr.pvalue)) + geom_histogram(bins=20) + 
  labs(title=paste0("Histogram of p-adjusted: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

if(output_save==TRUE){ png(file = paste0("histogram_p-adjusted_50bins_", run_id, ".png")) }
ggplot(res_table, aes(x=fdr.pvalue)) + geom_histogram(bins=50) + 
  labs(title=paste0("Histogram of p-adjusted: ", run_id), x="p-adjusted")
if(output_save==TRUE){ dev.off() }

# sort the table by p-values
res_table_sorted <- res_table[order(all_pval),]

# sort the table by p-adjusted
res_table_sorted_bis <- res_table[order(all_adjusted),]

# save tables
write.csv2(res_table_sorted, file=paste0("table_sorted_pvalues_", run_id, ".csv"))
write.csv2(res_table_sorted_bis, file=paste0("table_sorted_padjusted_", run_id, ".csv"))

# filter out statistically insignificant results (on p-value and p-adjusted)
significant_pval <- length(which(res_table$p.value < 0.05))
significant_padjusted <- length(which(res_table$fdr.pvalue < 0.05))

# write a summary table with numbers of significant results
summary_significant = data.frame(V1=run_id, 
                                 V2=significant_pval,
                                 V3=significant_padjusted)

colnames(summary_significant) <- c("", "p-valiue < 0.05", "p-adjusted < 0.05")

if(output_save==TRUE){ write.csv2(summary_significant, file=paste0("table_summary_significant_", run_id, ".csv")) }

