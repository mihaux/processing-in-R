# this script checks the confounding influence of 'number_of_days_on_steroids'
# => using Mann-Whitney test (split the data in 2 approximately equal groups in terms of the number of patients to get comparison groups for Mann-Whitney test using gene expression data as input)
# => using spearman correlation coefficients between the gene expressions and the values of the number of days on steroids; use the cor.test() function in R to obtain p-values as well; can also try using Pearson (second method)

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
          paste0(main_dir, "/ANALYSES/oct20_confounding"))

# Example of usage: 
# Rscript test_confounding_steroids.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load data RAW | VST | rlog (run one at the time)
#load(paste0(args[1], "/Raw_DESeq_dataset_all.Rda"), verbose = TRUE); dds <- dds_all
#load(paste0(args[1], "/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
load(paste0(args[1], "/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# load clinical data 
df_meta <- read.csv(args[2], row.names = 1, header = TRUE)

# add "ID_" to all rownames
rownames(df_meta) <- paste0("ID_", rownames(df_meta))

# df_meta$number.of.days.on.steroids.at.TAB

# NOTE: there are 2 samples with 16 days, can consider them as too extreme and exclude

# split as:
# group_1 => 0 - 5 days
# group_2 => 6 - 16 days
# group_2_bis => 6 - 11 days (when both samples with 16 days are excluded)

gr_1 <- df_meta[which(df_meta$number.of.days.on.steroids.at.TAB < 6),]
gr_2 <- df_meta[which(df_meta$number.of.days.on.steroids.at.TAB >= 6),]

gr_2_bis <- gr_2[-which(gr_2$number.of.days.on.steroids.at.TAB == 16),]
  
# change data format to matrix and integer
dat <- as.matrix(assay(dds))
storage.mode(dat) <- "integer"          # class: matrix | type: integer

# add new column to define comparison groups 1 and 2
steroids <- c(1:41)
steroids[which(df_meta$number.of.days.on.steroids.at.TAB < 6)] <- 1
steroids[which(df_meta$number.of.days.on.steroids.at.TAB >= 6)] <- 2

df_meta$number.of.days.on.steroids.at.TAB <- steroids

# define running ID
run_id <- "rlog"

# define groups for comparison
group = df_meta$number.of.days.on.steroids.at.TAB

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

### histogram of p-values ###
# create a df with IDs and p-values
df_pval <- data.frame(ID=rownames(dat), pval=all_pval)

# make a histogram of pvalues
if(output_save==TRUE){ png(file = paste0("histogram_pvalues_", run_id, ".png")) }

ggplot(df_pval, aes(x=pval)) + 
  geom_histogram(binwidth=0.01) + 
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
ggplot(res_table_sorted[1:752,], aes(x=fdr.pvalue)) + 
  geom_histogram(binwidth=0.01) + 
  labs(title=paste0("Histogram of significant p-adjusted (fdr): ", run_id), x="p-ajdusted")
if(output_save==TRUE){ dev.off() }

# check correlation between p-val and p-adjusted











################################################################################################
### with 2 samples of 16 days excluded

# exclude from count datasets and clinical
dat_raw_bis <- dat_raw[,-which(gr_2$number.of.days.on.steroids.at.TAB == 16)]
dat_vst_bis <- dat_vst[,-which(gr_2$number.of.days.on.steroids.at.TAB == 16)]
dat_rlog_bis <- dat_rlog[,-which(gr_2$number.of.days.on.steroids.at.TAB == 16)]

df_meta_bis <- df_meta[-which(gr_2$number.of.days.on.steroids.at.TAB == 16),]

# add new column to define comparison groups 1 and 2
steroids_bis <- c(1:39)
steroids_bis[which(df_meta_bis$number.of.days.on.steroids.at.TAB < 6)] <- 1
steroids_bis[which(df_meta_bis$number.of.days.on.steroids.at.TAB >= 6)] <- 2

df_meta_bis$number.of.days.on.steroids.at.TAB <- steroids_bis

# define groups for comparison
group_bis = df_meta_bis$number.of.days.on.steroids.at.TAB

# initiate list for results
res_raw_bis <- list()
res_vst_bis <- list()
res_rlog_bis <- list()

for(j in 1:nrow(dat_raw_bis)){
  
  # get data frame for one gene at a time
  temp_raw = dat_raw_bis[j, ]
  temp_vst = dat_vst_bis[j, ]
  temp_rlog = dat_rlog_bis[j, ]
  
  # perform statistical testing
  res_raw_bis[[j]] = wilcox.test(temp_raw[group_bis==1], 
                             temp_raw[group_bis==2], 
                             alternative = "two.sided",
                             exact = FALSE)               # to suppress the warning message saying that “cannot compute exact p-value with tie”
  # it comes from the assumption of a Wilcoxon test that the responses are continuous. 
  
  res_vst_bis[[j]] = wilcox.test(temp_vst[group_bis==1], temp_vst[group_bis==2], alternative = "two.sided", exact = FALSE)              
  res_rlog_bis[[j]] = wilcox.test(temp_rlog[group_bis==1], temp_rlog[group_bis==2], alternative = "two.sided", exact = FALSE)              
  
}

# retrieve its corresponding p-value
all_pval_raw_bis = unlist(lapply(res_raw_bis, function(x) x$p.value))
all_pval_vst_bis = unlist(lapply(res_vst_bis, function(x) x$p.value))
all_pval_rlog_bis = unlist(lapply(res_rlog_bis, function(x) x$p.value))

# Multiplicity adjustment
# After obtaining the p-values, you can make multiplicity correction to the pvalue. 
all_adj.pval_raw_bis <- p.adjust(all_pval_raw_bis, "fdr")    # Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
all_adj.pval_vst_bis <- p.adjust(all_pval_vst_bis, "fdr") 
all_adj.pval_rlog_bis <- p.adjust(all_pval_rlog_bis, "fdr") 

# result table
# create a table as R data frame
result.table2_raw_bis = data.frame(ID=rownames(dat_raw_bis), pvalue=all_pval_raw_bis, fdr.pvalue=all_adj.pval_raw_bis)
result.table2_vst_bis = data.frame(ID=rownames(dat_vst_bis), pvalue=all_pval_vst_bis, fdr.pvalue=all_adj.pval_vst_bis)
result.table2_rlog_bis = data.frame(ID=rownames(dat_rlog_bis), pvalue=all_pval_rlog_bis, fdr.pvalue=all_adj.pval_rlog_bis)

# sort the table based on the order of the (adjusted) p-values
result.table2.sorted_raw_bis = result.table2_raw_bis[order(all_adj.pval_raw_bis),]
result.table2.sorted_vst_bis = result.table2_vst_bis[order(all_adj.pval_vst_bis),]
result.table2.sorted_rlog_bis = result.table2_rlog_bis[order(all_adj.pval_rlog_bis),]


result.table2.sorted_raw_bis_pval = result.table2_raw_bis[order(all_pval_raw_bis),]
result.table2.sorted_vst_bis_pval = result.table2_vst_bis[order(all_pval_vst_bis),]
result.table2.sorted_rlog_bis_pval = result.table2_rlog_bis[order(all_pval_rlog_bis),]


significant.cutoff_pval_raw_bis <- length(which(result.table2.sorted_raw_bis$pvalue < 0.05))
significant.cutoff_pval_vst_bis <- length(which(result.table2.sorted_vst_bis$pvalue < 0.05))
significant.cutoff_pval_rlog_bis <- length(which(result.table2.sorted_rlog_bis$pvalue < 0.05))


