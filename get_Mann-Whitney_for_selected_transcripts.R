# extract Mann-Whitney results for selected transcripts
# data types: raw

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
output_save <- FALSE

# define directory with data (INPUT)
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/statistical_results")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_results/EXTRACTED_RESULTS")
setwd(dir_out)

# load data all results files at once (both vst and rlog)
nm <- list.files(path=data_dir, recursive = TRUE, pattern = "_raw.csv", full.names = TRUE)
all_results <- lapply(nm, function(x) read.csv(x, sep = ";"))

names(all_results) <- list.files(path=data_dir, recursive = TRUE, pattern = "\\_raw.csv$")
# define running ID (either "raw", "vst" pr "rlog")
#run_id <- "rlog"

# extract all transcripts of interest
# => in df, look for:
# SSTR1: "ENST00000267377.2"
# SSTR2: "ENST00000357585.3", "ENST00000579323.5"
# SSTR3: "ENST00000610913.1", "ENST00000617123.1"
# SSTR4: "ENST00000255008.4"
# SSTR5: "ENST00000293897.5"

idx <- lapply(all_results, function(x) which(x$ID %in% c("ENST00000267377.2", "ENST00000357585.3", "ENST00000579323.5", "ENST00000610913.1", "ENST00000617123.1", "ENST00000255008.4", "ENST00000293897.5")))

results_extracted <- list()
for (i in 1:length(all_results)) {
  results_extracted[[i]] <- all_results[[i]][idx[[i]],]
}

feat_names <- substr(names(all_results), 25, 50)
feat_names_fin <- gsub('.{0,8}$', '', feat_names)

# create a table with all the results (fdr.pvalue)
res.fdr.pvalue <- data.frame(feature=feat_names_fin,
                             sstr1_t1=unlist(lapply(results_extracted, function(x) x[1,4])),
                             sstr2_t1=unlist(lapply(results_extracted, function(x) x[2,4])),
                             sstr2_t2=unlist(lapply(results_extracted, function(x) x[3,4])),
                             sstr3_t1=unlist(lapply(results_extracted, function(x) x[4,4])),
                             sstr3_t2=unlist(lapply(results_extracted, function(x) x[5,4])),
                             sstr4_t1=unlist(lapply(results_extracted, function(x) x[5,4])),
                             sstr5_t1=unlist(lapply(results_extracted, function(x) x[5,4])))
                  

res.p.value <- data.frame(feature=feat_names_fin,
                          sstr1_t1=unlist(lapply(results_extracted, function(x) x[1,3])),
                          sstr2_t1=unlist(lapply(results_extracted, function(x) x[2,3])),
                          sstr2_t2=unlist(lapply(results_extracted, function(x) x[3,3])),
                          sstr3_t1=unlist(lapply(results_extracted, function(x) x[4,3])),
                          sstr3_t2=unlist(lapply(results_extracted, function(x) x[5,3])),
                          sstr4_t1=unlist(lapply(results_extracted, function(x) x[5,3])),
                          sstr5_t1=unlist(lapply(results_extracted, function(x) x[5,3])))

# save results table
if(output_save==TRUE){ 
  write.csv(res.fdr.pvalue, file = paste0("extracted_p-adjusted.csv")) 
  write.csv(res.p.value, file = paste0("extracted_p-values.csv"))
}

# separate clinical and histological features
res.fdr.pvalue_clinical <- res.fdr.pvalue[c(2,4,5,9),]
row.names(res.fdr.pvalue_clinical) <- res.fdr.pvalue_clinical$feature
res.fdr.pvalue_clinical$feature <- NULL

res.fdr.pvalue_histological <- res.fdr.pvalue[c(1,3,6,7,8),]
row.names(res.fdr.pvalue_histological) <- res.fdr.pvalue_histological$feature
res.fdr.pvalue_histological$feature <- NULL
  
res.p.value_clinical <- res.p.value[c(2,4,5,9),]
row.names(res.p.value_clinical) <- res.p.value_clinical$feature
res.p.value_clinical$feature <- NULL

res.p.value_histological <- res.p.value[c(1,3,6,7,8),]
row.names(res.p.value_histological) <- res.p.value_histological$feature
res.p.value_histological$feature <- NULL

library('plot.matrix')

png('clinical_raw_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_clinical[,1:7]), digits=6, las=2, main="clinical features - raw - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_raw_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_histological[,1:7]), digits=6, las=2, main="histological features - raw - fdr.pval", ylab="", xlab="")
dev.off()

png('clinical_raw_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_clinical[,1:7]), digits=6, las=2, main="clinical features - raw - p.value", ylab="", xlab="")
dev.off()

png('histological_raw_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_histological[,1:7]), digits=6, las=2, main="histological features - raw - p.value", ylab="", xlab="")
dev.off()





