# extract Mann-Whitney results for selected genes
# data types: Normalised_rlog | Normalised_vst 

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

# define directory with data (INPUT)
data_dir <- paste0(main_dir,"/ANALYSES/Dec20_rerun/RESULTS")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES/Dec20_rerun/EXTRACTED_RESULTS")
setwd(dir_out)

# load data all results files at once (both vst and rlog)
nm <- list.files(path=data_dir, recursive = TRUE, pattern = ".csv", full.names = TRUE)
all_results <- lapply(nm, function(x) read.csv(x, sep = ";"))

names(all_results) <- list.files(path=data_dir, recursive = TRUE, pattern = "\\_rlog.csv$|\\_vst.csv$")
# define running ID (either "raw", "vst" pr "rlog")
#run_id <- "rlog"

# extract SSTR1, SSTR2, SSTR3, SSTR4, SSTR5
idx <- lapply(all_results, function(x) which(x$ID %in% c("SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5")))

results_extracted <- list()
for (i in 1:length(all_results)) {
  results_extracted[[i]] <- all_results[[i]][idx[[i]],]
}


# create a table with all the results (fdr.pvalue)
res.fdr.pvalue <- data.frame(feature=gsub('.{4}$', '', dirname(names(all_results))),
                             sstr1=unlist(lapply(results_extracted, function(x) x[1,4])),
                             sstr2=unlist(lapply(results_extracted, function(x) x[2,4])),
                             sstr3=unlist(lapply(results_extracted, function(x) x[3,4])),
                             sstr4=unlist(lapply(results_extracted, function(x) x[4,4])),
                             sstr5=unlist(lapply(results_extracted, function(x) x[5,4])))
                  
res.p.value <- data.frame(feature=gsub('.{4}$', '', dirname(names(all_results))),
                          sstr1=unlist(lapply(results_extracted, function(x) x[1,3])),
                          sstr2=unlist(lapply(results_extracted, function(x) x[2,3])),
                          sstr3=unlist(lapply(results_extracted, function(x) x[3,3])),
                          sstr4=unlist(lapply(results_extracted, function(x) x[4,3])),
                          sstr5=unlist(lapply(results_extracted, function(x) x[5,3])))

# save results table
if(output_save==TRUE){ 
  write.csv(res.fdr.pvalue, file = paste0("extracted_p-adjusted.csv")) 
  write.csv(res.p.value, file = paste0("extracted_p-values.csv"))
}

# separate rlog and vst and clinical and histological
padj_rlog <- res.fdr.pvalue[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27),]

padj_rlog_clinic <- padj_rlog[c(1,3,7,8,14),]
row.names(padj_rlog_clinic) <- padj_rlog_clinic$feature
padj_rlog_clinic$feature <- NULL

padj_rlog_histo <- padj_rlog[c(2,4,5,6,9,10,11,12,13),]
row.names(padj_rlog_histo) <- padj_rlog_histo$feature
padj_rlog_histo$feature <- NULL

padj_vst <- res.fdr.pvalue[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28),]

padj_vst_clinic <- padj_vst[c(1,3,7,8,14),]
row.names(padj_vst_clinic) <- padj_vst_clinic$feature
padj_vst_clinic$feature <- NULL

padj_vst_histo <- padj_vst[c(2,4,5,6,9,10,11,12,13),]
row.names(padj_vst_histo) <- padj_vst_histo$feature
padj_vst_histo$feature <- NULL

pval_rlog <- res.p.value[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27),]

pval_rlog_clinic <- pval_rlog[c(1,3,7,8,14),]
row.names(pval_rlog_clinic) <- pval_rlog_clinic$feature
pval_rlog_clinic$feature <- NULL

pval_rlog_histo <- pval_rlog[c(2,4,5,6,9,10,11,12,13),]
row.names(pval_rlog_histo) <- pval_rlog_histo$feature
pval_rlog_histo$feature <- NULL

pval_vst <- res.p.value[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28),]

pval_vst_clinic <- pval_vst[c(1,3,7,8,14),]
row.names(pval_vst_clinic) <- pval_vst_clinic$feature
pval_vst_clinic$feature <- NULL

pval_vst_histo <- pval_vst[c(2,4,5,6,9,10,11,12,13),]
row.names(pval_vst_histo) <- pval_vst_histo$feature
pval_vst_histo$feature <- NULL


library('plot.matrix')

png('clinical_rlog_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(padj_rlog_clinic[,1:5]), digits=6, las=2, main="clinical features - rlog - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_rlog_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(padj_rlog_histo[,1:5]), digits=6, las=2, main="histological features - rlog - fdr.pval", ylab="", xlab="")
dev.off()

png('clinical_vst_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(padj_vst_clinic[,1:5]), digits=6, las=2, main="clinical features - vst - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_vst_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(padj_vst_histo[,1:5]), digits=6, las=2, main="histological features - vst - fdr.pval", ylab="", xlab="")
dev.off()

png('clinical_rlog_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(pval_rlog_clinic[,1:5]), digits=6, las=2, main="clinical features - rlog - p.value", ylab="", xlab="")
dev.off()

png('histological_rlog_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(pval_rlog_histo[,1:5]), digits=6, las=2, main="histological features - rlog - p.value", ylab="", xlab="")
dev.off()

png('clinical_vst_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(pval_vst_clinic[,1:5]), digits=6, las=2, main="clinical features - vst - p.value", ylab="", xlab="")
dev.off()

png('histological_vst_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(pval_vst_histo[,1:5]), digits=6, las=2, main="histological features - vst - p.value", ylab="", xlab="")
dev.off()




