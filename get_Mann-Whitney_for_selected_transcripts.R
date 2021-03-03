# extract Mann-Whitney results for selected transcripts
# data types: RAW & VST 

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(ggplot2)
library(qdapRegex)

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
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/statistical_testing/output/")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/output/extracted/transcript-level")
setwd(dir_out)

# load data all results files at once (both RAW and VST)
nm_clinical_raw <- list.files(path=paste0(data_dir, "clinical"), recursive = TRUE, pattern = "raw_transcript-level.csv", full.names = TRUE)
nm_clinical_VST <- list.files(path=paste0(data_dir, "clinical"), recursive = TRUE, pattern = "VST_transcript-level.csv", full.names = TRUE)

nm_histological_raw <- list.files(path=paste0(data_dir, "histological"), recursive = TRUE, pattern = "raw_transcript-level.csv", full.names = TRUE)
nm_histological_VST <- list.files(path=paste0(data_dir, "histological"), recursive = TRUE, pattern = "VST_transcript-level.csv", full.names = TRUE)

# keep only "table_sorted_by_pvalues_" (half of each list)
nm_clinical_raw_fin <- nm_clinical_raw[seq(2, length(nm_clinical_raw), by = 2)]
nm_clinical_VST_fin <- nm_clinical_VST[seq(2, length(nm_clinical_VST), by = 2)]

nm_histological_raw_fin <- nm_histological_raw[seq(2, length(nm_histological_raw), by = 2)]
nm_histological_VST_fin <- nm_histological_VST[seq(2, length(nm_histological_VST), by = 2)]

# read all table results and put them on a list
all_results_clinical_raw_fin <- lapply(nm_clinical_raw_fin, function(x) read.csv(x, sep = ";"))
all_results_clinical_VST_fin <- lapply(nm_clinical_VST_fin, function(x) read.csv(x, sep = ";"))

all_results_histological_raw_fin <- lapply(nm_histological_raw_fin, function(x) read.csv(x, sep = ";"))
all_results_histological_VST_fin <- lapply(nm_histological_VST_fin, function(x) read.csv(x, sep = ";"))

# give it proper names
names(all_results_clinical_raw_fin) <- unlist(qdapRegex::ex_between(nm_clinical_raw_fin, "clinical/", "/table"))
names(all_results_clinical_VST_fin) <- unlist(qdapRegex::ex_between(nm_clinical_VST_fin, "clinical/", "/table"))

names(all_results_histological_raw_fin) <- unlist(qdapRegex::ex_between(nm_histological_raw_fin, "histological/", "/table"))
names(all_results_histological_VST_fin) <- unlist(qdapRegex::ex_between(nm_histological_VST_fin, "histological/", "/table"))

# extract all transcripts of interest
# => in df, look for:
# SSTR1: "ENST00000267377.2"
# SSTR2: "ENST00000357585.3", "ENST00000579323.5"
# SSTR3: "ENST00000610913.1", "ENST00000617123.1"
# SSTR4: "ENST00000255008.4"
# SSTR5: "ENST00000293897.5"

idx_clinical_raw <- lapply(all_results_clinical_raw_fin, function(x) which(x$ID %in% c("ENST00000267377.2", "ENST00000357585.3", "ENST00000579323.5", "ENST00000610913.1", "ENST00000617123.1", "ENST00000255008.4", "ENST00000293897.5")))
idx_clinical_VST <- lapply(all_results_clinical_VST_fin, function(x) which(x$ID %in% c("ENST00000267377.2", "ENST00000357585.3", "ENST00000579323.5", "ENST00000610913.1", "ENST00000617123.1", "ENST00000255008.4", "ENST00000293897.5")))

idx_histological_raw <- lapply(all_results_histological_raw_fin, function(x) which(x$ID %in% c("ENST00000267377.2", "ENST00000357585.3", "ENST00000579323.5", "ENST00000610913.1", "ENST00000617123.1", "ENST00000255008.4", "ENST00000293897.5")))
idx_histological_VST <- lapply(all_results_histological_VST_fin, function(x) which(x$ID %in% c("ENST00000267377.2", "ENST00000357585.3", "ENST00000579323.5", "ENST00000610913.1", "ENST00000617123.1", "ENST00000255008.4", "ENST00000293897.5")))

# extract the results
results_extracted_clinical_raw <- list()
for (i in 1:length(all_results_clinical_raw_fin)) {
  results_extracted_clinical_raw[[i]] <- all_results_clinical_raw_fin[[i]][idx_clinical_raw[[i]],]
}

results_extracted_clinical_VST <- list()
for (i in 1:length(all_results_clinical_VST_fin)) {
  results_extracted_clinical_VST[[i]] <- all_results_clinical_VST_fin[[i]][idx_clinical_VST[[i]],]
}

results_extracted_histological_raw <- list()
for (i in 1:length(all_results_histological_raw_fin)) {
  results_extracted_histological_raw[[i]] <- all_results_histological_raw_fin[[i]][idx_histological_raw[[i]],]
}

results_extracted_histological_VST <- list()
for (i in 1:length(all_results_histological_VST_fin)) {
  results_extracted_histological_VST[[i]] <- all_results_histological_VST_fin[[i]][idx_histological_VST[[i]],]
}

# give them proper names
feat_names_fin_clinical_raw <- names(all_results_clinical_raw_fin)
feat_names_fin_clinical_VST <- names(all_results_clinical_VST_fin)

feat_names_fin_histological_raw <- names(all_results_histological_raw_fin)
feat_names_fin_histological_VST <- names(all_results_histological_VST_fin)

# create a table with all the results (fdr.pvalue) and (p.value)
### => clinical_raw
res.fdr.pvalue_clinical_raw <- data.frame(feature=feat_names_fin_clinical_raw,
                                          sstr1_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[1,4])),
                                          sstr2_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[2,4])),
                                          sstr2_t2=unlist(lapply(results_extracted_clinical_raw, function(x) x[3,4])),
                                          sstr3_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[4,4])),
                                          sstr3_t2=unlist(lapply(results_extracted_clinical_raw, function(x) x[5,4])),
                                          sstr4_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[6,4])),
                                          sstr5_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[7,4])))
                  

res.p.value_clinical_raw <- data.frame(feature=feat_names_fin_clinical_raw,
                                          sstr1_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[1,3])),
                                          sstr2_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[2,3])),
                                          sstr2_t2=unlist(lapply(results_extracted_clinical_raw, function(x) x[3,3])),
                                          sstr3_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[4,3])),
                                          sstr3_t2=unlist(lapply(results_extracted_clinical_raw, function(x) x[5,3])),
                                          sstr4_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[6,3])),
                                          sstr5_t1=unlist(lapply(results_extracted_clinical_raw, function(x) x[7,3])))

### => clinical_VST
res.fdr.pvalue_clinical_VST <- data.frame(feature=feat_names_fin_clinical_VST,
                                          sstr1_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[1,4])),
                                          sstr2_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[2,4])),
                                          sstr2_t2=unlist(lapply(results_extracted_clinical_VST, function(x) x[3,4])),
                                          sstr3_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[4,4])),
                                          sstr3_t2=unlist(lapply(results_extracted_clinical_VST, function(x) x[5,4])),
                                          sstr4_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[6,4])),
                                          sstr5_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[7,4])))


res.p.value_clinical_VST <- data.frame(feature=feat_names_fin_clinical_VST,
                                       sstr1_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[1,3])),
                                       sstr2_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[2,3])),
                                       sstr2_t2=unlist(lapply(results_extracted_clinical_VST, function(x) x[3,3])),
                                       sstr3_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[4,3])),
                                       sstr3_t2=unlist(lapply(results_extracted_clinical_VST, function(x) x[5,3])),
                                       sstr4_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[6,3])),
                                       sstr5_t1=unlist(lapply(results_extracted_clinical_VST, function(x) x[7,3])))

### => histological_raw
res.fdr.pvalue_histological_raw <- data.frame(feature=feat_names_fin_histological_raw,
                                              sstr1_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[1,4])),
                                              sstr2_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[2,4])),
                                              sstr2_t2=unlist(lapply(results_extracted_histological_raw, function(x) x[3,4])),
                                              sstr3_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[4,4])),
                                              sstr3_t2=unlist(lapply(results_extracted_histological_raw, function(x) x[5,4])),
                                              sstr4_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[6,4])),
                                              sstr5_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[7,4])))


res.p.value_histological_raw <- data.frame(feature=feat_names_fin_histological_raw,
                                           sstr1_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[1,3])),
                                           sstr2_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[2,3])),
                                           sstr2_t2=unlist(lapply(results_extracted_histological_raw, function(x) x[3,3])),
                                           sstr3_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[4,3])),
                                           sstr3_t2=unlist(lapply(results_extracted_histological_raw, function(x) x[5,3])),
                                           sstr4_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[6,3])),
                                           sstr5_t1=unlist(lapply(results_extracted_histological_raw, function(x) x[7,3])))

### => histological_VST
res.fdr.pvalue_histological_VST <- data.frame(feature=feat_names_fin_histological_VST,
                                              sstr1_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[1,4])),
                                              sstr2_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[2,4])),
                                              sstr2_t2=unlist(lapply(results_extracted_histological_VST, function(x) x[3,4])),
                                              sstr3_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[4,4])),
                                              sstr3_t2=unlist(lapply(results_extracted_histological_VST, function(x) x[5,4])),
                                              sstr4_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[6,4])),
                                              sstr5_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[7,4])))


res.p.value_histological_VST <- data.frame(feature=feat_names_fin_histological_VST,
                                           sstr1_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[1,3])),
                                           sstr2_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[2,3])),
                                           sstr2_t2=unlist(lapply(results_extracted_histological_VST, function(x) x[3,3])),
                                           sstr3_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[4,3])),
                                           sstr3_t2=unlist(lapply(results_extracted_histological_VST, function(x) x[5,3])),
                                           sstr4_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[6,3])),
                                           sstr5_t1=unlist(lapply(results_extracted_histological_VST, function(x) x[7,3])))


# save results table
if(output_save==TRUE){ 
  write.csv(res.fdr.pvalue_clinical_raw, file = paste0("extracted_p-adjusted_clinical_raw.csv")) 
  write.csv(res.p.value_clinical_raw, file = paste0("extracted_p-values_clinical_raw.csv"))

  write.csv(res.fdr.pvalue_clinical_VST, file = paste0("extracted_p-adjusted_clinical_VST.csv")) 
  write.csv(res.p.value_clinical_VST, file = paste0("extracted_p-values_clinical_VST.csv"))
  
  write.csv(res.fdr.pvalue_histological_raw, file = paste0("extracted_p-adjusted_histological_raw.csv")) 
  write.csv(res.p.value_histological_raw, file = paste0("extracted_p-values_histological_raw.csv"))
  
  write.csv(res.fdr.pvalue_histological_VST, file = paste0("extracted_p-adjusted_histological_VST.csv")) 
  write.csv(res.p.value_histological_VST, file = paste0("extracted_p-values_histological_VST.csv"))
  
}

if(FALSE){ # not needed anymore
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
}

library('plot.matrix')

#### => p-value

# add correct rownames
rownames(res.p.value_clinical_raw) <- res.p.value_clinical_raw$feature
rownames(res.p.value_clinical_VST) <- res.p.value_clinical_VST$feature

rownames(res.p.value_histological_raw) <- res.p.value_histological_raw$feature
rownames(res.p.value_histological_VST) <- res.p.value_histological_VST$feature

png('clinical_raw_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_clinical_raw[,2:8]), digits=6, las=2, main="clinical features - raw - p.value", ylab="", xlab="")
dev.off()

png('clinical_VST_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_clinical_VST[,2:8]), digits=6, las=2, main="clinical features - VST - p.value", ylab="", xlab="")
dev.off()

png('histological_raw_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_histological_raw[,2:8]), digits=6, las=2, main="histological features - raw - p.value", ylab="", xlab="")
dev.off()

png('histological_VST_pvalue.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_histological_VST[,2:8]), digits=6, las=2, main="histological features - VST - p.value", ylab="", xlab="")
dev.off()

#### => p-adjusted

# add correct rownames
rownames(res.fdr.pvalue_clinical_raw) <- res.fdr.pvalue_clinical_raw$feature
rownames(res.fdr.pvalue_clinical_VST) <- res.fdr.pvalue_clinical_VST$feature
rownames(res.fdr.pvalue_histological_raw) <- res.fdr.pvalue_histological_raw$feature
rownames(res.fdr.pvalue_histological_VST) <- res.fdr.pvalue_histological_VST$feature

png('clinical_raw_padjusted_clinical.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_clinical_raw[,2:8]), digits=6, las=2, main="clinical features - raw - fdr.pval", ylab="", xlab="")
dev.off()

png('clinical_VST_padjusted_clinical.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_clinical_VST[,2:8]), digits=6, las=2, main="clinical features - VST - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_raw_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_histological_raw[,2:8]), digits=6, las=2, main="histological features - raw - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_VST_padjusted.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_histological_VST[,2:8]), digits=6, las=2, main="histological features - VST - fdr.pval", ylab="", xlab="")
dev.off()


# SAME AS ABOVE BUT WITHOUT "Infiltrate around vasa vasorum", "Intima pattern", "Media pattern"
exclude <- c(which(rownames(res.p.value_histological_raw) == "Infiltrate_around_vasa_vasorum"),
             which(rownames(res.p.value_histological_raw) == "Intima_pattern"),
             which(rownames(res.p.value_histological_raw) == "Media_pattern"))

png('histological_raw_pvalue_MOD.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_histological_raw[-exclude,2:8]), digits=6, las=2, main="histological features - raw - p.value", ylab="", xlab="")
dev.off()

png('histological_VST_pvalue_MOD.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.p.value_histological_VST[-exclude,2:8]), digits=6, las=2, main="histological features - VST - p.value", ylab="", xlab="")
dev.off()

png('histological_raw_padjusted_MOD.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_histological_raw[-exclude,2:8]), digits=6, las=2, main="histological features - raw - fdr.pval", ylab="", xlab="")
dev.off()

png('histological_VST_padjusted_MOD.png', width = 1024, height = 768)
par(mar=c(11.1, 9.1, 11.1, 9.1))
plot(as.matrix(res.fdr.pvalue_histological_VST[-exclude,2:8]), digits=6, las=2, main="histological features - VST - fdr.pval", ylab="", xlab="")
dev.off()
