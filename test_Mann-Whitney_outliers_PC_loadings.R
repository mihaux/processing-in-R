# script to compare results from Mann-Whitney for outliers vs. non-outliers and PCA loadings

# install (if necessary) and load package
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

setwd(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/comparison_with_PC_loadings"))

# load tables with results (p-values, p-adjusted)
DE_res_raw <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/raw/table_sorted_padjusted_raw_mod.csv"), row.names = 1)
DE_res_vst <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/vst/table_sorted_padjusted_vst_mod.csv"), row.names = 1)
DE_res_rlog <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/rlog/table_sorted_padjusted_rlog_mod.csv"), row.names = 1)

# get gene names of the results with p-adj < 0.05
DE_genes_raw <- rownames(DE_res_raw)[which(DE_res_raw$fdr.pvalue < 0.05)]
DE_genes_vst <- rownames(DE_res_vst)[which(DE_res_vst$fdr.pvalue < 0.05)]
DE_genes_rlog <- rownames(DE_res_rlog)[which(DE_res_rlog$fdr.pvalue < 0.05)]

# run PCA using the built-in method and data with gene symbols
# define directory with data (INPUT)
data_dir_loadings <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts")

# load data RAW | VST | rlog (run for one data type at a time)
load(paste0(data_dir_loadings, "/Raw_DESeq_dataset_all.Rda")); dds <- dds_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# run PCA using the built-in function
pca_dds <- prcomp(t(assay(dds)))

# retrieve PCA loadings
loadings_dds <- as.data.frame(pca_dds$rotation)

# get PC1 and PC2 values for DE_genes_raw | DE_genes_vst | DE_genes_rlog
#loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)]
#loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)]

loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_vst)]
loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_vst)]

#loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_rlog)]
#loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_rlog)]

# make histograms of PC1 and PC2
#hist(loadings_dds$PC1, breaks = 5, main = "Histogram of PC1 loadings - raw")
#hist(loadings_dds$PC2, breaks = 5, main = "Histogram of PC2 loadings - raw")

hist(loadings_dds$PC1, breaks = 5, main = "Histogram of PC1 loadings - vst")
hist(loadings_dds$PC2, breaks = 5, main = "Histogram of PC2 loadings - vst")

#hist(loadings_dds$PC1, breaks = 5, main = "Histogram of PC1 loadings - rlog")
#hist(loadings_dds$PC2, breaks = 5, main = "Histogram of PC2 loadings - rlog")

# min and max value of PC1
#min(loadings_dds$PC1); max(loadings_dds$PC1)

# plot loadings of PC1
plot(loadings_dds$PC1, ylab="PC1 loadings")

# check the "extreme values"
# raw
#plot(loadings_dds$PC1[which(loadings_dds$PC1>0.05)], main="Extreme values of PC1 loadings", ylab="PC1 loadings > 0.05")
#text(loadings_dds$PC1[which(loadings_dds$PC1>0.05)], labels=rownames(loadings_dds)[which(loadings_dds$PC1>0.05)], cex=0.6, pos=3, col="red")   # print labels

# vst
# greater than 0.025 => which(loadings_dds$PC1 > 0.025)
plot(loadings_dds$PC1[which(loadings_dds$PC1>0.025)], main="Extreme values of PC1 loadings (positive)", ylab="PC1 loadings > 0.025")
text(loadings_dds$PC1[which(loadings_dds$PC1>0.025)], labels=rownames(loadings_dds)[which(loadings_dds$PC1>0.025)], cex=0.6, pos=3, col="red")   # print labels

# smaller than -0.025 => which(loadings_dds$PC1 < -0.025)
plot(loadings_dds$PC1[which(loadings_dds$PC1<(-0.025))], main="Extreme values of PC1 loadings (negative)", ylab="PC1 loadings < -0.025")
text(loadings_dds$PC1[which(loadings_dds$PC1<(-0.025))], labels=rownames(loadings_dds)[which(loadings_dds$PC1<=-0.025)], cex=0.6, pos=3, col="red")   # print labels

# rlog

# plot loadings of PC2
plot(loadings_dds$PC2, ylab="PC2 loadings")

# raw
# plot(loadings_dds$PC2[which(loadings_dds$PC2>0.015)], main="Extreme values of PC2 loadings (positive)", ylab="PC2 loadings > 0.015")
# text(loadings_dds$PC2[which(loadings_dds$PC2>0.015)], labels=rownames(loadings_dds)[which(loadings_dds$PC2>0.015)], cex=0.6, pos=3, col="red")   # print labels

# plot(loadings_dds$PC2[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], main="Extreme values of PC2 loadings (negative)", ylab="PC2 loadings < -0.05 & > -0.08")
# text(loadings_dds$PC2[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], labels=rownames(loadings_dds)[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], cex=0.6, pos=3, col="red")   # print labels

# vst

# rlog

# get loadings of DE genes from Mann-Whitney
which(rownames(loadings_dds) %in% DE_genes_raw)

rownames(loadings_dds)[which(rownames(loadings_dds) %in% DE_genes_raw)]

plot(loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)], main="PC1 loadings of 30 DE genes - raw data", ylab="PC1 loadings")
plot(loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)], main="PC2 loadings of 30 DE genes - raw data", ylab="PC2 loadings")

df_raw_final <- data.frame(genes=rownames(loadings_dds)[which(rownames(loadings_dds) %in% DE_genes_raw)],
                           PC1=loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)],
                           PC2=loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)])

write.csv(df_raw_final, file="df_raw_final.csv")

#####################################################################################
# average loading per each chromosome (for PC1 only) -> handle neg and pos values separately
#####################################################################################

# change warking directory
setwd(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/results_chromosomal_positions"))

# load positions on chromsomes
chr_positions <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_manhattan_plot/sorted.csv"), row.names = 1)


# get positions of genes on chr
names_chr1 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 1)]) # 2356
idx_chr1 <- which(rownames(loadings_dds) %in% names_chr1)
chr1_pos <- loadings_dds$PC1[idx_chr1]

names_chr2 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 2)]) # 1404
names_chr3 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 3)]) # 1237
names_chr4 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 4)]) # 825
names_chr5 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 5)]) # 984
names_chr6 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 6)]) # 1173
names_chr7 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 7)]) # 1071
names_chr8 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 8)]) # 763
names_chr9 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 9)]) # 900
names_chr10 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 10)]) # 858
names_chr11 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 11)]) # 1427
names_chr12 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 12)]) # 1154
names_chr13 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 13)]) # 398
names_chr14 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 14)]) # 772
names_chr15 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 15)]) # 810
names_chr16 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 16)]) # 939
names_chr17 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 17)]) # 1344
names_chr18 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 18)]) # 298
names_chr19 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 19)]) # 1573
names_chr20 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 20)]) # 621
names_chr21 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 21)]) # 296
names_chr22 <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == 22)]) # 538
names_chrX <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == "X")]) # 914
names_chrY <- as.character(chr_positions$Approved.symbol[which(chr_positions$Chromosome == "Y")]) # 70

names_all <- list(names_chr1, names_chr2, names_chr3, names_chr4, names_chr5, names_chr6, names_chr7,
                  names_chr8, names_chr9, names_chr10, names_chr11, names_chr12, names_chr13, names_chr14,
                  names_chr15, names_chr16, names_chr17, names_chr18, names_chr19, names_chr20, names_chr21,
                  names_chr22, names_chrX, names_chrY)

idx_all <- lapply(names_all, function(x) which(rownames(loadings_dds) %in% x))
chr_all_pos <- lapply(idx_all, function(x) loadings_dds$PC1[x])

ave_all <- lapply(chr_all_pos, function(x) mean(x))

min_all <- lapply(chr_all_pos, function(x) min(x))    # min(unlist(min_all)) [1] -0.004601971
max_all <- lapply(chr_all_pos, function(x) max(x))    # max(unlist(max_all)) [1] 0.7285054

#idx_chr1 <- which(rownames(loadings_dds) %in% names_chr1)
#chr1_pos <- loadings_dds$PC1[idx_chr1]
  
# Basic box plot
png("boxplot_average_raw_1.png", width = 780, height = 780)
boxplot(chr_all_pos[[1]], chr_all_pos[[2]], chr_all_pos[[3]], chr_all_pos[[4]], chr_all_pos[[5]],
        chr_all_pos[[6]], chr_all_pos[[7]], chr_all_pos[[8]], chr_all_pos[[9]], chr_all_pos[[10]],
        chr_all_pos[[11]], chr_all_pos[[12]], chr_all_pos[[13]], chr_all_pos[[14]], chr_all_pos[[15]],
        chr_all_pos[[16]], chr_all_pos[[17]], chr_all_pos[[18]], chr_all_pos[[19]], chr_all_pos[[20]],
        chr_all_pos[[21]], chr_all_pos[[22]], chr_all_pos[[23]], chr_all_pos[[24]],
        #main = "Average values of PC1 loadings per chromomsome",
        xlab = "Chromosomes",
        ylab = "PC1 loading",
        names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
                  "14", "15", "16","17", "18", "19", "20", "21", "22", "X", "Y"))
dev.off()

# specyfic range: from -0.004601971 to 0.001
png("boxplot_average_raw_2.png", width = 780, height = 780)
boxplot(chr_all_pos[[1]], chr_all_pos[[2]], chr_all_pos[[3]], chr_all_pos[[4]], chr_all_pos[[5]],
        chr_all_pos[[6]], chr_all_pos[[7]], chr_all_pos[[8]], chr_all_pos[[9]], chr_all_pos[[10]],
        chr_all_pos[[11]], chr_all_pos[[12]], chr_all_pos[[13]], chr_all_pos[[14]], chr_all_pos[[15]],
        chr_all_pos[[16]], chr_all_pos[[17]], chr_all_pos[[18]], chr_all_pos[[19]], chr_all_pos[[20]],
        chr_all_pos[[21]], chr_all_pos[[22]], chr_all_pos[[23]], chr_all_pos[[24]],
        #main = "Average values of PC1 loadings per chromomsome",
        xlab = "Chromosomes",
        ylab = "PC1 loading",
        col = "grey",
        border = "black",
        ylim = c(-0.004601971, 0.001),
        names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
                 "14", "15", "16","17", "18", "19", "20", "21", "22", "X", "Y"))
dev.off()


# specyfic range: from -0.001 to 0.001
png("boxplot_average_raw_3.png", width = 780, height = 780)
boxplot(chr_all_pos[[1]], chr_all_pos[[2]], chr_all_pos[[3]], chr_all_pos[[4]], chr_all_pos[[5]],
        chr_all_pos[[6]], chr_all_pos[[7]], chr_all_pos[[8]], chr_all_pos[[9]], chr_all_pos[[10]],
        chr_all_pos[[11]], chr_all_pos[[12]], chr_all_pos[[13]], chr_all_pos[[14]], chr_all_pos[[15]],
        chr_all_pos[[16]], chr_all_pos[[17]], chr_all_pos[[18]], chr_all_pos[[19]], chr_all_pos[[20]],
        chr_all_pos[[21]], chr_all_pos[[22]], chr_all_pos[[23]], chr_all_pos[[24]],
        #main = "Average values of PC1 loadings per chromomsome",
        xlab = "Chromosomes",
        ylab = "PC1 loading",
        col = "grey",
        border = "black",
        ylim = c(-0.001, 0.001),
        names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
                  "14", "15", "16","17", "18", "19", "20", "21", "22", "X", "Y"))
dev.off()



#####################################################################################
# chromosomal locations of the results
#####################################################################################

# check chromosomal positions of the differentially expressed transcripts from Mann-Whitney analysis
# how many of them are on X and Y chromosomes

#---------- raw ----------#
pos_raw <- which(chr_positions$Approved.symbol %in% DE_genes_raw)     # 2 missing
setdiff(DE_genes_raw, as.character(chr_positions$Approved.symbol[pos_raw]))
# "LOC105369937" - chr12 | source: https://www.genecards.org
# "LOC107984931" - chr1
df_plot_raw <- data.frame(table(as.character(chr_positions$Chromosome[pos_raw])))

# need to drop levels before ordering by column
df_plot_raw_final <- data.frame(chr=as.numeric(levels(droplevels(df_plot_raw$Var1))),
                              Freq=df_plot_raw$Freq)

df_plot_raw_sorted <- df_plot_raw_final[order(df_plot_raw_final$chr),]

# add one in Freq for chr 12 and one for chr1
df_plot_raw_sorted[df_plot_raw_sorted$chr==12,]$Freq <- df_plot_raw_sorted[df_plot_raw_sorted$chr==12,]$Freq + 1
df_plot_raw_sorted[df_plot_raw_sorted$chr==1,]$Freq <- df_plot_raw_sorted[df_plot_raw_sorted$chr==1,]$Freq + 1

png("chr_pos_freq_raw.png")
ggplot(data=df_plot_raw_sorted, aes(x=chr, y=Freq)) + labs(title="raw data", x ="Chromosomal location") +
  geom_bar(stat="identity") + scale_x_discrete(limits=c(c(seq(1,22)), "X", "Y"))
dev.off()

#---------- vst ----------#
pos_vst <- which(chr_positions$Approved.symbol %in% DE_genes_vst)     # 5 missing
setdiff(DE_genes_vst, as.character(chr_positions$Approved.symbol[pos_vst]))
# "LOC105378360" - chr10
# "LOC105369937" - chr12
# "LOC101927055" - chr2
# "LOC112267886" - chrX
# "LOC107985290" - chr19

# replace X with 23 and Y with 24 (if occur)
which(as.character(chr_positions$Chromosome[pos_vst])== "X")
pos_vst_new <- as.character(chr_positions$Chromosome[pos_vst])
pos_vst_new[which(pos_vst_new == "X")] <- 23

df_plot_vst <- data.frame(table(pos_vst_new))

# need to drop levels before ordering by column
df_plot_vst_final <- data.frame(chr=as.numeric(levels(droplevels(df_plot_vst$pos_vst_new))),
                                Freq=df_plot_vst$Freq)

df_plot_vst_sorted <- df_plot_vst_final[order(df_plot_vst_final$chr),]

# "X" is replaced with 23

# add one in Freq for chr2, chr10, chr12, chr19, chrX (23)
df_plot_vst_sorted[df_plot_vst_sorted$chr==2,]$Freq <- df_plot_vst_sorted[df_plot_vst_sorted$chr==2,]$Freq + 1
df_plot_vst_sorted[df_plot_vst_sorted$chr==10,]$Freq <- df_plot_vst_sorted[df_plot_vst_sorted$chr==10,]$Freq + 1
df_plot_vst_sorted[df_plot_vst_sorted$chr==12,]$Freq <- df_plot_vst_sorted[df_plot_vst_sorted$chr==12,]$Freq + 1
df_plot_vst_sorted[df_plot_vst_sorted$chr==19,]$Freq <- df_plot_vst_sorted[df_plot_vst_sorted$chr==19,]$Freq + 1
df_plot_vst_sorted[df_plot_vst_sorted$chr==23,]$Freq <- df_plot_vst_sorted[df_plot_vst_sorted$chr==23,]$Freq + 1

png("chr_pos_freq_vst.png")
ggplot(data=df_plot_vst_sorted, aes(x=chr, y=Freq)) + labs(title="vst data", x ="Chromosomal location") +
  geom_bar(stat="identity") + scale_x_discrete(limits=c(c(seq(1,22)), "X", "Y"))
dev.off()

#---------- rlog ----------#
pos_rlog <- which(chr_positions$Approved.symbol %in% DE_genes_rlog)   # 8 missing
setdiff(DE_genes_rlog, as.character(chr_positions$Approved.symbol[pos_rlog]))
# "LOC107984931" - chr1
# "LOC101927401" - chr1
# "LOC105379404" - chr4
# "LOC105369937" - chr12
# "LOC107985290" - chr19
# "LOC101927055" - chr2
# "LOC653513"    - chr1
# "LOC105374421" - chr4
df_plot_rlog <- data.frame(table(as.character(chr_positions$Chromosome[pos_rlog])))

# replace X with 23 and Y with 24 (if occur)
which(as.character(chr_positions$Chromosome[pos_rlog])== "X")
pos_rlog_new <- as.character(chr_positions$Chromosome[pos_rlog])
pos_rlog_new[which(pos_rlog_new == "X")] <- 23


df_plot_rlog <- data.frame(table(pos_rlog_new))

# need to drop levels before ordering by column
df_plot_rlog_final <- data.frame(chr=as.numeric(levels(droplevels(df_plot_rlog$pos_rlog_new))),
                                Freq=df_plot_rlog$Freq)

df_plot_rlog_sorted <- df_plot_rlog_final[order(df_plot_rlog_final$chr),]

# "X" is replaced with 23

# add three in Freq for chr1, one for2, two for chr4, chr12, chr19,
df_plot_rlog_sorted[df_plot_rlog_sorted$chr==1,]$Freq <- df_plot_rlog_sorted[df_plot_rlog_sorted$chr==1,]$Freq + 3
df_plot_rlog_sorted[df_plot_rlog_sorted$chr==2,]$Freq <- df_plot_rlog_sorted[df_plot_rlog_sorted$chr==2,]$Freq + 1
df_plot_rlog_sorted[df_plot_rlog_sorted$chr==4,]$Freq <- df_plot_rlog_sorted[df_plot_rlog_sorted$chr==4,]$Freq + 2
df_plot_rlog_sorted[df_plot_rlog_sorted$chr==12,]$Freq <- df_plot_rlog_sorted[df_plot_rlog_sorted$chr==12,]$Freq + 1
df_plot_rlog_sorted[df_plot_rlog_sorted$chr==19,]$Freq <- df_plot_rlog_sorted[df_plot_rlog_sorted$chr==19,]$Freq + 1

png("chr_pos_freq_rlog.png")
ggplot(data=df_plot_rlog_sorted, aes(x=chr, y=Freq)) + labs(title="rlog data", x ="Chromosomal location") +
  geom_bar(stat="identity") + scale_x_discrete(limits=c(c(seq(1,22)), "X", "Y"))
dev.off()

#---------- 29 common ----------#
# there are 29 common between all 3 data types
common <- c("CA3", "CKM", "COX6A2", "MYOT", "TNNI1", "MYL1", "MYH7", "KLHL41", "ACTN2", "TNNC1",
            "NEB", "CMYA5", "NRAP", "LMOD2", "MYBPC1", "SLN", "CA3-AS1", "MB", "XIRP2", "CSRP3", 
            "MYL2", "TCAP", "TNNT1", "RYR1", "C10orf71", "LOC105369937", "ANKRD2", "EEF1A2", "MHRT")

pos_common <- which(chr_positions$Approved.symbol %in% common) 
# LOC105369937 is missing | chr12 

df_common <- data.frame(table(as.character(chr_positions$Chromosome[pos_common])))

# need to drop levels before ordering by column
df_common_final <- data.frame(chr=as.numeric(levels(droplevels(df_common$Var1))),
                              Freq=df_common$Freq)

df_common_sorted <- df_common_final[order(df_common_final$chr),]

# add one in Freq for chr 12 (because of the missing one)
df_common_sorted[df_common_sorted$chr==12,]$Freq <- 3

png("chr_pos_freq_common.png")
ggplot(data=df_common_sorted, aes(x=chr, y=Freq)) + labs(title="29 common genes between all data types", x ="Chromosomal location") +
  geom_bar(stat="identity") + scale_x_discrete(limits=c(c(seq(1,22)), "X", "Y"))
dev.off()




