
# split df in 
# => chr 1-22
# => chr X and Y
# => all chr (no chenges needed)

# load lists with positions on chrX and chrY
chrX_IDs <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/chrX_positions.csv"), row.names = 1)
chrY_IDs <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/chrY_positions.csv"), row.names = 1)

# get union of chrX (l=2276) and chrY (l=303)
both_chr_IDs <- union(chrX_IDs$x, chrY_IDs$x) # l=2516

# subset df for chrX and chrY only
df_chrXY <- df[which(rownames(df) %in% both_chr_IDs),]
counts_mat_chrXY <- as.matrix(df_chrXY)                       
storage.mode(counts_mat_chrXY) <- "integer"                   # type: integer | dim: 2516    41

# subset df for chr1-22 only
df_chr1_22 <- df[-which(rownames(df) %in% both_chr_IDs),]
counts_mat_chr1_22 <- as.matrix(df_chr1_22)                   
storage.mode(counts_mat_chr1_22) <- "integer"                 # type: integer | dim: 56092    41

# df includes all chromosomes

dds_chrXY <- DESeqDataSetFromMatrix(countData = counts_mat_chrXY, colData = coldata, design = ~ gender)
dds_chr1_22 <- DESeqDataSetFromMatrix(countData = counts_mat_chr1_22, colData = coldata, design = ~ gender)


keep_chrXY <- rowSums(counts(dds_chrXY)) >= 10
dds_chrXY <- dds_chrXY[keep_chrXY,]                     # 2516 - 2249 = 267 rows dropped

keep_chr1_22 <- rowSums(counts(dds_chr1_22)) >= 10
dds_chr1_22 <- dds_chr1_22[keep_chr1_22,]               # 56092 - 49990 = 6102 rows dropped

dds_chrXY <- DESeq(dds_chrXY)
res_chrXY <- results(dds_chrXY)       # dim = 52239     6

dds_chr1_22 <- DESeq(dds_chr1_22)
res_chr1_22 <- results(dds_chr1_22)       # dim = 52239     6

counts_norm_chrXY <- counts(dds_chrXY, normalized = TRUE)
#write.csv(counts_norm_chrXY, file="norm_counts_chrXY.csv")

counts_norm_chr1_22 <- counts(dds_chr1_22, normalized = TRUE)
#write.csv(counts_norm_chr1_22, file="norm_counts_chr1_22.csv")

# unnormalised
dds_chr1_22
dds_chrXY

# it takes too long, especially for rlog, so just load saved data (below)
# vst (variance stabilizing transformation)
vst_chr1_22 <- vst(dds_chr1_22, blind=FALSE)
vst_chrXY   <- vst(dds_chrXY, blind=FALSE)

# rlog (it takes quite long time to run)
rlog_chr1_22   <- rlog(dds_chr1_22, blind=FALSE)
rlog_chrXY     <- rlog(dds_chrXY, blind=FALSE)

# save transformed data and unnormalised as well

#save(dds_chr1_22, file="Raw_DESeq_dataset_chr1_22.Rda")
#save(dds_chrXY, file="Raw_DESeq_dataset_chrXY.Rda")

#save(vst_chr1_22, file="Normalised_DESeq_vst_dataset_chr1_22.Rda")
#save(vst_chrXY, file="Normalised_DESeq_vst_dataset_chrXY.Rda")

#save(rlog_chr1_22, file= "Normalised_DESeq_rlog_dataset_chr1_22.Rda")
#save(rlog_chrXY, file="Normalised_DESeq_rlog_dataset_chrXY.Rda")

#----------------------------------------------------------------------------------------------------

#load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_vst_dataset_chr1_22.Rda"))
#load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_vst_dataset_chrXY.Rda"))

#load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_rlog_dataset_chr1_22.Rda"))
#load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_rlog_dataset_chrXY.Rda"))

#dds_chrXY_trans     <- DESeqTransform(dds_chrXY)
#dds_chr1_22_trans   <- DESeqTransform(dds_chr1_22)

#vst_chr1_22_trans   <- DESeqTransform(vst_chr1_22)
#vst_chrXY_trans     <- DESeqTransform(vst_chrXY)

#rlog_chr1_22_trans   <- DESeqTransform(rlog_chr1_22)
#rlog_chrXY_trans     <- DESeqTransform(rlog_chrXY)

#p2 <- plotPCA(dds_chrXY_trans, intgroup=c("gender")) + ggtitle("raw: chrX & chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
#p3 <- plotPCA(dds_chr1_22_trans, intgroup=c("gender")) + ggtitle("raw: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)

#p5 <- plotPCA(vst_chrXY_trans, intgroup=c("gender")) + ggtitle("vst: chrX & chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
#p6 <- plotPCA(vst_chr1_22_trans, intgroup=c("gender")) + ggtitle("vst: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)

#p8 <- plotPCA(rlog_chrXY_trans, intgroup=c("gender")) + ggtitle("rlog: chrX - chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)
#p9 <- plotPCA(rlog_chr1_22_trans, intgroup=c("gender")) + ggtitle("rlog: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)

# plots without titles

#p2_bis <- plotPCA(dds_chrXY_trans, intgroup=c("gender")) + theme(legend.position = "none")
#p3_bis <- plotPCA(dds_chr1_22_trans, intgroup=c("gender")) + theme(legend.position = "none")

# vst => variance stabilizing transformation
#p5_bis <- plotPCA(vst_chrXY_trans, intgroup=c("gender")) + theme(legend.position = "none")
#p6_bis <- plotPCA(vst_chr1_22_trans, intgroup=c("gender")) + theme(legend.position = "none")

# save plots
#pdf(file="PCA_plot_gender_raw_chrXY.pdf")   ; p2; dev.off()
#pdf(file="PCA_plot_gender_raw_chr1_22.pdf") ; p3; dev.off()

#pdf(file="PCA_plot_gender_vst_chrXY.pdf")   ; p5; dev.off()
#pdf(file="PCA_plot_gender_vst_chr1_22.pdf") ; p6; dev.off()

#pdf(file="PCA_plot_gender_rlog_chrXY.pdf")   ; p8; dev.off()
#pdf(file="PCA_plot_gender_rlog_chr1_22.pdf") ; p9; dev.off()


