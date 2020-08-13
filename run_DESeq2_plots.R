# this script performs plots on DESeq2 analysis output
# data types: Raw | Normalised_rlog | Normalised_vst

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(grid)
library(gridExtra)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

# define directory with data
data_dir <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/no_chrXY")

# load data RAW | VST | rlog
load(paste0(data_dir, "/Raw_DESeq_dataset_noXY.Rda"))
load(paste0(data_dir,"/Normalised_DESeq_vst_dataset_noXY.Rda"))
load(paste0(data_dir,"/Normalised_DESeq_rlog_dataset_noXY.Rda"))

# all the objects need to be passed through DESeqTransform() function before PCAplot()
dds_all_trans       <- DESeqTransform(dds_all)
vst_all_trans       <- DESeqTransform(vst_all)
rlog_all_trans      <- DESeqTransform(rlog_all)

#colnames(dds_all@colData) => labels for PCA plots:

#--- "visual_loss"
p1_visual_loss <- plotPCA(dds_all_trans, intgroup=c("visual_loss")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 2)
p1_visual_loss$labels$colour <- "visual_loss"

p2_visual_loss <- plotPCA(vst_all_trans, intgroup=c("visual_loss")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 2)
p2_visual_loss$labels$colour <- "visual_loss"

p3_visual_loss <- plotPCA(rlog_all_trans, intgroup=c("visual_loss")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 2) 
p3_visual_loss$labels$colour <- "visual_loss"

#--- "jaw_claudication"
p1_jaw_claudication <- plotPCA(dds_all_trans, intgroup=c("jaw_claudication")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 6)
p1_jaw_claudication$labels$colour <- "jaw_claudication"

p2_jaw_claudication <- plotPCA(vst_all_trans, intgroup=c("jaw_claudication")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 6)
p2_jaw_claudication$labels$colour <- "jaw_claudication"

p3_jaw_claudication <- plotPCA(rlog_all_trans, intgroup=c("jaw_claudication")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 6) 
p3_jaw_claudication$labels$colour <- "jaw_claudication"

#--- "ischaemic_features"
p1_ischaemic_features <- plotPCA(dds_all_trans, intgroup=c("ischaemic_features")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 7)
p1_ischaemic_features$labels$colour <- "ischaemic_features"

p2_ischaemic_features <- plotPCA(vst_all_trans, intgroup=c("ischaemic_features")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 7)
p2_ischaemic_features$labels$colour <- "ischaemic_features"

p3_ischaemic_features <- plotPCA(rlog_all_trans, intgroup=c("ischaemic_features")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) + scale_color_brewer(type = 'qual', palette = 7) 
p3_ischaemic_features$labels$colour <- "ischaemic_features"

#--- "gender"
p1_gender <- plotPCA(dds_all_trans, intgroup=c("gender")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender$labels$colour <- "gender"

p2_gender <- plotPCA(vst_all_trans, intgroup=c("gender")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender$labels$colour <- "gender"

p3_gender <- plotPCA(rlog_all_trans, intgroup=c("gender")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) 
p3_gender$labels$colour <- "gender"

#--- "year_TAB_collected"
p1_year_TAB_collected <- plotPCA(dds_all_trans, intgroup=c("year_TAB_collected")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) #+ scale_color_brewer(type = 'qual', palette = 1)
p1_year_TAB_collected$labels$colour <- "year_TAB_collected"

p2_year_TAB_collected <- plotPCA(vst_all_trans, intgroup=c("year_TAB_collected")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) #+ scale_color_brewer(type = 'qual', palette = 1)
p2_year_TAB_collected$labels$colour <- "year_TAB_collected"

p3_year_TAB_collected <- plotPCA(rlog_all_trans, intgroup=c("year_TAB_collected")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) #+ scale_color_brewer(type = 'qual', palette = 1) 
p3_year_TAB_collected$labels$colour <- "year_TAB_collected"

#--- "days_on_steroids"
p1_days_on_steroids <- plotPCA(dds_all_trans, intgroup=c("days_on_steroids")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_days_on_steroids$labels$colour <- "days_on_steroids"

p2_days_on_steroids <- plotPCA(vst_all_trans, intgroup=c("days_on_steroids")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) 
p2_days_on_steroids$labels$colour <- "days_on_steroids"

p3_days_on_steroids <- plotPCA(rlog_all_trans, intgroup=c("days_on_steroids")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p3_days_on_steroids$labels$colour <- "days_on_steroids"

#--- "age"   
p1_age <- plotPCA(dds_all_trans, intgroup=c("age")) + ggtitle("raw: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_age$labels$colour <- "age"

p2_age <- plotPCA(vst_all_trans, intgroup=c("age")) + ggtitle("vst: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) 
p2_age$labels$colour <- "age"

p3_age <- plotPCA(rlog_all_trans, intgroup=c("age")) + ggtitle("rlog: no chrX - chrY") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p3_age$labels$colour <- "age"


#--- "visual_loss"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_visual_loss$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_visual_loss$labels$colour))

# save plots
png(file=paste0(p1_visual_loss$labels$colour, "_raw_noXY_chr.png")); p1_visual_loss; dev.off()
png(file=paste0(p1_visual_loss$labels$colour, "_raw_noXY_chr_label.png"))
p1_visual_loss + geom_point(size = 3) + geom_text(label=p1_visual_loss$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_visual_loss$labels$colour, "_VST_noXY_chr.png")); p2_visual_loss; dev.off()
png(file=paste0(p2_visual_loss$labels$colour, "_VST_noXY_chr_label.png"))
p2_visual_loss + geom_point(size = 3) + geom_text(label=p2_visual_loss$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_visual_loss$labels$colour, "_rlog_noXY_chr.png")); p3_visual_loss; dev.off()
png(file=paste0(p3_visual_loss$labels$colour, "_rlog_noXY_chr_label.png"))
p3_visual_loss + geom_point(size = 3) + geom_text(label=p3_visual_loss$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

#--- "jaw_claudication"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_jaw_claudication$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_jaw_claudication$labels$colour))

# save plots
png(file=paste0(p1_jaw_claudication$labels$colour, "_raw_noXY_chr.png")); p1_jaw_claudication; dev.off()
png(file=paste0(p1_jaw_claudication$labels$colour, "_raw_noXY_chr_label.png"))
p1_jaw_claudication + geom_point(size = 3) + geom_text(label=p1_jaw_claudication$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_jaw_claudication$labels$colour, "_VST_noXY_chr.png")); p2_jaw_claudication; dev.off()
png(file=paste0(p2_jaw_claudication$labels$colour, "_VST_noXY_chr_label.png"))
p2_jaw_claudication + geom_point(size = 3) + geom_text(label=p2_jaw_claudication$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_jaw_claudication$labels$colour, "_rlog_noXY_chr.png")); p3_jaw_claudication; dev.off()
png(file=paste0(p3_jaw_claudication$labels$colour, "_rlog_noXY_chr_label.png"))
p3_jaw_claudication + geom_point(size = 3) + geom_text(label=p3_jaw_claudication$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

#--- "ischaemic_features"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_ischaemic_features$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_ischaemic_features$labels$colour))

# save plots
png(file=paste0(p1_ischaemic_features$labels$colour, "_raw_noXY_chr.png")); p1_ischaemic_features; dev.off()
png(file=paste0(p1_ischaemic_features$labels$colour, "_raw_noXY_chr_label.png"))
p1_ischaemic_features + geom_point(size = 3) + geom_text(label=p1_ischaemic_features$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_ischaemic_features$labels$colour, "_VST_noXY_chr.png")); p2_ischaemic_features; dev.off()
png(file=paste0(p2_ischaemic_features$labels$colour, "_VST_noXY_chr_label.png"))
p2_ischaemic_features + geom_point(size = 3) + geom_text(label=p2_ischaemic_features$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_ischaemic_features$labels$colour, "_rlog_noXY_chr.png")); p3_ischaemic_features; dev.off()
png(file=paste0(p3_ischaemic_features$labels$colour, "_rlog_noXY_chr_label.png"))
p3_ischaemic_features + geom_point(size = 3) + geom_text(label=p3_ischaemic_features$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

#--- "gender"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_gender$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_gender$labels$colour))

# save plots
png(file=paste0(p1_gender$labels$colour, "_raw_noXY_chr.png")); p1_gender; dev.off()
png(file=paste0(p1_gender$labels$colour, "_raw_noXY_chr_label.png"))
p1_gender + geom_point(size = 3) + geom_text(label=p1_gender$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_gender$labels$colour, "_VST_noXY_chr.png")); p2_gender; dev.off()
png(file=paste0(p2_gender$labels$colour, "_VST_noXY_chr_label.png"))
p2_gender + geom_point(size = 3) + geom_text(label=p2_gender$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_gender$labels$colour, "_rlog_noXY_chr.png")); p3_gender; dev.off()
png(file=paste0(p3_gender$labels$colour, "_rlog_noXY_chr_label.png"))
p3_gender + geom_point(size = 3) + geom_text(label=p3_gender$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

#--- "year_TAB_collected"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_year_TAB_collected$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_year_TAB_collected$labels$colour))

# save plots
png(file=paste0(p1_year_TAB_collected$labels$colour, "_raw_noXY_chr.png")); p1_year_TAB_collected; dev.off()
png(file=paste0(p1_year_TAB_collected$labels$colour, "_raw_noXY_chr_label.png"))
p1_year_TAB_collected + geom_point(size = 3) + geom_text(label=p1_year_TAB_collected$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_year_TAB_collected$labels$colour, "_VST_noXY_chr.png")); p2_year_TAB_collected; dev.off()
png(file=paste0(p2_year_TAB_collected$labels$colour, "_VST_noXY_chr_label.png"))
p2_year_TAB_collected + geom_point(size = 3) + geom_text(label=p2_days_on_steroids$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_year_TAB_collected$labels$colour, "_rlog_noXY_chr.png")); p3_year_TAB_collected; dev.off()
png(file=paste0(p3_year_TAB_collected$labels$colour, "_rlog_noXY_chr_label.png"))
p3_year_TAB_collected + geom_point(size = 3) + geom_text(label=p3_year_TAB_collected$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

#--- "days_on_steroids"
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_days_on_steroids$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_days_on_steroids$labels$colour))

# save plots
png(file=paste0(p1_days_on_steroids$labels$colour, "_raw_noXY_chr.png")); p1_days_on_steroids; dev.off()
png(file=paste0(p1_days_on_steroids$labels$colour, "_raw_noXY_chr_label.png"))
p1_days_on_steroids + geom_point(size = 3) + geom_text(label=p1_days_on_steroids$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_days_on_steroids$labels$colour, "_VST_noXY_chr.png")); p2_days_on_steroids; dev.off()
png(file=paste0(p2_days_on_steroids$labels$colour, "_VST_noXY_chr_label.png"))
p2_days_on_steroids + geom_point(size = 3) + geom_text(label=p2_days_on_steroids$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_days_on_steroids$labels$colour, "_rlog_noXY_chr.png")); p3_days_on_steroids; dev.off()
png(file=paste0(p3_days_on_steroids$labels$colour, "_rlog_noXY_chr_label.png"))
p3_days_on_steroids + geom_point(size = 3) + geom_text(label=p3_days_on_steroids$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()


#--- "age" 
# create diectory for output plots
dir.create(paste0(data_dir,"/PCA_plots/", p1_age$labels$colour))
setwd(paste0(data_dir,"/PCA_plots/", p1_age$labels$colour))

# save plots
png(file=paste0(p1_age$labels$colour, "_raw_noXY_chr.png")); p1_age; dev.off()
png(file=paste0(p1_age$labels$colour, "_raw_noXY_chr_label.png"))
p1_age + geom_point(size = 3) + geom_text(label=p1_age$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p2_age$labels$colour, "_VST_noXY_chr.png")); p2_age; dev.off()
png(file=paste0(p2_age$labels$colour, "_VST_noXY_chr_label.png"))
p2_age + geom_point(size = 3) + geom_text(label=p2_age$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

png(file=paste0(p3_age$labels$colour, "_rlog_noXY_chr.png")); p3_age; dev.off()
png(file=paste0(p3_age$labels$colour, "_rlog_noXY_chr_label.png"))
p3_age + geom_point(size = 3) + geom_text(label=p3_age$data$name, check_overlap=F, nudge_x=5, nudge_y=0, size=3)
dev.off()

# extract the legend and plot it separately
#p0 <- plotPCA(dds_all_trans, intgroup=c("gender"))
#leg <- get_legend(p0)
#pdf(file="PCA_legend.pdf", width=3, height=3)
#plot(leg)
#dev.off()

# there's more code for adding sample names as labels on the plots, etc









##################################################################################################

# extract transformed values (Variance stabilizing transformation)
vst <- vst(dds, blind=FALSE)
rlog <- rlog(dds, blind=FALSE)
#head(assay(vst), 3)

ntd <- normTransform(dds)

library(vsn)

meanSdPlot(assay(ntd))    # Normalized counts transformation
meanSdPlot(assay(vst))    # Variance stabilizing transformation
meanSdPlot(assay(rlog))    # Regularized log transformation

# save other count matrices as well
write.csv(assay(ntd), file="normalized_counts_transformation_counts.csv")
write.csv(assay(vst), file="variance_stabilizing_transformation_counts.csv")
write.csv(assay(rlog), file="regularized_log_transformation_counts.csv")

#-----------------------------------------------------------------------------------------------#
# copied from 'run_Distribution_plot.R'

# raw counts (as they came from featureCounts())
data.raw <- counts_mat

# raw filtered counts (pre-filtering applied to keep only rows that have at least 10 reads total)
data.raw.filtered <- counts(dds)

# logarithm transformation => it will get rid of some extreme values. 
data.log2.on.raw <- log2(data.raw + 1)

# logarithm transformation => it will get rid of some extreme values. 
data.log2.on.filtered <- log2(data.raw.filtered + 1)

# variance-stabilizing transformation (VST), implemented in the DESeq package (Anders and Huber, 2010)
data.vst.on.raw <- vst(data.raw)

# variance-stabilizing transformation (VST), implemented in the DESeq package (Anders and Huber, 2010)
data.vst.on.filtered <- vst(data.raw.filtered)

# normalised counts (by DESeq)
data.norm <- counts_norm

# other DESeq2 transformations
# assay(ntd)
# assay(vst)
# assay(rlog)

# visualise 4 genes: 1st, 2nd, 14th and 18th; can't just take random genes as there are many with counts around 0
# data.raw            |     data.log2.on.raw        |     data.vst.on.raw         | dim= 58608    41
# data.raw.filtered   |     data.log2.on.filtered   |     data.vst.on.filtered    | dim= 52239    41

# data.norm     | dim= 52239    41

### FIRST PLOT
pdf(file="density_plots_DESeq2-p1.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.raw[1,])), main=paste0("raw - ", rownames(data.raw)[1]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[1,])), main=paste0("log2 - ", rownames(data.log2.on.raw)[1]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[1,])), main=paste0("VST - ", rownames(data.vst.on.raw)[1]), cex.main=1)

plot(density(as.numeric(data.raw[14,])), main=paste0("raw - ", rownames(data.raw)[14]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[14,])), main=paste0("log2 - ", rownames(log2.on.raw)[14]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[14,])), main=paste0("VST - ", rownames(data.vst.on.raw)[14]), cex.main=1)

plot(density(as.numeric(data.raw[18,])), main=paste0("raw - ", rownames(data.raw)[18]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[18,])), main=paste0("log2 - ", rownames(log2.on.raw)[18]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[18,])), main=paste0("VST - ", rownames(data.raw)[18]), cex.main=1)
dev.off()

### SECOND PLOT
pdf(file="density_plots_DESeq2-p2.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.raw.filtered[1,])), main=paste0("raw (filtered) - ", rownames(data.raw.filtered)[1]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[1,])), main=paste0("log2 (filtered) - ", rownames(data.log2.on.filtered)[1]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[1,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[1]), cex.main=1)

plot(density(as.numeric(data.raw[14,])), main=paste0("raw (filtered) - ", rownames(data.raw)[14]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[14,])), main=paste0("log2 (filtered) - ", rownames(data.log2.on.filtered)[14]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[14,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[14]), cex.main=1)

plot(density(as.numeric(data.raw[18,])), main=paste0("raw (filtered) - ", rownames(data.raw)[18]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[18,])), main=paste0("log2 (filtered - ", rownames(data.log2.on.filtered)[18]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[18,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[18]), cex.main=1)
dev.off()

### THIRD PLOT
pdf(file="density_plots_DESeq2-p3.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.norm[1,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[1]), cex.main=1,)
plot(density(as.numeric(assay(vst)[1,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[1]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[1,])), main=paste0("transformed R-log (DESeq2) - ", rownames(assay(rlog))[1]), cex.main=1,)

plot(density(as.numeric(data.norm[14,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[14]), cex.main=1,)
plot(density(as.numeric(assay(vst)[14,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[14]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[14,])), main=paste0("transformed R-log (DESeq2) - ", rownames(assay(rlog))[14]), cex.main=1,)

plot(density(as.numeric(data.norm[18,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[18]), cex.main=1,)
plot(density(as.numeric(assay(vst)[18,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[18]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[18,])), main=paste0("transformed Rilog (DESeq2) - ", rownames(assay(rlog))[18]), cex.main=1,)
dev.off()

cat("Finished!")
cat("Created: density_plots_DESeq2-p1.pdf")
cat("Created: density_plots_DESeq2-p2.pdf")
cat("Created: density_plots_DESeq2-p3.pdf")

#-----------------------------------------------------------------------------------------------#

