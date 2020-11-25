# script to plot QC metrics from SE run against QC metrics from PE run
# "uniquely mapped reads" (STAR)
# "Assigned" (featureCounts)
# "Unassigned - no feature" (featureCounts)

# install (if necessary) and load package
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

# load table with QC metrics
# for STAR (_mod => removed "%" from each column manually)
df_star_se <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/4_alignment/SE/STAR_stats_all_SE_mod.csv"), row.names = 1)
df_star_pe <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/4_alignment/PE/STAR_stats_all_PE_mod.csv"), row.names = 1)

# for featureCounts()
df_fcounts_se <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/6_featCounts/SE/all_stats_nodups_run1_SE.csv"), row.names = 1)
df_fcounts_pe <- read.csv(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/6_featCounts/PE/all_stats_nodups_run1_PE.csv"), row.names = 1)

# check if colnames match between df_X_se and df_X_pe
any(colnames(df_star_se) == colnames(df_star_pe))
any(colnames(df_fcounts_se) == colnames(df_fcounts_se))

setwd(paste0(main_dir, "/ANALYSES_archived/QC_results_RNA-seq_pipeline/comparison_plots"))

######### STAR #########

# get vectors with "Uniquely mapped reads number" => NOTE: need to use percentages
#vec_star_se <- df_star_se[which(rownames(df_star_se) == "Uniquely mapped reads number"),]
#vec_star_pe <- df_star_pe[which(rownames(df_star_pe) == "Uniquely mapped reads number"),]
#df_star <- data.frame(se=as.numeric(as.matrix(vec_star_se)),
#                      pe=as.numeric(as.matrix(vec_star_pe)))

vec_star_se_mapped <- df_star_se[which(rownames(df_star_se) == "Uniquely mapped reads %"),]
vec_star_pe_mapped <- df_star_pe[which(rownames(df_star_pe) == "Uniquely mapped reads %"),]

vec_star_se_mismatch <- df_star_se[which(rownames(df_star_se) == "Mismatch rate per base, %" ),]
vec_star_pe_mismatch <- df_star_pe[which(rownames(df_star_pe) == "Mismatch rate per base, %" ),]

vec_star_se_multiple_loci <- df_star_se[which(rownames(df_star_se) == "% of reads mapped to multiple loci"),]
vec_star_pe_multiple_loci <- df_star_pe[which(rownames(df_star_pe) == "% of reads mapped to multiple loci"),]

vec_star_se_too_many_loci <- df_star_se[which(rownames(df_star_se) == "% of reads mapped to too many loci"),]
vec_star_pe_too_many_loci <- df_star_pe[which(rownames(df_star_pe) == "% of reads mapped to too many loci"),]

# vec_star_se_too_many_mismatch <- df_star_se[which(rownames(df_star_se) == "% of reads unmapped: too many mismatches"),] # zeros everywhere
# vec_star_pe_too_many_mismatch <- df_star_pe[which(rownames(df_star_pe) == "% of reads unmapped: too many mismatches"),] # zeros everywhere

vec_star_se_unmapped_too_short <- df_star_se[which(rownames(df_star_se) == "% of reads unmapped: too short"),]
vec_star_pe_unmapped_too_short <- df_star_pe[which(rownames(df_star_pe) == "% of reads unmapped: too short"),]

vec_star_se_unmapped_other <- df_star_se[which(rownames(df_star_se) == "% of reads unmapped: other"),]
vec_star_pe_unmapped_other <- df_star_pe[which(rownames(df_star_pe) == "% of reads unmapped: other"),]

# vec_star_se_chimeric <- df_star_se[which(rownames(df_star_se) == "% of chimeric reads" ),] # zeros everywhere
# vec_star_pe_chimeric <- df_star_pe[which(rownames(df_star_pe) == "% of chimeric reads" ),] # zeros everywhere

df_star_all <- data.frame(se_mapped=as.numeric(as.matrix(vec_star_se_mapped)),
                          pe_mapped=as.numeric(as.matrix(vec_star_pe_mapped)),
                          se_mismatch=as.numeric(as.matrix(vec_star_se_mismatch)),
                          pe_mismatch=as.numeric(as.matrix(vec_star_pe_mismatch)),
                          se_multiple_loci=as.numeric(as.matrix(vec_star_se_multiple_loci)),
                          pe_multiple_loci=as.numeric(as.matrix(vec_star_pe_multiple_loci)),
                          se_too_many_loci=as.numeric(as.matrix(vec_star_se_too_many_loci)),
                          pe_too_many_loci=as.numeric(as.matrix(vec_star_pe_too_many_loci)),
                          #se_too_many_mismatch=as.numeric(as.matrix(vec_star_se_too_many_mismatch)),
                          #pe_too_many_mismatch=as.numeric(as.matrix(vec_star_pe_too_many_mismatch)),
                          se_unmapped_too_short=as.numeric(as.matrix(vec_star_se_unmapped_too_short)),
                          pe_unmapped_too_short=as.numeric(as.matrix(vec_star_pe_unmapped_too_short)),
                          se_unmapped_other=as.numeric(as.matrix(vec_star_se_unmapped_other)),
                          pe_unmapped_other=as.numeric(as.matrix(vec_star_pe_unmapped_other)))
                          #se_chimeric=as.numeric(as.matrix(vec_star_se_chimeric)),
                          #pe_chimeric=as.numeric(as.matrix(vec_star_pe_chimeric)))

# plot one against another
png("STAR_SE-vs-PE_percentages.png")
ggplot(df_star_all, aes(se_mapped, pe_mapped)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of uniquely mapped reads during alignment - SE vs. PE")
dev.off()

png("STAR_SE-vs-PE_mismatch.png")
ggplot(df_star_all, aes(se_mismatch, pe_mismatch)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of mismatch rate per base during alignment - SE vs. PE")
dev.off()

png("STAR_SE-vs-PE_multiple_loci.png")
ggplot(df_star_all, aes(se_multiple_loci, pe_multiple_loci)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of reads mapped to multiple loci during alignment - SE vs. PE")
dev.off()

png("STAR_SE-vs-PE_too_many_loci.png")
ggplot(df_star_all, aes(se_too_many_loci, pe_too_many_loci)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of reads mapped to too many loci during alignment - SE vs. PE")
dev.off()

png("STAR_SE-vs-PE_unmapped_too_short.png")
ggplot(df_star_all, aes(se_unmapped_too_short, pe_unmapped_too_short)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of unmapped too short reads during alignment - SE vs. PE")
dev.off()

png("STAR_SE-vs-PE_unmapped_other.png")
ggplot(df_star_all, aes(se_unmapped_other, pe_unmapped_other)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of unmapped 'other' reads during alignment - SE vs. PE")
dev.off()

##### featureCounts ####

# get vectors with "Assigned" and "Unassigned - no feature"
vec_fcounts_se_1 <- df_fcounts_se[which(df_fcounts_se$Status == "Assigned"),]
vec_fcounts_pe_1 <- df_fcounts_pe[which(df_fcounts_se$Status == "Assigned"),]

vec_fcounts_se_2 <- df_fcounts_se[which(df_fcounts_se$Status == "Unassigned_NoFeatures"),]
vec_fcounts_pe_2 <- df_fcounts_pe[which(df_fcounts_se$Status == "Unassigned_NoFeatures"),]

df_fcounts_all <- data.frame(se_assigned=as.numeric(as.matrix(vec_fcounts_se_1[2:42])),
                             pe_assigned=as.numeric(as.matrix(vec_fcounts_pe_1[2:42])),
                             se_unassigned=as.numeric(as.matrix(vec_fcounts_se_2[2:42])),
                             pe_unassigned=as.numeric(as.matrix(vec_fcounts_pe_2[2:42])))

#png("featureCounts_SE-vs-PE_assigned.png")
#ggplot(df_fcounts_all, aes(se_assigned, pe_assigned)) + geom_point() + geom_abline(colour = "brown") + 
#  xlab("single-end") + ylab("paired-end") + ggtitle("Number of assigned reads during quantification - SE vs. PE")
#dev.off()

#png("featureCounts_SE-vs-PE_unassigned.png")
#ggplot(df_fcounts_all, aes(se_unassigned, pe_unassigned)) + geom_point() + geom_abline(colour = "brown") + 
#  xlab("single-end") + ylab("paired-end") + ggtitle("Number of unassigned reads during quantification - SE vs. PE")
#dev.off()

# get percentages of
# => unassigned / assigned [SE & PE]
# => assigned / total [SE & PE]
# => unassigned / total [SE & PE]

df_fcounts_all_percentages <- data.frame(se_unassigned_assigned=(df_fcounts_all$se_unassigned/df_fcounts_all$se_assigned)*100,
                                         pe_unassigned_assigned=(df_fcounts_all$pe_unassigned/df_fcounts_all$pe_assigned)*100,
                                         se_assigned_total=(df_fcounts_all$se_assigned/(df_fcounts_all$se_assigned+df_fcounts_all$se_unassigned))*100,
                                         pe_assigned_total=(df_fcounts_all$pe_assigned/(df_fcounts_all$pe_assigned+df_fcounts_all$pe_unassigned))*100,
                                         se_unassigned_total=(df_fcounts_all$se_unassigned/(df_fcounts_all$se_assigned+df_fcounts_all$se_unassigned))*100,
                                         pe_unassigned_total=(df_fcounts_all$se_unassigned/(df_fcounts_all$pe_assigned+df_fcounts_all$pe_unassigned))*100)
                                           
png("featureCounts_SE-vs-PE_percentage_unassigned_assigned.png")
ggplot(df_fcounts_all_percentages, aes(se_unassigned_assigned, pe_unassigned_assigned)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of unassigned / assigned - SE vs. PE")
dev.off()

png("featureCounts_SE-vs-PE_percentage_assigned_total.png")
ggplot(df_fcounts_all_percentages, aes(se_assigned_total, pe_assigned_total)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of assigned / total - SE vs. PE")
dev.off()                              

png("featureCounts_SE-vs-PE_percentage_unassigned_total.png")
ggplot(df_fcounts_all_percentages, aes(se_unassigned_total, pe_unassigned_total)) + geom_point() + geom_abline(colour = "brown") + 
  xlab("single-end") + ylab("paired-end") + ggtitle("Percentage of unassigned / total - SE vs. PE")
dev.off()       



