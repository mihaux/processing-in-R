# get expression profiles for a particular molecules (genes or transcripts)

# plots:
# (1) heatmap to show raw expression level for a particular molecule [gene/transcript]
# (2) heatmap to show log-transformed expression level for a particular molecule [gene/transcript]

# (3) boxplot to compare mean expression for one feature (TO BE VERIFIED)
# (4) heatmap to show expression level for one feature (TO BE VERIFIED)

library(ggplot2)

# load the parameter file
#source("/Users/ummz/Documents/OneDrive - University of Leeds/data/plots_for_Sarah/parameter-Sarah.R")
source("/Users/ummz/Documents/OneDrive - University of Leeds/colaborations/colaboration_with_Jason/parameter-Jason.R")

# define wheather output files should be saved or not (NOT IN USE FOR NOW)
output_save <- FALSE

###-------------- load counts data ---------------###

# load "all" counts on transcript- and gene-level
cts_gene <- read.csv(dir_gene, row.names = 1, stringsAsFactors = FALSE, check.names = F)                #  37,788 x 41 
cts_transcript <- read.csv(dir_transcript, row.names = 1, stringsAsFactors = FALSE, check.names = F)    # 175,775 x 41

# create log-transformed counts
cts_gene_log <- log2(cts_gene + 1)
cts_transcript_log <- log2(cts_transcript + 1)

###-------- genes/transcripts of interest --------###

# get the list of all gene names
all_gene_names <- unique(df_info$Gene)

# get unique list of gene and transcript IDs
gene_list_id <- unique(df_info$Ensembl_gene_ID)
transcript_list_id <- unique(df_info$Ensembl_transcript_ID)

# >>> check if there are any genes or transcripts that are not present in the dataset <<< #

if ( file.exists(file.path(dir_out, "genes-transcripts_not_found")) ){
  cat("Output directory already exists.")
} else {
  dir.create( file.path(dir_out, "genes-transcripts_not_found") )
}

genes_found <- rownames(cts_gene)[which(rownames(cts_gene) %in% gene_list_id)]
genes_NOT_found <- setdiff(gene_list_id, genes_found)

if ( length(genes_NOT_found) > 0 ) {
  
  # get their corresponding transcript names
  df_info_subset_genes_not_found <- df_info[which(gene_list_id %in% genes_NOT_found),]
  
  # save a list of genes that are not present in the dataset, with their corresponding transcript names
  write.csv(df_info_subset_genes_not_found[,c("Gene", "Ensembl_gene_ID", "Ensembl_transcript_ID")], file = file.path(dir_out, "genes-transcripts_not_found", "genes_NOT_found.csv"))
  
}

transcripts_found <- rownames(cts_transcript)[which(rownames(cts_transcript) %in% transcript_list_id)]
transcripts_NOT_found <- setdiff(transcript_list_id, transcripts_found)

if ( length(transcripts_NOT_found) > 0 ) {

  # get their corresponding transcript names
  df_info_subset_transcripts_not_found <- df_info[which(transcript_list_id %in% transcripts_NOT_found),]

  # save a list of transcripts that are not present in the dataset, with their corresponding gene names
  write.csv(df_info_subset_transcripts_not_found[,c("Gene", "Ensembl_gene_ID", "Ensembl_transcript_ID")], file = file.path(dir_out, "genes-transcripts_not_found", "transcripts_NOT_found.csv"))

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# extract all genes of interest
gene_extracted <- cts_gene[which(rownames(cts_gene) %in% gene_list_id),]
gene_extracted_log <- cts_gene_log[which(rownames(cts_gene) %in% gene_list_id),]

# check which genes / transcripts have no expression across all the samples
# rowSums(gene_extracted);        rowSums(gene_extracted_log)
# rowSums(transcript_extracted);  rowSums(transcript_extracted_log)

if ( file.exists(file.path(dir_out, "genes-transcripts_no_expression")) ){
  cat("Output directory already exists.")
} else {
  dir.create( file.path(dir_out, "genes-transcripts_no_expression") )
}

if( length(which(rowSums(gene_extracted) == 0)) ) {
  
  # get their corresponding transcript names
  df_info_subset_zero_exp_genes <- df_info[which(gene_list_id %in% names(which(rowSums(gene_extracted) == 0))),]
  
  # save a list of genes that are not expressed at all, with their corresponding transcript names
  write.csv(df_info_subset_zero_exp_genes[,c("Gene", "Ensembl_gene_ID", "Ensembl_transcript_ID")], file = file.path(dir_out, "genes-transcripts_no_expression", "transcripts_zero_expression.csv"))
  
}

# extract all transcripts of interest 
transcript_extracted <- cts_transcript[which(rownames(cts_transcript) %in% transcript_list_id),]
transcript_extracted_log <- cts_transcript_log[which(rownames(cts_transcript) %in% transcript_list_id),]

if( length(which(rowSums(transcript_extracted) == 0)) ) {

  # get their corresponding gene names
  df_info_subset_zero_exp_transcripts <- df_info[which(transcript_list_id %in% names(which(rowSums(transcript_extracted) == 0))),]
  
  # save a list of transcripts that are not present in the dataset, with their corresponding gene names
  write.csv(df_info_subset_zero_exp_transcripts[,c("Gene", "Ensembl_gene_ID", "Ensembl_transcript_ID")], file = file.path(dir_out, "genes-transcripts_no_expression", "transcripts_zero_expression.csv"))
  
}
  
### HEATMAPS ###

# create data.frames for heatmaps: 
# N - number of rows in the extracted dataset (i.e. number of molecules of interest)
# M - number of samples in the dataset
# A -> list of samples (N times)
# B -> list of all molecules IDs (N times, but each ID repeated M times one after another)
# C -> expression values (N x M)

create_DF_for_heatmap <- function(counts_matrix){
  df <- data.frame(A=rep(colnames(counts_matrix), nrow(counts_matrix)), 
                   B=unlist(lapply(rownames(counts_matrix), function(x) rep(x, length(colnames(counts_matrix))))),
                   C=c(apply(counts_matrix, 1, function(x) as.numeric(c(x)))),
                   row.names = NULL)
  return(df)
}

###---------- need to change names for plots ----------###

# for genes, use gene_names instead of ENSEMBL IDs
matched_ind_gene <- match(rownames(gene_extracted), unique(df_info$Ensembl_gene_ID))
matched_names_gene <- unique(df_info$Gene)[matched_ind_gene]

gene_extracted_with_names <- gene_extracted
rownames(gene_extracted_with_names) <- matched_names_gene

gene_extracted_log_with_names <- gene_extracted_log
rownames(gene_extracted_log_with_names) <- matched_names_gene

# for transcripts, merge gene_names with ENSEMBL IDs 
merged_g_t_names <- paste0(df_info$Gene, "-", df_info$Ensembl_transcript_ID)

matched_ind_transcript <- match(rownames(transcript_extracted), unique(df_info$Ensembl_transcript_ID))

transcript_extracted_with_names <- transcript_extracted
rownames(transcript_extracted_with_names) <- merged_g_t_names[matched_ind_transcript]

transcript_extracted_log_with_names <- transcript_extracted_log
rownames(transcript_extracted_log_with_names) <- merged_g_t_names[matched_ind_transcript]

### create dataframe for heatmaps on gene-level ###
df_heatmap_gene <- create_DF_for_heatmap(gene_extracted_with_names)
df_heatmap_gene_log <- create_DF_for_heatmap(gene_extracted_log_with_names)

### create dataframe for heatmaps on transcript-level ###
df_heatmap_transcript <- create_DF_for_heatmap(transcript_extracted_with_names)
df_heatmap_transcript_log <- create_DF_for_heatmap(transcript_extracted_log_with_names)

#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------- make heatmaps --------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

if ( file.exists(file.path(dir_out, "plots")) ){
  cat("Output directory already exists.")
} else {
  dir.create( file.path(dir_out, "plots", "heatmaps") ) 
}

png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_gene-level_unnormalised.png")), width = 800, height = 1000)

ggplot(df_heatmap_gene, aes(A, B, fill=C)) + geom_tile() + ggtitle("Expression profiles of selected genes - unnormalised counts") + 
  scale_fill_gradient(name = "gene expression", low = "#FFFFFF", high = "#012345") + geom_text(aes(label = round(C,1)), angle = 90, size=3, colour="grey47") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20)) + xlab("Samples")

dev.off()

png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_gene-level_log-transformed.png")), width = 800, height = 1000)

ggplot(df_heatmap_gene_log, aes(A, B, fill=C)) + geom_tile() + scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = "gene expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  xlab("Samples") + ggtitle("Expression profiles of selected genes - log-transformed counts")

dev.off()

png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_transcript-level_unnormalised.png")), width = 1000, height = 1200)

ggplot(df_heatmap_transcript, aes(A, B, fill=C)) + geom_tile() + ggtitle("Expression profiles of selected transcripts - unnormalised counts") + 
  scale_fill_gradient(name = "transcript expression", low = "#FFFFFF", high = "#012345") + geom_text(aes(label = round(C,1)), angle = 90, size=2, colour="grey47") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20)) + xlab("Samples")

dev.off()

png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_transcript-level_log-transformed.png")), width = 1000, height = 1200)

ggplot(df_heatmap_transcript_log, aes(A, B, fill=C)) + geom_tile() + scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = "transcript expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  xlab("Samples") + ggtitle("Expression profiles of selected transcripts - log-transformed counts")

dev.off()

##################################################################################################
########### heatmaps with samples clustered in two groups; e.g. “normal” vs “inflamed” ###########
##################################################################################################

# add additional column "depth" to df_heatmap_gene and df_heatmap_gene_log  
GCA_status <- df_histo$GCA_present
idx_0 <- which(GCA_status == 0) ; idx_1 <- which(GCA_status == 1)
GCA_status[idx_0] <- "normal"   ; GCA_status[idx_1] <- "inflamed"

df_heatmap_gene_mod <- cbind(df_heatmap_gene, rep(GCA_status, nrow(gene_extracted_with_names)))
colnames(df_heatmap_gene_mod)[4] <- "D"

df_heatmap_gene_log_mod <- cbind(df_heatmap_gene_log, rep(GCA_status, nrow(gene_extracted_with_names)))
colnames(df_heatmap_gene_log_mod)[4] <- "D"


png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_gene-level_unnormalised_grouped.png")), width = 800, height = 1000)

ggplot(df_heatmap_gene_mod, aes(A, B, fill=C)) + geom_tile() + facet_grid(~ D, scales = "free_x", space = "free_x") +
  ggtitle("Expression profiles of selected genes - unnormalised counts - grouped") + 
  scale_fill_gradient(name = "gene expression", low = "#FFFFFF", high = "#012345") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20)) + xlab("Samples")

dev.off()

png(file.path(dir_out, "plots", paste0("heatmap_expression_profiles_gene-level_log-transformed_grouped.png")), width = 800, height = 1000)

ggplot(df_heatmap_gene_log_mod, aes(A, B, fill=C)) + geom_tile() + facet_grid(~ D, scales = "free_x", space = "free_x") +
  ggtitle("Expression profiles of selected genes - log-transformed counts - grouped") +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, name = "gene expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + xlab("Samples") 

dev.off()

#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------- make boxplots --------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

#### ONLY for gene-level ####

# NOTE: df_spreadsheet has 39 elements and gene_extracted_with_names has 41
# remove 2 samples from gene_extracted_with_names
# diference <- setdiff(colnames(gene_extracted_with_names), df_histo$Database_number)
# gene_extracted_with_names_final <- gene_extracted_with_names[,-which(colnames(gene_extracted_with_names) %in% diference)]
# gene_extracted_log_with_names_final <- gene_extracted_log_with_names[,-which(colnames(gene_extracted_log_with_names) %in% diference)]

#------------------------------------ PLOT 1 ------------------------------------#
# create boxplots to compare mean gene expression for SSTR2 between "no" vs "yes"
#--------------------------------------------------------------------------------#

# define data type: "raw" | "log"
data_type <- "raw"

# specify feature for running: 
# "GCA_present" | "Giant_cells" | "Infiltrate_around_vasa_vasorum" | "Media_destruction"
# "Adventitia_pattern" | "Media_pattern" | "Intima_pattern" | "Occlusion_grade"
feat <- "GCA_present"

#transcript_extracted_with_names

# create dataframe with all genes and the feature of interest
df_for_boxplot <- data.frame(cbind(t(gene_extracted_with_names), df_histo[,feat]))
colnames(df_for_boxplot)[ncol(df_for_boxplot)] <- "feature"

df_for_boxplot_log <- data.frame(cbind(t(gene_extracted_log_with_names), df_histo[,feat]))
colnames(df_for_boxplot_log)[ncol(df_for_boxplot_log)] <- "feature"
  
makeBoxplotPerGene <- function( df , column_number ){   # column_number => correspond to the feature
  
  # create a subset
  df_sub <- df[,c(column_number, ncol(df))]
  colnames(df_sub)[1] <- "expression"
  
  # set 'feature' as factor
  df_sub$feature <- as.factor(df_sub$feature)

  # make a boxplot
  pl <- ggplot(df_sub, aes(x=feature, y=expression, color=feature)) + 
          geom_boxplot(outlier.shape=NA) +   theme(text = element_text(size=14)) + 
          geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
          labs(title=paste0("gene ", colnames(df_for_boxplot)[column_number]), x=feat, y="gene expression") + 
          theme(plot.title = element_text(hjust = 0.5))

  return(pl)
}

if ( file.exists(file.path(dir_out, "plots", "boxplots")) ){
  cat("Output directory already exists.")
} else {
  dir.create(file.path(dir_out, "plots", "boxplots"))
}

### ---- overate over all genes ---- ###
for (i in 1:5) {
  print(i)
  png(file.path(dir_out, "plots", "boxplots", paste0("boxplot_gene-level_", colnames(df_for_boxplot)[i], "_unnormalised.png")))
  print(makeBoxplotPerGene(df_for_boxplot, i)); dev.off()
  
  png(file.path(dir_out, "plots", "boxplots", paste0("boxplot_gene-level_", colnames(df_for_boxplot)[i], "_log-transformed.png")))
  print(makeBoxplotPerGene(df_for_boxplot_log, i)); dev.off()
}

################################################################################################
#library("gridExtra")
#empty_pl <- ggplot() + theme_void()

#grid.arrange(makeBoxplotPerGene(df_for_boxplot, 15), makeBoxplotPerGene(df_for_boxplot, 12), 
#             makeBoxplotPerGene(df_for_boxplot, 10), makeBoxplotPerGene(df_for_boxplot, 13), 
#             makeBoxplotPerGene(df_for_boxplot, 8), makeBoxplotPerGene(df_for_boxplot, 1),
#             makeBoxplotPerGene(df_for_boxplot, 9), makeBoxplotPerGene(df_for_boxplot, 3), 
#             makeBoxplotPerGene(df_for_boxplot, 4), ncol = 3, nrow = 3)

#grid.arrange(makeBoxplotPerGene(df_for_boxplot, 16), makeBoxplotPerGene(df_for_boxplot, 7), 
#             makeBoxplotPerGene(df_for_boxplot, 14), makeBoxplotPerGene(df_for_boxplot, 2),  
#             makeBoxplotPerGene(df_for_boxplot, 6), empty_pl, makeBoxplotPerGene(df_for_boxplot, 11), 
#             makeBoxplotPerGene(df_for_boxplot, 5),  ncol = 3, nrow = 3)
################################################################################################


# TO BE CONTINUES


### ---- feat_1 = "GCA_present" ---- ###
df_sstr2_boxplot_gene_feat_1 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_1])

df_sstr2_boxplot_gene_feat_1$feature <- as.factor(df_sstr2_boxplot_gene_feat_1$feature)
run_feat_1 <- colnames(df_spreadsheet)[feat_1]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_1, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_1, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_1), x =run_feat_1, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_2 = "Giant_cells" ---- ###
df_sstr2_boxplot_gene_feat_2 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_2])

df_sstr2_boxplot_gene_feat_2$feature <- as.factor(df_sstr2_boxplot_gene_feat_2$feature)
run_feat_2 <- colnames(df_spreadsheet)[feat_2]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_2, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_2, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_2), x =run_feat_2, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_3 = "Infiltrate_around_vasa_vasorum" ---- ###
df_sstr2_boxplot_gene_feat_3 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_3])

df_sstr2_boxplot_gene_feat_3$feature <- as.factor(df_sstr2_boxplot_gene_feat_3$feature)
run_feat_3 <- colnames(df_spreadsheet)[feat_3]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_3, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_3, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_3), x =run_feat_3, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_4 = "Media_destruction" ---- ###
df_sstr2_boxplot_gene_feat_4 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_4])

df_sstr2_boxplot_gene_feat_4$feature <- as.factor(df_sstr2_boxplot_gene_feat_4$feature)
run_feat_4 <- colnames(df_spreadsheet)[feat_4]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_4, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_4, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_4), x =run_feat_4, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_5 = "Adventitia_pattern" ---- ###
df_sstr2_boxplot_gene_feat_5 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_5])

df_sstr2_boxplot_gene_feat_5$feature <- as.factor(df_sstr2_boxplot_gene_feat_5$feature)
run_feat_5 <- colnames(df_spreadsheet)[feat_5]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_5, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_5, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_5), x =run_feat_5, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_6 = "Media_pattern" ---- ###
df_sstr2_boxplot_gene_feat_6 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_6])

df_sstr2_boxplot_gene_feat_6$feature <- as.factor(df_sstr2_boxplot_gene_feat_6$feature)
run_feat_6 <- colnames(df_spreadsheet)[feat_6]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_6, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_6, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_6), x =run_feat_6, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_7 = "Intima_pattern" ---- ####
df_sstr2_boxplot_gene_feat_7 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_7])

df_sstr2_boxplot_gene_feat_7$feature <- as.factor(df_sstr2_boxplot_gene_feat_7$feature)
run_feat_7 <- colnames(df_spreadsheet)[feat_7]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_7, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_7, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_7), x =run_feat_7, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### ---- feat_8 = "Occlusion_grade" ---- ####
df_sstr2_boxplot_gene_feat_8 <- data.frame(expression=as.numeric(gene_extracted_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_8])

df_sstr2_boxplot_gene_feat_8$feature <- as.factor(df_sstr2_boxplot_gene_feat_8$feature)
run_feat_8 <- colnames(df_spreadsheet)[feat_8]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_8, "_", data_type, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_8, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_8), x =run_feat_8, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
 

############# log-transformed data #############
data_type_log <- "log"

#feat_1 = "GCA_present" 
df_sstr2_boxplot_gene_feat_1 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_1])

df_sstr2_boxplot_gene_feat_1$feature <- as.factor(df_sstr2_boxplot_gene_feat_1$feature)
run_feat_1 <- colnames(df_spreadsheet)[feat_1]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_1, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_1, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_1), x =run_feat_1, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_2 = "Giant_cells"
df_sstr2_boxplot_gene_feat_2 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_2])

df_sstr2_boxplot_gene_feat_2$feature <- as.factor(df_sstr2_boxplot_gene_feat_2$feature)
run_feat_2 <- colnames(df_spreadsheet)[feat_2]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_2, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_2, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_2), x =run_feat_2, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_3 = "Infiltrate_around_vasa_vasorum"
df_sstr2_boxplot_gene_feat_3 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_3])

df_sstr2_boxplot_gene_feat_3$feature <- as.factor(df_sstr2_boxplot_gene_feat_3$feature)
run_feat_3 <- colnames(df_spreadsheet)[feat_3]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_3, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_3, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_3), x =run_feat_3, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_4 = "Media_destruction"
df_sstr2_boxplot_gene_feat_4 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_4])

df_sstr2_boxplot_gene_feat_4$feature <- as.factor(df_sstr2_boxplot_gene_feat_4$feature)
run_feat_4 <- colnames(df_spreadsheet)[feat_4]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_4, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_4, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_4), x =run_feat_4, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_5 = "Adventitia_pattern"
df_sstr2_boxplot_gene_feat_5 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_5])

df_sstr2_boxplot_gene_feat_5$feature <- as.factor(df_sstr2_boxplot_gene_feat_5$feature)
run_feat_5 <- colnames(df_spreadsheet)[feat_5]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_5, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_5, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_5), x =run_feat_5, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_6 = "Media_pattern"
df_sstr2_boxplot_gene_feat_6 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_6])

df_sstr2_boxplot_gene_feat_6$feature <- as.factor(df_sstr2_boxplot_gene_feat_6$feature)
run_feat_6 <- colnames(df_spreadsheet)[feat_6]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_6, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_6, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_6), x =run_feat_6, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_7 = "Intima_pattern"
df_sstr2_boxplot_gene_feat_7 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_7])

df_sstr2_boxplot_gene_feat_7$feature <- as.factor(df_sstr2_boxplot_gene_feat_7$feature)
run_feat_7 <- colnames(df_spreadsheet)[feat_7]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_7, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_7, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_7), x =run_feat_7, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#feat_8 = "Occlusion_grade"
df_sstr2_boxplot_gene_feat_8 <- data.frame(expression=as.numeric(gene_extracted_log_with_names_final["SSTR2",]),
                                           feature=df_spreadsheet[,feat_8])

df_sstr2_boxplot_gene_feat_8$feature <- as.factor(df_sstr2_boxplot_gene_feat_8$feature)
run_feat_8 <- colnames(df_spreadsheet)[feat_8]

png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_", run_feat_8, "_", data_type_log, ".png")))
ggplot(df_sstr2_boxplot_gene_feat_8, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2 - gene level - ", run_feat_8), x =run_feat_8, y="gene expression") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

##########################################################################################
###################################### for old data ######################################
##########################################################################################

#------------------- RAW
# create a data.frame for ploting
#df_sstr2_boxplot_gene_old_raw <- data.frame(expression=as.numeric(cts_gene_subset_old_raw[2,-1]),
#                                            feature=df_histological[,run_feat])
#df_sstr2_boxplot_gene_old_raw$feature <- as.factor(df_sstr2_boxplot_gene_old_raw$feature)

#png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_OLD_raw_", run_feat, "_", data_type, ".png")))

#ggplot(df_sstr2_boxplot_gene_old_raw, aes(x=feature, y=expression, color=feature)) + 
#  geom_boxplot(outlier.shape=NA) + 
#  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
#  theme(text = element_text(size=14)) + 
#  labs(title=paste0("SSTR2 - gene level - 'old' data - raw ", run_feat), x =run_feat, y="gene expression") + 
#  theme(plot.title = element_text(hjust = 0.5))

#dev.off()

#------------------- RLOG
#df_sstr2_boxplot_gene_old_rlog <- data.frame(expression=as.numeric(cts_gene_subset_old_rlog[2,-1]),
#                                             feature=df_histological[,run_feat])
#df_sstr2_boxplot_gene_old_rlog$feature <- as.factor(df_sstr2_boxplot_gene_old_rlog$feature)

#png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_OLD_rlog_", run_feat, "_", data_type, ".png")))

#ggplot(df_sstr2_boxplot_gene_old_rlog, aes(x=feature, y=expression, color=feature)) + 
#  geom_boxplot(outlier.shape=NA) + 
#  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
#  theme(text = element_text(size=14)) + 
#  labs(title=paste0("SSTR2 - gene level - 'old' data - rlog ", run_feat), x =run_feat, y="gene expression") + 
#  theme(plot.title = element_text(hjust = 0.5))

#dev.off()

#------------------- VST

#df_sstr2_boxplot_gene_old_vst <- data.frame(expression=as.numeric(cts_gene_subset_old_vst[2,-1]),
#                                            feature=df_histological[,run_feat])
#df_sstr2_boxplot_gene_old_vst$feature <- as.factor(df_sstr2_boxplot_gene_old_vst$feature)

#png(file.path(dir_out, paste0("boxplot_SSTR2_gene-level_OLD_vst_", run_feat, "_", data_type, ".png")))

#ggplot(df_sstr2_boxplot_gene_old_vst, aes(x=feature, y=expression, color=feature)) + 
#  geom_boxplot(outlier.shape=NA) + 
#  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
#  theme(text = element_text(size=14)) + 
#  labs(title=paste0("SSTR2 - gene level - 'old' data - vst ", run_feat), x =run_feat, y="gene expression") + 
#  theme(plot.title = element_text(hjust = 0.5))

#dev.off()

##########################################################################################
##########################################################################################
##########################################################################################

#-------------------------------------- PLOT 2 --------------------------------------#
# create heatmap for one continous feature (e.g. number of days on steroids)
#------------------------------------------------------------------------------------#

# create a data.frame:
# X -> list of samples (5 times)
# Y -> list of all receptors IDs (5 times, but each ID repeated 40 times one after another)
# Z -> feature name (40 x 5)

#png(file.path(dir_out, paste0("heatmap_expression_profiles_gene-level_unnormalised.png")), width = 800, height = 1000)

#ggplot(df_heatmap_gene, aes(A, B, fill=C)) + geom_tile() + ggtitle("Expression profiles of selected genes - unnormalised counts") + 
#  scale_fill_gradient(name = "gene expression", low = "#FFFFFF", high = "#012345") + geom_text(aes(label = round(C,1)), angle = 90, size=3, colour="grey47") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
#        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 20)) + xlab("Samples")

#dev.off()

# concatenate nb of days on steroids with IDs
which(nchar(df_clinical_steroids) !=1)
# 2 => 11 
# 5 => 16
# 6 => 16

df_clinical_steroids_final <- paste0("0", df_clinical_steroids)
df_clinical_steroids_final[c(2, 5, 6)] <- c("11", "16", "16")

df_heatmap_gene_steroids <- data.frame(X=paste0(df_clinical_steroids_final, "_days_", colnames(gene_extracted_with_names)), 
                                       Y=c(rep("SSTR2", 41)),
                                       Z=as.numeric(c(gene_extracted_with_names["SSTR2",])))

df_heatmap_gene_steroids_log <- data.frame(X=paste0(df_clinical_steroids_final, "_days_", colnames(gene_extracted_with_names)), 
                                           Y=c(rep("SSTR2", 41)),
                                           Z=as.numeric(c(gene_extracted_log_with_names["SSTR2",])))

run_feat_9 <- "steroids"

png(file.path(dir_out, paste0("heatmap_SSTR2_gene-level_", run_feat_9, "_", data_type, ".png")))

ggplot(df_heatmap_gene_steroids, aes(X, Y, fill=Z)) + 
  geom_tile() + facet_grid(scales = "free_x", space = "free_x") + 
  scale_fill_gradient(name = "Gene expression", low = "#FFFFFF", high = "#012345") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("Samples") + ggtitle("Number of days on steroids - raw data")

dev.off()

png(file.path(dir_out, paste0("heatmap_SSTR2_gene-level_", run_feat_9, "_", data_type_log, ".png")))

ggplot(df_heatmap_gene_steroids_log, aes(X, Y, fill=Z)) + 
  geom_tile() + facet_grid(scales = "free_x", space = "free_x") + 
  scale_fill_gradient(name = "Gene expression", low = "#FFFFFF", high = "#012345") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("Samples") + ggtitle("Number of days on steroids - log-transformed data")

dev.off()






(3.1) create heatmaps with samples ordered in a way that they show two groups: “normal” vs “inflamed”
(3.2) create boxplots to compare the average expression levels between the 2 groups





# TODO: to be organised from here onwards


#-------------------------------------- PLOT 2 --------------------------------------#
# create heatmap to compare gene expression level for all SSTRs between "no" vs "yes"
#------------------------------------------------------------------------------------#

# create a data.frame:
# X -> list of samples (5 times)
# Y -> list of all receptors IDs (5 times, but each ID repeated 40 times one after another)
# Z -> feature name (40 x 5)

data_gene <- data.frame(X=rep(colnames(cts_gene_subset[2,-1]), 5), 
                        Y=c(rep("SSTR1", 40), rep("SSTR2", 40), rep("SSTR3", 40), rep("SSTR4", 40), rep("SSTR5", 40)),
                        Depth=rep(df_histological[,run_feat], 5), 
                        Z=as.numeric(c(cts_gene_subset[1,-1], cts_gene_subset[2,-1], cts_gene_subset[3,-1], cts_gene_subset[4,-1], cts_gene_subset[5,-1])))

# create a heatmap without separation in "no" vs. "yes"
#if(output_save==TRUE){ 

png(file.path(dir_out, run_feat, paste0("heatmap_SSTRs_gene-level_", run_feat, "_", data_type, ".png")))
ggplot(data_gene, aes(X, Y, fill=Z)) + 
  geom_tile() + facet_grid(~ Depth, scales = "free_x", space = "free_x") + 
  scale_fill_gradient(name = "Gene expression", low = "#FFFFFF", high = "#012345") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  xlab("Samples") + ggtitle(run_feat)
dev.off()
#}

# same as above but only for SSTR2 and ordered as feature="no", then feature="yes"
#data_gene_bis <- data.frame(X=colnames(cts_gene_subset[2,-1]),
#                            Y=c(rep("SSTR2", 40)),
#                            Depth=df_histological[,run_feat],
#                            Z=as.numeric(c(cts_gene_subset[2,-1])))

#ggplot(data_gene_bis, aes(X, Y, fill=Z)) + 
#  geom_tile() + facet_grid(~ Depth, scales = "free_x", space = "free_x") +
#  scale_fill_gradient2(name = "Gene expression", low = "red", high = "red") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#        axis.title.y = element_blank()) + xlab("Samples")

##########################################################
#################### transcript-level ####################
##########################################################

# NOTE: gr_no & gr_yes are the same regardless data type

# "ENST00000267377.2" | SSTR1-201
df_1_tr <- which(startsWith(cts_transcript$X, "ENST00000267377")); cts_transcript$X[df_1_tr]

# "ENST00000357585.3" | SSTR2-201
df_2.1_tr <- which(startsWith(cts_transcript$X, "ENST00000357585")); cts_transcript$X[df_2.1_tr]

# "ENST00000579323.5" | SSTR2-202
df_2.2_tr <- which(startsWith(cts_transcript$X, "ENST00000579323")); cts_transcript$X[df_2.2_tr]

# "ENST00000610913.1" | SSTR3-201
df_3.1_tr <- which(startsWith(cts_transcript$X, "ENST00000610913")); cts_transcript$X[df_3.1_tr] 

# "ENST00000617123.1" | SSTR3-202
df_3.2_tr <- which(startsWith(cts_transcript$X, "ENST00000617123")); cts_transcript$X[df_3.2_tr] 

# "ENST00000255008.4"  | SSTR4-201
df_4_tr <- which(startsWith(cts_transcript$X, "ENST00000255008")); cts_transcript$X[df_4_tr] 

# "ENST00000293897.5"  | SSTR5-201
df_5_tr <- which(startsWith(cts_transcript$X, "ENST00000293897")); cts_transcript$X[df_5_tr] 

# create subsets including only SSTRs transcripts
cts_transcript_subset <- cts_transcript[c(df_1_tr, df_2.1_tr, df_2.2_tr, df_3.1_tr, df_3.2_tr, df_4_tr, df_5_tr),]

# split in half feature=0 | feature=1
cts_transcript_subset_no <- cts_transcript_subset[,gr_no]
cts_transcript_subset_yes <- cts_transcript_subset[,gr_yes]

# get average transcript-level expression of each SSTR and compare "no" vs. "yes
transcript_mean_no <- apply(cts_transcript_subset_no, 1, mean)
transcript_mean_yes <- apply(cts_transcript_subset_yes, 1, mean)

names(transcript_mean_no) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201")
names(transcript_mean_yes) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201")

#--------------------------------------- PLOT 3 ---------------------------------------#
# create boxplots to compare mean transcript expression for SSTR2 between "no" vs "yes"
#--------------------------------------------------------------------------------------#

# create a data.frame for ploting
df_sstr2_boxplot_transcript <- data.frame(expression=as.numeric(cts_transcript_subset[2,-1]),
                                          feature=df_histological[,run_feat])

df_sstr2_boxplot_transcript$feature <- as.factor(df_sstr2_boxplot_transcript$feature)

#if(output_save==TRUE){ 
png(file.path(dir_out, run_feat, paste0("boxplot_SSTR2_transcript-level_", run_feat, "_", data_type, ".png")))

ggplot(df_sstr2_boxplot_transcript, aes(x=feature, y=expression, color=feature)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size=3) + 
  theme(text = element_text(size=14)) + 
  labs(title=paste0("SSTR2-201 - transcript level - ", run_feat), x=run_feat, y="transcript expression") + 
  theme(plot.title = element_text(hjust = 0.5))

dev.off()
#}

#----------------------------------------- PLOT 4 -----------------------------------------#
# create heatmap to compare transcript expression level for all SSTRs between "no" vs "yes"
#------------------------------------------------------------------------------------------#

# create a data.frame:
# X -> list of samples (7 times)
# Y -> list of all receptors IDs (5 times, but each ID repeated 40 times one after another)
# Z -> feature name (40 x 5)

data_transcript <- data.frame(A=rep(colnames(cts_transcript_subset[2,-1]), 7), 
                              B=c(rep("SSTR1-201", 40), rep("SSTR2-201", 40), rep("SSTR2-202", 40), rep("SSTR3-201", 40), rep("SSTR3-202", 40), rep("SSTR4-201", 40), rep("SSTR5-201", 40)),
                              Depth=rep(df_histological[,run_feat], 7), 
                              C=as.numeric(c(cts_transcript_subset[1,-1], cts_transcript_subset[2,-1], cts_transcript_subset[3,-1], cts_transcript_subset[4,-1], cts_transcript_subset[5,-1], cts_transcript_subset[6,-1], cts_transcript_subset[7,-1])))

# create a heatmap without separation in "no" vs. "yes"
#if(output_save==TRUE){ 
png(file.path(dir_out, run_feat, paste0("heatmap_SSTRs_transcript-level_", run_feat, "_", data_type, ".png")))

ggplot(data_transcript, aes(A, B, fill=C)) + 
  geom_tile() + facet_grid(~ Depth, scales = "free_x", space = "free_x") + 
  scale_fill_gradient(name = "Transcript expression", low = "#FFFFFF", high = "#012345") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("Samples") + ggtitle(run_feat)

dev.off()

#}


########################################################################################
################################### ADDITIONAL  CODE ###################################
########################################################################################

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
# load normalised data (NOT IN USE ANYMORE)
#dir_gene_norm <- 
#dir_transcript_norm <- 
#cts_gene_norm <- read.csv(dir_gene_norm, stringsAsFactors = FALSE, row.names = 1) 
#cts_transcript_norm <- read.csv(dir_transcript_norm, stringsAsFactors = FALSE, row.names = 1) 

# load "old" data (i.e. obtained with featureCounts()) (NOT IN USE ANYMORE)
#dir_cts_old_raw <- "/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_raw_no_outliers.csv"
#dir_cts_old_rlog <- "/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_rlog_no_outliers.csv"
#dir_cts_old_vst <- "/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_vst_no_outliers.csv"

#cts_old_raw <- read.csv(dir_cts_old_raw, stringsAsFactors = FALSE) 
#cts_old_rlog <- read.csv(dir_cts_old_rlog, stringsAsFactors = FALSE) 
#cts_old_vst <- read.csv(dir_cts_old_vst, stringsAsFactors = FALSE) 
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#

