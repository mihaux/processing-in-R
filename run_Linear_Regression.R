# script to perform Mann-Whitney U test (also known as Wilcoxon rank sum test or Mann-Whitney test)
# 2 methods: wilcox.test() [built-in R function] and wilcoxTest() [GSALightning package]

# very clear source: http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/RProgramming/NonParametricStats.html

# The unpaired two-samples Wilcoxon test (also known as Wilcoxon rank sum test or Mann-Whitney test) 
# => is a non-parametric alternative to the unpaired two-samples t-test, 
# => can be used to compare two independent groups of samples. 
# => is used when your data are not normally distributed.

# NOTE: if there are confounders, then you would need to rather use something like logistic regression, not Mann-Whitney test

# WHAT TO TEST: statistical test between 2 groups: visual loss vs no visual loss,
# to see if patients with visual loss are more likely to have GCA or not (if not, then there's no association)

# To perform two-samples Wilcoxon test comparing the means of two independent samples (x & y), 
# we perform a t-test on each gene (i.e. each row) by running the function wilcox.test() for each row, in a two-sample setting between group_1 and group_2.

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GSALightning", quietly = TRUE)) BiocManager::install("GSALightning"); library(GSALightning)
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt"); library(biomaRt)


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

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to _x.csv file with count data (or just the directory), 
       \n(2 - annotation) path to _x.csv annotation file and
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

# Raw_DESeq_dataset_all.Rda             | Raw_DESeq_dataset_chr1_22.Rda               | Raw_DESeq_dataset_chrXY.Rda
# Normalised_DESeq_rlog_dataset_all.Rda | Normalised_DESeq_rlog_dataset_chr1_22.Rda   | Normalised_DESeq_rlog_dataset_chrXY.Rda
# Normalised_DESeq_vst_dataset_all.Rda  | Normalised_DESeq_vst_dataset_chr1_22.Rda    | Normalised_DESeq_vst_dataset_chrXY.Rda

args <- c(paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/all_chr/"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/all_chr/mann_whitney/"))

# Example of usage: 
# Rscript test_Mann_Whitney.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load normalised counts from DESeq
load(paste0(args[1], "Raw_DESeq_dataset_all.Rda"))              # dds_all
load(paste0(args[1], "Normalised_DESeq_vst_dataset_all.Rda"))   # vst_all
load(paste0(args[1], "Normalised_DESeq_rlog_dataset_all.Rda"))  # rlog_all

# change data format to matrix and integer
#dat <- as.matrix(df)   
#storage.mode(dat) <- "integer"             # class: matrix | type: integer

dat_raw <- as.matrix(assay(dds_all))
storage.mode(dat_raw) <- "integer"             # class: matrix | type: integer

dat_vst <- as.matrix(assay(vst_all))
storage.mode(dat_vst) <- "integer"             # class: matrix | type: integer

dat_rlog <- as.matrix(assay(rlog_all))
storage.mode(dat_rlog) <- "integer"             # class: matrix | type: integer

# load annotation (clinical) data which will be used for coldata
anno_1 <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames for clinical data
rownames(anno_1) <- paste0("ID_", rownames(anno_1))

# load slide scores data as well
anno_2 <- read.csv(paste0(main_dir, "/data/metadata/slide_scores/slide_scores_v6.csv"), row.names = 1)

# NOTE: sample "ID_14058" is missing in 'slide_sc', needs to be removed from 

# add "ID_" to all rownames for slide scores data
rownames(anno_2) <- paste0("ID_", rownames(anno_2))

# features to be used:

# => clinical_feature (anno_1)
# gender..1.male..2.female.                 # => anno_1$gender                                    (1 vs. 2)
# visual.loss.at.BL..0.no..1.yes.           # => anno_1$visual.loss.at.BL..0.no..1.yes.           (0 vs. 1)
# jaw.claudication.at.BL...0.no..1.yes.     # => anno_1$jaw.claudication.at.BL...0.no..1.yes.     (0 vs. 1)
# ischaemic.features.at.BL...0.no..1.yes.   # => anno_1$ischaemic.features.at.BL...0.no..1.yes.   (0 vs. 1)

# => slide_scores (anno_2)
# GCA.present.                              # => anno_2$GCA.present.                              (0 vs. 1)
# Giant.cells.                              # => anno_2$Giant.cells.                              (0 vs. 1)
# Media.destruction.                        # => anno_2$Media.destruction.                        (0 vs. 1)
# Occlusion.grade.                          # => anno_2$Occlusion.grade. [Occlusion.grade.new]    (0 vs. 1)
# Neoangiogenesis.                          # => anno_2$Neoangiogenesis.                          (0 vs. 1)
# Intima.pattern.                           # => anno_2$Intima.pattern. [Intima.pattern.new]      (0 vs. 1)
# Infiltrate.around.vasa.vasorum.           # => anno_2$Infiltrate.around.vasa.vasorum.           (0 vs. 1)

# modify Occlusion.grade. and Intima.pattern. to make it 2 comparison groups
Occlusion.grade.new <- c(1:40)
Occlusion.grade.new[which(anno_2$Occlusion.grade. < 3)] <- 0      # less than 3           => replace with 0
Occlusion.grade.new[which(anno_2$Occlusion.grade. >= 3)] <- 1     # equal or more than 3  => replace with 1

Intima.pattern.new <- c(1:40)
Intima.pattern.new[which(anno_2$Intima.pattern. <= 1)] <- 0       # 0 and 1                => replace with 0
Intima.pattern.new[which(anno_2$Intima.pattern. > 1)] <- 1        # 2 and 3                => replace with 1

# define groups for comparison

# => from 'anno_1'
#group = anno_1$gender..1.male..2.female.                   # (1 vs. 2)
#group = anno_1$visual.loss.at.BL..0.no..1.yes.             # (0 vs. 1)
#group = anno_1$jaw.claudication.at.BL...0.no..1.yes.       # (0 vs. 1)
#group = anno_1$ischaemic.features.at.BL...0.no..1.yes.    # (0 vs. 1)

# => from 'anno_2'
group = anno_2$GCA.present.                               # (0 vs. 1)
#group = anno_2$Giant.cells.                               # (0 vs. 1)
# group = anno_2$Media.destruction.                         # (0 vs. 1)
# group = Occlusion.grade.new                               # (0 vs. 1)
# group = anno_2$Neoangiogenesis.                           # (0 vs. 1)
# group = Intima.pattern.new                                # (0 vs. 1)
# group = anno_2$Infiltrate.around.vasa.vasorum.            # (0 vs. 1)

# RUN STATISTICAL TESTING

running <- 'GCA_present'

# HYPOTHISES: we want to know if the mean of group 1 differs from the mean of group 2.

# NOTE: no metter the comparison group_1 vs group_2 or group_2 vs group_1, 
# the results are the same, just the 'statistic' component differs (but all is proportional)

# 'alternative' argument => "two,sided" (default)
# if you want to test whether the mean of group_1 is less than the mean of group_2        => alternative = "less"
# if you want to test whether the mean of group_1 is greater than the mean of group_2     => alternative = "greater"

# INTERPRETATION:
# if the p-value of the test is less than the significance level alpha = 0.05. 
# we can conclude that the mean of group 1 is significantly different from the mean of group 2 with a p-value of [p-value =].

# if running for slide scores (anno_2) => remove sample "ID_14058" from all three datasets
#which(colnames(dat_raw) == "ID_14058")
dat_raw <- dat_raw[,-which(colnames(dat_raw) == "ID_14058")]
dat_vst <- dat_vst[,-which(colnames(dat_vst) == "ID_14058")]
dat_rlog <- dat_rlog[,-which(colnames(dat_rlog) == "ID_14058")]

# initiate list for results
res_raw <- list()
res_vst <- list()
res_rlog <- list()

# iteration over all genes (as we perform statistical testing for each gene)
for(j in 1:nrow(dat_raw)){
  
  # get data frame for one gene at a time
  temp_raw = dat_raw[j, ]
  temp_vst = dat_vst[j, ]
  temp_rlog = dat_rlog[j, ]
  
  # perform statistical testing
  res_raw[[j]] = wilcox.test(temp_raw[group==0], 
                         temp_raw[group==1], 
                         alternative = "two.sided",
                         exact = FALSE)               # to suppress the warning message saying that “cannot compute exact p-value with tie”
                                                      # it comes from the assumption of a Wilcoxon test that the responses are continuous. 
  
  res_vst[[j]] = wilcox.test(temp_vst[group==0], temp_vst[group==1], alternative = "two.sided", exact = FALSE)              
  res_rlog[[j]] = wilcox.test(temp_rlog[group==0], temp_rlog[group==1], alternative = "two.sided", exact = FALSE)              
  
}
# might need to add as.numeric() if does not work

#---------------------------------------------------------------------------------------------------------
# running for 'gender' => group==1 AND group==2 
if(FALSE){
for(j in 1:nrow(dat_raw)){

  # get data frame for one gene at a time
  print(j)
  temp_raw = dat_raw[j, ]
  temp_vst = dat_vst[j, ]
  temp_rlog = dat_rlog[j, ]
  
  # perform statistical testing 
  res_raw[[j]] = wilcox.test(temp_raw[group==1], 
                             temp_raw[group==2], 
                             alternative = "two.sided",
                             exact = FALSE)               # to suppress the warning message saying that “cannot compute exact p-value with tie”
  # it comes from the assumption of a Wilcoxon test that the responses are continuous. 
  
  res_vst[[j]] = wilcox.test(temp_vst[group==1], temp_vst[group==2], alternative = "two.sided", exact = FALSE)              
  res_rlog[[j]] = wilcox.test(temp_rlog[group==1], temp_rlog[group==2], alternative = "two.sided", exact = FALSE)              
  
  }
}
#---------------------------------------------------------------------------------------------------------

# the object 'res' contains the result of the test | length = 52239
# res[[1]]
# Wilcoxon rank sum test with continuity correction
# data:  temp[group == 2] and temp[group == 1]
# W = 254.5, p-value = 0.07999
# alternative hypothesis: true location shift is not equal to 0

# in the OUTPUT: (assume there are 2 groups: A and B)
# => W is a rank sum subtracted by a constant. 
# => if group A is the reference level in the factor variable group, the rank sum is computed for the data in group A
# => the value of W is the rank sum for group A subtracted by nA(nA+1)/2. 
# => the value of W is not important, as we perform the Wilcoxon Mann-Whitney test to determine if there is significant difference between the two groups. 
# => we are primarily interested in the p-value
# => by default wilcox.test() calculates an exact p-value if the samples contain less than 50 finite values and there are no ties. 
# => Otherwise, a normal approximation is used. R uses normal approximation to calculate the p-value.
  
# retrieve t-statistic for each gene,
all_stat_raw = unlist(lapply(res_raw, function(x) x$statistic))
all_stat_vst = unlist(lapply(res_vst, function(x) x$statistic))
all_stat_rlog = unlist(lapply(res_rlog, function(x) x$statistic))

# retrieve its corresponding p-value
all_pval_raw = unlist(lapply(res_raw, function(x) x$p.value))
all_pval_vst = unlist(lapply(res_vst, function(x) x$p.value))
all_pval_rlog = unlist(lapply(res_rlog, function(x) x$p.value))

# Multiplicity adjustment
# After obtaining the p-values, you can make multiplicity correction to the pvalue. 
all_adj.pval_raw <- p.adjust(all_pval_raw, "fdr")    # Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
all_adj.pval_vst <- p.adjust(all_pval_vst, "fdr") 
all_adj.pval_rlog <- p.adjust(all_pval_rlog, "fdr") 

# result table
# create a table as R data frame
result.table2_raw = data.frame(ID=rownames(dat_raw), statistic=all_stat_raw, pvalue=all_pval_raw, fdr.pvalue=all_adj.pval_raw)
result.table2_vst = data.frame(ID=rownames(dat_vst), statistic=all_stat_vst, pvalue=all_pval_vst, fdr.pvalue=all_adj.pval_vst)
result.table2_rlog = data.frame(ID=rownames(dat_rlog), statistic=all_stat_rlog, pvalue=all_pval_rlog, fdr.pvalue=all_adj.pval_rlog)

# sort the table based on the order of the (adjusted) p-values
result.table2.sorted_raw = result.table2_raw[order(all_adj.pval_raw),]
result.table2.sorted_vst = result.table2_vst[order(all_adj.pval_vst),]
result.table2.sorted_rlog = result.table2_rlog[order(all_adj.pval_rlog),]

# the 'statistic' column is not necessary
result.table2.sorted_final_raw <- result.table2.sorted_raw[,-2]
result.table2.sorted_final_vst <- result.table2.sorted_vst[,-2]
result.table2.sorted_final_rlog <- result.table2.sorted_rlog[,-2]

# and then show the top 10 (most) significant genes.
#result.table2.sorted[1:10,]       # listing the top 10 genes

# select only significant genes
significant.cutoff_raw <- length(which(result.table2.sorted_raw$fdr.pvalue < 0.05))
significant.cutoff_vst <- length(which(result.table2.sorted_vst$fdr.pvalue < 0.05))
significant.cutoff_rlog <- length(which(result.table2.sorted_rlog$fdr.pvalue < 0.05))

result.table2.significant_raw <- result.table2.sorted_raw[1:significant.cutoff_raw,]
result.table2.significant_vst <- result.table2.sorted_vst[1:significant.cutoff_vst,]
result.table2.significant_rlog <- result.table2.sorted_rlog[1:significant.cutoff_rlog,]

# get the number of genes with pval < 5 % and padj < 5 %
length(which(result.table2_raw$pvalue < 0.05))
length(which(result.table2_vst$pvalue < 0.05))
length(which(result.table2_rlog$pvalue < 0.05))

length(which(result.table2_raw$fdr.pvalue < 0.05))
length(which(result.table2_vst$fdr.pvalue < 0.05))
length(which(result.table2_rlog$fdr.pvalue < 0.05))

# create a data frame with summary of genes with pval < 5 % and padj < 5 %
sum_table_all <- matrix(c("number of genes with pval < 0.05", 
                          length(which(result.table2_raw$pvalue < 0.05)), length(which(result.table2_vst$pvalue < 0.05)), length(which(result.table2_rlog$pvalue < 0.05)),
                          "number of genes with padj (FDR) < 0.05 ", 
                          length(which(result.table2_raw$fdr.pvalue < 0.05)), length(which(result.table2_vst$fdr.pvalue < 0.05)), length(which(result.table2_rlog$fdr.pvalue < 0.05))), 
                        nrow = 2, ncol = 4, byrow = TRUE)

colnames(sum_table_all) <- c("metrics", "raw", "vst", "rlog")

# create output directory
dir.create(paste0(args[3], running))
setwd(paste0(args[3], running))

# save result tables
write.csv(sum_table_all, file = paste0("results_", running, "_summary.csv"))

# save list with pval and padj
write.csv(result.table2.significant_raw, file = paste0("results_raw_", running, "_list.csv"))
write.csv(result.table2.significant_vst, file = paste0("results_vst_", running, "_list.csv"))
write.csv(result.table2.significant_rlog, file = paste0("results_rlog_", running, "_list.csv"))


if(FALSE){
# switch to gene names and check chromosome location for obtained genes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#refseq <- c()

g_names_raw <- getBM(filters="refseq_mrna", attributes=c("hgnc_symbol", "chromosome_name"), values=result.table2.significant_raw$ID, mart=mart)
g_names_vst <- getBM(filters="refseq_mrna", attributes=c("hgnc_symbol", "chromosome_name"), values=result.table2.significant_vst$ID, mart=mart)
g_names_rlog <- getBM(filters="refseq_mrna", attributes=c("hgnc_symbol", "chromosome_name"), values=result.table2.significant_rlog$ID, mart=mart)

# save results
write.csv(g_names_raw, file = paste0("gene_table_", running, "_raw.csv"))
write.csv(g_names_vst, file = paste0("gene_table_", running, "_vst.csv"))
write.csv(g_names_rlog, file = paste0("gene_table_", running, "_rlog.csv"))


sink(paste0("gene_list_", running, "_raw.txt"))
g_names_raw$hgnc_symbol
sink()

sink(paste0("gene_list_", running, "_vst.txt"))
for (i in 1:length(g_names_vst$hgnc_symbol)) {
  print(g_names_vst$hgnc_symbol[i])
}
sink()

sink(paste0("gene_list_", running, "_rlog.txt"))
for (i in 1:length(g_names_rlog$hgnc_symbol)) {
  print(g_names_rlog$hgnc_symbol[i])
}
sink()
}

# use gene_list to create venn diagram:
# http://bioinformatics.psb.ugent.be/cgi-bin/liste/Venn/calculate_venn.htpl


# create venn diagram to see how many genes are in common between these 3 methods
#library(VennDiagram)
#venn.diagram(list(g_names_raw$hgnc_symbol, g_names_raw$hgnc_symbol, g_names_raw$hgnc_symbol),
#             category.names = c("Raw", "VST", "rlog"),
#             filename = 'venn_diagramm.png')
  

# save results
#write.csv(result.table2.sorted_final, file="results_table_Mann-Whitney_gender_norm.csv")
#write.csv(sum_table, file="summary_table_Mann-Whitney_gender_norm.csv")
#write.csv(result.table2.sorted_final, file="results_table_Mann-Whitney_gender_VST.csv")
#write.csv(sum_table, file="summary_table_Mann-Whitney_gender_VST.csv")

#write.csv(result.table2.sorted_final, file="results_table_Mann-Whitney_visual-loss_norm.csv")
#write.csv(sum_table, file="summary_table_Mann-Whitney_visual-loss_norm.csv")

#write.csv(result.table2.sorted_final, file="results_table_Mann-Whitney_visual-loss_VST.csv")
#write.csv(sum_table, file="summary_table_Mann-Whitney_visual-loss_VST.csv")

#---------------------------- use another method GSALightning package ----------------------------#

# source: https://www.bioconductor.org/packages/release/bioc/vignettes/GSALightning/inst/doc/vignette.html

# the GSALightning package offers the Mann-Whitney U test for single-gene testing. 
# Mann-Whitney U test is the non-parametric version of the independent t-test for two-sample problem. 
# To perform the Mann-Whitney U test, call the wilcoxTest() function:
  
#singleWilcox <- wilcoxTest(eset = dat, fac = factor(anno$gender..1.male..2.female.), tests = "unpaired")
singleWilcox <- wilcoxTest(eset = dat, fac = factor(anno$visual.loss.at.BL..0.no..1.yes.), tests = "unpaired")


# save results
#write.csv(singleWilcox, file="results_table_Mann-Whitney_gender_GSALightning_norm.csv")
write.csv(singleWilcox, file="results_table_Mann-Whitney_visual-loss_GSALightning_norm.csv")


