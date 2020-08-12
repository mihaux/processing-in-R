
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

# MArk's response:
# So one question is what difference the two different approaches make. You've noted that the number achieving p<0.05 (you should refer to it thus rather than as 5%) differs 
# I wouldn't expect them to be exactly the same, but I'd expect a reasonable correlation.

# TODO: could you plot the p-values (or, better, log10 p-values) against one another 

# gender using normalised counts against gender using vst transformed counts 
# then and calculate the correlation between them 

# visual loss using normalised counts against visual loss using vst transformed counts 
# then and calculate the correlation between them 

# load data
gender_norm <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/results_table_Mann-Whitney_gender_norm.csv"), row.names = 1)
gender_vst <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/results_table_Mann-Whitney_gender_VST.csv"), row.names = 1)

# sort by transcript names to make it match in terms of transcripts IDs
gender_norm_ord <- gender_norm[order(gender_norm$ID),] 
gender_vst_ord <- gender_vst[order(gender_vst$ID),]

visual_norm <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/results_table_Mann-Whitney_visual-loss_norm.csv"), row.names = 1)
visual_vst <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/results_table_Mann-Whitney_visual-loss_VST.csv"), row.names = 1)

visual_norm_ord <- visual_norm[order(visual_norm$ID),] 
visual_vst_ord <- visual_vst[order(visual_vst$ID),]

# transform to log10 p-values and plot one again another
# NOTE: the tables are ordered and do not match in transcript IDs
pdf(file="plots_log10pval_gender.pdf", width=12, height=12)
plot(log10(gender_norm_ord$pvalue), log10(gender_vst_ord$pvalue), main = "Gender - log10(pvalue)")
#plot(log10(gender_vst_ord$pvalue), log10(gender_norm_ord$pvalue), main = "Gender - log10(pvalue)")

dev.off()

pdf(file="plots_log10pval_visual_loss.pdf", width=12, height=12)
plot(log10(visual_norm_ord$pvalue), log10(visual_vst_ord$pvalue), main = "Visual loss - log10(pvalue)")
#plot(log10(visual_vst_ord$pvalue), log10(visual_norm_ord$pvalue), main = "Visual loss - log10(pvalue)")
dev.off()

pearson_gender <- cor.test(log10(gender_norm_ord$pvalue), log10(gender_vst_ord$pvalue), method = c("pearson"))
write.csv(unlist(pearson_gender), file="cor_test_pearson_log10pval_gender.csv")


kot <- cor.test(log10(gender_norm_ord$pvalue), log10(gender_vst_ord$pvalue), method = c("pearson"))

cor.test(log10(visual_norm_ord$pvalue), log10(visual_vst_ord$pvalue), method = c("pearson"))

#cor.test(log10(gender_norm_ord$pvalue), log10(gender_vst_ord$pvalue), method = c("kendall"))
#cor.test(log10(gender_norm_ord$pvalue), log10(gender_vst_ord$pvalue), method = c("spearman"))


#cor.test(log10(visual_norm_ord$pvalue), log10(visual_vst_ord$pvalue), method = c("kendall"))
#cor.test(log10(visual_norm_ord$pvalue), log10(visual_vst_ord$pvalue), method = c("spearman"))




