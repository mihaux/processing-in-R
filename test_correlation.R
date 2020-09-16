# script to perform the analysis on slide scores 

# => use Mann-Whitney test for for continuous variables (e.g. number of days on steroids or age) 
# => use Chi-squared or Fisher's exact test for categorical variables (Chi-squared if they fulfill the condition of size [i.e. at least 5 instances in the expected category], otherwise Fisher's exact)

# sources about Fisher's exact test: 
# => https://towardsdatascience.com/fishers-exact-test-in-r-independence-test-for-a-small-sample-56965db48e87

# sources about Chi-squared test: 
# => https://www.statsandr.com/blog/chi-square-test-of-independence-by-hand/

# Chi-square tests of independence test whether two qualitative variables are independent, 
# that is, whether there exists a relationship between two categorical variables. 
# In other words, this test is used to determine whether the values of one of the 2 qualitative variables depend on the values of the other qualitative variable.

# HYPOTHESES: null AND alternative
# (H0): the variables are independent, there is no relationship between the two categorical variables. Knowing the value of one variable does not help to predict the value of the other variable
# (H1): the variables are dependent, there is a relationship between the two categorical variables. Knowing the value of one variable helps to predict the value of the other variable

# How the test works?
# it works by comparing the observed frequencies (so the frequencies observed in your sample) 
# to the expected frequencies if there was no relationship between the two categorical variables 
# (so the expected frequencies if the null hypothesis was true).

# if the difference between the observed frequencies and the expected frequencies is small, 
# we cannot reject the null hypothesis of independence and thus we cannot reject the fact that the two variables are not related. 

# if the difference between the observed frequencies and the expected frequencies is large, 
# we can reject the null hypothesis of independence and thus we can conclude that the two variables are related.

# The threshold comes from the Chi-square distribution. This value, referred as the critical value, depends on the significance level 

# RESULTS of chi-square testing:
# If the test shows no association between the two variables (i.e., the variables are independent), 
# it means that knowing the value of one variable gives no information about the value of the other variable. 

# If the test shows a relationship between the variables (i.e., the variables are dependent), 
# it means that knowing the value of one variable provides information about the value of the other variable.

# sources about Mann-Whitney test: 
# => see another script (test_Mann_Whitney.R)

# from Wikipedia: 
# Fisher's exact test is a statistical significance test used in the analysis of contingency tables
# TBC
# TBC

# WHAT TO TEST: 
# => we want to check if there's correlation between the features (from slide score table) and the outcome (mostly visual loss, but could check others as well)
# outcome: visual_loss (the main one), jaw_claudication, ischaemic_features, can also test for gender

# => comparison groups labaling visual_loss vs. no visual_loss
# => data to be used: binary values from slide score tables

# EXAMPLE:
# => var_1: visual_loss         -> (outcome)
# => var_2: Media_destruction   -> (feature)
# goal: we want to see if there's correlation between Media_destruction and visual_loss

# create contingency table
#                             visual_loss
#                         |   0   |   1   |   TOTAL   |
# Media_destruction   0   |   11  |   8   |     19    |
#                     1   |   12  |   9   |     21    |
#                   TOTAL |   23  |   17  |     40    |

# compute the expected counts in the case the variables were independent.	
# the expected frequencies are computed for each subgroup one by one with the following formula:

# total of the row * total of the column / total number of observations

# expected counts table
#                             visual_loss
#                         |   0       |   1       |   TOTAL   |
# Media_destruction   0   |   10.925  |   8.075   |     19    |
#                     1   |   12.075  |   8.925   |     21    |
#                   TOTAL |   23      |   17      |     40    |

# NOTE: the Chi-square test of independence should only be done when 
# the expected frequencies in all groups are equal to or greater than 5. 

# If the condition is not met, the Fisher’s exact test is preferred.

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt"); library(biomaRt)
suppressMessages(library(DESeq2))

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

args <- commandArgs(trailingOnly = TRUE)

# let's run for visual_loss only

if (length(args)!=4) {
  stop("4 arguments must be supplied: 
       \n(1 - clinical data) path to .csv file with clinical data, 
       \n(2 - histological data) path to .csv file with slide score data,
       \n(3 - outcome) name of the clinical feature (outcome) for running and
       \n(4 - output) path where output files should be stored", call.=FALSE)
} # NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 4

args <- c(paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/data/metadata/slide_scores/slide_scores_v6.csv"),
          "gender", #"visual_loss", # "jaw_claudication", "ischaemic_features", "gender"
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/correlation_tests"))

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[4], sep="\n")
setwd(args[4])

# load clinical and histological data
df_clinical <- read.csv(args[1], row.names = 1)
df_histological <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames for clinical data
rownames(df_clinical) <- paste0("ID_", rownames(df_clinical))
rownames(df_histological) <- paste0("ID_", rownames(df_histological))

# NOTE: sample "ID_14058" is missing in df_histological, so it needs to be removed from df_clinical to make it equl length

# subset the column with outcome from clinical data
df_clinical_mod <- as.data.frame(df_clinical[,args[3]][-which(rownames(df_clinical) == "ID_14058")])
rownames(df_clinical_mod) <- rownames(df_clinical)[-which(rownames(df_clinical) == "ID_14058")]
colnames(df_clinical_mod) <- args[3]

# create contingency tables for each histological feature and visual_loss
# to make sure that there are at least 5 instances in each expected group

# features selected for statistical testing
tb_GCA_present <- table(cbind(df_clinical_mod, df_histological$GCA_present))
tb_Giant_cells <- table(cbind(df_clinical_mod, df_histological$Giant_cells))
tb_Infiltrate_around_vasa_vasorum <- table(cbind(df_clinical_mod, df_histological$Infiltrate_around_vasa_vasorum))
tb_Media_destruction <- table(cbind(df_clinical_mod, df_histological$Media_destruction))

cat("Running for: ", args[3])

chi_GCA_present                     <- chisq.test(tb_GCA_present)
chi_Giant_cells                     <- chisq.test(tb_Giant_cells)
chi_Infiltrate_around_vasa_vasorum  <- chisq.test(tb_Infiltrate_around_vasa_vasorum)
chi_Media_destruction               <- chisq.test(tb_Media_destruction)

# other features
tb_Any.granulomatous.infiltrate           <- table(cbind(df_clinical_mod, df_histological$Any.granulomatous.infiltrate))
tb_Granulomatous.infiltrate.in.adventitia <- table(cbind(df_clinical_mod, df_histological$Granulomatous.infiltrate.in.adventitia))
tb_Granulomatous.infiltrate.in.media      <- table(cbind(df_clinical_mod, df_histological$Granulomatous.infiltrate.in.media))
tb_Granulomatous.infiltrate.in.intima     <- table(cbind(df_clinical_mod, df_histological$Granulomatous.infiltrate.in.intima))
tb_Any.lymphocytic.infiltrate             <- table(cbind(df_clinical_mod, df_histological$Any.lymphocytic.infiltrate))
tb_Lymphocytic.infiltrate.in.adventitia   <- table(cbind(df_clinical_mod, df_histological$Lymphocytic.infiltrate.in.adventitia))
tb_Lymphocytic.infiltrate.in.media        <- table(cbind(df_clinical_mod, df_histological$Lymphocytic.infiltrate.in.media))
tb_Lymphocytic.infiltrate.in.intima       <- table(cbind(df_clinical_mod, df_histological$Lymphocytic.infiltrate.in.intima))
tb_Aggregates                             <- table(cbind(df_clinical_mod, df_histological$Aggregates))
tb_PALI                                   <- table(cbind(df_clinical_mod, df_histological$PALI))
tb_Neoangiogenesis                        <- table(cbind(df_clinical_mod, df_histological$Neoangiogenesis))
tb_Hyperplasia                            <- table(cbind(df_clinical_mod, df_histological$Hyperplasia))
tb_Fibrosis                               <- table(cbind(df_clinical_mod, df_histological$Fibrosis))
tb_Oedema                                 <- table(cbind(df_clinical_mod, df_histological$Oedema))

chi_Any.granulomatous.infiltrate            <- chisq.test(tb_Any.granulomatous.infiltrate)     
chi_Granulomatous.infiltrate.in.adventitia  <- chisq.test(tb_Granulomatous.infiltrate.in.adventitia)
chi_Granulomatous.infiltrate.in.media       <- chisq.test(tb_Granulomatous.infiltrate.in.media)     
chi_Granulomatous.infiltrate.in.intima      <- chisq.test(tb_Granulomatous.infiltrate.in.intima)
chi_Any.lymphocytic.infiltrate              <- chisq.test(tb_Any.lymphocytic.infiltrate)  
chi_Lymphocytic.infiltrate.in.adventitia    <- chisq.test(tb_Lymphocytic.infiltrate.in.adventitia)
chi_Lymphocytic.infiltrate.in.media         <- chisq.test(tb_Lymphocytic.infiltrate.in.media)
chi_Lymphocytic.infiltrate.in.intima        <- chisq.test(tb_Lymphocytic.infiltrate.in.intima) 
chi_Aggregates                              <- chisq.test(tb_Aggregates)                            
chi_PALI                                    <- chisq.test(tb_PALI)                            
chi_Neoangiogenesis                         <- chisq.test(tb_Neoangiogenesis)                 
chi_Hyperplasia                             <- chisq.test(tb_Hyperplasia)                 
chi_Fibrosis                                <- chisq.test(tb_Fibrosis)                          
chi_Oedema                                  <- chisq.test(tb_Oedema)

# create a list with all tests
chi_all <- list(chi_GCA_present, 
                chi_Giant_cells, 
                chi_Infiltrate_around_vasa_vasorum, 
                chi_Media_destruction, 
                chi_Any.granulomatous.infiltrate, 
                chi_Granulomatous.infiltrate.in.adventitia, 
                chi_Granulomatous.infiltrate.in.media, 
                chi_Granulomatous.infiltrate.in.intima, 
                chi_Any.lymphocytic.infiltrate, 
                chi_Lymphocytic.infiltrate.in.adventitia, 
                chi_Lymphocytic.infiltrate.in.media, 
                chi_Lymphocytic.infiltrate.in.intima, 
                chi_Aggregates, 
                chi_PALI, 
                chi_Neoangiogenesis, 
                chi_Hyperplasia, 
                chi_Fibrosis, 
                chi_Oedema)

# extract data.names
unlist(lapply(chi_all, function(x) x$data.name))

# extract p-values
unlist(lapply(chi_all, function(x) x$p.value))

# check if the expected frequencies in all groups are equal to or greater than 5
all(as.vector(chi_Oedema$expected) >= 5)

chi_results <- as.data.frame(cbind(unlist(lapply(chi_all, function(x) x$data.name)), 
                                   unlist(lapply(chi_all, function(x) x$p.value)),
                                   unlist(lapply(chi_all, function(x) all(as.vector(x$expected) >= 5)))))

colnames(chi_results) <- c("features", "p.value", "group_frequencies>=5")

# save chi_squared results
#write.csv(chi_results, file=paste0("results_Chi_squared_", args[3], ".csv"))

### Fisher's Exact Test

fish_GCA_present                      <- fisher.test(tb_GCA_present)
fish_Giant_cells                      <- fisher.test(tb_Giant_cells)
fish_Infiltrate_around_vasa_vasorum   <- fisher.test(tb_Infiltrate_around_vasa_vasorum)
fish_Media_destruction                <- fisher.test(tb_Media_destruction)

fish_Any.granulomatous.infiltrate             <- fisher.test(tb_Any.granulomatous.infiltrate)		
fish_Granulomatous.infiltrate.in.adventitia   <- fisher.test(tb_Granulomatous.infiltrate.in.adventitia)
#fish_Granulomatous.infiltrate.in.media        <- fisher.test(tb_Granulomatous.infiltrate.in.media)
fish_Granulomatous.infiltrate.in.intima       <- fisher.test(tb_Granulomatous.infiltrate.in.intima)	
fish_Any.lymphocytic.infiltrate               <- fisher.test(tb_Any.lymphocytic.infiltrate)				
fish_Lymphocytic.infiltrate.in.adventitia     <- fisher.test(tb_Lymphocytic.infiltrate.in.adventitia)
fish_Lymphocytic.infiltrate.in.media          <- fisher.test(tb_Lymphocytic.infiltrate.in.media)
fish_Lymphocytic.infiltrate.in.intima         <- fisher.test(tb_Lymphocytic.infiltrate.in.intima)		
fish_Aggregates                               <- fisher.test(tb_Aggregates)								
fish_PALI                                     <- fisher.test(tb_PALI)										
fish_Neoangiogenesis                          <- fisher.test(tb_Neoangiogenesis)							
fish_Hyperplasia                              <- fisher.test(tb_Hyperplasia)								
fish_Fibrosis                                 <- fisher.test(tb_Fibrosis)									
fish_Oedema                                   <- fisher.test(tb_Oedema)

fish_all <- list(fish_GCA_present, 
                 fish_Giant_cells, 
                 fish_Infiltrate_around_vasa_vasorum, 
                 fish_Media_destruction, 
                 fish_Any.granulomatous.infiltrate, 
                 fish_Granulomatous.infiltrate.in.adventitia, 
                 #fish_Granulomatous.infiltrate.in.media, 
                 fish_Granulomatous.infiltrate.in.intima, 
                 fish_Any.lymphocytic.infiltrate, 
                 fish_Lymphocytic.infiltrate.in.adventitia, 
                 fish_Lymphocytic.infiltrate.in.media, 
                 fish_Lymphocytic.infiltrate.in.intima, 
                 fish_Aggregates, 
                 fish_PALI, 
                 fish_Neoangiogenesis, 
                 fish_Hyperplasia, 
                 fish_Fibrosis, 
                 fish_Oedema)

# extract data.names
unlist(lapply(fish_all, function(x) x$data.name))

# extract p-values
unlist(lapply(fish_all, function(x) x$p.value))

# extract odds ratios
unlist(lapply(fish_all, function(x) x$estimate))

# ectract 95% confidence interval [CI] 
unlist(lapply(fish_all, function(x) x$conf.int[1]))
unlist(lapply(fish_all, function(x) x$conf.int[2]))

fish_results <- as.data.frame(cbind(unlist(lapply(fish_all, function(x) x$data.name)), 
                                    unlist(lapply(fish_all, function(x) x$p.value)),
                                    as.vector(unlist(lapply(fish_all, function(x) x$estimate))),
                                    unlist(lapply(fish_all, function(x) x$conf.int[1])),
                                    unlist(lapply(fish_all, function(x) x$conf.int[2]))))

colnames(fish_results) <- c("features", "p.value", "odds.ratio", "conf.int_1", "conf.int_2")

# save chi_squared results
#write.csv(fish_results, file=paste0("results_Fishers_exact_", args[3], ".csv"))

################################################################################################
### run for histological features with multiple categories
################################################################################################

tb_Barcelona.score    <- table(cbind(df_clinical_mod, df_histological$Barcelona.score))
tb_Adventitia.pattern <- table(cbind(df_clinical_mod, df_histological$Adventitia.pattern))
tb_Media_pattern      <- table(cbind(df_clinical_mod, df_histological$Media_pattern))
tb_Intima_pattern     <- table(cbind(df_clinical_mod, df_histological$Intima_pattern))
tb_Occlusion_grade    <- table(cbind(df_clinical_mod, df_histological$Occlusion_grade))

chi_Barcelona.score     <- chisq.test(tb_Barcelona.score)
chi_Adventitia.pattern  <- chisq.test(tb_Adventitia.pattern)
chi_Media_pattern       <- chisq.test(tb_Media_pattern)
chi_Intima_pattern      <- chisq.test(tb_Intima_pattern)
chi_Occlusion_grade     <- chisq.test(tb_Occlusion_grade)

# create a list with all tests
chi_all_multiple <- list(chi_Barcelona.score,
                        chi_Adventitia.pattern,
                        chi_Media_pattern, 
                        chi_Intima_pattern,
                        chi_Occlusion_grade)

# extract data.names
unlist(lapply(chi_all_multiple, function(x) x$data.name))

# extract p-values
unlist(lapply(chi_all_multiple, function(x) x$p.value))

# check if the expected frequencies in all groups are equal to or greater than 5
#all(as.vector(chi_Oedema$expected) >= 5)

chi_results_multiple <- as.data.frame(cbind(unlist(lapply(chi_all_multiple, function(x) x$data.name)), 
                                            unlist(lapply(chi_all_multiple, function(x) x$p.value)),
                                            unlist(lapply(chi_all_multiple, function(x) all(as.vector(x$expected) >= 5)))))

colnames(chi_results_multiple) <- c("features", "p.value", "group_frequencies>=5")

# save chi_squared results
write.csv(chi_results_multiple, file=paste0("results_Chi_squared_multiple_", args[3], ".csv"))



fish_Barcelona.score     <- fisher.test(tb_Barcelona.score)
fish_Adventitia.pattern  <- fisher.test(tb_Adventitia.pattern)
fish_Media_pattern       <- fisher.test(tb_Media_pattern)
fish_Intima_pattern      <- fisher.test(tb_Intima_pattern)
fish_Occlusion_grade     <- fisher.test(tb_Occlusion_grade)

# transform to binary and run the same analysis, to make it 2 comparison groups
Occlusion_grade_new <- c(1:40)
Occlusion_grade_new[which(df_histological$Occlusion_grade < 3)] <- 0      # less than 3           => replace with 0
Occlusion_grade_new[which(df_histological$Occlusion_grade >= 3)] <- 1     # equal or more than 3  => replace with 1

Intima_pattern_new <- c(1:40)
Intima_pattern_new[which(df_histological$Intima_pattern <= 1)] <- 0       # 0 and 1                => replace with 0
Intima_pattern_new[which(df_histological$Intima_pattern > 1)] <- 1        # 2 and 3                => replace with 1

Adventitia_pattern_new <- c(1:40)
Adventitia_pattern_new[which(df_histological$Adventitia.pattern <= 1)] <- 0       # 0 and 1                => replace with 0
Adventitia_pattern_new[which(df_histological$Adventitia.pattern > 1)] <- 1        # 2 and 3                => replace with 1

Media_pattern_new <- c(1:40)
Media_pattern_new[which(df_histological$Media_pattern <= 1)] <- 0       # 0 and 1                => replace with 0
Media_pattern_new[which(df_histological$Media_pattern > 1)] <- 1        # 2 and 3                => replace with 1

tb_Adventitia.pattern_trans <- table(cbind(df_clinical_mod, Adventitia_pattern_new))
tb_Media_pattern_trans      <- table(cbind(df_clinical_mod, Media_pattern_new))
tb_Intima_pattern_trans     <- table(cbind(df_clinical_mod, Intima_pattern_new))
tb_Occlusion_grade_trans    <- table(cbind(df_clinical_mod, Occlusion_grade_new))

chi_Adventitia.pattern_trans  <- chisq.test(tb_Adventitia.pattern_trans)
chi_Media_pattern_trans       <- chisq.test(tb_Media_pattern_trans)
chi_Intima_pattern_trans      <- chisq.test(tb_Intima_pattern_trans)
chi_Occlusion_grade_trans     <- chisq.test(tb_Occlusion_grade_trans)

# create a list with all tests
chi_all_multiple_trans <- list(chi_Adventitia.pattern_trans,
                                chi_Media_pattern_trans, 
                                chi_Intima_pattern_trans,
                                chi_Occlusion_grade_trans)

# extract data.names
unlist(lapply(chi_all_multiple_trans, function(x) x$data.name))

# extract p-values
unlist(lapply(chi_all_multiple_trans, function(x) x$p.value))

# check if the expected frequencies in all groups are equal to or greater than 5
#all(as.vector(chi_Oedema$expected) >= 5)

chi_results_multiple_trans <- as.data.frame(cbind(unlist(lapply(chi_all_multiple_trans, function(x) x$data.name)), 
                                            unlist(lapply(chi_all_multiple_trans, function(x) x$p.value)),
                                            unlist(lapply(chi_all_multiple_trans, function(x) all(as.vector(x$expected) >= 5)))))

colnames(chi_results_multiple_trans) <- c("features", "p.value", "group_frequencies>=5")

# save chi_squared results
write.csv(chi_results_multiple_trans, file=paste0("results_Chi_squared_multiple_trans_", args[3], ".csv"))


fish_Adventitia.pattern_trans  <- fisher.test(tb_Adventitia.pattern_trans)
fish_Media_pattern_trans       <- fisher.test(tb_Media_pattern_trans)
fish_Intima_pattern_trans      <- fisher.test(tb_Intima_pattern_trans)
fish_Occlusion_grade_trans     <- fisher.test(tb_Occlusion_grade_trans)

################################################################################################
### run Mann-Whitney for continous features 
################################################################################################
# => "year.TAB.sample.was.collected"                    
# => "number.of.days.on.steroids.at.TAB"                
# => "number.of.days.between.TAB.and.BL.blood.sample"   
# => "age.at.BL"

# visual loss vs. no visual loss 

# create a subset with these 4 features of interest
df_visual_loss <- as.data.frame(cbind(df_clinical$visual_loss, 
                                      df_clinical$year.TAB.sample.was.collected,
                                      df_clinical$number.of.days.on.steroids.at.TAB,
                                      df_clinical$number.of.days.between.TAB.and.BL.blood.sample,
                                      df_clinical$age.at.BL))

colnames(df_visual_loss) <- c("visual_loss", "year_TAB", "nb_days_steroids", "nb_days_TAB_blood", "age")

# get vectors with age for patients with visual_loos and no_visual_loss
cond_1 <- which(df_visual_loss$visual_loss == 0)
cond_2 <- which(df_visual_loss$visual_loss == 1)


mann_whitney_year_TAB_visual_loss <- wilcox.test(as.numeric(df_visual_loss$year_TAB[cond_1]), 
                                     as.numeric(df_visual_loss$year_TAB[cond_2]), 
                                     alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_steroids_visual_loss <- wilcox.test(as.numeric(df_visual_loss$nb_days_steroids[cond_1]), 
                                             as.numeric(df_visual_loss$nb_days_steroids[cond_2]), 
                                             alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_TAB_blood_visual_loss <- wilcox.test(as.numeric(df_visual_loss$nb_days_TAB_blood[cond_1]), 
                                             as.numeric(df_visual_loss$nb_days_TAB_blood[cond_2]), 
                                             alternative = "two.sided",exact = FALSE) 
# NOTE there's one missing values in nb_days_TAB_blood

mann_whitney_age_visual_loss <- wilcox.test(as.numeric(df_visual_loss$age[cond_1]), 
                                as.numeric(df_visual_loss$age[cond_2]), 
                                alternative = "two.sided",exact = FALSE) 

# jaw_claudication vs. no_jaw_claudication

# create a subset with these 4 features of interest
df_jaw_claudication <- as.data.frame(cbind(df_clinical$jaw_claudication, 
                                           df_clinical$year.TAB.sample.was.collected,
                                           df_clinical$number.of.days.on.steroids.at.TAB,
                                           df_clinical$number.of.days.between.TAB.and.BL.blood.sample,
                                           df_clinical$age.at.BL))

colnames(df_jaw_claudication) <- c("jaw_claudication", "year_TAB", "nb_days_steroids", "nb_days_TAB_blood", "age")

# get vectors with age for patients with visual_loos and no_jaw_claudication
cond_1 <- which(df_jaw_claudication$jaw_claudication == 0)
cond_2 <- which(df_jaw_claudication$jaw_claudication == 1)


mann_whitney_year_TAB_jaw_claudication <- wilcox.test(as.numeric(df_jaw_claudication$year_TAB[cond_1]), 
                                                      as.numeric(df_jaw_claudication$year_TAB[cond_2]), 
                                                      alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_steroids_jaw_claudication <- wilcox.test(as.numeric(df_jaw_claudication$nb_days_steroids[cond_1]), 
                                                              as.numeric(df_jaw_claudication$nb_days_steroids[cond_2]), 
                                                              alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_TAB_blood_jaw_claudication <- wilcox.test(as.numeric(df_jaw_claudication$nb_days_TAB_blood[cond_1]), 
                                                               as.numeric(df_jaw_claudication$nb_days_TAB_blood[cond_2]), 
                                                               alternative = "two.sided",exact = FALSE) 
# NOTE there's one missing values in nb_days_TAB_blood

mann_whitney_age_jaw_claudication <- wilcox.test(as.numeric(df_jaw_claudication$age[cond_1]), 
                                                 as.numeric(df_jaw_claudication$age[cond_2]), 
                                                 alternative = "two.sided",exact = FALSE) 

# ischaemic_features vs. no_ischaemic_features

# create a subset with these 4 features of interest
df_ischaemic_features <- as.data.frame(cbind(df_clinical$ischaemic_features, 
                                             df_clinical$year.TAB.sample.was.collected,
                                             df_clinical$number.of.days.on.steroids.at.TAB,
                                             df_clinical$number.of.days.between.TAB.and.BL.blood.sample,
                                             df_clinical$age.at.BL))

colnames(df_ischaemic_features) <- c("ischaemic_features", "year_TAB", "nb_days_steroids", "nb_days_TAB_blood", "age")

# get vectors with age for patients with visual_loos and no_ischaemic_features
cond_1 <- which(df_ischaemic_features$ischaemic_features == 0)
cond_2 <- which(df_ischaemic_features$ischaemic_features == 1)


mann_whitney_year_TAB_ischaemic_features <- wilcox.test(as.numeric(df_ischaemic_features$year_TAB[cond_1]), 
                                                        as.numeric(df_ischaemic_features$year_TAB[cond_2]), 
                                                        alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_steroids_ischaemic_features <- wilcox.test(as.numeric(df_ischaemic_features$nb_days_steroids[cond_1]), 
                                                                as.numeric(df_ischaemic_features$nb_days_steroids[cond_2]), 
                                                                alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_TAB_blood_ischaemic_features <- wilcox.test(as.numeric(df_ischaemic_features$nb_days_TAB_blood[cond_1]), 
                                                                 as.numeric(df_ischaemic_features$nb_days_TAB_blood[cond_2]), 
                                                                 alternative = "two.sided",exact = FALSE) 
# NOTE there's one missing values in nb_days_TAB_blood

mann_whitney_age_ischaemic_features <- wilcox.test(as.numeric(df_ischaemic_features$age[cond_1]), 
                                                   as.numeric(df_ischaemic_features$age[cond_2]), 
                                                   alternative = "two.sided",exact = FALSE) 

# gender: male vs. female


# create a subset with these 4 features of interest
df_gender <- as.data.frame(cbind(df_clinical$gender, 
                                 df_clinical$year.TAB.sample.was.collected,
                                 df_clinical$number.of.days.on.steroids.at.TAB,
                                 df_clinical$number.of.days.between.TAB.and.BL.blood.sample,
                                 df_clinical$age.at.BL))

colnames(df_gender) <- c("gender", "year_TAB", "nb_days_steroids", "nb_days_TAB_blood", "age")

# get vectors with age for patients with visual_loos and no_gender
cond_1 <- which(df_gender$gender == 1)
cond_2 <- which(df_gender$gender == 2)


mann_whitney_year_TAB_gender <- wilcox.test(as.numeric(df_gender$year_TAB[cond_1]), 
                                            as.numeric(df_gender$year_TAB[cond_2]), 
                                            alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_steroids_gender <- wilcox.test(as.numeric(df_gender$nb_days_steroids[cond_1]), 
                                                    as.numeric(df_gender$nb_days_steroids[cond_2]), 
                                                    alternative = "two.sided",exact = FALSE) 

mann_whitney_nb_days_TAB_blood_gender <- wilcox.test(as.numeric(df_gender$nb_days_TAB_blood[cond_1]), 
                                                     as.numeric(df_gender$nb_days_TAB_blood[cond_2]), 
                                                     alternative = "two.sided",exact = FALSE) 
# NOTE there's one missing values in nb_days_TAB_blood

mann_whitney_age_gender <- wilcox.test(as.numeric(df_gender$age[cond_1]), 
                                       as.numeric(df_gender$age[cond_2]), 
                                       alternative = "two.sided",exact = FALSE) 

# create table with results:
mann_whitney_year_TAB_visual_loss 
mann_whitney_nb_days_steroids_visual_loss 
mann_whitney_nb_days_TAB_blood_visual_loss 
mann_whitney_age_visual_loss

mann_whitney_visual_loss_all <- c(mann_whitney_year_TAB_visual_loss$p.value,
                                  mann_whitney_nb_days_steroids_visual_loss$p.value, 
                                  mann_whitney_nb_days_TAB_blood_visual_loss$p.value, 
                                  mann_whitney_age_visual_loss$p.value)

mann_whitney_jaw_claudication_all <- c(mann_whitney_year_TAB_jaw_claudication$p.value,
                                        mann_whitney_nb_days_steroids_jaw_claudication$p.value, 
                                        mann_whitney_nb_days_TAB_blood_jaw_claudication$p.value, 
                                        mann_whitney_age_jaw_claudication$p.value)

mann_whitney_ischaemic_features_all <- c(mann_whitney_year_TAB_ischaemic_features$p.value,
                                          mann_whitney_nb_days_steroids_ischaemic_features$p.value, 
                                          mann_whitney_nb_days_TAB_blood_ischaemic_features$p.value, 
                                          mann_whitney_age_ischaemic_features$p.value)

mann_whitney_gender_all <- c(mann_whitney_year_TAB_gender$p.value,
                              mann_whitney_nb_days_steroids_gender$p.value, 
                              mann_whitney_nb_days_TAB_blood_gender$p.value, 
                              mann_whitney_age_gender$p.value)

results_mann_whitney <- cbind(mann_whitney_visual_loss_all, 
                            mann_whitney_jaw_claudication_all,
                            mann_whitney_ischaemic_features_all,
                            mann_whitney_gender_all)

rownames(results_mann_whitney) <- c("year.TAB.sample.was.collected", 
                                  "number.of.days.on.steroids.at.TAB", 
                                  "number.of.days.between.TAB.and.BL.blood.sample", 
                                  "age.at.BL")

# save Mann-Whitney results
write.csv(results_mann_whitney, file=paste0("results_Mann_Whitney.csv"))

