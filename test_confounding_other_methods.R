# script to assess confounding based on Cancer_Epidemiology_Principles_and_Methods.pdf (Chapter 14)
# stored at: /Users/michal/Documents/OneDrive - University of Leeds/x_other_useful_pdfs/

# features suspected of being confounders: age, gender, presents of steroid

# 2 main methods: stratification and regression modelling

# other possible ways: 
# => check if any of the confounder can predict the outcome using some regression analysis
# => run a paired test to see if each of the suspected confounders is associated with visual loss 

# resources about paired test:
# https://www.statisticssolutions.com/manova-analysis-paired-sample-t-test/
# http://www.sthda.com/english/wiki/paired-samples-t-test-in-r => perform in R

# other resources about confounding:
# http://bayes.cs.ucla.edu/BOOK-2K/ch6-2.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6302489/ 

library(epitools)

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
       \n(1 - input) path to _x.csv file with count data, 
       \n(2 - annotation) path to _x.csv annotation file and
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/DESeq2/norm_counts.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/"))

# Example of usage: 
# Rscript test_confounding.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load normalised counts from DESeq
df <- read.csv(args[1], row.names = 1)

# change data format to matrix and integer
dat <- as.matrix(df)   
storage.mode(dat) <- "integer"             # class: matrix | type: integer

# load annotation (clinical) data which will be used for coldata
anno <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# create table for odds ratio for gender 
# NOTE: this method can be used for 'number.of.days.on.steroids.at.TAB' or 'age' as they are not binary
vl_male <- length(which(anno$visual.loss.at.BL..0.no..1.yes. == 1 & anno$gender..1.male..2.female. == 1))       # visual loss: male
vl_female <- length(which(anno$visual.loss.at.BL..0.no..1.yes. == 1 & anno$gender..1.male..2.female. == 2))     # visual loss: female
ctrl_male <- length(which(anno$visual.loss.at.BL..0.no..1.yes. == 0 & anno$gender..1.male..2.female. == 1))     # control (no visual loss): male
ctrl_female <- length(which(anno$visual.loss.at.BL..0.no..1.yes. == 0 & anno$gender..1.male..2.female. == 2))   # control (no visual loss): female

#                   male    | female    | total
# visual_loss       7       |   10      | 17              7/17 = 41.18 %
# no_visual_loss    9       |   15      | 24              9/24 = 37.5 %
# total             16      |   25      | 41

# ( 7/10 ) / ( 9/15 ) = 1.166667


genders <- c("male", "female")
conds <- c("visual_loss", "no_visual_loss")

dat <- matrix(c(vl_male, vl_female, ctrl_male, ctrl_female), nrow = 2, ncol = 2, byrow = TRUE)
dimnames(dat) <- list("Conditions" = conds, "Gender" = genders)

# compute the odds ratio
or_fit <- oddsratio(dat)

# another method
#library(questionr)
#odds.ratio(dat)

### use Simple Linear Regression 
# source: https://www.r-bloggers.com/r-tutorial-series-simple-linear-regression/


# lm() function => “linear model” 
# => the function can be used to create a simple regression model. 

# The following list explains the two most commonly used parameters.
# -> formula: describes the model
# NOTE: the formula argument follows a specific format. For simple linear regression, this is “YVAR ~ XVAR” where 
# => YVAR is the dependent, or predicted, variable and 
# => XVAR is the independent, or predictor, variable.

# -> data: the variable that contains the dataset
# NOTE: It is recommended that you save a newly created linear model into a variable. 
# By doing so, the model can be used in subsequent calculations and analyses without having to retype the entire lm() function each time. 

# create a linear model using 
lm(FORMULA, DATAVAR)

#predict the fall enrollment (ROLL) using the unemployment rate (UNEM)
linearModelVar <- lm(ROLL ~ UNEM, datavar)
#display linear model
linearModelVar










