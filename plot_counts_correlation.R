# script to plot counts from one dataset against counts from another dataset (to see if there's correlation)

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

args <- commandArgs(trailingOnly = TRUE)

# let's run for visual_loss only

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with dataset_1, 
       \n(2 - input) path to .csv file with dataset_2 and
       \n(3 - output) path where output files should be stored", call.=FALSE)
} # NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/run_12_Aug20/5_counting/SE_all/all_counts_nodups_SE_all_mod.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/5_counting/PE_all/all_counts_nodups_PE_all_mod.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/other_plots/"))

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load both datasets
df_1 <- read.csv(args[1], row.names = 1)
df_2 <- read.csv(args[2], row.names = 1)

# check if colnames and rownames match between df_1 and df_2
any(rownames(df_1) == rownames(df_2))
any(colnames(df_1) == colnames(df_2))

# calculate corellation coefficients for each sample between df_1 and df_2
spears <- c()
pears <- c()
#ken <- c()
for (i in 1:ncol(df_1)) {
  pears[i] <- cor(df_1[,i], df_2[,i], method = 'pearson')
  spears[i] <- cor(df_1[,i], df_2[,i], method = 'spearman')
#  ken[i] <- cor(df_1[,i], df_2[,i], method = 'kendall')      # takes too long to process
}

df_coefficients <- as.data.frame(cbind(colnames(df_1) ,pears, spears))
colnames(df_coefficients) <- c("IDs", "pearson", "spearman")

write.csv(df_coefficients, file = "coefficients.csv")

# for each sample, plot SE agains PE
df_temp <- list()
for (i in 1:ncol(df_1)) {
  df_temp[[i]] <- as.data.frame(cbind(log(df_1[,i]+1), log(df_2[,i]+1)))
  colnames(df_temp[[i]]) <- c(paste0(colnames(df_1)[i], "_SE"), paste0(colnames(df_1)[i], "_PE"))
}

for (i in 1:41) {
  print(names(df_temp[[i]]))
}

png("correlation_SE_PE_ID_11026.png"); ggplot(df_temp[[1]], aes(ID_11026_SE, ID_11026_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off() 
png("correlation_SE_PE_ID_11028.png"); ggplot(df_temp[[2]], aes(ID_11028_SE, ID_11028_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_11037.png"); ggplot(df_temp[[3]], aes(ID_11037_SE, ID_11037_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12223.png"); ggplot(df_temp[[4]], aes(ID_12223_SE, ID_12223_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12330.png"); ggplot(df_temp[[5]], aes(ID_12330_SE, ID_12330_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12331.png"); ggplot(df_temp[[6]], aes(ID_12331_SE, ID_12331_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12388.png"); ggplot(df_temp[[7]], aes(ID_12388_SE, ID_12388_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12389.png"); ggplot(df_temp[[8]], aes(ID_12389_SE, ID_12389_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12849.png"); ggplot(df_temp[[9]], aes(ID_12849_SE, ID_12849_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12855.png"); ggplot(df_temp[[10]], aes(ID_12855_SE, ID_12855_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12898.png"); ggplot(df_temp[[11]], aes(ID_12898_SE, ID_12898_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_12954.png"); ggplot(df_temp[[12]], aes(ID_12954_SE, ID_12954_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14058.png"); ggplot(df_temp[[13]], aes(ID_14058_SE, ID_14058_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14455.png"); ggplot(df_temp[[14]], aes(ID_14455_SE, ID_14455_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14912.png"); ggplot(df_temp[[15]], aes(ID_14912_SE, ID_14912_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14914.png"); ggplot(df_temp[[16]], aes(ID_14914_SE, ID_14914_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14915.png"); ggplot(df_temp[[17]], aes(ID_14915_SE, ID_14915_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14916.png"); ggplot(df_temp[[18]], aes(ID_14916_SE, ID_14916_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14919.png"); ggplot(df_temp[[19]], aes(ID_14919_SE, ID_14919_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14920.png"); ggplot(df_temp[[20]], aes(ID_14920_SE, ID_14920_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14922.png"); ggplot(df_temp[[21]], aes(ID_14922_SE, ID_14922_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14925.png"); ggplot(df_temp[[22]], aes(ID_14925_SE, ID_14925_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14926.png"); ggplot(df_temp[[23]], aes(ID_14926_SE, ID_14926_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14927.png"); ggplot(df_temp[[24]], aes(ID_14927_SE, ID_14927_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14929.png"); ggplot(df_temp[[25]], aes(ID_14929_SE, ID_14929_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14930.png"); ggplot(df_temp[[26]], aes(ID_14930_SE, ID_14930_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_14931.png"); ggplot(df_temp[[27]], aes(ID_14931_SE, ID_14931_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16082.png"); ggplot(df_temp[[28]], aes(ID_16082_SE, ID_16082_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16098.png"); ggplot(df_temp[[29]], aes(ID_16098_SE, ID_16098_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16109.png"); ggplot(df_temp[[30]], aes(ID_16109_SE, ID_16109_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16111.png"); ggplot(df_temp[[31]], aes(ID_16111_SE, ID_16111_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16125.png"); ggplot(df_temp[[32]], aes(ID_16125_SE, ID_16125_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16137.png"); ggplot(df_temp[[33]], aes(ID_16137_SE, ID_16137_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16142.png"); ggplot(df_temp[[34]], aes(ID_16142_SE, ID_16142_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16153.png"); ggplot(df_temp[[35]], aes(ID_16153_SE, ID_16153_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16154.png"); ggplot(df_temp[[36]], aes(ID_16154_SE, ID_16154_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16167.png"); ggplot(df_temp[[37]], aes(ID_16167_SE, ID_16167_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16187.png"); ggplot(df_temp[[38]], aes(ID_16187_SE, ID_16187_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16648.png"); ggplot(df_temp[[39]], aes(ID_16648_SE, ID_16648_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_16649.png"); ggplot(df_temp[[40]], aes(ID_16649_SE, ID_16649_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()
png("correlation_SE_PE_ID_8546.png"); ggplot(df_temp[[41]], aes(ID_8546_SE, ID_8546_PE)) + geom_point() + geom_abline(colour = "brown"); dev.off()

