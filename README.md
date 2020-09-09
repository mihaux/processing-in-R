# processing-in-R

LAST UPDATED: 2020-07-09

# TODO: no need to place usage and example here; move them to each script (at the beginning) and here just brief description about each script

# current list of files:
# => tutorials
# PCA_Tutorial.R


# compare_norm_vs_vst.R		=>

# featureCounts_merge.R		=>

# get_FastQCreport.R		=>

# get_Histogram_for_features.R	=> TO 

# get_Nb_of_cases_per_group.R	=> 

# get_PicardMetrics.sh		=> TO BE CHECKED

# get_STARstats.R		=> TO BE CHEKED

# get_TranscriptByChromosome.R	=> TO BE CHECKED

# get_TrimStats_PE.sh		=> TO BE MERGE WITH _SE

# get_TrimStats_SE.sh		=> TO BE MARGED WITH _PE

# get_Venn_diagrams.R		=> NOT WORKING, TO BE CHECKED

# get_reads-distribution.R	=> TO BE CHECKED

# run_DESeq2_Ian.R		=> TO BE CHEKED, THEN REMOVED

# run_DESeq2_Michal.R		=> TO BE CHECKED

# run_DESeq2_plots.R		=> TO BE CHECKED

# run_DESeq2_split_chr.R	=> TO BE CHECKED AND POSSIBLY REMOVED IF NOT NEEDED

# run_DE_analysis_TAI.R		=> TO BE REMOVED

# run_Distribution_plot.R	=> TO BE CHECKED

# run_Filter_NormaliseTMM.R	=> TO BE CHECKED AND POSSIBLY REMOVED IF NOT NEEDED ANYMORE

# run_Hierarchical_clust.R	=> TO BE COMPLETED (empty file)

# run_Linear_Regression.R	=> TO BE COMPLETED (empty file)

# run_Modify_colnames.R		=> TO BE CHECKED

# run_PCA_DESeq2.R		=> PRACTICALLY FINISHED (to be cjecked)

# run_PCA_built-in.R		=> TO BE COMPLETED

# run_PCA_built-in_another.R	=> TO BE CHECKED AND POSSIBLY REMOVED IF NOT NEEDED

# run_Process_metadata.R	=> to be chcecked

# run_Switch_ENTREZID-SYMBOL.R	=> TO BE CHECKED IF IT'S WORKING CORRECTLY (possibly to be extended for other types of annotation as well)

# test_Fisher_s-exact.R		=> TO BE COMPLETED

# test_Mann-Whitney.R		=> PRACTICALLY FINSIHED (to be checked)

# test_confounding.R		=> TO BE COMPLETED


## 1) get_FastQCreport.R => 

USAGE: `Rscript get_FastQCreport.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript get_FastQCreport.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/report /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/postprocessed
```

## 2) get_STARstats.R 

USAGE: `Rscript get_STARstats.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript get_STARstats.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_I_Nov19/4_alignement
```

## 3) featureCounts_merge.R

USAGE: `Rscript featureCounts_merge.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript featureCounts_merge.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I
```

## 4) runPCA_built-in.R and runPCA_built-in_another.R

USAGE: `Rscript runPCA_built-in.R [counts.csv] [annotation.csv] [directory-for-output]`

EXAMPLE: 
```
Rscript runPCA_built-in.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end
```

## 5) run_Switch_ENTREZID-SYMBOL.R

USAGE: `Rscript run_Switch_ENTREZID-SYMBOL.R [counts.csv] [directory-for-output]`

EXAMPLE:  
**for mode_I**
```
Rscript run_Switch_ENTREZID-SYMBOL.R mode_I /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_I/counts_merged.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_I/
```
and
```
Rscript run_Switch_ENTREZID-SYMBOL.R mode_I /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I/counts_merged.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I/
```
**for mode_II**
```
Rscript changeAnno.R mode_II /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/
```
and
```
Rscript changeAnno.R mode_II /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/all_counts_PE.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/
