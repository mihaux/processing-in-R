# processing-in-R

LAST UPDATED: 2020-07-09

# TODO: no need to place usage and example here; move them to each script (at the beginning) and here just brief description about each script

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
```

featureCounts_merge.R
get_PicardMetrics.sh
get_TrimStats_PE.sh
get_TrimStats_SE.sh
report_p1.Rmd
report_p2.Rmd
runPCA_DESeq2.R
run_Normalise_TMM.R
run_Process_metadata.R
test_Mann-Whitney.R
test_confounding.R
