# make PCA plots with and without chrX and chrY

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
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

# define wheather output files should be saved or not [TRUE / FALSE]
output_save <- FALSE

# define directory with data (INPUT)
args <- c(paste0(main_dir,"/data/count_matrices/outliers_excluded/counts_vst_no_outliers.csv"))

run_id <- "vst"

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES/Dec20_pca/")
setwd(dir_out)

# define directory with metadata
dir_metadata <- paste0(main_dir, "/data/metadata/outliers_excluded/cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv")

# define directory for chromosome info
dir_chromosome <- paste0(main_dir, "/ANALYSES/Dec20_pca/chr_info/")

# load data RAW | VST | rlog (run for one data type at a time)
dds <- read.csv(args[1], row.names = 1, header = TRUE)
dat <- as.matrix(dds)

# load clinical data 
df_meta <- read.csv(dir_metadata, row.names = 1, header = TRUE)

# load chromosome info
# source: https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart&attributes=hgnc_gene__hgnc_gene_id_1010%2Chgnc_gene__approved_symbol_1010%2Chgnc_gene__chromosome_1010%2Chgnc_gene__chromosome_location_1010&hgnc_gene__chromosome_1010=mitochondria%2Cc10_B%2Creserved%2CNA
#chrX <- read.csv(paste0(dir_chromosome, "results_X.csv"), header = TRUE, sep = "\t")
#chrY <- read.csv(paste0(dir_chromosome, "results_Y.csv"), header = TRUE, sep = "\t")
#chr_XY <- read.csv(paste0(dir_chromosome, "results_X_and_Y.csv"), header = TRUE, sep = "\t")
#chr_all <- read.csv(paste0(dir_chromosome, "results_all.csv"), header = TRUE, sep = "\t")
#chr_1_22 <- read.csv(paste0(dir_chromosome, "results_1-22.csv"), header = TRUE, sep = "\t")
#chr_other <- read.csv(paste0(dir_chromosome, "results_other.csv"), header = TRUE, sep = "\t")

# find intersections
#length(intersect(rownames(dat), chrX$Approved.symbol))          # => 913, out of 2016
#length(intersect(rownames(dat), chrY$Approved.symbol))          # => 70, out of 502
#length(intersect(rownames(dat), chr_1_22$Approved.symbol))      # => 21 740, out of 41 054
#length(intersect(rownames(dat), chr_all$Approved.symbol))       # => 22 747, out of 43 612
#length(intersect(rownames(dat), chr_XY$Approved.symbol))        # => 24, out of 40

if(FALSE){
chr_anno <- read.csv(paste0(main_dir, "/ANALYSES/run_12_Aug20/5_counting/PE_all/all_annotations_dups_PE_all.csv"), sep=",")

dim(dat) # => 26 486

chr_info_all <- str_split_fixed(chr_anno$Chr, ";", n=2)
chr_info_first_col <- chr_info_all[,1]

gene_names <- chr_anno$GeneID

# create a dataframe with gene_names and chr_info
df_chr <- data.frame(genes=gene_names,
                     chr=chr_info_first_col)

df_chr_final <- df_chr[which(gene_names %in% rownames(dat)),]

write.csv(df_chr_final, "chromosome_info.csv")
}

df_chr_final <- read.csv(paste0(dir_out, "chromosome_info.csv"))
# length(intersect(rownames(dat), gene_names))  # => all of them are included

# extract gene names for chrX and chrY
idx_chrX <- which(df_chr_final$chr == "chrX")
idx_chrY <- which(df_chr_final$chr == "chrY")

# chr 1 - 22
dat_chr1_22 <- dat[-c(idx_chrX, idx_chrY),]
  
# chr 1 - 22 + chrX
dat_chr1_22_X <- dat[-idx_chrY,]
  
# chr 1 - 22 + chrY
dat_chr1_22_Y <- dat[-idx_chrX,]

# create DESeqDataSetFromMatrix (in order to use plotPCA() function)
gender <- c(rep("female", 40))
gender[which(df_meta$gender == 1)] <- "male"

extra <- c(rep("example_1", 20), rep("example_1", 20))
coldata_1 <- cbind(gender, extra)
                   
dds_all_chr <- DESeqDataSetFromMatrix(countData = round(dat), colData = coldata_1, design = ~ gender)
tr_all_chr <- DESeqTransform(dds_all_chr)

dds_chr1_22 <- DESeqDataSetFromMatrix(countData = round(dat_chr1_22), colData = coldata_1, design = ~ gender)
tr_chr1_22 <- DESeqTransform(dds_chr1_22)

dds_chr1_22_X <- DESeqDataSetFromMatrix(countData = round(dat_chr1_22_X), colData = coldata_1, design = ~ gender)
tr_chr1_22_X <- DESeqTransform(dds_chr1_22_X)

dds_chr1_22_Y <- DESeqDataSetFromMatrix(countData = round(dat_chr1_22_Y), colData = coldata_1, design = ~ gender)
tr_chr1_22_Y <- DESeqTransform(dds_chr1_22_Y)

# make pca plots using DESeq2 package command

# all chromosomes
if(output_save==TRUE){ png(file = paste0("PCA_plot_all_chr_", run_id, ".png")) }
plotPCA(tr_all_chr, intgroup=c("gender")) + geom_point(size = 3)
if(output_save==TRUE){ dev.off() }
# plotPCA(tr_all_chr, intgroup=c("gender")) + geom_text(aes(label=colnames(dds_all_chr)))

# chr 1 - 22
if(output_save==TRUE){ png(file = paste0("PCA_plot_chr_1-22_", run_id, ".png")) }
plotPCA(tr_chr1_22, intgroup=c("gender")) + geom_point(size = 3)
if(output_save==TRUE){ dev.off() }
#plotPCA(tr_chr1_22, intgroup=c("gender")) + geom_text(aes(label=colnames(dds_chr1_22)))

# chr 1 - 22 + chrX
if(output_save==TRUE){ png(file = paste0("PCA_plot_chr_1-22-X_", run_id, ".png")) }
plotPCA(tr_chr1_22_X, intgroup=c("gender")) + geom_point(size = 3)
if(output_save==TRUE){ dev.off() }
#plotPCA(tr_chr1_22_X, intgroup=c("gender")) + geom_text(aes(label=colnames(dds_chr1_22_X)))

# chr 1 - 22 + chrY
if(output_save==TRUE){ png(file = paste0("PCA_plot_chr_1-22-Y_", run_id, ".png")) }
plotPCA(tr_chr1_22_Y, intgroup=c("gender")) + geom_point(size = 3)
if(output_save==TRUE){ dev.off() }
#plotPCA(tr_chr1_22_Y, intgroup=c("gender")) + geom_text(aes(label=colnames(dds_chr1_22_Y)))


##########################################################################################
# use spearman coefficient to compute correlation between
########################################################################################## 

# create vectors with PC1 for all three sets
chr_all <- plotPCA(tr_all_chr, intgroup=c("gender"))
chr_1_22 <- plotPCA(tr_chr1_22, intgroup=c("gender"))
chr_1_22_X <- plotPCA(tr_chr1_22_X, intgroup=c("gender"))
chr_1_22_Y <- plotPCA(tr_chr1_22_Y, intgroup=c("gender"))

# => PC2 for all_chromosomes vs. PC2 for chr1-22
res_spearman_all_1_22 <- cor.test(chr_all$data$PC2, chr_1_22$data$PC2, alternative = c("two.sided"), method = c("spearman"), conf.level = 0.95, exact=FALSE)

# => PC2 for all_chromosomes vs. PC2 for chr1-22+chrX
res_spearman_all_1_22_X <- cor.test(chr_all$data$PC2, chr_1_22_X$data$PC2, alternative = c("two.sided"), method = c("spearman"), conf.level = 0.95, exact=FALSE)

# => PC2 for all_chromosomes vs. PC2 for chr1-22+chrY
res_spearman_all_1_22_Y <- cor.test(chr_all$data$PC2, chr_1_22_Y$data$PC2, alternative = c("two.sided"), method = c("spearman"), conf.level = 0.95, exact=FALSE)


# Multiple testing correction
p_adj_pearson <- p.adjust(unlist(lapply(res_pearson, function(x) x$p.value)), "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")
p_adj_spearman <- p.adjust(unlist(lapply(res_spearman, function(x) x$p.value)), "fdr")    # Benjamini & Hochberg ("BH" or its alias "fdr")

# create table with results
res_table_pearson <- data.frame(ID=rownames(dat),
                                p.value=unlist(lapply(res_pearson, function(x) x$p.value)),
                                p.adj=p_adj_pearson,
                                pearson_coef=unlist(lapply(res_pearson, function(x) x$estimate)),
                                CI_1=unlist(lapply(res_pearson, function(x) x$conf.int[1])),
                                CI_2=unlist(lapply(res_pearson, function(x) x$conf.int[2])))


res_table_spearman <- data.frame(ID=rownames(dat),
                                 p.value=unlist(lapply(res_spearman, function(x) x$p.value)),
                                 p.adj=p_adj_spearman,
                                 spearman_coef=unlist(lapply(res_spearman, function(x) x$estimate)))

# write a summary table with numbers of significant results
summary_significant = data.frame(V1=run_id, 
                                 V2=length(which(res_table_pearson$p.value < 0.05)),
                                 V3=length(which(res_table_pearson$p.adj < 0.05)),
                                 V4=length(which(res_table_spearman$p.value < 0.05)),
                                 V5=length(which(res_table_spearman$p.adj < 0.05)))

colnames(summary_significant) <- c("", "p-value < 0.05 (Pearson)", "p-adjusted < 0.05 (Pearson)",
                                   "p-value < 0.05 (Spearman)", "p-adjusted < 0.05 (Spearman)")

if(output_save==TRUE){  
  write.csv2(summary_significant, file=paste0("table_summary_significant_", run_id, ".csv"))
  write.csv2(res_table_pearson, file=paste0("table_pearson_", run_id, ".csv"))
  write.csv2(res_table_spearman, file=paste0("table_spearman_", run_id, ".csv"))
}





# run PCA using the built-in function
#pca_dds <- prcomp(t(assay(dds)))
pca_dds <- prcomp(t(dat))

# retrieve PCA loadings
loadings_dds <- as.data.frame(pca_dds$rotation)
# dim(loadings_dds)   => [1] 26486    41
# loadings_dds[1:5,1:5]
#                   PC1           PC2           PC3           PC4          PC5
#DDX11L1   -0.011953693 -0.0053474555 -0.0045972092 -0.0001899563  0.004498413
#WASH7P     0.001688707  0.0005548291 -0.0062353401 -0.0132954071  0.007916622
#MIR6859-1  0.002086512  0.0003553188 -0.0003207413 -0.0070682739  0.002156764
#MIR1302-2 -0.012560021 -0.0092091513  0.0035144072  0.0013492019 -0.005964795
#FAM138A   -0.017429927 -0.0118937876  0.0068809726  0.0042712606 -0.007199880

# save transcript IDs 
sink("rownames.txt")
for (i in 1:length(rownames(loadings_dds))) {
  cat(rownames(loadings_dds)[i], "\n")
}
sink()

# load chromosome information
chr_info <- read.csv(paste0(dir_out, "mart_export.csv"))

# get only needed information: chr | position (start) | ID
chr_info_final <- chr_info[,c(5,6,8)]

library(stringr)
# get PATCHes and exclud them
patches_id <- which(str_detect(chr_info_final$Chromosome.scaffold.name, "PATCH"))
chr_info_final_no_patches <- chr_info_final[-c(patches_id),]

# remove other unneccesary stuff:
# CHR_HG151_NOVEL_TEST 
to_be_ex_1 <- which(chr_info_final_no_patches$Chromosome.scaffold.name == "CHR_HG151_NOVEL_TEST")

# CHR_HG142_HG150_NOVEL_TEST 
to_be_ex_2 <- which(chr_info_final_no_patches$Chromosome.scaffold.name == "CHR_HG142_HG150_NOVEL_TEST")

# GL000213.1 
to_be_ex_3 <- which(chr_info_final_no_patches$Chromosome.scaffold.name == "GL000213.1")

to_be_ex_123 <- c(to_be_ex_1, to_be_ex_2, to_be_ex_3)

chr_info_final_no_patches_no_other <- chr_info_final_no_patches[-c(to_be_ex_123),]

# NOTE: some IDs are missing, REMOVE THEM
to_be_ex_id <- which(chr_info_final_no_patches_no_other$RefSeq.match.transcript == "")

chr_info_final_END <- chr_info_final_no_patches_no_other[-c(to_be_ex_id),]

write.csv(chr_info_final_END, "chr_info_final_END.csv")

chr_info_final_END_mod <- read.csv(paste0(dir_out, "chr_info_final_END_mod.csv"), row.names = 1)

#write.csv(chr_info_final_no_patches_no_other, "chr_pos.csv")
#kot <- unique(chr_info_final_no_patches_no_other$Chromosome.scaffold.name)

#chr_info_final_no_patches_no_other$Chromosome.scaffold.name[which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR1_"))]
#vec_chr_1 <- c(rep(1, 109))

# replace with chromosome names
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR1_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR2_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR3_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR4_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR5_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR6_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR7_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR8_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR9_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR10"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR11_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR12_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR13_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR14_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR15_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR16_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR17_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR18_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR19"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR20_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR21_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHR22_"))
which(str_detect(chr_info_final_no_patches_no_other$Chromosome.scaffold.name, "CHR_HSCHRX_"))


library(lattice)

manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab="PC1 loadings",
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.1,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(pvalue,thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=sig.level, lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

#FAKE SAMPLE DATA
createSampleGwasData<-function(chr.count=10, include.X=F) {
  chr<-c(); pos<-c()
  for(i in 1:chr.count) {
    chr <- c(chr,rep(i, 1000))
    pos <- c(pos,ceiling(runif(1000)*(chr.count-i+1)*25*1e3))
  }
  if(include.X) {
    chr <- c(chr,rep("X", 1000))
    pos <- c(pos,ceiling(runif(1000)*5*25*1e3))
  }
  pvalue <- runif(length(pos))
  return(data.frame(chr, pos,pvalue))
}

dd<-createSampleGwasData()
dd$pvalue[3000] <- 1e-7 #include a significant result

manhattan.plot(dd$chr, dd$pos, dd$pvalue)


# create a dataframe with loadings for PC1
load_PC1 <- data.frame(ID=rownames(loadings_dds),
                       PC1=loadings_dds[,1])
  
# match load_PC1 with chr_info_final_END_mod
load_PC1_sub <- load_PC1[which(load_PC1$ID %in% chr_info_final_END_mod$RefSeq.match.transcript),]

load_PC1_sub_ordered <- load_PC1_sub[order(load_PC1_sub$ID),]

chr_info_final_END_mod_sub <- chr_info_final_END_mod[unique(chr_info_final_END_mod$RefSeq.match.transcript),]

write.csv(chr_info_final_END_mod_sub, "chr_info_final_END_mod_sub.csv", row.names = FALSE)

chr_info_final_END_mod_sub_mod <- read.csv(paste0(dir_out, "chr_info_final_END_mod_sub_mod.csv"))

id_t_or_f <- chr_info_final_END_mod_sub_mod$RefSeq.match.transcript %in% load_PC1_sub_ordered$ID

# it's still too long, exclude these rows whicha re not included in load_PC1_sub
setdiff(chr_info_final_END_mod_sub$RefSeq.match.transcript, load_PC1_sub_ordered$ID)

# create final df for manhattan plot
df_manhattan <- data.frame(chr=chr_info_final_END_mod_sub_mod$Chromosome.scaffold.name,
                           pos=chr_info_final_END_mod_sub_mod$Gene.start..bp.,
                           pvalue=chr_info_final_END_mod_sub_mod$RefSeq.match.transcript)


load_PC1_sub[1:15071,]

df_manhattan_FINAL <- df_manhattan[order(df_manhattan$chr),]

manhattan.plot(c(df_manhattan_FINAL$chr), c(df_manhattan_FINAL$pos), load_PC1_sub[1:15071,2])






plotPCA.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
}

# plotPCA(dds_tr, intgroup=c("gender"))    => object: dds_tr (or dds, depands what input data needed)

# retrieve variance for each gene
rv <- rowVars(assay(dds))

# create a dataframe with rownames and rv values
df_rv_new <- as.data.frame(cbind(rownames(dds), rv))

newdata <- df_rv_new[order(-rv),] 
newdata_bis <- newdata[1:500,]

# load positions on chromosmes for each gene
tbch <- read.csv2(paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/table_TranscriptsByChromosome_modified.csv"), sep = ",", row.names = 1)

# match genes from the vector of 500 most variable with their chromosome names 
chr_pos_temp <- list()
for (i in 1:500){
  chr_pos_temp[[i]] <- unique(tbch[which(tbch$IDs == as.vector(newdata_bis$V1[i])),])[1,1]
}

# NOTE: can also try for the whole datadet, but it might take long
chr_pos <- as.vector(unlist(chr_pos_temp))

# add chr_pos to newdata_bis
newdata_tris <- data.frame(ID=newdata_bis$V1,
                           gene_var=newdata_bis$rv,
                           gene_posistion=chr_pos,
                           PC1=loadings_dds$PC1[newdata_bis$V1],
                           PC2=loadings_dds$PC2[newdata_bis$V1])



barplot(as.numeric(as.vector(newdata_bis$rv)), main = "500 top most variable genes", names.arg=chr_pos, las=2, xlab="Chromomsomes where the gene is located")

barplot(as.numeric(as.vector(newdata_tris$PC1)[1:20]), main = "500 top most variable genes", names.arg=chr_pos[1:20], las=2, xlab="Chromomsomes where the gene is located")

barplot(as.numeric(as.vector(newdata_tris$PC2)[1:20]), main = "500 top most variable genes", names.arg=chr_pos[1:20], las=2, xlab="Chromomsomes where the gene is located")


# make a scater scatterplot with either
# (a) the variance values of each gene or 
# (b) use the PC loadings to see how the values are spread, in regard to the chromosomes


chr_1 <- which(newdata_tris$gene_posistion == "chr1")

barplot(newdata_tris$gene_var[chr_1])

# the genes need to be ordered by chromosome
plot(newdata_tris$gene_posistion[1:10], newdata_tris$gene_var[1:10])

library(ggplot2) 
ggplot(newdata_tris) + geom_histogram(aes(x=gene_var),binwidth=1e6) + facet_grid(gene_posistion ~ID)

# You need a data frame with colums Chr, Start and Sample. If you just want to see points, change geom_histogram(...) for geom_point(aes(x=start, y=0))


# unlist and get how man per each chromosomes










###########################################################################################
# => run for different numbers of retained genes (TO BE MODIFIED)
###########################################################################################

setwd(paste0(dir_out,"/PCA_plots/", p1_gender$labels$colour))

p1_gender_100 <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=100) + ggtitle("raw: all chr | ntop=100") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_100$labels$colour <- "gender"

p1_gender_250 <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=250) + ggtitle("raw: all chr | ntop=250") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_250$labels$colour <- "gender"

p1_gender_500 <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=500) + ggtitle("raw: all chr | ntop=500") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_500$labels$colour <- "gender"

p1_gender_1000 <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=1000) + ggtitle("raw: all chr | ntop=1 000") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_1000$labels$colour <- "gender"

p1_gender_10000 <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=10000) + ggtitle("raw: all chr | ntop=10 000") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_10000$labels$colour <- "gender"

p1_gender_all <- plotPCA(dds_all_trans, intgroup=c("gender"), ntop=26486) + ggtitle("raw: all chr | ntop=26 486") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p1_gender_all$labels$colour <- "gender"

# save plots
png(file=paste0(p1_gender_100$labels$colour, "_raw_all_chr_100.png")); p1_gender_100; dev.off()
png(file=paste0(p1_gender_250$labels$colour, "_raw_all_chr_250.png")); p1_gender_250; dev.off()
png(file=paste0(p1_gender_500$labels$colour, "_raw_all_chr_500.png")); p1_gender_500; dev.off()
png(file=paste0(p1_gender_1000$labels$colour, "_raw_all_chr_1000.png")); p1_gender_1000; dev.off()
png(file=paste0(p1_gender_10000$labels$colour, "_raw_all_chr_10000.png")); p1_gender_10000; dev.off()
png(file=paste0(p1_gender_all$labels$colour, "_raw_all_chr_26486.png")); p1_gender_all; dev.off()


p2_gender_100 <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=100) + ggtitle("vst: all chr | ntop=100") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_100$labels$colour <- "gender"

p2_gender_250 <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=250) + ggtitle("vst: all chr | ntop=250") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_250$labels$colour <- "gender"

p2_gender_500 <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=500) + ggtitle("vst: all chr | ntop=500") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_500$labels$colour <- "gender"

p2_gender_1000 <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=1000) + ggtitle("vst: all chr | ntop=1 000") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_1000$labels$colour <- "gender"

p2_gender_10000 <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=10000) + ggtitle("vst: all chr | ntop=10 000") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_10000$labels$colour <- "gender"

p2_gender_all <- plotPCA(vst_all_trans, intgroup=c("gender"), ntop=26486) + ggtitle("vst: all chr | ntop=26 486") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1)
p2_gender_all$labels$colour <- "gender"

png(file=paste0(p2_gender_100$labels$colour, "_VST_all_chr_100.png")); p2_gender_100; dev.off()
png(file=paste0(p2_gender_250$labels$colour, "_VST_all_chr_250.png")); p2_gender_250; dev.off()
png(file=paste0(p2_gender_500$labels$colour, "_VST_all_chr_500.png")); p2_gender_500; dev.off()
png(file=paste0(p2_gender_1000$labels$colour, "_VST_all_chr_1000.png")); p2_gender_1000; dev.off()
png(file=paste0(p2_gender_10000$labels$colour, "_VST_all_chr_10000.png")); p2_gender_10000; dev.off()
png(file=paste0(p2_gender_all$labels$colour, "_VST_all_chr_26486.png")); p2_gender_all; dev.off()



# source: 
# => https://www.biostars.org/p/237730/
# => https://support.bioconductor.org/p/51270/ 

# run PCA using the built-in function
pca_dds <- prcomp(t(assay(dds_all_trans)))
pca_rlog <- prcomp(t(assay(rlog_all_trans)))
pca_vst <- prcomp(t(assay(vst_all_trans)))

# get PCA loadings
loadings_dds <- as.data.frame(pca_dds$rotation)
loadings_rlog <- as.data.frame(pca_rlog$rotation)
loadings_vst <- as.data.frame(pca_vst$rotation)

# load modified subset of gtf file

gtf_new <- read.csv(file = paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/", "table_TranscriptsByChromosome_modified.csv"), row.names = 1)

# => need to look at the difference between chrX and chrY (check thair contribution in PCA loadings)

dim(gtf_new)                  # [1] 3968352       2
length(unique(gtf_new$IDs))

# TODO: match rownames(loadings_dds) [l=26 486] with unique(gtf_new$IDs) [l=38 214]


kot <- c()
for (i in 1:length(rownames(loadings_dds))) {
  print(which(rownames(loadings_dds)[i] == unique(gtf_new$IDs)))
  kot[[i]] <- which(rownames(loadings_dds)[i] == unique(gtf_new$IDs))
}

# sort obtained vector and check if it has the same order as unsorted
kot_sorted <- sort(kot)

any(kot_sorted == kot)


my_chr_fin_1 <- gsub("chr", "", my_chr)

# need to replace "chrM"  "chrX" and "chrY" with numerics: 23, 24, 25
my_chr_fin_2 <- gsub("M", "23", my_chr_fin_1)
my_chr_fin_3 <- gsub("X", "24", my_chr_fin_2)
my_chr_fin_4 <- gsub("Y", "25", my_chr_fin_3)
my_chr_final <- as.numeric(as.vector(my_chr_fin_4))


mygwas_list <- list(my_snp, as.integer(my_chr_final), as.integer(my_bp), as.numeric(my_p1[-unnecessary_indices])*(-1))
names(mygwas_list) <- c("SNP", "CHR", "BP", "P")

mygwas_list_bis <- list(as.character(as.vector(pc1_ordered$ID)), as.integer(my_chr_final), as.integer(my_bp), pc1_ordered$PC1*(-1))

#mygwas <- cbind(as.character(my_snp), as.vector(my_chr_final), as.integer(my_bp), as.numeric(my_p1[-unnecessary_indices]))
mygwas_fin <- as.data.frame(mygwas_list)
mygwas_fin_bis <- as.data.frame(mygwas_list_bis)

names(mygwas_fin) <- c("SNP", "CHR", "BP", "P")
names(mygwas_fin_bis) <- c("SNP", "CHR", "BP", "P")

manhattan(mygwas_fin, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylim = c(-0.05, 0.45))

manhattan(mygwas_fin_bis, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylim = c(-0.05, 0.45))

manhattan(mygwas_fin_bis, main = "Manhattan Plot - PC1", ylim = c(-0.05, 0.45), cex = 0.6, logp=FALSE,
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:22, "M", "X", "Y"))

# => what exactly they are 









####################################################################
####################################################################
####################################################################




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

