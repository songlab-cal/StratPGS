#############################################
########## Perturbation Experiment ##########
########## Using GWAS cutoff 1e-5 ###########
####### Sign flips and eSize Shuffles #######
#############################################
# This script compares sensitivity to shuffles vs
# sensitivity to signflips. 

## Log file / libraries --------------------------------------------------------
sink('/illumina/scratch/deep_learning/aaw/051723/logs/gwas1e-5_cutoff1e-6_compare_perturb_sensitivity.log')
sink(stdout(), type = "message")
library(data.table)
library(dplyr)

## Variables -------------------------------------------------------------------
R.workbench <- FALSE
gwas.thres <- '1e-5'
sig.cutoff <- '1e-6'
message('gwas.thres = ', gwas.thres)
message('sig.cutoff = ', sig.cutoff)

if (R.workbench) {
  out.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/comparison/'
  perturb.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
} else {
  out.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/comparison/'
  perturb.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
}

## Main Body -------------------------------------------------------------------
# Read perturbation sensitivity metrics
prs_shuffle_max <- readr::read_csv(paste0(perturb.dir,'gwas',
                                          gwas.thres,'_cutoff',
                                          sig.cutoff,'_shuffle_metrics_df.csv'),
                                   show_col_types = FALSE) %>% 
  as.data.frame()
prs_signflip_max <- readr::read_csv(paste0(perturb.dir,'gwas',
                                           gwas.thres,'_cutoff',
                                           sig.cutoff,'_signflip_metrics_df.csv'),
                                    show_col_types = FALSE) %>% 
  as.data.frame()
prs_shuffle_usual <- readr::read_csv(paste0(perturb.dir,'gwas',
                                          gwas.thres,'_cutoff',
                                          sig.cutoff,'_shuffle_usual_metrics_df.csv'),
                                     show_col_types = FALSE) %>% 
  as.data.frame()
prs_signflip_usual <- readr::read_csv(paste0(perturb.dir,'gwas',
                                           gwas.thres,'_cutoff',
                                           sig.cutoff,'_signflip_usual_metrics_df.csv'),
                                      show_col_types = FALSE) %>% 
  as.data.frame()

# Combine 
metric.vector <- colnames(prs_shuffle_max)[-1]
n.tests <- length(metric.vector)
message(date(), ": No. tests performed = ", n.tests)
fwer.thres.cutoff <- 0.05/n.tests
message(date(), ": FWER corrected p-value threshold = ", fwer.thres.cutoff)
metrics.df <- data.frame(METRIC=character(),
                         SHUFFLE_GREATER_SIGNFLIP_PVAL=character(),
                         SIGNFLIP_GREATER_SHUFFLE_PVAL=character())
for (metric in metric.vector) {
  message(date(), ": Comparing ", metric)
  max_shuffle_counts <- prs_shuffle_max %>% 
    select(c('PHENO',metric)) %>% 
    `colnames<-`(c('Pheno','Shuffle'))
  max_signflip_counts <- prs_signflip_max %>% 
    select(c('PHENO',metric)) %>% 
    `colnames<-`(c('Pheno','Signflip'))
  max_combined_df <- merge(max_shuffle_counts, max_signflip_counts,by='Pheno')
  max_shuffle_beat_signflip <- wilcox.test(x=max_combined_df$Shuffle,
                                           y=max_combined_df$Signflip,
                                           paired=TRUE,
                                           alternative='greater')$p.value
  max_signflip_beat_shuffle <- wilcox.test(x=max_combined_df$Shuffle,
                                           y=max_combined_df$Signflip,
                                           paired=TRUE,
                                           alternative='less')$p.value
  usual_shuffle_counts <- prs_shuffle_usual %>% 
    select(c('PHENO',metric)) %>% 
    `colnames<-`(c('Pheno','Shuffle'))
  usual_signflip_counts <- prs_signflip_usual %>% 
    select(c('PHENO',metric)) %>% 
    `colnames<-`(c('Pheno','Signflip'))
  usual_combined_df <- merge(usual_shuffle_counts, usual_signflip_counts,by='Pheno')
  usual_shuffle_beat_signflip <- wilcox.test(x=usual_combined_df$Shuffle,
                                             y=usual_combined_df$Signflip,
                                             paired=TRUE,
                                             alternative='greater')$p.value
  usual_signflip_beat_shuffle <- wilcox.test(x=usual_combined_df$Shuffle,
                                             y=usual_combined_df$Signflip,
                                             paired=TRUE,
                                             alternative='less')$p.value
  
  if (max_shuffle_beat_signflip < fwer.thres.cutoff) {
    message('MAX: Phenotypes more sensitive to Shuffle than to Signflip (one-sided p-value = ', 
            max_shuffle_beat_signflip, ')')
  }
  if (max_signflip_beat_shuffle < fwer.thres.cutoff) {
    message('MAX: Phenotypes more sensitive to Signflip than to Shuffle (one-sided p-value = ', 
            max_signflip_beat_shuffle, ')')
  }
  
  if (usual_shuffle_beat_signflip < fwer.thres.cutoff) {
    message('USUAL: Phenotypes more sensitive to Shuffle than to Signflip (one-sided p-value = ', 
            usual_shuffle_beat_signflip, ')')
  }
  if (usual_signflip_beat_shuffle < fwer.thres.cutoff) {
    message('USUAL: Phenotypes more sensitive to Signflip than to Shuffle (one-sided p-value = ', 
            usual_signflip_beat_shuffle, ')')
  }
  
  metrics.df <- rbind(metrics.df,
                      data.frame(METRIC=metric,
                                 MAX_SHUFFLE_GREATER_SIGNFLIP_PVAL=max_shuffle_beat_signflip,
                                 MAX_SIGNFLIP_GREATER_SHUFFLE_PVAL=max_signflip_beat_shuffle,
                                 USUAL_SHUFFLE_GREATER_SIGNFLIP_PVAL=usual_shuffle_beat_signflip,
                                 USUAL_SIGNFLIP_GREATER_SHUFFLE_PVAL=usual_signflip_beat_shuffle))
  
}
readr::write_csv(metrics.df,
                 file = paste0(out.dir,
                               "gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_signflip_comparisons.csv"))
sink()

## Plotting --------------------------------------------------------------------