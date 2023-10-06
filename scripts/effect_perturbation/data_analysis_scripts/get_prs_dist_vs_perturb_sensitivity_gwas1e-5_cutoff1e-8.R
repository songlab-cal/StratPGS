#############################################
########## Perturbation Experiment ##########
########## Using GWAS cutoff 1e-5 ###########
####### Sign flips and eSize Shuffles #######
#############################################
# This script combines PRS distributional metrics (beta sizes, etc.)
# and sensitivity of PRS to perturbations on a training set
# Created on 6/2/23: This version takes max(f(y,PRS^+),f(y,PRS^-))
#                     Also computes the f(y,PRS) version (aka 'usual approach')
# objects with 'usual' appended in front correspond to the f(y,PRS) version of things

## Log file / libraries --------------------------------------------------------
sink('/illumina/scratch/deep_learning/aaw/051723/logs/gwas1e-5_cutoff1e-8_prs_dist_vs_perturb_sensitivity.log')
sink(stdout(), type = "message")
library(data.table)
library(dplyr)
#library(bigsnpr)

## Variables -------------------------------------------------------------------
R.workbench <- FALSE
gwas.thres <- '1e-5'
sig.cutoff <- '1e-8'
message('gwas.thres = ', gwas.thres)
message('sig.cutoff = ', sig.cutoff)

if (R.workbench) {
  phenos.shortlist <- data.table::fread('/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
                                        header = FALSE)$V1
  out.dir <- '/deep_learning/aaw/051723/results/prs_dist_vs_perturbation_sensitivity/'
  prs.dist.dir <- '/deep_learning/aaw/051723/results/prs_dist_summaries/'
  #prs.strat.dir <- '/deep_learning/aaw/051723/results/prs_pc_stratification/'
  #pheno.strat.dir <- '/deep_learning/aaw/051723/results/pheno_pc_stratification/'
  perturb.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
} else {
  phenos.shortlist <- data.table::fread('/illumina/scratch/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
                                        header = FALSE)$V1
  out.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_dist_vs_perturbation_sensitivity/'
  prs.dist.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_dist_summaries/'
  #prs.strat.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_pc_stratification/'
  #pheno.strat.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/pheno_pc_stratification/'
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

# Read PRS distribution metrics
prs_dist <- readr::read_csv(paste0(prs.dist.dir,
                                    'prs_dist_summaries_gwas',
                                    gwas.thres,
                                    '_cutoff',
                                    as.numeric(sig.cutoff),
                                    '.csv'),
                             show_col_types = FALSE) %>% 
  as.data.frame()
colnames(prs_dist)[1] <- 'PHENO'
# Combine 
shuffle_max_prs_dist <- merge(prs_shuffle_max,prs_dist,by='PHENO')
signflip_max_prs_dist <- merge(prs_signflip_max,prs_dist,by='PHENO')
shuffle_usual_prs_dist <- merge(prs_shuffle_usual,prs_dist,by='PHENO')
signflip_usual_prs_dist <- merge(prs_signflip_usual,prs_dist,by='PHENO')

readr::write_csv(shuffle_max_prs_dist,
                 file = paste0(out.dir,"/raw_tables/gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_max_combined.csv"))
readr::write_csv(signflip_max_prs_dist,
                 file = paste0(out.dir,"/raw_tables/gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_max_combined.csv"))
readr::write_csv(shuffle_usual_prs_dist,
                 file = paste0(out.dir,"/raw_tables/gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_usual_combined.csv"))
readr::write_csv(signflip_usual_prs_dist,
                 file = paste0(out.dir,"/raw_tables/gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_usual_combined.csv"))

# Compute p-values and correlations
shuffle_max_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_shuffle_max)-1,ncol=ncol(prs_dist)-1)
signflip_max_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_signflip_max)-1,ncol=ncol(prs_dist)-1)
shuffle_usual_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_shuffle_usual)-1,ncol=ncol(prs_dist)-1)
signflip_usual_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_signflip_usual)-1,ncol=ncol(prs_dist)-1)

shuffle_max_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_shuffle_max)-1,ncol=ncol(prs_dist)-1)
signflip_max_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_signflip_max)-1,ncol=ncol(prs_dist)-1)
shuffle_usual_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_shuffle_usual)-1,ncol=ncol(prs_dist)-1)
signflip_usual_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_signflip_usual)-1,ncol=ncol(prs_dist)-1)

colnames(shuffle_max_vs_strat_corr_array) <- colnames(prs_dist)[-1]
rownames(shuffle_max_vs_strat_corr_array) <- colnames(prs_shuffle_max)[-1]
colnames(signflip_max_vs_strat_corr_array) <- colnames(prs_dist)[-1]
rownames(signflip_max_vs_strat_corr_array) <- colnames(prs_signflip_max)[-1]
colnames(shuffle_usual_vs_strat_corr_array) <- colnames(prs_dist)[-1]
rownames(shuffle_usual_vs_strat_corr_array) <- colnames(prs_shuffle_usual)[-1]
colnames(signflip_usual_vs_strat_corr_array) <- colnames(prs_dist)[-1]
rownames(signflip_usual_vs_strat_corr_array) <- colnames(prs_signflip_usual)[-1]

colnames(shuffle_max_vs_strat_pvals_array) <- colnames(prs_dist)[-1]
rownames(shuffle_max_vs_strat_pvals_array) <- colnames(prs_shuffle_max)[-1]
colnames(signflip_max_vs_strat_pvals_array) <- colnames(prs_dist)[-1]
rownames(signflip_max_vs_strat_pvals_array) <- colnames(prs_signflip_max)[-1]
colnames(shuffle_usual_vs_strat_pvals_array) <- colnames(prs_dist)[-1]
rownames(shuffle_usual_vs_strat_pvals_array) <- colnames(prs_shuffle_usual)[-1]
colnames(signflip_usual_vs_strat_pvals_array) <- colnames(prs_dist)[-1]
rownames(signflip_usual_vs_strat_pvals_array) <- colnames(prs_signflip_usual)[-1]

n.tests <- nrow(shuffle_max_vs_strat_corr_array)*ncol(shuffle_max_vs_strat_corr_array)
message(date(), ": No. tests performed = ", n.tests)
fwer.thres.cutoff <- 0.05/n.tests
message(date(), ": FWER corrected p-value threshold = ", fwer.thres.cutoff)
for (i in 1:nrow(shuffle_max_vs_strat_corr_array)) {
  for (j in 1:ncol(shuffle_max_vs_strat_corr_array)) {
    # Shuffle, max
    shuffle_max_vs_strat_corr_array[i,j] <- cor(shuffle_max_prs_dist[[rownames(shuffle_max_vs_strat_corr_array)[i]]],
                                                shuffle_max_prs_dist[[colnames(shuffle_max_vs_strat_corr_array)[j]]])
    p_val <- cor.test(shuffle_max_prs_dist[[rownames(shuffle_max_vs_strat_corr_array)[i]]],
                      shuffle_max_prs_dist[[colnames(shuffle_max_vs_strat_corr_array)[j]]])$p.value
    shuffle_max_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SHUFFLE, MAX) Constant value observed for ', rownames(shuffle_max_vs_strat_corr_array)[i])
      message('(SHUFFLE, MAX) Value = ', mean(shuffle_max_prs_dist[[rownames(shuffle_max_vs_strat_corr_array)[i]]]))
    } else if (p_val < fwer.thres.cutoff) {
      message('(SHUFFLE, MAX): Significant correlation (p = ',p_val,') between ', 
              rownames(shuffle_max_vs_strat_corr_array)[i], ' and ', 
              colnames(shuffle_max_vs_strat_corr_array)[j])
      message('(SHUFFLE, MAX): Pearson r = ', shuffle_max_vs_strat_corr_array[i,j])
    }
    
    # Signflip, max
    signflip_max_vs_strat_corr_array[i,j] <- cor(signflip_max_prs_dist[[rownames(signflip_max_vs_strat_corr_array)[i]]],
                                                signflip_max_prs_dist[[colnames(signflip_max_vs_strat_corr_array)[j]]])
    p_val <- cor.test(signflip_max_prs_dist[[rownames(signflip_max_vs_strat_corr_array)[i]]],
                      signflip_max_prs_dist[[colnames(signflip_max_vs_strat_corr_array)[j]]])$p.value
    signflip_max_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SIGNFLIP, MAX) Constant value observed for ', rownames(signflip_max_vs_strat_corr_array)[i])
      message('(SIGNFLIP, MAX) Value = ', mean(signflip_max_prs_dist[[rownames(signflip_max_vs_strat_corr_array)[i]]]))
    } else if (p_val < fwer.thres.cutoff) {
      message('(SIGNFLIP, MAX): Significant correlation (p = ',p_val,') between ', 
              rownames(signflip_max_vs_strat_corr_array)[i], ' and ', 
              colnames(signflip_max_vs_strat_corr_array)[j])
      message('(SIGNFLIP, MAX): Pearson r = ', signflip_max_vs_strat_corr_array[i,j])
    }
    
    # Shuffle, usual
    shuffle_usual_vs_strat_corr_array[i,j] <- cor(shuffle_usual_prs_dist[[rownames(shuffle_usual_vs_strat_corr_array)[i]]],
                                                shuffle_usual_prs_dist[[colnames(shuffle_usual_vs_strat_corr_array)[j]]])
    p_val <- cor.test(shuffle_usual_prs_dist[[rownames(shuffle_usual_vs_strat_corr_array)[i]]],
                      shuffle_usual_prs_dist[[colnames(shuffle_usual_vs_strat_corr_array)[j]]])$p.value
    shuffle_usual_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SHUFFLE, USUAL) Constant value observed for ', rownames(shuffle_usual_vs_strat_corr_array)[i])
      message('(SHUFFLE, USUAL) Value = ', mean(shuffle_usual_prs_dist[[rownames(shuffle_usual_vs_strat_corr_array)[i]]]))
    } else if (p_val < fwer.thres.cutoff) {
      message('(SHUFFLE, USUAL): Significant correlation (p = ',p_val,') between ', 
              rownames(shuffle_usual_vs_strat_corr_array)[i], ' and ', 
              colnames(shuffle_usual_vs_strat_corr_array)[j])
      message('(SHUFFLE, USUAL): Pearson r = ', shuffle_usual_vs_strat_corr_array[i,j])
    }
    
    # Signflip, usual
    signflip_usual_vs_strat_corr_array[i,j] <- cor(signflip_usual_prs_dist[[rownames(signflip_usual_vs_strat_corr_array)[i]]],
                                                   signflip_usual_prs_dist[[colnames(signflip_usual_vs_strat_corr_array)[j]]])
    p_val <- cor.test(signflip_usual_prs_dist[[rownames(signflip_usual_vs_strat_corr_array)[i]]],
                      signflip_usual_prs_dist[[colnames(signflip_usual_vs_strat_corr_array)[j]]])$p.value
    signflip_usual_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SIGNFLIP, USUAL) Constant value observed for ', rownames(signflip_usual_vs_strat_corr_array)[i])
      message('(SIGNFLIP, USUAL) Value = ', mean(signflip_usual_prs_dist[[rownames(signflip_usual_vs_strat_corr_array)[i]]]))
    } else if (p_val < fwer.thres.cutoff) {
      message('(SIGNFLIP, USUAL): Significant correlation (p = ',p_val,') between ', 
              rownames(signflip_usual_vs_strat_corr_array)[i], ' and ', 
              colnames(signflip_usual_vs_strat_corr_array)[j])
      message('(SIGNFLIP, USUAL): Pearson r = ', signflip_usual_vs_strat_corr_array[i,j])
    }
  }
}

# Convert to dataframe and include sensitivity metrics as a new feature column
shuffle_max_vs_strat_corr_array <- as.data.frame(shuffle_max_vs_strat_corr_array)
shuffle_max_vs_strat_corr_array$SENSITIVITY_METRIC <- rownames(shuffle_max_vs_strat_corr_array)
shuffle_max_vs_strat_pvals_array <- as.data.frame(shuffle_max_vs_strat_pvals_array)
shuffle_max_vs_strat_pvals_array$SENSITIVITY_METRIC <- rownames(shuffle_max_vs_strat_pvals_array)

signflip_max_vs_strat_corr_array <- as.data.frame(signflip_max_vs_strat_corr_array)
signflip_max_vs_strat_corr_array$SENSITIVITY_METRIC <- rownames(signflip_max_vs_strat_corr_array)
signflip_max_vs_strat_pvals_array <- as.data.frame(signflip_max_vs_strat_pvals_array)
signflip_max_vs_strat_pvals_array$SENSITIVITY_METRIC <- rownames(signflip_max_vs_strat_pvals_array)

shuffle_usual_vs_strat_corr_array <- as.data.frame(shuffle_usual_vs_strat_corr_array)
shuffle_usual_vs_strat_corr_array$SENSITIVITY_METRIC <- rownames(shuffle_usual_vs_strat_corr_array)
shuffle_usual_vs_strat_pvals_array <- as.data.frame(shuffle_usual_vs_strat_pvals_array)
shuffle_usual_vs_strat_pvals_array$SENSITIVITY_METRIC <- rownames(shuffle_usual_vs_strat_pvals_array)

signflip_usual_vs_strat_corr_array <- as.data.frame(signflip_usual_vs_strat_corr_array)
signflip_usual_vs_strat_corr_array$SENSITIVITY_METRIC <- rownames(signflip_usual_vs_strat_corr_array)
signflip_usual_vs_strat_pvals_array <- as.data.frame(signflip_usual_vs_strat_pvals_array)
signflip_usual_vs_strat_pvals_array$SENSITIVITY_METRIC <- rownames(signflip_usual_vs_strat_pvals_array)

# Save files
readr::write_csv(shuffle_max_vs_strat_corr_array,
                 file = paste0(out.dir,"corr_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_max_vs_prs_dist.csv"))
readr::write_csv(shuffle_max_vs_strat_pvals_array,
                 file = paste0(out.dir,"pvals_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_max_vs_prs_dist.csv"))

readr::write_csv(signflip_max_vs_strat_corr_array,
                 file = paste0(out.dir,"corr_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_max_vs_prs_dist.csv"))
readr::write_csv(signflip_max_vs_strat_pvals_array,
                 file = paste0(out.dir,"pvals_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_max_vs_prs_dist.csv"))

readr::write_csv(shuffle_usual_vs_strat_corr_array,
                 file = paste0(out.dir,"corr_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_usual_vs_prs_dist.csv"))
readr::write_csv(shuffle_usual_vs_strat_pvals_array,
                 file = paste0(out.dir,"pvals_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_usual_vs_prs_dist.csv"))

readr::write_csv(signflip_usual_vs_strat_corr_array,
                 file = paste0(out.dir,"corr_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_usual_vs_prs_dist.csv"))
readr::write_csv(signflip_usual_vs_strat_pvals_array,
                 file = paste0(out.dir,"pvals_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_usual_vs_prs_dist.csv"))
sink()

## Manual plotting -------------------------------------------------------------
# prefix.dir <- '/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/raw_tables/'
# shuffle_max_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_shuffle_max_combined.csv'))
# signflip_max_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_signflip_max_combined.csv'))
# shuffle_usual_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_shuffle_usual_combined.csv'))
# signflip_usual_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_signflip_usual_combined.csv'))
# 