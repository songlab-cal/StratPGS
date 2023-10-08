#############################################
########## Perturbation Experiment ##########
########## Using GWAS cutoff 1e-5 ###########
####### Sign flips and eSize Shuffles #######
#############################################
# This script combines stratification of PRS metrics by PCs
# and sensitivity of PRS to perturbations on a training set
# Created on 5/30/23: This version takes max(f(y,PRS^+),f(y,PRS^-))
#                     Also computes the f(y,PRS) version (aka 'usual approach')
# objects with 'usual' appended in front correspond to the f(y,PRS) version of things
# Modified on 6/1/23: for gwas.thres = 1e-5
# Modified on 10/6/23: perform Spearman correlation tests, change directories 
# to ensure Github-compatibility

## Log file / libraries --------------------------------------------------------
sink('~/Documents/GitHub/StratPRS/logs/gwas1e-5_cutoff1e-6_strat_vs_perturb_sensitivity.log')
sink(stdout(), type = "message")
library(data.table)
library(dplyr)
#library(bigsnpr)

## Variables -------------------------------------------------------------------
#R.workbench <- FALSE
gwas.thres <- '1e-5'
sig.cutoff <- '1e-6'
message('gwas.thres = ', gwas.thres)
message('sig.cutoff = ', sig.cutoff)

# if (R.workbench) {
#   phenos.shortlist <- data.table::fread('/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
#                                         header = FALSE)$V1
#   out.dir <- '/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/'
#   prs.strat.dir <- '/deep_learning/aaw/051723/results/prs_pc_stratification/'
#   pheno.strat.dir <- '/deep_learning/aaw/051723/results/pheno_pc_stratification/'
#   perturb.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
# } else {
#   phenos.shortlist <- data.table::fread('/illumina/scratch/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
#                                         header = FALSE)$V1
#   out.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/'
#   prs.strat.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_pc_stratification/'
#   pheno.strat.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/pheno_pc_stratification/'
#   perturb.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
# }
# 
## Helper function for computing Spearman correlation --------------------------
# The function below is based on cor.test, and does the following.
# It performs a Spearman correlation test without the t-distribution approximation
# (i.e., the parameter exact=TRUE) in case there are no ties in ranks,
# but switches to the t-distribution approximation (i.e., the parameter exact=FALSE)
# when ties in ranks are detected.
spearman_cor_test <- function(x, y, ...) {
  # Check for ties in ranks
  ties_check <- length(unique(rank(x))) < length(x) || length(unique(rank(y))) < length(y)

  if (ties_check) {
    message("Ties detected, using exact = FALSE...")
    # Ties detected, use t-distribution approximation
    result <- cor.test(x, y, method = "spearman", exact = FALSE, ...)
  } else {
    # No ties, use exact calculation
    result <- cor.test(x, y, method = "spearman", exact = TRUE, ...)
  }

  return(result)
}

# ## Main Body -------------------------------------------------------------------
# # Read perturbation sensitivity metrics
# prs_shuffle_max <- readr::read_csv(paste0(perturb.dir,'gwas',
#                                           gwas.thres,'_cutoff',
#                                           sig.cutoff,'_shuffle_metrics_df.csv'),
#                                    show_col_types = FALSE) %>% 
#   as.data.frame()
# prs_signflip_max <- readr::read_csv(paste0(perturb.dir,'gwas',
#                                            gwas.thres,'_cutoff',
#                                            sig.cutoff,'_signflip_metrics_df.csv'),
#                                     show_col_types = FALSE) %>% 
#   as.data.frame()
# prs_shuffle_usual <- readr::read_csv(paste0(perturb.dir,'gwas',
#                                           gwas.thres,'_cutoff',
#                                           sig.cutoff,'_shuffle_usual_metrics_df.csv'),
#                                      show_col_types = FALSE) %>% 
#   as.data.frame()
# prs_signflip_usual <- readr::read_csv(paste0(perturb.dir,'gwas',
#                                            gwas.thres,'_cutoff',
#                                            sig.cutoff,'_signflip_usual_metrics_df.csv'),
#                                       show_col_types = FALSE) %>% 
#   as.data.frame()
# 
# # Read PRS stratification metrics
# prs_strat <- readr::read_csv(paste0(prs.strat.dir,
#                                     'train_n288728_prs_gwas_',
#                                     gwas.thres,
#                                     '_gPC_metrics.csv'),
#                              show_col_types = FALSE) %>% 
#   as.data.frame()
# prs_strat$PC1_COSSIM <- NULL;prs_strat$PC1_PEARSON <- NULL;prs_strat$PC1_SPEARMAN <- NULL
# 
# # Combine 
# shuffle_max_prs_strat <- merge(prs_shuffle_max,prs_strat,by='PHENO')
# signflip_max_prs_strat <- merge(prs_signflip_max,prs_strat,by='PHENO')
# shuffle_usual_prs_strat <- merge(prs_shuffle_usual,prs_strat,by='PHENO')
# signflip_usual_prs_strat <- merge(prs_signflip_usual,prs_strat,by='PHENO')
# 
# readr::write_csv(shuffle_max_prs_strat,
#                  file = paste0(out.dir,"/raw_tables/gwas",
#                                gwas.thres,"_cutoff",
#                                sig.cutoff,"_shuffle_max_combined.csv"))
# readr::write_csv(signflip_max_prs_strat,
#                  file = paste0(out.dir,"/raw_tables/gwas",
#                                gwas.thres,"_cutoff",
#                                sig.cutoff,"_signflip_max_combined.csv"))
# 
# readr::write_csv(shuffle_usual_prs_strat,
#                  file = paste0(out.dir,"/raw_tables/gwas",
#                                gwas.thres,"_cutoff",
#                                sig.cutoff,"_shuffle_usual_combined.csv"))
# readr::write_csv(signflip_usual_prs_strat,
#                  file = paste0(out.dir,"/raw_tables/gwas",
#                                gwas.thres,"_cutoff",
#                                sig.cutoff,"_signflip_usual_combined.csv"))

# Read raw tables 
user_defined_results_dir <- "~/Documents/GitHub/StratPRS/results/effect_perturbation/"
prs_strat <- readr::read_csv(paste0(user_defined_results_dir,
                                    'prs_pc_stratification/train_n288728_prs_gwas_',
                                    gwas.thres,
                                    '_gPC_metrics.csv'),
                             show_col_types = FALSE) %>% 
  as.data.frame()
prs_strat$PC1_COSSIM <- NULL;prs_strat$PC1_PEARSON <- NULL;prs_strat$PC1_SPEARMAN <- NULL

prs_shuffle_usual <- readr::read_csv(paste0(user_defined_results_dir,
                                            'prs_perturbation_sensitivity/gwas',
                                            gwas.thres,'_cutoff',
                                            sig.cutoff,'_shuffle_usual_metrics_df.csv'),
                                     show_col_types = FALSE) %>% 
  as.data.frame()
prs_shuffle_max <- readr::read_csv(paste0(user_defined_results_dir,
                                            'prs_perturbation_sensitivity/gwas',
                                            gwas.thres,'_cutoff',
                                            sig.cutoff,'_shuffle_metrics_df.csv'),
                                     show_col_types = FALSE) %>% 
  as.data.frame()
prs_signflip_usual <- readr::read_csv(paste0(user_defined_results_dir,
                                            'prs_perturbation_sensitivity/gwas',
                                            gwas.thres,'_cutoff',
                                            sig.cutoff,'_signflip_usual_metrics_df.csv'),
                                     show_col_types = FALSE) %>% 
  as.data.frame()
prs_signflip_max <- readr::read_csv(paste0(user_defined_results_dir,
                                          'prs_perturbation_sensitivity/gwas',
                                          gwas.thres,'_cutoff',
                                          sig.cutoff,'_signflip_metrics_df.csv'),
                                   show_col_types = FALSE) %>% 
  as.data.frame()

signflip_max_prs_strat <- readr::read_csv(paste0(user_defined_results_dir,
                                                 "/stratification_vs_perturbation_sensitivity_raw_tables/gwas",
                                                 gwas.thres,"_cutoff",
                                                 sig.cutoff,"_signflip_max_combined.csv")) %>%
  as.data.frame()
shuffle_max_prs_strat <- readr::read_csv(paste0(user_defined_results_dir,
                                                "/stratification_vs_perturbation_sensitivity_raw_tables/gwas",
                                                gwas.thres,"_cutoff",
                                                sig.cutoff,"_shuffle_max_combined.csv"))%>%
  as.data.frame()
signflip_usual_prs_strat <- readr::read_csv(paste0(user_defined_results_dir,
                                                   "/stratification_vs_perturbation_sensitivity_raw_tables/gwas",
                                                   gwas.thres,"_cutoff",
                                                   sig.cutoff,"_signflip_usual_combined.csv"))%>%
  as.data.frame()
shuffle_usual_prs_strat <- readr::read_csv(paste0(user_defined_results_dir,
                                                  "/stratification_vs_perturbation_sensitivity_raw_tables/gwas",
                                                  gwas.thres,"_cutoff",
                                                  sig.cutoff,"_shuffle_usual_combined.csv"))%>%
  as.data.frame()

# Compute p-values and correlations
shuffle_max_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_shuffle_max)-1,ncol=ncol(prs_strat)-1)
signflip_max_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_signflip_max)-1,ncol=ncol(prs_strat)-1)
shuffle_usual_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_shuffle_usual)-1,ncol=ncol(prs_strat)-1)
signflip_usual_vs_strat_corr_array <- matrix(NA,nrow=ncol(prs_signflip_usual)-1,ncol=ncol(prs_strat)-1)

shuffle_max_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_shuffle_max)-1,ncol=ncol(prs_strat)-1)
signflip_max_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_signflip_max)-1,ncol=ncol(prs_strat)-1)
shuffle_usual_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_shuffle_usual)-1,ncol=ncol(prs_strat)-1)
signflip_usual_vs_strat_pvals_array <- matrix(NA,nrow=ncol(prs_signflip_usual)-1,ncol=ncol(prs_strat)-1)

colnames(shuffle_max_vs_strat_corr_array) <- colnames(prs_strat)[-1]
rownames(shuffle_max_vs_strat_corr_array) <- colnames(prs_shuffle_max)[-1]
colnames(signflip_max_vs_strat_corr_array) <- colnames(prs_strat)[-1]
rownames(signflip_max_vs_strat_corr_array) <- colnames(prs_signflip_max)[-1]
colnames(shuffle_usual_vs_strat_corr_array) <- colnames(prs_strat)[-1]
rownames(shuffle_usual_vs_strat_corr_array) <- colnames(prs_shuffle_usual)[-1]
colnames(signflip_usual_vs_strat_corr_array) <- colnames(prs_strat)[-1]
rownames(signflip_usual_vs_strat_corr_array) <- colnames(prs_signflip_usual)[-1]

colnames(shuffle_max_vs_strat_pvals_array) <- colnames(prs_strat)[-1]
rownames(shuffle_max_vs_strat_pvals_array) <- colnames(prs_shuffle_max)[-1]
colnames(signflip_max_vs_strat_pvals_array) <- colnames(prs_strat)[-1]
rownames(signflip_max_vs_strat_pvals_array) <- colnames(prs_signflip_max)[-1]
colnames(shuffle_usual_vs_strat_pvals_array) <- colnames(prs_strat)[-1]
rownames(shuffle_usual_vs_strat_pvals_array) <- colnames(prs_shuffle_usual)[-1]
colnames(signflip_usual_vs_strat_pvals_array) <- colnames(prs_strat)[-1]
rownames(signflip_usual_vs_strat_pvals_array) <- colnames(prs_signflip_usual)[-1]

n.tests <- nrow(shuffle_max_vs_strat_corr_array)*ncol(shuffle_max_vs_strat_corr_array)
message(date(), ": No. tests performed = ", n.tests)
fwer.thres.cutoff <- 0.05/n.tests
message(date(), ": FWER corrected p-value threshold = ", fwer.thres.cutoff)
for (i in 1:nrow(shuffle_max_vs_strat_corr_array)) {
  for (j in 1:ncol(shuffle_max_vs_strat_corr_array)) {
    # Print pair 
    message('----- Working on (', rownames(shuffle_max_vs_strat_corr_array)[i],
            ', ',colnames(shuffle_max_vs_strat_corr_array)[j],') -----')
    # Shuffle, max
    shuffle_max_res <- spearman_cor_test(shuffle_max_prs_strat[[rownames(shuffle_max_vs_strat_corr_array)[i]]],
                                         shuffle_max_prs_strat[[colnames(shuffle_max_vs_strat_corr_array)[j]]])
    shuffle_max_vs_strat_corr_array[i,j] <- shuffle_max_res$estimate
    p_val <- shuffle_max_res$p.value
    shuffle_max_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SHUFFLE, MAX): NaN p-value observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(SHUFFLE, MAX): Significant correlation, p = ',p_val)
      message('(SHUFFLE, MAX): Spearman rho = ', shuffle_max_vs_strat_corr_array[i,j])
    }
    
    # Signflip, max
    signflip_max_res <- spearman_cor_test(signflip_max_prs_strat[[rownames(signflip_max_vs_strat_corr_array)[i]]],
                                          signflip_max_prs_strat[[colnames(signflip_max_vs_strat_corr_array)[j]]])
    signflip_max_vs_strat_corr_array[i,j] <- signflip_max_res$estimate
    p_val <- signflip_max_res$p.value
    signflip_max_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SIGNFLIP, MAX): NaN p-value observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(SIGNFLIP, MAX): Significant correlation, p = ',p_val)
      message('(SIGNFLIP, MAX): Spearman rho = ', signflip_max_vs_strat_corr_array[i,j])
    }
    
    # Shuffle, usual
    shuffle_usual_res <- spearman_cor_test(shuffle_usual_prs_strat[[rownames(shuffle_usual_vs_strat_corr_array)[i]]],
                                           shuffle_usual_prs_strat[[colnames(shuffle_usual_vs_strat_corr_array)[j]]])
    shuffle_usual_vs_strat_corr_array[i,j] <- shuffle_usual_res$estimate
    p_val <- shuffle_usual_res$p.value
    shuffle_usual_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SHUFFLE, USUAL): NaN p-value observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(SHUFFLE, USUAL): Significant correlation, p = ',p_val)
      message('(SHUFFLE, USUAL): Spearman rho = ', shuffle_usual_vs_strat_corr_array[i,j])
    }
    
    # Signflip, usual
    signflip_usual_res <- spearman_cor_test(signflip_usual_prs_strat[[rownames(signflip_usual_vs_strat_corr_array)[i]]],
                                            signflip_usual_prs_strat[[colnames(signflip_usual_vs_strat_corr_array)[j]]])
    signflip_usual_vs_strat_corr_array[i,j] <- signflip_usual_res$estimate
    p_val <- signflip_usual_res$p.value
    signflip_usual_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(SIGNFLIP, USUAL): NaN observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(SIGNFLIP, USUAL): Significant correlation, p = ',p_val)
      message('(SIGNFLIP, USUAL): Spearman rho = ', signflip_usual_vs_strat_corr_array[i,j])
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
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/corr_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_max_vs_strat.csv"))
readr::write_csv(shuffle_max_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/pvals_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_max_vs_strat.csv"))

readr::write_csv(signflip_max_vs_strat_corr_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/corr_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_max_vs_strat.csv"))
readr::write_csv(signflip_max_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/pvals_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_max_vs_strat.csv"))

readr::write_csv(shuffle_usual_vs_strat_corr_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/corr_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_usual_vs_strat.csv"))
readr::write_csv(shuffle_usual_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/pvals_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_shuffle_usual_vs_strat.csv"))

readr::write_csv(signflip_usual_vs_strat_corr_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/corr_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_usual_vs_strat.csv"))
readr::write_csv(signflip_usual_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "stratification_vs_perturbation_sensitivity_pvals_and_cors/pvals_spearman_gwas",
                               gwas.thres,"_cutoff",
                               sig.cutoff,"_signflip_usual_vs_strat.csv"))
sink()