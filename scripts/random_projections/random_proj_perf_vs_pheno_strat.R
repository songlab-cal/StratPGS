#####################################################
########## Analysis of Random Projections ###########
#####################################################
# Created on 10/6/23: Perform Spearman correlation tests, define directories 
# to ensure Github-compatibility

## Log file / libraries --------------------------------------------------------
sink('~/Documents/GitProjects/StratPRS/logs/pheno_strat_vs_random_proj_perf.log')
sink(stdout(), type = "message")
library(data.table)
library(dplyr)

user_defined_results_dir <- "/Users/alanaw/Documents/GitProjects/StratPRS/results/random_projections/"
raw_orig_1stvisit <- readr::read_csv(paste0(user_defined_results_dir,
                                            "stratification_vs_predictability/ORIG_1stVisit_phenos_metrics.csv"))
IRNT_orig <- readr::read_csv(paste0(user_defined_results_dir,
                                    "stratification_vs_predictability/IRNT_phenos_metrics.csv"))

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

## Main Body -------------------------------------------------------------------
# Get p-values and correlations
# I follow code in analysis_perturb_gwas1e-5_cutoff1e-8.R 
# to generate matrices of p-values as well as raw correlations

# For original 1st visit phenos (141) -- OLD VERSION, UNCOMMENT FOR USE
#strat_metrics <- c('PC1_ABS_COS_SIM','PC1_ABS_PEARSON','RANK_PC1_COS_SIM','RANK_PC1_PEARSON','EVEN_COS_SIM','EVEN_PEARSON')
#non_random_perf_metrics <- c('VAR_PP_CORS','MAX_ABS_PP_CORS','MEAN_ABS_PP_CORS','VAR_BETAS','MAX_ABS_BETAS','MEAN_ABS_BETAS',
#                             'MEAN.DELTA.R2','MEAN.DELTA.ADJ.R2')
# For original 1st visit phenos (141) -- MODIFIED 6/22/23
strat_metrics <- c('PC1_ABS_COS_SIM','PC1_ABS_PEARSON',
                   'RANK_PC1_COS_SIM','RANK_PC1_PEARSON',
                   'EVEN_COS_SIM','EVEN_PEARSON')
non_random_perf_metrics <- c('MEAN_ABS_PP_CORS','MEAN.DELTA.R2')

# Original raw, 1st visit results
orig_rp_perf_vs_strat_corr_array <- matrix(NA,
                                      nrow=length(non_random_perf_metrics),
                                      ncol=length(strat_metrics))
orig_rp_perf_vs_strat_pvals_array <- matrix(NA,
                                       nrow=length(non_random_perf_metrics),
                                       ncol=length(strat_metrics))

colnames(orig_rp_perf_vs_strat_corr_array) <- strat_metrics
rownames(orig_rp_perf_vs_strat_corr_array) <- non_random_perf_metrics
colnames(orig_rp_perf_vs_strat_pvals_array) <- strat_metrics
rownames(orig_rp_perf_vs_strat_pvals_array) <- non_random_perf_metrics

# IRNT raw results
IRNT_rp_perf_vs_strat_corr_array <- matrix(NA,
                                           nrow=length(non_random_perf_metrics),
                                           ncol=length(strat_metrics))
IRNT_rp_perf_vs_strat_pvals_array <- matrix(NA,
                                            nrow=length(non_random_perf_metrics),
                                            ncol=length(strat_metrics))

colnames(IRNT_rp_perf_vs_strat_corr_array) <- strat_metrics
rownames(IRNT_rp_perf_vs_strat_corr_array) <- non_random_perf_metrics
colnames(IRNT_rp_perf_vs_strat_pvals_array) <- strat_metrics
rownames(IRNT_rp_perf_vs_strat_pvals_array) <- non_random_perf_metrics

n.tests <- nrow(orig_rp_perf_vs_strat_corr_array)*ncol(orig_rp_perf_vs_strat_corr_array)
message(date(), ": No. tests performed = ", n.tests)
fwer.thres.cutoff <- 0.05/n.tests
message(date(), ": FWER corrected p-value threshold = ", fwer.thres.cutoff)
for (i in 1:nrow(orig_rp_perf_vs_strat_corr_array)) {
  for (j in 1:ncol(orig_rp_perf_vs_strat_corr_array)) {
    # Print pair 
    message('----- Working on (', rownames(orig_rp_perf_vs_strat_corr_array)[i],
            ', ',colnames(orig_rp_perf_vs_strat_corr_array)[j],') -----')
    # Raw, original 1st visit
    raw_orig_1stvisit_res <- spearman_cor_test(raw_orig_1stvisit[[rownames(orig_rp_perf_vs_strat_corr_array)[i]]],
                                               raw_orig_1stvisit[[colnames(orig_rp_perf_vs_strat_corr_array)[j]]])
    orig_rp_perf_vs_strat_corr_array[i,j] <- raw_orig_1stvisit_res$estimate
    p_val <- raw_orig_1stvisit_res$p.value
    orig_rp_perf_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(Raw, original 1st visit): NaN p-value observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(Raw, original 1st visit): Significant correlation, p = ',p_val)
      message('(Raw, original 1st visit): Spearman rho = ', orig_rp_perf_vs_strat_corr_array[i,j])
    }
    
    # IRNT, original
    IRNT_orig_res <- spearman_cor_test(IRNT_orig[[rownames(IRNT_rp_perf_vs_strat_corr_array)[i]]],
                                       IRNT_orig[[colnames(IRNT_rp_perf_vs_strat_corr_array)[j]]])
    IRNT_rp_perf_vs_strat_corr_array[i,j] <- IRNT_orig_res$estimate
    p_val <- IRNT_orig_res$p.value
    IRNT_rp_perf_vs_strat_pvals_array[i,j] <- p_val
    if (is.na(p_val)) {
      message('(IRNT, original): NaN observed')
    } else if (p_val < fwer.thres.cutoff) {
      message('(IRNT, original): Significant correlation, p = ',p_val)
      message('(IRNT, original): Spearman rho = ', IRNT_rp_perf_vs_strat_corr_array[i,j])
    }
  }
}

# Convert to dataframe and include sensitivity metrics as a new feature column
orig_rp_perf_vs_strat_corr_array <- as.data.frame(orig_rp_perf_vs_strat_corr_array)
orig_rp_perf_vs_strat_corr_array$PERFORMANCE_METRIC <- rownames(orig_rp_perf_vs_strat_corr_array)
orig_rp_perf_vs_strat_pvals_array <- as.data.frame(orig_rp_perf_vs_strat_pvals_array)
orig_rp_perf_vs_strat_pvals_array$PERFORMANCE_METRIC <- rownames(orig_rp_perf_vs_strat_pvals_array)

IRNT_rp_perf_vs_strat_corr_array <- as.data.frame(IRNT_rp_perf_vs_strat_corr_array)
IRNT_rp_perf_vs_strat_corr_array$PERFORMANCE_METRIC <- rownames(IRNT_rp_perf_vs_strat_corr_array)
IRNT_rp_perf_vs_strat_pvals_array <- as.data.frame(IRNT_rp_perf_vs_strat_pvals_array)
IRNT_rp_perf_vs_strat_pvals_array$PERFORMANCE_METRIC <- rownames(IRNT_rp_perf_vs_strat_pvals_array)

# Save files
readr::write_csv(orig_rp_perf_vs_strat_corr_array,
                 file = paste0(user_defined_results_dir,
                               "pvals_and_cors/corr_spearman_raw_orig_1stvisit.csv"))
readr::write_csv(orig_rp_perf_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "pvals_and_cors/pvals_spearman_raw_orig_1stvisit.csv"))

readr::write_csv(IRNT_rp_perf_vs_strat_corr_array,
                 file = paste0(user_defined_results_dir,
                               "pvals_and_cors/corr_spearman_IRNT_orig.csv"))
readr::write_csv(IRNT_rp_perf_vs_strat_pvals_array,
                 file = paste0(user_defined_results_dir,
                               "pvals_and_cors/pvals_spearman_IRNT_orig.csv"))

sink()
