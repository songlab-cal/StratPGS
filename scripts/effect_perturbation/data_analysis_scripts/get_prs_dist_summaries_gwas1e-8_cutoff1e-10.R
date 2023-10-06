#############################################
########### Shuffling Experiment ############
########## Using GWAS cutoff 1e-8 ###########
### Computation of variance and magnitudes ##
#############################################
# This script computes variances and magnitudes
# of effect sizes obtained from statistical 
# model.

# Modified 2/28/23 to include two-sample test of stochastic dominance 
# Modified 4/28/23 to compute for GWAS 1e-8 PRS obtained from C&T
# of background variant effect sizes over target variant effect sizes
# Modified 5/21/23 to compute summary stats on the new PRS files (sans 
# duplicate variants)

#################
## Directories ##
#################
# Create output file 
sink('/deep_learning/aaw/051723/logs/get_prs_dist_summaries_gwas1e-8_cutoff1e-10.log')
sink(stdout(), type = "message")

library(data.table)
library(dplyr)
library(bigsnpr)

R.workbench <- TRUE
gwas.thres <- '1e-8'
sig.cutoff <- 1e-10
message('gwas.thres = ', gwas.thres)
message('sig.cutoff = ', sig.cutoff)

if (R.workbench) {
  phenos.shortlist <- data.table::fread('/deep_learning/aaw/051723/results/pheno_names_gwas1e8.txt',
                                        header = FALSE)$V1
  pheno.gwas.dir <- '/deep_learning/aaw/051723/prs/'
  out.dir <- '/deep_learning/aaw/051723/results/prs_dist_summaries/'
} else {
  phenos.shortlist <- data.table::fread('/illumina/scratch/deep_learning/aaw/051723/results/pheno_names_gwas1e8.txt',
                                        header = FALSE)$V1
  pheno.gwas.dir <- '/illumina/scratch/deep_learning/aaw/051723/prs/'
  out.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_dist_summaries/'
}

# Save intermediate results to a connection just in case
# fileConn <- file(paste0(out.dir,"output_031423.txt"))

##########
## Main ##
##########
message(date(), ": Computing magnitudes and variances...")
mag.var.df <- data.frame(PHENOTYPE = character(),
                         N_SNPS = numeric(),
                         N_TARGET_SNPS = numeric(),
                         MAX_MAG_BKGRD = numeric(),
                         MEDIAN_MAG_BKGRD = numeric(),
                         AVE_MAG_BKGRD = numeric(),
                         BETA_WILCOX_PVAL = numeric(),
                         VAR_BKGRD = numeric(),
                         SUM_BKGRD = numeric(),
                         MAX_MAG_TARGET = numeric(),
                         MEDIAN_MAG_TARGET = numeric(),
                         AVE_MAG_TARGET = numeric(),
                         VAR_TARGET = numeric(),
                         SUM_TARGET = numeric(),
                         SUM_RATIO_TARGET_BKGRD = numeric())
message(date(),": ", length(phenos.shortlist), " phenotypes detected.")
for (pheno in phenos.shortlist) {
  # Print message
  message(date(), ": Working on ", pheno)
  
  # Load GWAS results 
  message("Loading GWAS results...")
  pheno.gwas <- data.table::fread(paste0(pheno.gwas.dir,
                                         pheno,
                                         ".loci.common.",gwas.thres,"_no_dups.txt"))
  
  # Compute variance and magnitudes of effect sizes
  # for participating SNPs and background SNPs
  message("Collecting metrics...")
  rel.inds <- which(pheno.gwas$gwas_p_value >= sig.cutoff)
  n.snps <- length(pheno.gwas$gwas_p_value)
  n.target.snps <- length(rel.inds)
  
  var.target <- var(pheno.gwas$beta[rel.inds])
  var.bkgrd <- var(pheno.gwas$beta[-rel.inds])
  
  abs.betas <- abs(pheno.gwas$beta)
  max.mag.target <- max(abs.betas[rel.inds])
  mean.mag.target <- mean(abs.betas[rel.inds])
  median.mag.target <- median(abs.betas[rel.inds])
  
  max.mag.bkgrd <- max(abs.betas[-rel.inds])
  mean.mag.bkgrd <- mean(abs.betas[-rel.inds])
  median.mag.bkgrd <- median(abs.betas[-rel.inds])
  
  # Are bkgrd variant esizes > target variant esizes?
  beta.wilcox.p <- wilcox.test(abs.betas[rel.inds],
                               abs.betas[-rel.inds],
                               alternative = 'greater')$p.value
  
  # What's the sum(|beta_target|)/sum(|beta_causal|)
  sum.target <- sum(abs.betas[rel.inds])
  sum.bkgrd <- sum(abs.betas[-rel.inds])
    
  new.row <- data.frame(PHENOTYPE = pheno,
                        N_SNPS = n.snps,
                        N_TARGET_SNPS = n.target.snps,
                        MAX_MAG_BKGRD = max.mag.bkgrd,
                        MEDIAN_MAG_BKGRD = median.mag.bkgrd,
                        AVE_MAG_BKGRD = mean.mag.bkgrd,
                        BETA_WILCOX_PVAL = beta.wilcox.p,
                        VAR_BKGRD = var.bkgrd,
                        SUM_BKGRD = sum.bkgrd,
                        MAX_MAG_TARGET = max.mag.target,
                        MEDIAN_MAG_TARGET = median.mag.target,
                        AVE_MAG_TARGET = mean.mag.target,
                        VAR_TARGET = var.target,
                        SUM_TARGET = sum.target,
                        SUM_RATIO_TARGET_BKGRD = sum.target/sum.bkgrd)
  
  # writeLines(as.vector(new.row) %>% `colnames<-`(NULL) %>% as.character(), 
  #            fileConn)

  # Add new row
  mag.var.df <- rbind(mag.var.df,
                      new.row)
  
}

readr::write_csv(mag.var.df, 
                 file = paste0(out.dir, "prs_dist_summaries_gwas",gwas.thres,"_cutoff",sig.cutoff,".csv"))
# close(fileConn)
sink()


#################
#### Analysis ###
#################
prs_dist_metrics <- readr::read_csv("/deep_learning/aaw/051723/results/prs_dist_summaries/prs_dist_summaries_gwas1e-8_cutoff1e-10.csv")
skimr::skim(prs_dist_metrics[,-1])

old_prs_dist_metrics <- readr::read_csv("/deep_learning/aaw/042023/tables/vars_and_mags_df_gwas1e-8_cutoff1e-10.csv")
skimr::skim(old_prs_dist_metrics[,-1])

library(dplyr);library(ggplot2)
new.df <- prs_dist_metrics %>% 
  select(c('PHENOTYPE', 'SUM_RATIO_TARGET_BKGRD')) 
old.df <- old_prs_dist_metrics %>% 
  select(c('PHENOTYPE', 'SUM_RATIO_TARGET_BKGRD'))
colnames(old.df)[2] <- 'Old'; colnames(new.df)[2] <- 'New'
new.vs.old.df <- left_join(x=new.df, y=old.df, by='PHENOTYPE')


plot4 <- ggplot(new.vs.old.df,aes(x=Old, y=New)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1,intercept=0,lty='dashed')+
  xlim(c(0,0.25)) + ylim(c(0,0.25)) +
  ggtitle('(GWAS,cutoff)=(1e-8,1e-10)')

library(gridExtra)

combined.plots <- gridExtra::grid.arrange(plot1,plot2,plot3,plot4,
                        ncol=2,nrow=2)
ggsave(combined.plots,
       filename = paste0('/deep_learning/aaw/051723/plots/beta_sum_ratios_old_dist_vs_new_dist.jpg'),
       width = 9.2, height = 9.2,
       dpi = 400)
# colnames(vars_and_mags_df)[1]<-'PHENO'
# # library(GGally)
# # pairwise.plot <- ggpairs(vars_and_mags_df[,-1])
# # ggsave(pairwise.plot,
# #        filename = "/deep_learning/aaw/022823/vars_and_mags_pairwise.jpg",
# #        width = 15, height = 15,
# #        dpi = 400)
# # 
# # wilcox.test(vars_and_mags_df$MAX_MAG_BKGRD,
# #             vars_and_mags_df$MAX_MAG_TARGET,
# #             paired = TRUE,
# #             alternative='greater') # p-value = 1.322e-15
# # wilcox.test(vars_and_mags_df$MEDIAN_MAG_BKGRD,
# #             vars_and_mags_df$MEDIAN_MAG_TARGET,
# #             paired = TRUE,
# #             alternative='less') # p-value = 1.275e-15
# # wilcox.test(vars_and_mags_df$AVE_MAG_BKGRD,
# #             vars_and_mags_df$AVE_MAG_TARGET,
# #             paired = TRUE,
# #             alternative='less') # p-value = 1.674e-12
# # wilcox.test(vars_and_mags_df$VAR_BKGRD,
# #             vars_and_mags_df$VAR_TARGET,
# #             paired = TRUE,
# #             alternative='less') # p-value = 0.1617
# 
# ## Shuffle ---------------------------------------------------------------------
# shuffle_metrics_df <- readr::read_csv("/deep_learning/aaw/042023/tables/gwas_1e-8_cutoff1e-10_shuffle_metrics_df.csv")
# merged.df <- merge(vars_and_mags_df,shuffle_metrics_df, by='PHENO')
# colnames(merged.df)[16:20] <- c('PEARSON_CORR', 'SPEARMAN_RHO', 'COSINE_SIM',
#                                 'TOP10PCT_PREV', 'TOP1PCT_PREV')
# 
# corr.df <- data.frame(EMP_DIST_METRIC = colnames(merged.df)[2:15],
#                       PEARSON_CORR = rep(0, length(colnames(merged.df)[2:15])),
#                       SPEARMAN_RHO = rep(0, length(colnames(merged.df)[2:15])),
#                       COSINE_SIM = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP10PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP1PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP10PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP1PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PREV_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PREV_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PHENO_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])))
# 
# corr.pval.df <- data.frame(EMP_DIST_METRIC = colnames(merged.df)[2:15],
#                            PEARSON_CORR = rep(0, length(colnames(merged.df)[2:15])),
#                            SPEARMAN_RHO = rep(0, length(colnames(merged.df)[2:15])),
#                            COSINE_SIM = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP10PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP1PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP10PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP1PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PREV_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PREV_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PHENO_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])))
# 
# for (i in 1:nrow(corr.df)) {
#   for (j in 1:(ncol(corr.df)-1)) {
#     corr.df[i,][[colnames(corr.df)[j+1]]] <- cor(merged.df[[corr.df$EMP_DIST_METRIC[i]]],merged.df[[colnames(corr.df)[j+1]]]) 
#     corr.pval.df[i,][[colnames(corr.pval.df)[j+1]]] <- cor.test(merged.df[[corr.pval.df$EMP_DIST_METRIC[i]]],merged.df[[colnames(corr.pval.df)[j+1]]])$p.value 
#   }
# }
# 
# 
# apply(corr.pval.df[,-1],2,function(x) { (x < 0.05/156)}) # 1 significant p-value
# 
# readr::write_csv(corr.pval.df, 
#                  file = '/deep_learning/aaw/042023/tables/polygenic_vec_vs_insensitivity_pvals_gwas1e-8_cutoff1e-10_shuffle.csv')
# readr::write_csv(corr.df, 
#                  file = '/deep_learning/aaw/042023/tables/polygenic_vec_vs_insensitivity_corrs_gwas1e-8_cutoff1e-10_shuffle.csv')
# 
# ggplot(merged.df,aes(x=N_TARGET_SNPS,y=SPEARMAN_RHO)) +
#   geom_point() +
#   #geom_abline(slope=1,intercept=0) +
#   #guides(colour=guide_colourbar(title="Insensitivity\nDegree")) +
#   #scale_colour_gradient(low = "#fdbb84", high = "#990000") +
#   #ylim(c(0,0.007)) +
#   #xlim(c(0,0.007)) +
#   theme_bw()
# 
# # corr(SUM_RATIO_TARGET_BKGRD,PEARSON_CORR) p-value = 0.0001850996
# 
# ## Signflip --------------------------------------------------------------------
# signflip_metrics_df <- readr::read_csv("/deep_learning/aaw/042023/tables/gwas_1e-8_cutoff1e-10_signflip_metrics_df.csv")
# merged.df <- merge(vars_and_mags_df,signflip_metrics_df, by='PHENO')
# colnames(merged.df)[16:20] <- c('PEARSON_CORR', 'SPEARMAN_RHO', 'COSINE_SIM',
#                                 'TOP10PCT_PREV', 'TOP1PCT_PREV')
# 
# corr.df <- data.frame(EMP_DIST_METRIC = colnames(merged.df)[2:15],
#                       PEARSON_CORR = rep(0, length(colnames(merged.df)[2:15])),
#                       SPEARMAN_RHO = rep(0, length(colnames(merged.df)[2:15])),
#                       COSINE_SIM = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP10PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP1PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP10PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                       TOP1PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PREV_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PREV_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PHENO_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                       PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])))
# 
# corr.pval.df <- data.frame(EMP_DIST_METRIC = colnames(merged.df)[2:15],
#                            PEARSON_CORR = rep(0, length(colnames(merged.df)[2:15])),
#                            SPEARMAN_RHO = rep(0, length(colnames(merged.df)[2:15])),
#                            COSINE_SIM = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP10PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP1PCT_PREV = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP10PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                            TOP1PCT_OR = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PREV_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PREV_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PHENO_SLOPE = rep(0, length(colnames(merged.df)[2:15])),
#                            PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, length(colnames(merged.df)[2:15])))
# 
# for (i in 1:nrow(corr.df)) {
#   for (j in 1:(ncol(corr.df)-1)) {
#     corr.df[i,][[colnames(corr.df)[j+1]]] <- cor(merged.df[[corr.df$EMP_DIST_METRIC[i]]],merged.df[[colnames(corr.df)[j+1]]]) 
#     corr.pval.df[i,][[colnames(corr.pval.df)[j+1]]] <- cor.test(merged.df[[corr.pval.df$EMP_DIST_METRIC[i]]],merged.df[[colnames(corr.pval.df)[j+1]]])$p.value 
#   }
# }
# 
# 
# apply(corr.pval.df[,-1],2,function(x) { (x < 0.05/156)}) # No significant p-values
# 
# readr::write_csv(corr.pval.df, 
#                  file = '/deep_learning/aaw/042023/tables/polygenic_vec_vs_insensitivity_pvals_gwas1e-8_cutoff1e-10_signflip.csv')
# readr::write_csv(corr.df, 
#                  file = '/deep_learning/aaw/042023/tables/polygenic_vec_vs_insensitivity_corrs_gwas1e-8_cutoff1e-10_signflip.csv')
