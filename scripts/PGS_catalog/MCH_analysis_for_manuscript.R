####################################################
########### Mean Corpuscular Hemoglobin ############
############### Shuffling Experiment ###############
############# Shuffle and Sign flips  ##############
####################################################
# This script analyzes original PRS performance 
# and sensitivity of perturbations of each PRS

#################
## Directories ##
#################
library(data.table)
library(dplyr)
library(bigsnpr)
sig.cutoff<-'1e-8'
R.workbench <- TRUE
message('sig.cutoff = ', sig.cutoff)

if (R.workbench) {
  prs.dist.dir <- '/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/tables/'
  c_and_t.prs.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/phenotypes/'
  c_and_t.dist.dir <- '/deep_learning/aaw/051723/results/prs_dist_summaries/'
  shuffled.prs.dir <- '/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/prs/'
  pheno.gwas.dir <- '/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/gwas_matched_prs/'
  out.dir <- '/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/'
} else {
  prs.dist.dir <- '/illumina/scratch/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/tables/'
  c_and_t.prs.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/phenotypes/'
  c_and_t.dist.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_dist_summaries/'
  shuffled.prs.dir <- '/illumina/scratch/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/prs/'
  pheno.gwas.dir <- '/illumina/scratch/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/gwas_matched_prs/'
  out.dir <- '/illumina/scratch/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/'
}

## Main Body -------------------------------------------------------------------
message(date(), ": Computing original performances of Mean Corpuscular Haemoglobin PRSs")
# FULL IDS (UNCOMMENT TO RUN FOR ALL PRSs) -------------
prs.ids <- c("PGS000099","PGS000174","PGS001219","PGS001989","PGS002206","PGS002339",
             "PGS002371","PGS002411","PGS002460","PGS002509","PGS002558","PGS002607",
             "PGS002656","PGS002705","PGS003560")
# PARTIAL IDS (SINCE SOME TAKE A WHILE TO RUN) ---------
# prs.ids <- c("PGS000099","PGS000174","PGS001219","PGS001989",
#              "PGS002206","PGS002411","PGS002460","PGS002509",
#              "PGS002558","PGS002607","PGS002656","PGS003560")
message(date(), ": No. PRSes = ", length(prs.ids))

## Analysis of original PRS ----------------------------------------------------
orig_prs_df <- readr::read_csv('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/tables/orig_prs_performance.csv',
                               show_col_types = FALSE)
orig_prs_df_minus_BOLTLMM_BBJ <- orig_prs_df %>% subset(PRS_ID!='PGS002371')
apply(orig_prs_df_minus_BOLTLMM_BBJ,2,min)
orig_prs_ranked_df <- readr::read_csv('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/tables/orig_prs_ranked_performance.csv',
                                      show_col_types = FALSE)
prs_dist_cutoff <- readr::read_csv(paste0(prs.dist.dir,
                                          'prs_dist_summaries_cutoff1e-10.csv'),
                                   show_col_types = FALSE) %>% 
  as.data.frame()
colnames(prs_dist_cutoff)[1] <- 'PRS_ID'
to_join <- prs_dist_cutoff %>% select(c('PRS_ID','N_SNPS'))
ranks_vs_no_snps <- left_join(ranked.PRS.scores.df %>% select(c('PRS_ID','AGGREGATED')),
                              to_join, by='PRS_ID')
ranks_vs_no_snps$N_SNPS[16] <- 372
ranks_vs_no_snps$N_SNPS[17] <- 320

# PLOT A: Original PRS performance vs no. variants
plot_A <- ggplot(ranks_vs_no_snps,aes(y=AGGREGATED,
                            x=log10(N_SNPS),
                            label=PRS_ID)) +
  geom_point() +
  #geom_text(hjust=0,vjust=0) +
  geom_label_repel(aes(label = PRS_ID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   alpha=0.5,
                   segment.color = 'grey50') +
  theme_bw() +
  xlab(expression(paste(log[10],'(No. Variants)'))) +
  xlim(c(2,7)) +
  ylab('Aggregate Rank Across Metrics') +
  ggtitle('A. Performance of MCH PRSs') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16))

## Analysis of shuffle sensitivity
metric_names <- c("PEARSON_CORR","SPEARMAN_RHO","COSINE_SIM",
                  "TOP10PCT_PREV","TOP1PCT_PREV","TOP10PCT_OR",
                  "TOP1PCT_OR","PERCENTILE_AVE_PREV_SLOPE","PERCENTILE_AVE_PREV_SPEARMAN",
                  "PERCENTILE_AVE_PHENO_SLOPE","PERCENTILE_AVE_PHENO_SPEARMAN")
getSensitivityScores <- function(score_df) {
  reduced_score_df <- score_df %>% select(c('PRS_ID',all_of(metric_names)))
  sensitivity_scores <- rowSums(reduced_score_df[,-1])
  return(data.frame(PRS_ID=reduced_score_df[['PRS_ID']],
                    SENSITIVITY=sensitivity_scores))
}

getSpecificScores <- function(score_df, choice) {
  reduced_score_df <- score_df %>% select(c('PRS_ID',all_of(metric_names)))
  assertthat::assert_that(choice %in% c('min','max',metric_names),
                          msg='Specify a valid choice of eval metric, or min/max')
  if (choice=='min') {
    sensitivity_scores <- apply(reduced_score_df[,-1],1,min)
  } else if (choice=='max') {
    sensitivity_scores <- apply(reduced_score_df[,-1],1,max)
  } else {
    sensitivity_scores <- reduced_score_df[[choice]]
  }
  return(data.frame(PRS_ID=reduced_score_df[['PRS_ID']],
                    SENSITIVITY=sensitivity_scores))
}

max_shuffle_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
max_signflip_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
usual_shuffle_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
usual_signflip_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))

max_shuffle_min_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
max_signflip_min_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
usual_shuffle_min_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))
usual_signflip_min_sensitivity_score <- data.frame(PRS_ID=c(prs.ids,'Lenient','Stringent'))

for (sig.cutoff in c('1e-6','1e-7','1e-8','1e-10')) {
  # Load files
  shuffle_max_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
                                                     sig.cutoff,'_shuffle_max_combined.csv'),
                                              show_col_types = FALSE)
  signflip_max_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
                                                 sig.cutoff,'_signflip_max_combined.csv'),
                                          show_col_types = FALSE)
  shuffle_usual_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
                                                 sig.cutoff,'_shuffle_usual_combined.csv'),
                                          show_col_types = FALSE)
  signflip_usual_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
                                                 sig.cutoff,'_signflip_usual_combined.csv'),
                                          show_col_types = FALSE)
  if (sig.cutoff=='1e-10') {
    shuffle_max_score_df$PRS_ID[16] <- 'Stringent'
    signflip_max_score_df$PRS_ID[16] <- 'Stringent'
    shuffle_usual_score_df$PRS_ID[16] <- 'Stringent'
    signflip_usual_score_df$PRS_ID[16] <- 'Stringent'
  } else {
    shuffle_max_score_df$PRS_ID[16] <- 'Lenient'
    signflip_max_score_df$PRS_ID[16] <- 'Lenient'
    shuffle_usual_score_df$PRS_ID[16] <- 'Lenient'
    signflip_usual_score_df$PRS_ID[16] <- 'Lenient'
  }
  # Add up sensitivities across all metrics
  shuffle_max_nc <- getSensitivityScores(shuffle_max_score_df); colnames(shuffle_max_nc)[2] <- paste0('CUTOFF_',sig.cutoff)
  signflip_max_nc <- getSensitivityScores(signflip_max_score_df); colnames(signflip_max_nc)[2] <- paste0('CUTOFF_',sig.cutoff)
  shuffle_usual_nc <- getSensitivityScores(shuffle_usual_score_df); colnames(shuffle_usual_nc)[2] <- paste0('CUTOFF_',sig.cutoff)
  signflip_usual_nc <- getSensitivityScores(signflip_usual_score_df); colnames(signflip_usual_nc)[2] <- paste0('CUTOFF_',sig.cutoff)
  
  # Collect scores for specific metric
  shuffle_max_min_score <- getSpecificScores(shuffle_max_score_df, choice='min'); colnames(shuffle_max_min_score)[2] <- paste0('CUTOFF_',sig.cutoff)
  signflip_max_min_score <- getSpecificScores(signflip_max_score_df, choice='min'); colnames(signflip_max_min_score)[2] <- paste0('CUTOFF_',sig.cutoff)
  shuffle_usual_min_score <- getSpecificScores(shuffle_usual_score_df, choice='min'); colnames(shuffle_usual_min_score)[2] <- paste0('CUTOFF_',sig.cutoff)
  signflip_usual_min_score <- getSpecificScores(signflip_usual_score_df, choice='min'); colnames(signflip_usual_min_score)[2] <- paste0('CUTOFF_',sig.cutoff)
  
  # Merge
  max_shuffle_sensitivity_score <- left_join(max_shuffle_sensitivity_score, shuffle_max_nc,by='PRS_ID')
  max_signflip_sensitivity_score <- left_join(max_signflip_sensitivity_score, signflip_max_nc,by='PRS_ID')
  usual_shuffle_sensitivity_score <- left_join(usual_shuffle_sensitivity_score, shuffle_usual_nc,by='PRS_ID')
  usual_signflip_sensitivity_score <- left_join(usual_signflip_sensitivity_score, signflip_usual_nc,by='PRS_ID')
  
  max_shuffle_min_sensitivity_score <- left_join(max_shuffle_min_sensitivity_score, shuffle_max_min_score,by='PRS_ID')
  max_signflip_min_sensitivity_score <- left_join(max_signflip_min_sensitivity_score, signflip_max_min_score,by='PRS_ID')
  usual_shuffle_min_sensitivity_score <- left_join(usual_shuffle_min_sensitivity_score, shuffle_usual_min_score,by='PRS_ID')
  usual_signflip_min_sensitivity_score <- left_join(usual_signflip_min_sensitivity_score, signflip_usual_min_score,by='PRS_ID')
}

readr::write_csv(max_shuffle_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_max_shuffle_sensitivity_sum.csv'))
readr::write_csv(max_signflip_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_max_signflip_sensitivity_sum.csv'))
readr::write_csv(usual_shuffle_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_usual_shuffle_sensitivity_sum.csv'))
readr::write_csv(usual_signflip_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_usual_signflip_sensitivity_sum.csv'))

readr::write_csv(max_shuffle_min_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_max_shuffle_sensitivity_min.csv'))
readr::write_csv(max_signflip_min_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_max_signflip_sensitivity_min.csv'))
readr::write_csv(usual_shuffle_min_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_usual_shuffle_sensitivity_min.csv'))
readr::write_csv(usual_signflip_min_sensitivity_score,
                 file = paste0(out.dir,'tables/prs_usual_signflip_sensitivity_min.csv'))

to_join <- rbind(to_join, data.frame(PRS_ID='Lenient',N_SNPS=372))

# Plot B: Sensitivity Score vs no. variants
sensitivity_vs_no_snps <- merge(max_shuffle_min_sensitivity_score %>% select(c('PRS_ID','CUTOFF_1e-6')),
                                to_join,
                                by='PRS_ID')
colnames(sensitivity_vs_no_snps)[2] <- 'AGGREGATED'
library(ggrepel)
plot_B <- ggplot(sensitivity_vs_no_snps,aes(y=AGGREGATED,
                                      x=log10(N_SNPS),
                                      label=PRS_ID)) +
  geom_point() +
  #geom_text(hjust=0,vjust=0) +
  geom_label_repel(aes(label = PRS_ID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   alpha=0.5,
                   segment.color = 'grey50') +
  theme_bw() +
  xlab(expression(paste(log[10],'(No. Variants)'))) +
  xlim(c(2,7)) +
  ylim(c(0.75,1)) +
  ylab('Minimum Sensitivity Degree Across Metrics') +
  ggtitle('B. Sensitivity of MCH PRSs') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16))

plots_A_and_B <- gridExtra::grid.arrange(plot_A,plot_B,ncol=2)
# Plot sensitivity curves
melted_max_shuffle_score <- max_shuffle_sensitivity_score[1:15,] %>% reshape2::melt(id.vars='PRS_ID')
melted_max_signflip_score <- max_signflip_sensitivity_score[1:15,] %>% reshape2::melt(id.vars='PRS_ID')
#install.packages("ggsci")
library(ggsci)
# Plot C: Shuffle Sensitivity Curves 
plot_C <- ggplot(melted_max_shuffle_score,
       aes(x=variable,
           y=value,
           colour=PRS_ID,
           group=PRS_ID)) +
  geom_point(alpha=0.8) +
  geom_line(alpha=0.7) +
  theme_bw() +
  ylab('Aggregate Sensitivity Degree Across Metrics') +
  ylim(c(7,11)) +
  xlab('Perturbation p-value Cutoff') +
  scale_colour_igv() +
  labs(colour="PRS") +
  scale_x_discrete(labels=c(expression(10^-6),expression(10^-7),expression(10^-8),expression(10^-10))) +
  ggtitle('C. Shuffle Sensitivity Diagnostic Curve') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16))

# Plot D: Signflip Sensitivity Curves 
plot_D <- ggplot(melted_max_signflip_score,
                 aes(x=variable,
                     y=value,
                     colour=PRS_ID,
                     group=PRS_ID)) +
  geom_point(alpha=0.8) +
  geom_line(alpha=0.7) +
  theme_bw() +
  ylab('Aggregate Sensitivity Degree Across Metrics') +
  ylim(c(7,11)) +
  xlab('Perturbation p-value Cutoff') +
  scale_colour_igv() +
  labs(colour="PRS") +
  scale_x_discrete(labels=c(expression(10^-6),expression(10^-7),expression(10^-8),expression(10^-10))) +
  ggtitle('D. Sign Flipping Sensitivity Diagnostic Curve') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14.5))


# Code for grabbing legend of plot
getLegend <- function(plot) {

  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot))

  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")

  # extract legend
  legend <- plot_table$grobs[[legend_plot]]

  # return legend
  return(legend)
}

plot.legend <- getLegend(plot_D)

two.plots <- gridExtra::grid.arrange(plot_C + theme(legend.position="None"), 
                                     plot_D + theme(legend.position="None")+ylab(''), ncol=2, widths = c(5,4.9))
two.plots.with.legend <- gridExtra::grid.arrange(two.plots,plot.legend,ncol=2,widths=c(9.9,1))

final_plot <- gridExtra::grid.arrange(plots_A_and_B,
                                      two.plots.with.legend,
                                      nrow=2)

ggsave(final_plot,
       filename = "/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/figures/MainText_Fig5_1e-6B.jpg",
       width = 16, height = 15,
       dpi = 300)
# Plot metrics variation 
#sig.cutoff <- '1e-8'
#shuffle_max_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
#                                               sig.cutoff,'_shuffle_max_combined.csv'),
#                                        show_col_types = FALSE)
#signflip_max_score_df <- readr::read_csv(paste0('/deep_learning/aaw/051723/Mean_Corpuscular_Hemoglobin/results/strat_vs_perturb_sensitivity/raw_tables/cutoff',
#                                                sig.cutoff,'_signflip_max_combined.csv'),
#                                         show_col_types = FALSE)
#shuffle_max_rel_metrics <- shuffle_max_score_df %>% select(c('PRS_ID',all_of(metric_names)))
#signflip_max_rel_metrics <- signflip_max_score_df %>% select(c('PRS_ID',all_of(metric_names)))
