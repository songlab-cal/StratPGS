## Manual plotting -------------------------------------------------------------
# prefix.dir <- '/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/raw_tables/'
# shuffle_max_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_shuffle_max_combined.csv'))
# signflip_max_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_signflip_max_combined.csv'))
# shuffle_usual_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_shuffle_usual_combined.csv'))
# signflip_usual_prs_strat <- readr::read_csv(paste0(prefix.dir,'gwas',gwas.thres,'_cutoff',sig.cutoff,'_signflip_usual_combined.csv'))
# 

# Plot A: Perturbation Sensitivity Rate by PRS Performance Metric
library(reshape2)
gwas.thres <- '1e-5'; sig.cutoff <- '1e-8'
prs_shuffle_max_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas",
                                             gwas.thres,"_cutoff",
                                             sig.cutoff,"_shuffle_max_summary_plotting.csv"))
prs_shuffle_max_df_melted <- prs_shuffle_max_df %>% melt(id.vars = 'Metric')
new_plot_A <- ggplot(prs_shuffle_max_df_melted,
                     aes(x = Metric, y = value)) +
  geom_bar(stat = 'identity',
           aes(fill = variable),
           position = "dodge2") + 
  theme_bw() +
  scale_fill_discrete("Sensitivity\nCutoff",
                      labels = c("> 90","> 95","> 99")) +
  scale_x_discrete(name ="", 
                   labels=c("CosineSim(y,PRS)", "Cor(y,PRS)", "Cor(Ave. y,Perc)",
                            "ρ(Ave. y,Perc)", "Cor(Prev,Perc)",
                            "ρ(Prev,Perc)", "ρ(y,PRS)",
                            "Top 10% OR", "Top 10% Prev", "Top 1% OR", "Top 1% Prev")) +
  ylab('No. Phenotypes Classified Sensitive') +
  ggtitle("A. Shuffle Sensitivity Rate by PRS Performance Metric")+ 
  #geom_hline(yintercept = 103, lty = 'dashed') +
  theme(plot.title = element_text(face='bold',size = 16),
        legend.position = c(0.9, 0.7),
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14,angle = 45, 
                                   vjust = 1, hjust=1)) 

# Plot B: Example phenotype with sensitivity shown w.r.t multiple performance metrics 
gwas.thres<-'1e-5'; sig.cutoff<-'1e-8'
pheno <- 'Gamma_glutamyltransferase'
pheno.perturbed.metrics.df <- readr::read_csv(paste0('/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/phenotypes/',
                                                     pheno,'_gwas_',gwas.thres,'_cutoff',sig.cutoff,'_metrics_df.csv'))
shuffle.and.orig.only.df <- pheno.perturbed.metrics.df %>% subset(TYPE != 'signflip') 
#shuffle.and.orig.only.df$TYPE <- c('original',rep('shuffled',100))
shuffle.and.orig.only.df <- shuffle.and.orig.only.df %>% dplyr::select(c('TYPE','SPEARMAN_RHO', 'TOP10PCT_OR', 'PERCENTILE_AVE_PHENO_SPEARMAN'))
colnames(shuffle.and.orig.only.df ) <- c('Type','PRS-Phenotype Spearman Rho', 'Top 10% Odds Ratio', 'Pctile-Average Phenotype Spearman Rho')

shuffle.and.orig.melted.df <- reshape2::melt(shuffle.and.orig.only.df, id.vars = 'Type')

new_plot_B <- ggplot(shuffle.and.orig.melted.df, aes(x=value,fill=factor(Type))) +
  geom_histogram(position = "dodge") +
  facet_wrap(.~variable,
             scales = "free") +
  theme_bw() +
  #xlab('Prevalence') +
  guides(fill=guide_legend("PRS Type")) +
  ggtitle('B. Relative Performance of Gamma Glutamyltransferase PRS (Original vs Shuffled)') +
  theme(legend.position = 'right',
        plot.title = element_text(face = 'bold',size=16),
        strip.text.x = element_text(size = 14.5),
        axis.title = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

# Plots C, D and E: Perturbation Sensitivity Rate by PRS Performance Metric
# new plot C
new_plot_C_df <- shuffle_max_prs_strat %>% 
  mutate(SIGNIFICANT = (COSINE_SIM > 0.95)) %>%
  select(c('SIGNIFICANT', 'PC1_ABS_COSSIM'))

absCossim.vs.cossim.shuffle.max <- ggplot(new_plot_C_df, aes(x=SIGNIFICANT,y=PC1_ABS_COSSIM)) +
  geom_violin(aes(fill = SIGNIFICANT), width=1,alpha=0.4) +
  geom_boxplot(width=0.05, alpha=0.8) +
  #geom_jitter(colour="black", size=0.4, alpha=0.9) +
  theme_bw() +
  ylab('Unsigned Cosine Similarity of PRS\nwith PC1 (Training Cohort)') +
  scale_x_discrete(name="Fraction of Shuffled PRSs\nBeaten by Original PRS (Test Cohort)",
                   labels=c('At most 95', '> 95')) +
  ggtitle('C. PC1 Stratification vs\n  PRS Sensitivity') +
  theme(plot.title = element_text(face='bold',size=16),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        legend.position='none')

# old plot C
# absCossim.vs.cossim.shuffle.max <- ggplot(shuffle_max_prs_strat %>% mutate(SIGNIFICANT = (COSINE_SIM >= 0.95)),
#                                           aes(x=PC1_ABS_COSSIM,
#                                               y=COSINE_SIM)) +
#   geom_point(aes(colour = SIGNIFICANT))+
#   scale_color_manual("Sensitive",values=c("black", "red")) +
#   theme_bw() +
#   stat_smooth(method = "lm",
#               formula = y ~ x + I(x^2) + I(x^3),
#               se = FALSE,
#               alpha=0.5) +
#   xlab('Unsigned Cosine Similarity of PRS\nwith PC1 (Training Cohort)') +
#   ylab('Fraction of Cosine Simlarities for Shuffled PRS\nBeaten by Original PRS (Test Cohort)') +
#   ggtitle('B. PRS Sensitivity vs\nPC1 Stratification') +
#   scale_y_continuous(breaks = seq(0,1.1,by=0.1),
#                      labels = seq(0,1.1,by=0.1)) +
#   scale_x_continuous(labels = c('0','0.25','0.5','0.75','1')) +
#   theme(plot.title = element_text(face='bold',size=16),
#         axis.title.x = element_text(size = 15),
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.y = element_text(size = 15),
#         legend.position='none')
#
# new plot D
new_plot_D_df <- shuffle_max_prs_strat %>% 
  mutate(SIGNIFICANT = (COSINE_SIM > 0.95)) %>%
  select(c('SIGNIFICANT', 'EVENNESS_COSSIM'))

cossimEven.vs.cossim.shuffle.max <- ggplot(new_plot_D_df, aes(x=SIGNIFICANT,y=EVENNESS_COSSIM)) +
  geom_violin(aes(fill = SIGNIFICANT), width=1,alpha=0.4) +
  geom_boxplot(width=0.05, alpha=0.8) +
  #geom_jitter(colour="black", size=0.4, alpha=0.9) +
  theme_bw() +
  ylab('Entropy of Unsigned Cosine Similarity\n of PRS with 40 PCs (Training Cohort)') +
  scale_x_discrete(name="Fraction of Shuffled PRSs\nBeaten by Original PRS (Test Cohort)",
                   labels=c('At most 95', '> 95')) +
  ggtitle('D. PC-Similarity Evenness vs\n  PRS Sensitivity') +
  theme(plot.title = element_text(face='bold',size=16),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        legend.position='none')
# old plot D
# cossimEven.vs.cossim.shuffle.max <- ggplot(shuffle_max_prs_strat %>% mutate(SIGNIFICANT = (COSINE_SIM >= 0.95)),
#                                            aes(x=EVENNESS_COSSIM,
#                                                y=COSINE_SIM)) +
#   geom_point(aes(colour = SIGNIFICANT))+
#   scale_color_manual("Sensitive",values=c("black", "red")) +
#   theme_bw() +
#   stat_smooth(method = "lm",
#               formula = y ~ x + I(x^2) + I(x^3),
#               se = FALSE,
#               alpha=0.5) +
#   xlab('Entropy of Unsigned Cosine Similarity\n of PRS with 40 PCs (Training Cohort)') +
#   ylab('Fraction of Cosine Simlarities for Shuffled PRS\nBeaten by Original PRS (Test Cohort)') +
#   ggtitle('C. PRS Sensitivity vs\nPC-Similarity Evenness') +
#   scale_y_continuous(breaks = seq(0,1.1,by=0.1),
#                      labels = seq(0,1.1,by=0.1)) +
#   theme(plot.title = element_text(face='bold',size=16),
#         legend.position = c(0.8, 0.4),
#         axis.title.x = element_text(size = 15),
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.y = element_text(size = 15))

#
# getLegend <- function(plot) {
#   
#   # get tabular interpretation of plot
#   plot_table <- ggplot_gtable(ggplot_build(plot)) 
#   
#   #  Mark only legend in plot
#   legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
#   
#   # extract legend
#   legend <- plot_table$grobs[[legend_plot]]
#   
#   # return legend
#   return(legend) 
# }
#
#plot.legend <- getLegend(absCossim.vs.cossim.shuffle.max)
two.plots <- gridExtra::grid.arrange(absCossim.vs.cossim.shuffle.max, 
                                     cossimEven.vs.cossim.shuffle.max, ncol=2, widths = c(5,5))

evenness.vs.pc1strat <- ggplot(shuffle_max_prs_strat,
                               aes(x=PC1_ABS_COSSIM,
                                   y=EVENNESS_COSSIM)) +
  geom_point(aes(colour = COSINE_SIM)) +
  theme_bw() +
  xlab('Unsigned Cosine Similarity of PRS with PC1') +
  ylab('Entropy of Unsigned Cosine Similarity\nof PRS with 40 PCs') +
  scale_color_continuous("Fraction Beaten\nby Original PRS") +
  scale_x_continuous(labels = c('0','0.25','0.5','0.75','1')) +
  ggtitle('E. Evenness vs PC1 Stratification') +
  theme(plot.title = element_text(face='bold',size = 16),
        legend.position = c(0.8, 0.7),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15))

combined.plot2 <- gridExtra::grid.arrange(two.plots, evenness.vs.pc1strat, ncol=2, widths = c(10, 5))
combined.plot3 <- gridExtra::grid.arrange(new_plot_A, new_plot_B, combined.plot2,nrow=3,heights=c(6,4.75,6))

ggsave(combined.plot3,
       filename = "/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/MainText_Fig4.jpg",
       width = 14.5, height = 14.5,
       dpi = 300)
# ggsave(combined.plot2,
#        filename = "/deep_learning/aaw/051723/results/stratification_vs_perturbation_sensitivity/gwas1e5_cutoff1e8_max_shuffle_2.jpg",
#        width = 14.5, height = 4.8,
#        dpi = 300)
# 
# cossimEven.vs.cossim.signflip.max  <- ggplot(signflip_max_prs_strat, 
#                               aes(x=EVENNESS_COSSIM,
#                                   y=COSINE_SIM)) +
#   geom_point()+
#   theme_bw() +
#   stat_smooth(method = "lm",
#               formula = y ~ x + I(x^2) + I(x^3),
#               se = FALSE,
#               alpha=0.5) +
#   ylim(c(0,1.1)) +
#   xlab('Entropy of Unsigned Cosine Similarity of PRS with 40 PCs\n(Training Cohort)') +
#   ylab('Fraction of Cosine Simlarities for Sign-flipped PRS\nBeaten by Original PRS (Test Cohort)') +
#   ggtitle('PRS Sensitivity vs PC-Similarity Evenness\n(GWAS,cutoff)=(1e-5,1e-6)') +
#   scale_y_continuous(breaks = seq(0,1.1,by=0.1)) +
#   theme(plot.title = element_text(face='bold'))
# 

# Numerical comparison across sliding cutoff thresholds, for fixed GWAS threshold = 1e-5
#
shuffle_max_1e6_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-6_shuffle_max_summary_plotting.csv"))
shuffle_max_1e7_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-7_shuffle_max_summary_plotting.csv"))
shuffle_max_1e8_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-8_shuffle_max_summary_plotting.csv"))
shuffle_max_1e10_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-8_cutoff1e-10_shuffle_max_summary_plotting.csv"))

shuffle.max.cutoff95.df <- data.frame(Metric=shuffle_max_1e6_df$Metric,
                                      SigCutoff1e6=103-shuffle_max_1e6_df$Cutoff95,
                                      SigCutoff1e7=103-shuffle_max_1e7_df$Cutoff95,
                                      SigCutoff1e8=103-shuffle_max_1e8_df$Cutoff95,
                                      SigCutoff1e8gwas1e10=68-shuffle_max_1e10_df$Cutoff95)

#
signflip_max_1e6_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-6_signflip_max_summary_plotting.csv"))
signflip_max_1e7_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-7_signflip_max_summary_plotting.csv"))
signflip_max_1e8_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-5_cutoff1e-8_signflip_max_summary_plotting.csv"))
signflip_max_1e10_df <- readr::read_csv(paste0("/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/gwas1e-8_cutoff1e-10_signflip_max_summary_plotting.csv"))

signflip.max.cutoff95.df <- data.frame(Metric=signflip_max_1e6_df$Metric,
                                       SigCutoff1e6=103-signflip_max_1e6_df$Cutoff95,
                                       SigCutoff1e7=103-signflip_max_1e7_df$Cutoff95,
                                       SigCutoff1e8=103-signflip_max_1e8_df$Cutoff95,
                                       SigCutoff1e8gwas1e10=68-signflip_max_1e10_df$Cutoff95)
shuffle.max.cutoff95.df$Metric <- c('Cor(y,PRS)', 'rho(y,PRS)', 'CosineSim(y,PRS)', 
                                    'Top 10% Prevalence', 'Top 1% Prevalence',
                                    'Top 10% Odds Ratio', 'Top 1% Odds Ratio', 
                                    'Cor(Percentile,Prevalence)', 'rho(Percentile,Prevalence)',
                                    'Cor(Percentile,Ave. y)', 'rho(Percentile, Ave. y)')
kableExtra::kable(shuffle.max.cutoff95.df)
# 2 x 2 contingency table for each metric, focus on SigCutoff1e8 and shuffle max case


