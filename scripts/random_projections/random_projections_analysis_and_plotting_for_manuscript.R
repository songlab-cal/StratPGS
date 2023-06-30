#####################################################
########## Analysis of Random Projections ###########
#####################################################
## Analysis of Random Projections ----------------------------------------------
R.workbench <- FALSE

# Define directories 
if (R.workbench) {
  analysis.dir <- '/deep_learning/aaw/050123/random-projections/'
  save.folder <- '/deep_learning/aaw/050123/random-projections/results/'
} else {
  analysis.dir <- '/illumina/scratch/deep_learning/aaw/050123/random-projections/'
  save.folder <- '/illumina/scratch/deep_learning/aaw/050123/random-projections/results/'
}

## Load libraries, check memory ------------------------------------------------
# Report available RAM
message(date(), "Available memory: ", memuse::Sys.meminfo()$freeram)

library(dplyr)
library(tidyr)
library(reshape2)
library(readr)
library(doParallel)

## Main Body -------------------------------------------------------------------
# Plotting against pop struct bias and risk of spurious signal
base.df.2 <- readr::read_csv('/deep_learning/aaw/112122/R2_metrics_112822.csv')
gPC.metrics.df <- readr::read_csv('/deep_learning/aaw/110122/gpc_metrics_10pct.csv')
merged.df <- merge(gPC.metrics.df, base.df.2, by = 'PHENOTYPE')

raw.orig.ids <- as.vector(sapply(merged.df$PHENOTYPE, function(x) {grepl("original.RAW",x)}))
irnt.orig.ids <- as.vector(sapply(merged.df$PHENOTYPE, function(x) {grepl("original.IRNT",x)}))

merged.raw.orig.df <- merged.df[raw.orig.ids,]
merged.irnt.orig.df <- merged.df[irnt.orig.ids,]
merged.raw.orig.1stvisit.df <- merged.raw.orig.df[which(!sapply(merged.raw.orig.df$PHENOTYPE, 
                                                                function(x) grepl('2nd_visit',x))),]

readr::write_csv(merged.df,
                 file = "/deep_learning/aaw/050123/random-projections/results/all_phenos_metrics.csv")
readr::write_csv(merged.raw.orig.1stvisit.df %>%
                   select(c('PHENOTYPE',
                            'PC1_ABS_COS_SIM','PC1_ABS_PEARSON',
                            'RANK_PC1_COS_SIM','RANK_PC1_PEARSON',
                            'EVEN_COS_SIM','EVEN_PEARSON',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/ORIG_1stVisit_phenos_metrics.csv")
readr::write_csv(merged.irnt.orig.df %>%
                   select(c('PHENOTYPE',
                            'PC1_ABS_COS_SIM','PC1_ABS_PEARSON',
                            'RANK_PC1_COS_SIM','RANK_PC1_PEARSON',
                            'EVEN_COS_SIM','EVEN_PEARSON',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/IRNT_phenos_metrics.csv")

# Get p-values and correlations
# I follow code in analysis_perturb_gwas1e-5_cutoff1e-7.R (lines 510 onward)
# to generate matrices of p-values as well as raw correlations

# For original 1st visit phenos (141) -- OLD VERSION, UNCOMMENT FOR USE
#strat_metrics <- c('PC1_ABS_COS_SIM','PC1_ABS_PEARSON','RANK_PC1_COS_SIM','RANK_PC1_PEARSON','EVEN_COS_SIM','EVEN_PEARSON')
#non_random_perf_metrics <- c('VAR_PP_CORS','MAX_ABS_PP_CORS','MEAN_ABS_PP_CORS','VAR_BETAS','MAX_ABS_BETAS','MEAN_ABS_BETAS',
#                             'MEAN.DELTA.R2','MEAN.DELTA.ADJ.R2')
# For original 1st visit phenos (141) -- MODIFIED 6/22/23
strat_metrics <- c('PC1_ABS_COS_SIM','PC1_ABS_PEARSON',
                   'RANK_PC1_COS_SIM','RANK_PC1_PEARSON',
                   'EVEN_COS_SIM','EVEN_PEARSON')
non_random_perf_metrics <- c('VAR_PP_CORS','MAX_ABS_PP_CORS','MEAN_ABS_PP_CORS','MEAN.DELTA.R2')
stratMetric.vs.nonRandomPerf.array <- matrix(NA,nrow=length(strat_metrics),ncol=length(non_random_perf_metrics))
stratMetric.vs.nonRandomPerf.pval.array <- matrix(NA,nrow=length(strat_metrics),ncol=length(non_random_perf_metrics))

stratify.metrics.df4 <- merged.raw.orig.1stvisit.df[,strat_metrics]
stratify.metrics.df5 <- merged.raw.orig.1stvisit.df[,non_random_perf_metrics]

for (i in 1:nrow(stratMetric.vs.nonRandomPerf.array)) {
  for (j in 1:ncol(stratMetric.vs.nonRandomPerf.array)) {
    stratMetric.vs.nonRandomPerf.array[i,j] <- cor(stratify.metrics.df4[,i],stratify.metrics.df5[,j])
    stratMetric.vs.nonRandomPerf.pval.array[i,j] <- cor.test(stratify.metrics.df4[,i],stratify.metrics.df5[,j])$p.value
  }
}

stratMetric.vs.nonRandomPerf.array <- as.data.frame(stratMetric.vs.nonRandomPerf.array)
colnames(stratMetric.vs.nonRandomPerf.array) <- non_random_perf_metrics
stratMetric.vs.nonRandomPerf.array$METRIC <- strat_metrics
readr::write_csv(stratMetric.vs.nonRandomPerf.array %>%
                   select(c('METRIC','VAR_PP_CORS',
                            'MAX_ABS_PP_CORS','MEAN_ABS_PP_CORS',
                            'MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/pvals_and_cors/corr_ORIG_1stVisit_strat_vs_predictability.csv")

stratMetric.vs.nonRandomPerf.pval.array <- as.data.frame(stratMetric.vs.nonRandomPerf.pval.array)
colnames(stratMetric.vs.nonRandomPerf.pval.array) <- non_random_perf_metrics
stratMetric.vs.nonRandomPerf.pval.array$METRIC <- strat_metrics
readr::write_csv(stratMetric.vs.nonRandomPerf.pval.array %>% 
                   select(c('METRIC','VAR_PP_CORS',
                            'MAX_ABS_PP_CORS','MEAN_ABS_PP_CORS',
                            'MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/pvals_and_cors/pvals_ORIG_1stVisit_strat_vs_predictability.csv")

bool <- as.data.frame(stratMetric.vs.nonRandomPerf.pval.array < 0.05/24)
bool$METRIC <- strat_metrics

# For IRNT phenos (141)
stratMetric.vs.nonRandomPerf.array <- matrix(NA,nrow=length(strat_metrics),ncol=length(non_random_perf_metrics))
stratMetric.vs.nonRandomPerf.pval.array <- matrix(NA,nrow=length(strat_metrics),ncol=length(non_random_perf_metrics))

stratify.metrics.df6 <- merged.irnt.orig.df[,strat_metrics]
stratify.metrics.df7 <- merged.irnt.orig.df[,non_random_perf_metrics]

for (i in 1:nrow(stratMetric.vs.nonRandomPerf.array)) {
  for (j in 1:ncol(stratMetric.vs.nonRandomPerf.array)) {
    stratMetric.vs.nonRandomPerf.array[i,j] <- cor(stratify.metrics.df6[,i],stratify.metrics.df7[,j])
    stratMetric.vs.nonRandomPerf.pval.array[i,j] <- cor.test(stratify.metrics.df6[,i],stratify.metrics.df7[,j])$p.value
  }
}

stratMetric.vs.nonRandomPerf.array <- as.data.frame(stratMetric.vs.nonRandomPerf.array)
colnames(stratMetric.vs.nonRandomPerf.array) <- non_random_perf_metrics
stratMetric.vs.nonRandomPerf.array$METRIC <- strat_metrics
readr::write_csv(stratMetric.vs.nonRandomPerf.array %>%
                   select(c('METRIC',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/pvals_and_cors/corr_IRNT_strat_vs_predictability.csv")

stratMetric.vs.nonRandomPerf.pval.array <- as.data.frame(stratMetric.vs.nonRandomPerf.pval.array)
colnames(stratMetric.vs.nonRandomPerf.pval.array) <- non_random_perf_metrics
stratMetric.vs.nonRandomPerf.pval.array$METRIC <- strat_metrics
readr::write_csv(stratMetric.vs.nonRandomPerf.pval.array %>%
                   select(c('METRIC',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN.DELTA.R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/pvals_and_cors/pvals_IRNT_strat_vs_predictability.csv")

bool <- as.data.frame(stratMetric.vs.nonRandomPerf.pval.array < 0.05/24)
bool$METRIC <- strat_metrics
bool

# Get randomization results for each phenotype
randomization.betas.dir <- '/deep_learning/aaw/110122/'
randomization.R2.dir <- '/deep_learning/aaw/010623/'

# These dataframes show the number of times, out of 100, that the PRS relationship with
# original phenotype strictly beats the relationship with a randomized phenotype 
# w.r.t. a particular metric 
betas.randomization.df <- readr::read_csv(paste0(randomization.betas.dir,
                                                 'less_pval_100randomization_pvals.csv'),
                                          show_col_types = FALSE)
R2.randomization.df <- readr::read_csv(paste0(randomization.R2.dir,
                                              'less_pval_100randomization_pvals.csv'),
                                       show_col_types = FALSE)
fixed.betas.randomization <- apply(betas.randomization.df[,-1],2,function(x){
  (to.return <- (x*101 - 1)/100)
  ifelse(abs(to.return) < 1e-10,0,to.return)}) %>% as.data.frame()
fixed.betas.randomization$PHENOTYPE <- betas.randomization.df$PHENOTYPE

fixed.R2.randomization <- apply(R2.randomization.df[,-1],2,function(x){
  (to.return <- (x*101 - 1)/100)
  ifelse(abs(to.return) < 1e-10,0,to.return)}) %>% as.data.frame()
fixed.R2.randomization$PHENOTYPE <- R2.randomization.df$PHENOTYPE

merged.df <- merge(fixed.betas.randomization, fixed.R2.randomization,
                   by='PHENOTYPE')

raw.orig.ids <- as.vector(sapply(merged.df$PHENOTYPE, function(x) {grepl("original.RAW",x)}))
irnt.orig.ids <- as.vector(sapply(merged.df$PHENOTYPE, function(x) {grepl("original.IRNT",x)}))

merged.raw.orig.df <- merged.df[raw.orig.ids,]
merged.irnt.orig.df <- merged.df[irnt.orig.ids,]
merged.raw.orig.1stvisit.df <- merged.raw.orig.df[which(!sapply(merged.raw.orig.df$PHENOTYPE, 
                                                                function(x) grepl('2nd_visit',x))),]

readr::write_csv(merged.df,
                 file = "/deep_learning/aaw/050123/random-projections/results/predictability_inflation/all_phenos_inflation_metrics.csv")
readr::write_csv(merged.raw.orig.1stvisit.df %>%
                   select(c('PHENOTYPE',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN_DELTA_R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/predictability_inflation/ORIG_1stVisit_phenos_inflation_metrics.csv")
readr::write_csv(merged.irnt.orig.df %>%
                   select(c('PHENOTYPE',
                            'VAR_PP_CORS','MAX_ABS_PP_CORS',
                            'MEAN_ABS_PP_CORS','MEAN_DELTA_R2')),
                 file = "/deep_learning/aaw/050123/random-projections/results/predictability_inflation/IRNT_phenos_inflation_metrics.csv")

# Raw Original
raw.orig.ids <- as.vector(sapply(fixed.betas.randomization$PHENOTYPE, function(x) {grepl("original.RAW",x)}))
fixed.betas.randomization.raw.orig <- fixed.betas.randomization[raw.orig.ids,]
fixed.betas.randomization.raw.orig.1stvisit.df <- fixed.betas.randomization.raw.orig[which(!sapply(fixed.betas.randomization.raw.orig$PHENOTYPE, 
                                                                                                   function(x) grepl('2nd_visit',x))),]
for (i in 1:6) {
  metric.name <- colnames(fixed.betas.randomization.raw.orig.1stvisit.df)[i]
  message(date(), ": Results for ", metric.name)
  message("No. phenotypes with > 90 = ", sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.9))
  message("No. phenotypes with > 95 = ", sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.95))
  message("No. phenotypes with > 99 = ", sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.99))
}



raw.orig.ids <- as.vector(sapply(fixed.R2.randomization$PHENOTYPE, function(x) {grepl("original.RAW",x)}))
irnt.orig.ids <- as.vector(sapply(fixed.R2.randomization$PHENOTYPE, function(x) {grepl("original.IRNT",x)}))
fixed.R2.randomization.raw.orig <- fixed.R2.randomization[raw.orig.ids,];fixed.R2.randomization.irnt.orig <- fixed.R2.randomization[irnt.orig.ids,]
fixed.R2.randomization.raw.orig.1stvisit.df <- fixed.R2.randomization.raw.orig[which(!sapply(fixed.R2.randomization.raw.orig$PHENOTYPE, 
                                                                                             function(x) grepl('2nd_visit',x))),]
for (i in 1:3) {
  metric.name <- colnames(fixed.R2.randomization.raw.orig.1stvisit.df)[i]
  message(date(), ": Results for ", metric.name)
  message("No. phenotypes with > 90 = ", sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.9))
  message("No. phenotypes with > 95 = ", sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.95))
  message("No. phenotypes with > 99 = ", sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.99))
}

raw.orig.1stvisit.summary.df <- data.frame(METRIC = character(),
                                           BEAT_90 = numeric(),
                                           BEAT_95 = numeric(),
                                           BEAT_99 = numeric())
for (i in 1:9) {
  if (i <= 6) {
    metric.name <- colnames(fixed.betas.randomization.raw.orig.1stvisit.df)[i]
    raw.orig.1stvisit.summary.df <- rbind(raw.orig.1stvisit.summary.df,
                                          data.frame(METRIC = metric.name,
                                                     BEAT_90 = sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.9),
                                                     BEAT_95 = sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.95),
                                                     BEAT_99 = sum(fixed.betas.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.99)))
  } else {
    metric.name <- colnames(fixed.R2.randomization.raw.orig.1stvisit.df)[i-6]
    raw.orig.1stvisit.summary.df <- rbind(raw.orig.1stvisit.summary.df,
                                          data.frame(METRIC = metric.name,
                                                     BEAT_90 = sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.9),
                                                     BEAT_95 = sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.95),
                                                     BEAT_99 = sum(fixed.R2.randomization.raw.orig.1stvisit.df[[metric.name]] >= 0.99)))
  }
}

readr::write_csv(raw.orig.1stvisit.summary.df,
                 file = '/deep_learning/aaw/050123/random-projections/tables/raw_orig_1stvisit_randomization_summary_df.csv')

#[!]
fig1_plotA <- ggplot(fixed.R2.randomization.raw.orig.1stvisit.df,aes(x=MEAN_DELTA_R2*100)) +
  geom_histogram(fill='gray',colour='black') +
  theme_bw() +
  xlab('No. Randomized Phenotypes with Lower Predictability\nby Random Projections than Original') +
  ylab('count') +
  geom_vline(xintercept = 95, lty = 'dashed') +
  ggtitle('A. Distribution of Predictability Inflation\nAcross Phenotypes') +
  theme(plot.title = element_text(face='bold'))

plotB.df <- raw.orig.1stvisit.summary.df[4:7,] %>% 
  mutate(METRIC = c('Var(Cor(Prev,Perc))',
                    'Max(Cor(Prev,Perc))',
                    'Mean(Cor(Prev,Perc))',
                    'Mean(Incre. R2)')) %>% 
  melt(id.vars = 'METRIC')

#[!]
fig1_plotB <- ggplot(plotB.df, aes(x=METRIC,y=value)) +
  geom_bar(stat = 'identity', 
           aes(fill = variable),
           position = position_dodge2(preserve = "single")) +
  theme_bw() +
  #ylab('No. Phenotypes') +
  scale_y_continuous(name='No. Phenotypes',
                     breaks=seq(from=0,to=140,by=20),
                     labels=seq(from=0,to=140,by=20),
                     limits=c(0,141)) +
  xlab('') + 
  ggtitle('B. Inflation Rates by Predictability Measure') +
  scale_fill_discrete(name = "Inflation\nCutoff",
                      labels=c('> 90', '> 95', '> 99')) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5),
        plot.title = element_text(face='bold'))

figC_plotting_df <- merged.raw.orig.1stvisit.df %>% 
  select(c('PHENOTYPE',
           'PC1_ABS_COS_SIM',
           'EVEN_COS_SIM',
           'MEAN.DELTA.R2')) %>% 
  merge(x=fixed.R2.randomization.raw.orig.1stvisit.df,by='PHENOTYPE')

#[!]
fig1_plotC <- ggplot(figC_plotting_df,aes(x=PC1_ABS_COS_SIM,y=MEAN.DELTA.R2)) +
  geom_point(aes(colour = MEAN_DELTA_R2)) +
  theme_bw() +
  ylab(expression(paste('Mean Incremental ', R^2))) +
  xlab('Unsigned Cosine Similarity\nof Phenotype with PC1') +
  stat_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE,
              alpha=0.5) +
  ggtitle('C.') +
  guides(colour=guide_colourbar(title="Inflation\nDegree")) +
  scale_colour_gradient(low = "#fdbb84", high = "#990000") +
  theme(plot.title = element_text(face='bold'))

#[!]
fig1_plotD <- ggplot(figC_plotting_df,aes(x=EVEN_COS_SIM,y=MEAN.DELTA.R2)) +
  geom_point(aes(colour = MEAN_DELTA_R2)) +
  theme_bw() +
  ylab(expression(paste('Mean Incremental ', R^2))) +
  xlab('Entropy of Unsigned Cosine Similarity\nof Phenotype with 40 PCs') +
  stat_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE,
              alpha=0.5) +
  ggtitle('D.') +
  guides(colour=guide_colourbar(title="Inflation\nDegree")) +
  scale_colour_gradient(low = "#fdbb84", high = "#990000") +
  theme(plot.title = element_text(face='bold'))

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

fig1_bottom_shared_legend <- extract_legend(fig1_plotC)

library(gridExtra)
fig1_bottom <- grid.arrange(arrangeGrob(fig1_plotC+theme(legend.position='None'), 
                                        fig1_plotD+theme(legend.position='None')+ylab(''), 
                                        ncol = 2, widths = c(5.5,5.4)),
                            fig1_bottom_shared_legend, nrow = 1, widths = c(10,1))

fig1_top <- grid.arrange(fig1_plotA,fig1_plotB,ncol=2)
fig1_combined <- gridExtra::grid.arrange(fig1_top,fig1_bottom, nrow = 2)
ggsave(fig1_combined,
       filename = paste0('/deep_learning/aaw/050123/random-projections/plots/PRS_paper_fig1.jpg'),
       width = 9.2, height = 9.2,
       dpi = 300)

# IRNT Original
irnt.orig.ids <- as.vector(sapply(fixed.betas.randomization$PHENOTYPE, function(x) {grepl("original.IRNT",x)}))
fixed.betas.randomization.irnt.orig <- fixed.betas.randomization[irnt.orig.ids,]

for (i in 1:6) {
  metric.name <- colnames(fixed.betas.randomization.irnt.orig)[i]
  message(date(), ": Results for ", metric.name)
  message("No. phenotypes with > 90 = ", sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.9))
  message("No. phenotypes with > 95 = ", sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.95))
  message("No. phenotypes with > 99 = ", sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.99))
}

fixed.R2.randomization <- apply(R2.randomization.df[,-1],2,function(x){
  (to.return <- (x*101 - 1)/100)
  ifelse(abs(to.return) < 1e-10,0,to.return)}) %>% as.data.frame()
fixed.R2.randomization$PHENOTYPE <- R2.randomization.df$PHENOTYPE
raw.orig.ids <- as.vector(sapply(fixed.R2.randomization$PHENOTYPE, function(x) {grepl("original.RAW",x)}))
irnt.orig.ids <- as.vector(sapply(fixed.R2.randomization$PHENOTYPE, function(x) {grepl("original.IRNT",x)}))
fixed.R2.randomization.irnt.orig <- fixed.R2.randomization[irnt.orig.ids,]

for (i in 1:3) {
  metric.name <- colnames(fixed.R2.randomization.irnt.orig)[i]
  message(date(), ": Results for ", metric.name)
  message("No. phenotypes with > 90 = ", sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.9))
  message("No. phenotypes with > 95 = ", sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.95))
  message("No. phenotypes with > 99 = ", sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.99))
}

irnt.orig.1stvisit.summary.df <- data.frame(METRIC = character(),
                                           BEAT_90 = numeric(),
                                           BEAT_95 = numeric(),
                                           BEAT_99 = numeric())
for (i in 1:9) {
  if (i <= 6) {
    metric.name <- colnames(fixed.betas.randomization.irnt.orig)[i]
    irnt.orig.1stvisit.summary.df <- rbind(irnt.orig.1stvisit.summary.df,
                                          data.frame(METRIC = metric.name,
                                                     BEAT_90 = sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.9),
                                                     BEAT_95 = sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.95),
                                                     BEAT_99 = sum(fixed.betas.randomization.irnt.orig[[metric.name]] >= 0.99)))
  } else {
    metric.name <- colnames(fixed.R2.randomization.raw.orig.1stvisit.df)[i-6]
    irnt.orig.1stvisit.summary.df <- rbind(irnt.orig.1stvisit.summary.df,
                                          data.frame(METRIC = metric.name,
                                                     BEAT_90 = sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.9),
                                                     BEAT_95 = sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.95),
                                                     BEAT_99 = sum(fixed.R2.randomization.irnt.orig[[metric.name]] >= 0.99)))
  }
}
readr::write_csv(irnt.orig.1stvisit.summary.df,
                 file = '/deep_learning/aaw/050123/random-projections/tables/irnt_orig_1stvisit_randomization_summary_df.csv')

x.vec <- rnorm(20); y.vec <- rnorm(20)
cor(x.vec,y.vec)
lm(y.vec~x.vec-1)
