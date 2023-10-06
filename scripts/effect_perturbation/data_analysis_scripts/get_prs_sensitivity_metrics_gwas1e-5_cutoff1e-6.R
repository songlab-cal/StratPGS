#############################################
########### Shuffling Experiment ############
########## Using GWAS cutoff 1e-5 ###########
####### Sign flips and eSize Shuffles #######
#############################################
# This script analyses results across phenotypes
# and also within each phenotype. 
# Created on 4/22/23: This version takes max(f(y,PRS^+),f(y,PRS^-))
# Modified on 5/22/23: Also computes the f(y,PRS) version (aka 'usual approach')
# objects with 'usual' appended in front correspond to the f(y,PRS) version of things

## Log file / libraries --------------------------------------------------------
sink('/illumina/scratch/deep_learning/aaw/051723/logs/gwas1e-5_cutoff1e-6_sensitivity_metrics.log')
sink(stdout(), type = "message")
library(data.table)
library(dplyr)
library(bigsnpr)

## Variables -------------------------------------------------------------------
R.workbench <- FALSE
gwas.thres <- '1e-5'
sig.cutoff <- '1e-6'
message('gwas.thres = ', gwas.thres)
message('sig.cutoff = ', sig.cutoff)

if (R.workbench) {
  phenos.shortlist <- data.table::fread('/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
                                        header = FALSE)$V1
  out.dir <- '/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
  prs.orig.perturb.data.dir <- '/deep_learning/aaw/051723/results/prs/'
  phenos.df <- readr::read_csv('/deep_learning/aaw/101922/data_10pct/n487309_phenos.csv',
                               show_col_types = FALSE)
  autosome.metadata <- readr::read_csv('/deep_learning/aaw/022823/ind_autosome_metadata.csv',
                                       show_col_types = FALSE)
  phenos.subset.df <- readr::read_csv('/deep_learning/aaw/032323/phenos/phenos_subset_df.csv',
                                      show_col_types = FALSE)
} else {
  phenos.shortlist <- data.table::fread('/illumina/scratch/deep_learning/aaw/051723/results/pheno_names_gwas1e5.txt',
                                        header = FALSE)$V1
  out.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs_perturbation_sensitivity/'
  prs.orig.perturb.data.dir <- '/illumina/scratch/deep_learning/aaw/051723/results/prs/'
  phenos.df <- readr::read_csv('/illumina/scratch/deep_learning/aaw/101922/data_10pct/n487309_phenos.csv',
                               show_col_types = FALSE)
  autosome.metadata <- readr::read_csv('/illumina/scratch/deep_learning/aaw/022823/ind_autosome_metadata.csv',
                                       show_col_types = FALSE)
  phenos.subset.df <- readr::read_csv('/illumina/scratch/deep_learning/aaw/032323/phenos/phenos_subset_df.csv',
                                      show_col_types = FALSE)
}

## Main Body -------------------------------------------------------------------
test.ids <- which(autosome.metadata$TEST & autosome.metadata$EURO)
test.metadata <- autosome.metadata[test.ids,]; colnames(test.metadata)[1] <- 'sample_id'

# Left join to ensure ordering of phenotype matches ordering of genotypes
phenos.subset.df <- left_join(test.metadata,phenos.subset.df, by = 'sample_id')
leakage.ids <- which(phenos.subset.df$TEST & phenos.subset.df$TRAIN) 

# Create dataframes for saving
signflip.metrics.df <- data.frame(PHENO = character(),
                                  PEARSON_CORR = numeric(),
                                  SPEARMAN_RHO = numeric(),
                                  COSINE_SIM = numeric(),
                                  TOP10PCT_PREV = numeric(),
                                  TOP1PCT_PREV = numeric(),
                                  TOP10PCT_OR = numeric(),
                                  TOP1PCT_OR = numeric(),
                                  PERCENTILE_AVE_PREV_SLOPE = numeric(),
                                  PERCENTILE_AVE_PREV_SPEARMAN = numeric(),
                                  PERCENTILE_AVE_PHENO_SLOPE = numeric(),
                                  PERCENTILE_AVE_PHENO_SPEARMAN = numeric())
shuffle.metrics.df <- data.frame(PHENO = character(),
                                 PEARSON_CORR = numeric(),
                                 SPEARMAN_RHO = numeric(),
                                 COSINE_SIM = numeric(),
                                 TOP10PCT_PREV = numeric(),
                                 TOP1PCT_PREV = numeric(),
                                 TOP10PCT_OR = numeric(),
                                 TOP1PCT_OR = numeric(),
                                 PERCENTILE_AVE_PREV_SLOPE = numeric(),
                                 PERCENTILE_AVE_PREV_SPEARMAN = numeric(),
                                 PERCENTILE_AVE_PHENO_SLOPE = numeric(),
                                 PERCENTILE_AVE_PHENO_SPEARMAN = numeric())
usual.signflip.metrics.df <- data.frame(PHENO = character(),
                                        PEARSON_CORR = numeric(),
                                        SPEARMAN_RHO = numeric(),
                                        COSINE_SIM = numeric(),
                                        TOP10PCT_PREV = numeric(),
                                        TOP1PCT_PREV = numeric(),
                                        TOP10PCT_OR = numeric(),
                                        TOP1PCT_OR = numeric(),
                                        PERCENTILE_AVE_PREV_SLOPE = numeric(),
                                        PERCENTILE_AVE_PREV_SPEARMAN = numeric(),
                                        PERCENTILE_AVE_PHENO_SLOPE = numeric(),
                                        PERCENTILE_AVE_PHENO_SPEARMAN = numeric())
usual.shuffle.metrics.df <- data.frame(PHENO = character(),
                                       PEARSON_CORR = numeric(),
                                       SPEARMAN_RHO = numeric(),
                                       COSINE_SIM = numeric(),
                                       TOP10PCT_PREV = numeric(),
                                       TOP1PCT_PREV = numeric(),
                                       TOP10PCT_OR = numeric(),
                                       TOP1PCT_OR = numeric(),
                                       PERCENTILE_AVE_PREV_SLOPE = numeric(),
                                       PERCENTILE_AVE_PREV_SPEARMAN = numeric(),
                                       PERCENTILE_AVE_PHENO_SLOPE = numeric(),
                                       PERCENTILE_AVE_PHENO_SPEARMAN = numeric())

for (pheno in phenos.shortlist) {
  message(date(), ": Working on ", pheno)
  
  # Initialize empty pheno-specific metrics (for visualization)
  pheno.spec.df <- data.frame(TYPE = c('original',rep('shuffle',100), rep('signflip',100)),
                              PEARSON_CORR = rep(0, 201),
                              SPEARMAN_RHO = rep(0, 201),
                              COSINE_SIM = rep(0, 201),
                              TOP10PCT_PREV = rep(0, 201),
                              TOP1PCT_PREV = rep(0, 201),
                              TOP10PCT_OR = rep(0, 201),
                              TOP1PCT_OR = rep(0, 201),
                              PERCENTILE_AVE_PREV_SLOPE = rep(0, 201),
                              PERCENTILE_AVE_PREV_SPEARMAN = rep(0, 201),
                              PERCENTILE_AVE_PHENO_SLOPE = rep(0, 201),
                              PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, 201))
  usual.pheno.spec.df <- data.frame(TYPE = c('original',rep('shuffle',100), rep('signflip',100)),
                                    PEARSON_CORR = rep(0, 201),
                                    SPEARMAN_RHO = rep(0, 201),
                                    COSINE_SIM = rep(0, 201),
                                    TOP10PCT_PREV = rep(0, 201),
                                    TOP1PCT_PREV = rep(0, 201),
                                    TOP10PCT_OR = rep(0, 201),
                                    TOP1PCT_OR = rep(0, 201),
                                    PERCENTILE_AVE_PREV_SLOPE = rep(0, 201),
                                    PERCENTILE_AVE_PREV_SPEARMAN = rep(0, 201),
                                    PERCENTILE_AVE_PHENO_SLOPE = rep(0, 201),
                                    PERCENTILE_AVE_PHENO_SPEARMAN = rep(0, 201))
  
  # Load phenotype values
  pheno.name <- paste0(pheno,'.all_ethnicities.both.original.IRNT')
  pheno.vals <- phenos.subset.df[,c('sample_id',pheno.name)]
  pheno.vals <- pheno.vals[[pheno.name]]
  message(date(), ": Replacing leaked samples with NaNs...")
  pheno.vals[leakage.ids] <- NA
  message("No. NaNs in observed phenotypes = ", (is.na(pheno.vals) %>% sum()) - 5153)
  binarized.pheno.vals <- (pheno.vals > median(na.omit(pheno.vals))) 
  
  # Helper functions, depending on binarized.pheno.vals,
  # for computing prevalence 
  getPrevalence <- function(x) { # Top 10%
    cutoff <- as.numeric(quantile(na.omit(x), probs = 0.9))
    return(mean(na.omit(binarized.pheno.vals[which(x > cutoff)])))
  } 
  getPrevalenceC <- function(x) { # Top 10%
    cutoff <- as.numeric(quantile(na.omit(x), probs = 0.9))
    return(mean(na.omit(binarized.pheno.vals[which(x <= cutoff)])))
  } 
  getPrevalence2 <- function(x) { # Top 1%
    cutoff <- as.numeric(quantile(na.omit(x), probs = 0.99))
    return(mean(na.omit(binarized.pheno.vals[which(x > cutoff)])))
  } 
  getPrevalence2C <- function(x) { # Top 10%
    cutoff <- as.numeric(quantile(na.omit(x), probs = 0.99))
    return(mean(na.omit(binarized.pheno.vals[which(x <= cutoff)])))
  } 
  
  # for computing prevalence vector / average phenotype vector
  getPrevVector <- function(x) {
    cutoffs <- quantile(na.omit(x), probs = seq(0,1,by=0.01))
    cutoffs[101] <- cutoffs[101]+1 # make sure top scoring individual is included
    to.return <- c()
    for (i in 0:99) {
      to.return <- c(to.return,
                     mean(na.omit(binarized.pheno.vals[which(x >= cutoffs[[paste0(i,'%')]] & 
                                                               x < cutoffs[[paste0(i+1,'%')]])])))
    }
    return(to.return)
  }
  
  getAvePhenoVector <- function(x) {
    cutoffs <- quantile(na.omit(x), probs = seq(0,1,by=0.01))
    cutoffs[101] <- cutoffs[101]+1 # make sure top scoring individual is included
    to.return <- c()
    for (i in 0:99) {
      to.return <- c(to.return,
                     mean(na.omit(pheno.vals[which(x >= cutoffs[[paste0(i,'%')]] & 
                                                     x < cutoffs[[paste0(i+1,'%')]])])))
    }
    return(to.return)
  }
  
  # Load all perturbed results 
  # - shuffle, sign flip, inversion
  perturbed.res <- readRDS(paste0(prs.orig.perturb.data.dir, 
                                  pheno, 
                                  "_perturb_gwas",
                                  gwas.thres,
                                  "_cutoff",
                                  sig.cutoff,
                                  "_test_pred.rds"))
  
  # Load original PRS predictions (PRS_r + PRS_s)
  original.pred <- readRDS(paste0(prs.orig.perturb.data.dir, 
                                  pheno, 
                                  "_gwas_",
                                  gwas.thres,
                                  "_orig_test_train_pred.rds"))
  
  # Compute metrics
  shuffled.df <- perturbed.res$TEST_SHUFFLE; shuffled.df[leakage.ids,] <- NA
  signflip.df <- perturbed.res$TEST_SIGNFLIP; signflip.df[leakage.ids,] <- NA
  orig.test <- original.pred$TEST; orig.test[leakage.ids] <- NA
  invert.test <- perturbed.res$INVERSION; invert.test[leakage.ids] <- NA
  
  # Get PRS_s 
  PRS_s <- (orig.test + invert.test)/2

  orig.v.pheno.spearman <- cor(orig.test, pheno.vals, method = 'spearman',
                               use = "pairwise.complete.obs") 
  orig.v.pheno.pearson <- cor(orig.test, pheno.vals, method = 'pearson',
                              use = "pairwise.complete.obs") 
  orig.v.pheno.cossim <- sum(orig.test*pheno.vals,na.rm=TRUE)/(norm(na.omit(pheno.vals),'2')*norm(na.omit(orig.test),'2'))
  orig.v.pheno.prev <- getPrevalence(orig.test)
  orig.v.pheno.prev2 <- getPrevalence2(orig.test)
  orig.v.pheno.prevC <- getPrevalenceC(orig.test)
  orig.v.pheno.prev2C <- getPrevalence2C(orig.test)
  orig.top10pct.or <- orig.v.pheno.prev * (1-orig.v.pheno.prevC) / 
    (orig.v.pheno.prevC * (1-orig.v.pheno.prev))
  orig.top1pct.or <- orig.v.pheno.prev2 * (1-orig.v.pheno.prev2C) / 
    (orig.v.pheno.prev2C * (1-orig.v.pheno.prev2))
  orig.v.pheno.aveprev.pearson <- cor(seq(0,0.99,by=0.01),getPrevVector(orig.test))
  orig.v.pheno.aveprev.spearman <- cor(seq(0,0.99,by=0.01),getPrevVector(orig.test),method='spearman')
  orig.v.pheno.avepheno.pearson <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(orig.test))
  orig.v.pheno.avepheno.spearman <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(orig.test),method='spearman')
  
  # SHUFFLE --------------------------------------------------------------------
  # Need to compute max(f(PRS^-,y),f(PRS^+,y))
  # Original code computes just f(PRS^+,y)
  # To get PRS^- from PRS^+ for each perturbation df
  # PRS^- = PRS_s - hat(PRS_r) = 2*PRS_s - (PRS_s + hat(PRS_r)) = 2*PRS_s - PRS^+
  shuffled.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    minus <- cor(2*PRS_s-x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    return(max(plus,minus))
  })
  shuffled.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    minus <- cor(2*PRS_s-x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    return(max(plus,minus))
  })
  shuffled.cosine.sim.vals <- apply(shuffled.df,2,function(x) {
    plus <- sum(x*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(x),'2')))
    minus <- sum((2*PRS_s-x)*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(2*PRS_s-x),'2')))
    return(max(plus,minus))
  })
  # -----
  usual.shuffled.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    return(plus)
  })
  usual.shuffled.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    return(plus)
  })
  usual.shuffled.cosine.sim.vals <- apply(shuffled.df,2,function(x) {
    plus <- sum(x*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(x),'2')))
    return(plus)
  })
  # -----
  shuffled.prev.vals.plus <- apply(shuffled.df,2,function(x) {
    getPrevalence(x)
  })
  shuffled.prev2.vals.plus <- apply(shuffled.df,2,function(x) {
    getPrevalence2(x)
  })
  shuffled.prevC.vals.plus <- apply(shuffled.df,2,function(x) {
    getPrevalenceC(x)
  })
  shuffled.prev2C.vals.plus <- apply(shuffled.df,2,function(x) {
    getPrevalence2C(x)
  })
  shuffled.prev.vals.minus <- apply(shuffled.df,2,function(x) {
    getPrevalence(2*PRS_s-x)
  })
  shuffled.prev2.vals.minus <- apply(shuffled.df,2,function(x) {
    getPrevalence2(2*PRS_s-x)
  })
  shuffled.prevC.vals.minus <- apply(shuffled.df,2,function(x) {
    getPrevalenceC(2*PRS_s-x)
  })
  shuffled.prev2C.vals.minus <- apply(shuffled.df,2,function(x) {
    getPrevalence2C(2*PRS_s-x)
  })
  shuffled.prev.vals <- pmax(shuffled.prev.vals.plus,
                            shuffled.prev.vals.minus)
  shuffled.prev2.vals <- pmax(shuffled.prev2.vals.plus,
                             shuffled.prev2.vals.minus)
  # -----
  shuffled.top10pct.or.vals <- pmax(
    shuffled.prev.vals.plus * (1-shuffled.prevC.vals.plus) / 
    (shuffled.prevC.vals.plus * (1-shuffled.prev.vals.plus)),
    shuffled.prev.vals.minus * (1-shuffled.prevC.vals.minus) / 
      (shuffled.prevC.vals.minus * (1-shuffled.prev.vals.minus))
  )
  shuffled.top1pct.or.vals <- pmax(
    shuffled.prev2.vals.plus * (1-shuffled.prev2C.vals.plus) / 
    (shuffled.prev2C.vals.plus * (1-shuffled.prev2.vals.plus)),
    shuffled.prev2.vals.minus * (1-shuffled.prev2C.vals.minus) / 
      (shuffled.prev2C.vals.minus * (1-shuffled.prev2.vals.minus))
  )
  shuffled.aveprev.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x))
    minus <- cor(seq(0,0.99,by=0.01),getPrevVector(2*PRS_s-x))
    return(max(plus,minus))
  })
  shuffled.aveprev.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x),method='spearman')
    minus <- cor(seq(0,0.99,by=0.01),getPrevVector(2*PRS_s-x),method='spearman')
    return(max(plus,minus))
  })
  shuffled.avepheno.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x))
    minus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(2*PRS_s-x))
    return(max(plus,minus))
  })
  shuffled.avepheno.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x),method='spearman')
    minus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(2*PRS_s-x),method='spearman')
    return(max(plus,minus))
  })
  # -----
  usual.shuffled.top10pct.or.vals <- (
    shuffled.prev.vals.plus * (1-shuffled.prevC.vals.plus) / 
      (shuffled.prevC.vals.plus * (1-shuffled.prev.vals.plus))
  )
  usual.shuffled.top1pct.or.vals <- (
    shuffled.prev2.vals.plus * (1-shuffled.prev2C.vals.plus) / 
      (shuffled.prev2C.vals.plus * (1-shuffled.prev2.vals.plus))
  )
  usual.shuffled.aveprev.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x))
    return(plus)
  })
  usual.shuffled.aveprev.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x),method='spearman')
    return(plus)
  })
  usual.shuffled.avepheno.pearson.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x))
    return(plus)
  })
  usual.shuffled.avepheno.spearman.vals <- apply(shuffled.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x),method='spearman')
    return(plus)
  })
  
  # SIGNFLIP -------------------------------------------------------------------
  signflip.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    minus <- cor(2*PRS_s-x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    return(max(plus,minus))
  })
  signflip.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    minus <- cor(2*PRS_s-x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    return(max(plus,minus))
  })
  signflip.cosine.sim.vals <- apply(signflip.df,2,function(x) {
    plus <- sum(x*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(x),'2')))
    minus <- sum((2*PRS_s-x)*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(2*PRS_s-x),'2')))
    return(max(plus,minus))
  })
  # -----
  usual.signflip.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'spearman',use = "pairwise.complete.obs")
    return(plus)
  })
  usual.signflip.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(x,pheno.vals, method = 'pearson',use = "pairwise.complete.obs")
    return(plus)
  })
  usual.signflip.cosine.sim.vals <- apply(signflip.df,2,function(x) {
    plus <- sum(x*pheno.vals,na.rm=TRUE)/((norm(na.omit(pheno.vals),'2')*norm(na.omit(x),'2')))
    return(plus)
  })
  # -----
  signflip.prev.vals.plus <- apply(signflip.df,2,function(x) {
    getPrevalence(x)
  })
  signflip.prev2.vals.plus <- apply(signflip.df,2,function(x) {
    getPrevalence2(x)
  })
  signflip.prevC.vals.plus <- apply(signflip.df,2,function(x) {
    getPrevalenceC(x)
  })
  signflip.prev2C.vals.plus <- apply(signflip.df,2,function(x) {
    getPrevalence2C(x)
  })
  signflip.prev.vals.minus <- apply(signflip.df,2,function(x) {
    getPrevalence(2*PRS_s-x)
  })
  signflip.prev2.vals.minus <- apply(signflip.df,2,function(x) {
    getPrevalence2(2*PRS_s-x)
  })
  signflip.prevC.vals.minus <- apply(signflip.df,2,function(x) {
    getPrevalenceC(2*PRS_s-x)
  })
  signflip.prev2C.vals.minus <- apply(signflip.df,2,function(x) {
    getPrevalence2C(2*PRS_s-x)
  })
  signflip.prev.vals <- pmax(signflip.prev.vals.plus,
                             signflip.prev.vals.minus)
  signflip.prev2.vals <- pmax(signflip.prev2.vals.plus,
                              signflip.prev2.vals.minus)
  # -----
  signflip.top10pct.or.vals <- pmax(
    signflip.prev.vals.plus * (1-signflip.prevC.vals.plus) / 
      (signflip.prevC.vals.plus * (1-signflip.prev.vals.plus)),
    signflip.prev.vals.minus * (1-signflip.prevC.vals.minus) / 
      (signflip.prevC.vals.minus * (1-signflip.prev.vals.minus))
  )
  signflip.top1pct.or.vals <- pmax(
    signflip.prev2.vals.plus * (1-signflip.prev2C.vals.plus) / 
      (signflip.prev2C.vals.plus * (1-signflip.prev2.vals.plus)),
    signflip.prev2.vals.minus * (1-signflip.prev2C.vals.minus) / 
      (signflip.prev2C.vals.minus * (1-signflip.prev2.vals.minus))
  )
  signflip.aveprev.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x))
    minus <- cor(seq(0,0.99,by=0.01),getPrevVector(2*PRS_s-x))
    return(max(plus,minus))
  })
  signflip.aveprev.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x),method='spearman')
    minus <- cor(seq(0,0.99,by=0.01),getPrevVector(2*PRS_s-x),method='spearman')
    return(max(plus,minus))
  })
  signflip.avepheno.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x))
    minus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(2*PRS_s-x))
    return(max(plus,minus))
  })
  signflip.avepheno.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x),method='spearman')
    minus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(2*PRS_s-x),method='spearman')
    return(max(plus,minus))
  })
  # -----
  usual.signflip.top10pct.or.vals <- (
    signflip.prev.vals.plus * (1-signflip.prevC.vals.plus) / 
      (signflip.prevC.vals.plus * (1-signflip.prev.vals.plus))
  )
  usual.signflip.top1pct.or.vals <- (
    signflip.prev2.vals.plus * (1-signflip.prev2C.vals.plus) / 
      (signflip.prev2C.vals.plus * (1-signflip.prev2.vals.plus))
  )
  usual.signflip.aveprev.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x))
    return(plus)
  })
  usual.signflip.aveprev.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getPrevVector(x),method='spearman')
    return(plus)
  })
  usual.signflip.avepheno.pearson.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x))
    return(plus)
  })
  usual.signflip.avepheno.spearman.vals <- apply(signflip.df,2,function(x) {
    plus <- cor(seq(0,0.99,by=0.01),getAvePhenoVector(x),method='spearman')
    return(plus)
  })
  
  # Fill in the pheno-specific dataframe for plotting
  # -----
  pheno.spec.df$PEARSON_CORR <- c(orig.v.pheno.pearson, 
                                  shuffled.pearson.vals, 
                                  signflip.pearson.vals)
  pheno.spec.df$SPEARMAN_RHO <- c(orig.v.pheno.spearman, 
                                  shuffled.spearman.vals,
                                  signflip.spearman.vals)
  pheno.spec.df$COSINE_SIM <- c(orig.v.pheno.cossim, 
                                shuffled.cosine.sim.vals,
                                signflip.cosine.sim.vals)
  pheno.spec.df$TOP10PCT_PREV <- c(orig.v.pheno.prev,
                                   shuffled.prev.vals,
                                   signflip.prev.vals)
  pheno.spec.df$TOP1PCT_PREV <- c(orig.v.pheno.prev2,
                                  shuffled.prev2.vals,
                                  signflip.prev2.vals)
  pheno.spec.df$TOP10PCT_OR <- c(orig.top10pct.or,
                                 shuffled.top10pct.or.vals,
                                 signflip.top10pct.or.vals)
  pheno.spec.df$TOP1PCT_OR <- c(orig.top1pct.or,
                                shuffled.top1pct.or.vals,
                                signflip.top1pct.or.vals)
  pheno.spec.df$PERCENTILE_AVE_PREV_SLOPE <- c(orig.v.pheno.aveprev.pearson,
                                               shuffled.aveprev.pearson.vals,
                                               signflip.aveprev.pearson.vals)
  pheno.spec.df$PERCENTILE_AVE_PREV_SPEARMAN <- c(orig.v.pheno.aveprev.spearman,
                                                  shuffled.aveprev.spearman.vals,
                                                  signflip.aveprev.spearman.vals)
  pheno.spec.df$PERCENTILE_AVE_PHENO_SLOPE <- c(orig.v.pheno.avepheno.pearson,
                                                shuffled.avepheno.pearson.vals,
                                                signflip.avepheno.pearson.vals)
  pheno.spec.df$PERCENTILE_AVE_PHENO_SPEARMAN <- c(orig.v.pheno.avepheno.spearman,
                                                   shuffled.avepheno.spearman.vals,
                                                   signflip.avepheno.spearman.vals)
  # -----
  usual.pheno.spec.df$PEARSON_CORR <- c(orig.v.pheno.pearson, 
                                        usual.shuffled.pearson.vals, 
                                        usual.signflip.pearson.vals)
  usual.pheno.spec.df$SPEARMAN_RHO <- c(orig.v.pheno.spearman, 
                                        usual.shuffled.spearman.vals,
                                        usual.signflip.spearman.vals)
  usual.pheno.spec.df$COSINE_SIM <- c(orig.v.pheno.cossim, 
                                      usual.shuffled.cosine.sim.vals,
                                      usual.signflip.cosine.sim.vals)
  usual.pheno.spec.df$TOP10PCT_PREV <- c(orig.v.pheno.prev,
                                         shuffled.prev.vals.plus,
                                         signflip.prev.vals.plus)
  usual.pheno.spec.df$TOP1PCT_PREV <- c(orig.v.pheno.prev2,
                                        shuffled.prev2.vals.plus,
                                        signflip.prev2.vals.plus)
  usual.pheno.spec.df$TOP10PCT_OR <- c(orig.top10pct.or,
                                       usual.shuffled.top10pct.or.vals,
                                       usual.signflip.top10pct.or.vals)
  usual.pheno.spec.df$TOP1PCT_OR <- c(orig.top1pct.or,
                                      usual.shuffled.top1pct.or.vals,
                                      usual.signflip.top1pct.or.vals)
  usual.pheno.spec.df$PERCENTILE_AVE_PREV_SLOPE <- c(orig.v.pheno.aveprev.pearson,
                                                     usual.shuffled.aveprev.pearson.vals,
                                                     usual.signflip.aveprev.pearson.vals)
  usual.pheno.spec.df$PERCENTILE_AVE_PREV_SPEARMAN <- c(orig.v.pheno.aveprev.spearman,
                                                        usual.shuffled.aveprev.spearman.vals,
                                                        usual.signflip.aveprev.spearman.vals)
  usual.pheno.spec.df$PERCENTILE_AVE_PHENO_SLOPE <- c(orig.v.pheno.avepheno.pearson,
                                                      usual.shuffled.avepheno.pearson.vals,
                                                      usual.signflip.avepheno.pearson.vals)
  usual.pheno.spec.df$PERCENTILE_AVE_PHENO_SPEARMAN <- c(orig.v.pheno.avepheno.spearman,
                                                         usual.shuffled.avepheno.spearman.vals,
                                                         usual.signflip.avepheno.spearman.vals)
  # Saving metrics for individual phenotype
  message("Saving metrics for ", pheno)
  readr::write_csv(pheno.spec.df,
                   file = paste0(out.dir,'phenotypes/',pheno,'_gwas_',
                                 gwas.thres,
                                 '_cutoff',
                                 sig.cutoff,
                                 '_metrics_df.csv'))
  readr::write_csv(usual.pheno.spec.df,
                   file = paste0(out.dir,'phenotypes/',pheno,'_gwas_',
                                 gwas.thres,
                                 '_cutoff',
                                 sig.cutoff,
                                 '_usual_metrics_df.csv'))
  # (Note: these aren't p-values, but fractions of randomized PRS beaten by original!)
  #
  shuffled.spearman.corr.pval <- mean(orig.v.pheno.spearman > shuffled.spearman.vals)
  shuffled.pearson.corr.pval <- mean(orig.v.pheno.pearson > shuffled.pearson.vals)
  shuffled.cosine.sim.pval <- mean(orig.v.pheno.cossim > shuffled.cosine.sim.vals)
  shuffled.prev.pval <- mean(orig.v.pheno.prev > shuffled.prev.vals)
  shuffled.prev2.pval <- mean(orig.v.pheno.prev2 > shuffled.prev2.vals)
  shuffled.top10pct.or.pval <- mean(orig.top10pct.or > shuffled.top10pct.or.vals)
  shuffled.top1pct.or.pval <- mean(orig.top1pct.or > shuffled.top1pct.or.vals)
  shuffled.aveprev.pearson.pval <- mean(orig.v.pheno.aveprev.pearson > shuffled.aveprev.pearson.vals)
  shuffled.aveprev.spearman.pval <- mean(orig.v.pheno.aveprev.spearman > shuffled.aveprev.spearman.vals)
  shuffled.avepheno.pearson.pval <- mean(orig.v.pheno.avepheno.pearson > shuffled.avepheno.pearson.vals)
  shuffled.avepheno.spearman.pval <- mean(orig.v.pheno.avepheno.spearman > shuffled.avepheno.spearman.vals)

  #
  signflip.spearman.corr.pval <- mean(orig.v.pheno.spearman > signflip.spearman.vals)
  signflip.pearson.corr.pval <- mean(orig.v.pheno.pearson > signflip.pearson.vals)
  signflip.cosine.sim.pval <- mean(orig.v.pheno.cossim > signflip.cosine.sim.vals)
  signflip.prev.pval <- mean(orig.v.pheno.prev > signflip.prev.vals)
  signflip.prev2.pval <- mean(orig.v.pheno.prev2 > signflip.prev2.vals)
  signflip.top10pct.or.pval <- mean(orig.top10pct.or > signflip.top10pct.or.vals)
  signflip.top1pct.or.pval <- mean(orig.top1pct.or > signflip.top1pct.or.vals)
  signflip.aveprev.pearson.pval <- mean(orig.v.pheno.aveprev.pearson > signflip.aveprev.pearson.vals)
  signflip.aveprev.spearman.pval <- mean(orig.v.pheno.aveprev.spearman > signflip.aveprev.spearman.vals)
  signflip.avepheno.pearson.pval <- mean(orig.v.pheno.avepheno.pearson > signflip.avepheno.pearson.vals)
  signflip.avepheno.spearman.pval <- mean(orig.v.pheno.avepheno.spearman > signflip.avepheno.spearman.vals)
  
  # 
  usual.shuffled.spearman.corr.pval <- mean(orig.v.pheno.spearman > usual.shuffled.spearman.vals)
  usual.shuffled.pearson.corr.pval <- mean(orig.v.pheno.pearson > usual.shuffled.pearson.vals)
  usual.shuffled.cosine.sim.pval <- mean(orig.v.pheno.cossim > usual.shuffled.cosine.sim.vals)
  usual.shuffled.prev.pval <- mean(orig.v.pheno.prev > shuffled.prev.vals.plus)
  usual.shuffled.prev2.pval <- mean(orig.v.pheno.prev2 > shuffled.prev2.vals.plus)
  usual.shuffled.top10pct.or.pval <- mean(orig.top10pct.or > usual.shuffled.top10pct.or.vals)
  usual.shuffled.top1pct.or.pval <- mean(orig.top1pct.or > usual.shuffled.top1pct.or.vals)
  usual.shuffled.aveprev.pearson.pval <- mean(orig.v.pheno.aveprev.pearson > usual.shuffled.aveprev.pearson.vals)
  usual.shuffled.aveprev.spearman.pval <- mean(orig.v.pheno.aveprev.spearman > usual.shuffled.aveprev.spearman.vals)
  usual.shuffled.avepheno.pearson.pval <- mean(orig.v.pheno.avepheno.pearson > usual.shuffled.avepheno.pearson.vals)
  usual.shuffled.avepheno.spearman.pval <- mean(orig.v.pheno.avepheno.spearman > usual.shuffled.avepheno.spearman.vals)
  
  #
  usual.signflip.spearman.corr.pval <- mean(orig.v.pheno.spearman > usual.signflip.spearman.vals)
  usual.signflip.pearson.corr.pval <- mean(orig.v.pheno.pearson > usual.signflip.pearson.vals)
  usual.signflip.cosine.sim.pval <- mean(orig.v.pheno.cossim > usual.signflip.cosine.sim.vals)
  usual.signflip.prev.pval <- mean(orig.v.pheno.prev > signflip.prev.vals.plus)
  usual.signflip.prev2.pval <- mean(orig.v.pheno.prev2 > signflip.prev2.vals.plus)
  usual.signflip.top10pct.or.pval <- mean(orig.top10pct.or > usual.signflip.top10pct.or.vals)
  usual.signflip.top1pct.or.pval <- mean(orig.top1pct.or > usual.signflip.top1pct.or.vals)
  usual.signflip.aveprev.pearson.pval <- mean(orig.v.pheno.aveprev.pearson > usual.signflip.aveprev.pearson.vals)
  usual.signflip.aveprev.spearman.pval <- mean(orig.v.pheno.aveprev.spearman > usual.signflip.aveprev.spearman.vals)
  usual.signflip.avepheno.pearson.pval <- mean(orig.v.pheno.avepheno.pearson > usual.signflip.avepheno.pearson.vals)
  usual.signflip.avepheno.spearman.pval <- mean(orig.v.pheno.avepheno.spearman > usual.signflip.avepheno.spearman.vals)
  
  # Add to summary dataframes
  shuffle.metrics.df <- rbind(shuffle.metrics.df,
                              data.frame(PHENO = pheno,
                                         PEARSON_CORR = shuffled.pearson.corr.pval,
                                         SPEARMAN_RHO = shuffled.spearman.corr.pval,
                                         COSINE_SIM = shuffled.cosine.sim.pval,
                                         TOP10PCT_PREV = shuffled.prev.pval,
                                         TOP1PCT_PREV = shuffled.prev2.pval,
                                         TOP10PCT_OR = shuffled.top10pct.or.pval,
                                         TOP1PCT_OR = shuffled.top1pct.or.pval,
                                         PERCENTILE_AVE_PREV_SLOPE = shuffled.aveprev.pearson.pval,
                                         PERCENTILE_AVE_PREV_SPEARMAN = shuffled.aveprev.spearman.pval,
                                         PERCENTILE_AVE_PHENO_SLOPE = shuffled.avepheno.pearson.pval,
                                         PERCENTILE_AVE_PHENO_SPEARMAN = shuffled.avepheno.spearman.pval))
  signflip.metrics.df <- rbind(signflip.metrics.df,
                               data.frame(PHENO = pheno,
                                          PEARSON_CORR = signflip.pearson.corr.pval,
                                          SPEARMAN_RHO = signflip.spearman.corr.pval,
                                          COSINE_SIM = signflip.cosine.sim.pval,
                                          TOP10PCT_PREV = signflip.prev.pval,
                                          TOP1PCT_PREV = signflip.prev2.pval,
                                          TOP10PCT_OR = signflip.top10pct.or.pval,
                                          TOP1PCT_OR = signflip.top1pct.or.pval,
                                          PERCENTILE_AVE_PREV_SLOPE = signflip.aveprev.pearson.pval,
                                          PERCENTILE_AVE_PREV_SPEARMAN = signflip.aveprev.spearman.pval,
                                          PERCENTILE_AVE_PHENO_SLOPE = signflip.avepheno.pearson.pval,
                                          PERCENTILE_AVE_PHENO_SPEARMAN = signflip.avepheno.spearman.pval))
  usual.shuffle.metrics.df <- rbind(usual.shuffle.metrics.df,
                                    data.frame(PHENO = pheno,
                                               PEARSON_CORR = usual.shuffled.pearson.corr.pval,
                                               SPEARMAN_RHO = usual.shuffled.spearman.corr.pval,
                                               COSINE_SIM = usual.shuffled.cosine.sim.pval,
                                               TOP10PCT_PREV = usual.shuffled.prev.pval,
                                               TOP1PCT_PREV = usual.shuffled.prev2.pval,
                                               TOP10PCT_OR = usual.shuffled.top10pct.or.pval,
                                               TOP1PCT_OR = usual.shuffled.top1pct.or.pval,
                                               PERCENTILE_AVE_PREV_SLOPE = usual.shuffled.aveprev.pearson.pval,
                                               PERCENTILE_AVE_PREV_SPEARMAN = usual.shuffled.aveprev.spearman.pval,
                                               PERCENTILE_AVE_PHENO_SLOPE = usual.shuffled.avepheno.pearson.pval,
                                               PERCENTILE_AVE_PHENO_SPEARMAN = usual.shuffled.avepheno.spearman.pval))
  usual.signflip.metrics.df <- rbind(usual.signflip.metrics.df,
                                     data.frame(PHENO = pheno,
                                                PEARSON_CORR = usual.signflip.pearson.corr.pval,
                                                SPEARMAN_RHO = usual.signflip.spearman.corr.pval,
                                                COSINE_SIM = usual.signflip.cosine.sim.pval,
                                                TOP10PCT_PREV = usual.signflip.prev.pval,
                                                TOP1PCT_PREV = usual.signflip.prev2.pval,
                                                TOP10PCT_OR = usual.signflip.top10pct.or.pval,
                                                TOP1PCT_OR = usual.signflip.top1pct.or.pval,
                                                PERCENTILE_AVE_PREV_SLOPE = usual.signflip.aveprev.pearson.pval,
                                                PERCENTILE_AVE_PREV_SPEARMAN = usual.signflip.aveprev.spearman.pval,
                                                PERCENTILE_AVE_PHENO_SLOPE = usual.signflip.avepheno.pearson.pval,
                                                PERCENTILE_AVE_PHENO_SPEARMAN = usual.signflip.avepheno.spearman.pval))
}

# Save files
readr::write_csv(shuffle.metrics.df,
                 file = paste0(out.dir,
                               'gwas',
                               gwas.thres,'_cutoff',
                               sig.cutoff,'_shuffle_metrics_df.csv'))
readr::write_csv(signflip.metrics.df,
                 file = paste0(out.dir,
                               'gwas',
                               gwas.thres,'_cutoff',
                               sig.cutoff,'_signflip_metrics_df.csv'))
readr::write_csv(usual.shuffle.metrics.df,
                 file = paste0(out.dir,
                               'gwas',
                               gwas.thres,'_cutoff',
                               sig.cutoff,'_shuffle_usual_metrics_df.csv'))
readr::write_csv(usual.signflip.metrics.df,
                 file = paste0(out.dir,
                               'gwas',
                               gwas.thres,'_cutoff',
                               sig.cutoff,'_signflip_usual_metrics_df.csv'))
sink()

################################################################################
# esize.dist.metrics.df <- readr::read_csv('/deep_learning/aaw/022823/results/vars_and_mags_df.csv')
# colnames(esize.dist.metrics.df)[1]<-'PHENO'
# 
# old.shuffle.metrics.df <- readr::read_csv('/deep_learning/aaw/032323/tables/shuffle_metrics_df.csv')
# old.signflip.metrics.df <- readr::read_csv('/deep_learning/aaw/032323/tables/signflip_metrics_df.csv')
# old.shuffle.signflip.metrics.df <- readr::read_csv('/deep_learning/aaw/032323/tables/shuffle_signflip_metrics_df.csv')
# merged.old.shuffle.df <- merge(esize.dist.metrics.df %>% select(c('PHENO','SUM_BKGRD','SUM_TARGET','SUM_RATIO_TARGET_BKGRD')),
#                                old.shuffle.metrics.df %>% select(c('PHENO','PEARSON_CORR_PVAL','SPEARMAN_RHO_PVAL')),
#                                by='PHENO')
# merged.old.shuffle.df$PEARSON_CORR_PVAL<- 1-merged.old.shuffle.df$PEARSON_CORR_PVAL
# merged.old.shuffle.df$SPEARMAN_RHO_PVAL<- 1-merged.old.shuffle.df$SPEARMAN_RHO_PVAL
# #merged.old.shuffle.df$TOP10PCT_PREV_PVAL<- 1-merged.old.shuffle.df$TOP10PCT_PREV_PVAL
# #merged.old.shuffle.df$TOP1PCT_PREV_PVAL<- 1-merged.old.shuffle.df$TOP1PCT_PREV_PVAL
# 
# 
# combined.shuffle.old <- GGally::ggpairs(merged.old.shuffle.df[,-1],
#                                         upper = list(continuous = wrap("cor", method = "spearman")))
# ggsave(combined.shuffle.old,
#        filename = "/deep_learning/aaw/032923/plots/signflip_old_sum_ratio_target_bkgrd_spearman.jpg",
#        width = 12, height = 12, 
#        dpi = 300)
# 
# spearman.rho.old <- ggplot(merged.old.shuffle.df,aes(x=SUM_RATIO_TARGET_BKGRD,y=SPEARMAN_RHO_PVAL)) +
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method='loess',se=FALSE) +
#   stat_cor(label.x=6, label.y=0.45) +
#   stat_cor(label.x=6, label.y=0.5, method='spearman',cor.coef.name='rho') +
#   #xlab(expression(paste('Median |',beta['target'],'|'))) +
#   xlab(expression(paste(Sigma,'|',beta['target'],'|', '/', Sigma,'|',beta['bkgrd'],'|'))) +
#   ylab(expression(paste('Fraction of Perturbed PRSes with ', rho, ' as Large as Original'))) +
#   #ggtitle(expression(paste('Sensitivity of PRS Spearman ',rho,' to Sign-flipping vs Median Target Variant Effect Sizes |',beta,'|')))
#   ggtitle(expression(paste('Sensitivity of PRS Spearman ',rho,' to Shuffling vs Ratio of ', Sigma,'|',beta,'|')))
# 
# pearson.r.old <- ggplot(merged.old.shuffle.df,aes(x=SUM_RATIO_TARGET_BKGRD,y=PEARSON_CORR_PVAL)) +
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method='loess',se=FALSE) +
#   stat_cor(label.x=6, label.y=0.65) +
#   stat_cor(label.x=6, label.y=0.7, method='spearman',cor.coef.name='rho') +
#   #xlab(expression(paste('Median |',beta['target'],'|'))) +
#   xlab(expression(paste(Sigma,'|',beta['target'],'|', '/', Sigma,'|',beta['bkgrd'],'|'))) +
#   ylab(expression(paste('Fraction of Perturbed PRSes with ', r, ' as Large as Original'))) +
#   #ggtitle(expression(paste('Sensitivity of PRS Pearson ',r,' to Sign-flipping vs Median Target Variant Effect Sizes |',beta,'|')))
#   ggtitle(expression(paste('Sensitivity of PRS Pearson ',r,' to Shuffling vs Ratio of ', Sigma,'|',beta,'|')))
# 
# combined.shuffle.old <- gridExtra::grid.arrange(spearman.rho.old,pearson.r.old,ncol=2)
# ggsave(combined.shuffle.old,
#        filename = "/deep_learning/aaw/032923/plots/shuffle_old_combined_sum_ratio_target_bkgrd.jpg",
#        width = 13, height = 5, 
#        dpi = 300)
# 
# merged.new.shuffle.df <- merge(esize.dist.metrics.df %>% select(c('PHENO','SUM_BKGRD','SUM_TARGET','SUM_RATIO_TARGET_BKGRD')),
#                                shuffle.metrics.df %>% select(c('PHENO','PEARSON_CORR_PVAL','SPEARMAN_RHO_PVAL')),
#                                by='PHENO')
# merged.new.shuffle.df$PEARSON_CORR_PVAL<- 1-merged.new.shuffle.df$PEARSON_CORR_PVAL
# merged.new.shuffle.df$SPEARMAN_RHO_PVAL<- 1-merged.new.shuffle.df$SPEARMAN_RHO_PVAL
# 
# ## Plot boxplots for slide [4/10/23] -------------------------------------------
# merged.new.shuffle.df$SIGNIFICANT <- sapply(merged.new.shuffle.df$SPEARMAN_RHO_PVAL, function(x) {ifelse(x<= 0.05,
#                                                                                                          'Yes','No')})
# 
# boxplot.fig <- ggplot(merged.new.shuffle.df, aes(x = SIGNIFICANT, y = SUM_RATIO_TARGET_BKGRD)) +
#   geom_violin(aes(fill = SIGNIFICANT), width=1.4) +
#   geom_boxplot(width=0.05, color="black", alpha=0.2) +
#   theme_bw() +
#   ggtitle('Does original PRS have larger\nSpearman correlation than 95% of\nshuffled PRS?') +
#   xlab("") +
#   ylab(expression(paste('Ratio of Sum of Effect Sizes (',Sigma,'|',beta['perturb'],'|', '/', Sigma,'|',beta['fix'],'|)'))) +
#   theme(legend.position="none",
#         title = element_text(size=16,hjust=0.5),
#         axis.title = element_text(size = 17),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 17)) 
#   
# ggsave(boxplot.fig,
#        filename = "/deep_learning/aaw/032923/plots/boxplot_shuffle_new_combined_sum_ratio_target_bkgrd.jpg",
#        width = 6, height = 6, 
#        dpi = 300)
# ## End plot boxplots for slide [4/10/23] ---------------------------------------
# 
# GGally::ggpairs(merged.new.shuffle.df[,-1],
#                 upper = list(continuous = wrap("cor", method = "spearman")))
# 
# spearman.rho.new <- ggplot(merged.new.shuffle.df,aes(x=SUM_RATIO_TARGET_BKGRD,y=SPEARMAN_RHO_PVAL)) +
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method='loess',se=FALSE) +
#   stat_cor(label.x=6, label.y=0.7) +
#   stat_cor(label.x=6, label.y=0.75, method='spearman',cor.coef.name='rho') +
#   #xlab(expression(paste('Median |',beta['target'],'|'))) +
#   xlab(expression(paste(Sigma,'|',beta['target'],'|', '/', Sigma,'|',beta['bkgrd'],'|'))) +
#   ylab(expression(paste('Fraction of Perturbed PRSes with ', rho, ' as Large as Original'))) +
#   #ggtitle(expression(paste('Sensitivity of PRS Spearman ',rho,' to Sign-flipping vs Median Target Variant Effect Sizes |',beta,'|')))
#   ggtitle(expression(paste('Sensitivity of PRS Spearman ',rho,' to Shuffling vs Ratio of ', Sigma,'|',beta,'|')))
# 
# pearson.r.new <- ggplot(merged.new.shuffle.df,aes(x=SUM_RATIO_TARGET_BKGRD,y=PEARSON_CORR_PVAL)) +
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method='loess',se=FALSE) +
#   stat_cor(label.x=6, label.y=0.75) +
#   stat_cor(label.x=6, label.y=0.8, method='spearman',cor.coef.name='rho') +
#   #xlab(expression(paste('Median |',beta['target'],'|'))) +
#   xlab(expression(paste(Sigma,'|',beta['target'],'|', '/', Sigma,'|',beta['bkgrd'],'|'))) +
#   ylab(expression(paste('Fraction of Perturbed PRSes with ', r, ' as Large as Original'))) +
#   #ggtitle(expression(paste('Sensitivity of PRS Pearson ',r,' to Sign-flipping vs Median Target Variant Effect Sizes |',beta,'|')))
#   ggtitle(expression(paste('Sensitivity of PRS Pearson ',r,' to Shuffling vs Ratio of ', Sigma,'|',beta,'|')))
# 
# combined.shuffle.new <- gridExtra::grid.arrange(spearman.rho.new,pearson.r.new,ncol=2)
# ggsave(combined.shuffle.new,
#        filename = "/deep_learning/aaw/032923/plots/shuffle_new_combined_sum_ratio_target_bkgrd.jpg",
#        width = 13, height = 5, 
#        dpi = 300)
# 
# metric.names <- colnames(old.shuffle.metrics.df)[-1]
# 
# for (metric in metric.names) {
#   message(date(), ': Working on ', metric)
#   message('No. insignificant (old) = ', sum(old.signflip.metrics.df[[metric]] < 0.95))
#   message('No. insignificant (new) = ', sum(signflip.metrics.df[[metric]] < 0.95))
# }
# 
