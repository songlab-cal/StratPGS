############################################
########### Auxiliary functions ############
############################################
covar_corrs <- readr::read_csv('covar_corrs.csv')
covar_cossim <- readr::read_csv('covar_cossim.csv')
phenos_corrs <- readr::read_csv('phenos_corrs.csv')
phenos_cossim <- readr::read_csv('phenos_cossim.csv')
phenos_random_stats <- readr::read_csv('phenos_randomization_stats.csv')
ukbb_table <- readr::read_csv('ukbb_table.csv')

metric_names_df <- data.frame(english_name=c("Pearson r",
                                             "Spearman \u03c1",
                                             "Cosine Similarity",
                                             "Prevalence at Top 10%ile",
                                             "Prevalence at Top 1%ile",
                                             "Odds Ratio at Top 10%ile",
                                             "Odds Ratio at Top 1%ile",
                                             "Percentile-Prevalence Pearson r",
                                             "Percentile-Prevalence Spearman \u03c1",
                                             "Percentile-Average Phenotype Pearson r",
                                             "Percentile-Average Phenotype Spearman \u03c1"),
                              col_name=c("PEARSON_CORR","SPEARMAN_RHO","COSINE_SIM",
                                         "TOP10PCT_PREV","TOP1PCT_PREV","TOP10PCT_OR",
                                         "TOP1PCT_OR","PERCENTILE_AVE_PREV_SLOPE",
                                         "PERCENTILE_AVE_PREV_SPEARMAN",
                                         "PERCENTILE_AVE_PHENO_SLOPE",
                                         "PERCENTILE_AVE_PHENO_SPEARMAN"))

removeLeadingZero <- function(inputString) {
  if (substr(inputString, 1, 1) == "0") {
    return(substr(inputString, 2, nchar(inputString)))
  } else {
    return(inputString)
  }
}

getCorrPlot <- function(x) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating corrplot for ', x, '...')
  corr.with.gpcs <- covar_corrs %>% subset(covar == x)
  melted.corr.with.gpcs <- reshape2::melt(corr.with.gpcs,id.vars = 'covar')
  melted.corr.with.gpcs$value <- as.numeric(melted.corr.with.gpcs$value)
  melted.corr.with.gpcs$max_mag <- abs(melted.corr.with.gpcs$value) == max(abs(melted.corr.with.gpcs$value))
  
  corr.plot <- ggplot(melted.corr.with.gpcs, aes(x = variable, y = value)) +
    geom_bar(stat = 'identity',
             aes(fill = max_mag)) + 
    theme_bw() +
    ggtitle('Correlation with Genetic Principal Components') +
    ylab('Pearson Correlation (r)') +
    xlab('') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          text=element_text(family='DejaVu Sans')) +
    scale_fill_manual(values=c('#636363', '#e31a1c'))
  return(corr.plot)
}

getCossimPlot <- function(x) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating corrplot for ', x, '...')
  cossim.with.gpcs <- covar_cossim %>% subset(covar == x)
  melted.cossim.with.gpcs <- reshape2::melt(cossim.with.gpcs,id.vars = 'covar')
  melted.cossim.with.gpcs$value <- as.numeric(melted.cossim.with.gpcs$value)
  melted.cossim.with.gpcs$max_mag <- abs(melted.cossim.with.gpcs$value) == max(abs(melted.cossim.with.gpcs$value))
  
  cossim.plot <- ggplot(melted.cossim.with.gpcs, aes(x = variable, y = value)) +
    geom_bar(stat = 'identity',
             aes(fill = max_mag)) + 
    theme_bw() +
    ggtitle('Cosine Similarity with Genetic Principal Components') +
    ylab('Cosine Similarity') +
    xlab('') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          text=element_text(family='DejaVu Sans')) +
    scale_fill_manual(values=c('#636363', '#e31a1c'))
  
  return(cossim.plot)
}

#' Get Phenotype Pearson Correlation plot
#' 
getPhenoCorrPlot <- function(x,y) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating correlation plot for ', x, ' (',y,')...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  relevant_name <- ifelse(y=='Original', 
                          relevant_row$orig_name,
                          relevant_row$irnt_name)
  corr_with_gpcs <- phenos_corrs %>% subset(PHENOTYPE==relevant_name)
  melted_corr_with_gpcs <- reshape2::melt(corr_with_gpcs,id.vars='PHENOTYPE')
  melted_corr_with_gpcs$value <- as.numeric(melted_corr_with_gpcs$value)
  melted_corr_with_gpcs$max_mag <- abs(melted_corr_with_gpcs$value) == 
    max(abs(melted_corr_with_gpcs$value))
  
  corr_plot <- ggplot(melted_corr_with_gpcs, aes(x=variable,
                                                 y=value)) +
    geom_bar(stat = 'identity',
             aes(fill = max_mag)) + 
    theme_bw() +
    ggtitle('Correlation with Genetic Principal Components') +
    ylab('Pearson Correlation (r)') +
    xlab('') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          text=element_text(family='DejaVu Sans')) +
    scale_fill_manual(values=c('#636363', '#e31a1c'))
  return(corr_plot)
}

#' Get Phenotype Cosine Similarity plot
#' 
getPhenoCossimPlot <- function(x,y) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating cosine similarity plot for ', x, ' (',y,')...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  relevant_name <- ifelse(y=='Original', 
                          relevant_row$orig_name,
                          relevant_row$irnt_name)
  cossim_with_gpcs <- phenos_cossim %>% subset(PHENOTYPE==relevant_name)
  melted_cossim_with_gpcs <- reshape2::melt(cossim_with_gpcs,id.vars='PHENOTYPE')
  melted_cossim_with_gpcs$value <- as.numeric(melted_cossim_with_gpcs$value)
  melted_cossim_with_gpcs$max_mag <- abs(melted_cossim_with_gpcs$value)==
    max(abs(melted_cossim_with_gpcs$value))
  
  cossim_plot<- ggplot(melted_cossim_with_gpcs, aes(x=variable, 
                                                    y=value)) +
    geom_bar(stat = 'identity',
             aes(fill = max_mag)) + 
    theme_bw() +
    ggtitle('Cosine Similarity with Genetic Principal Components') +
    ylab('Cosine Similarity') +
    xlab('') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          text=element_text(family='DejaVu Sans')) +
    scale_fill_manual(values=c('#636363', '#e31a1c'))
  return(cossim_plot)
}

#' Get Phenotype stratification summary statistics
#' 
getPhenoStats <- function(x,y) {
  message(date(), ': Generating phenotype stratification statistics for ', x, ' (',y,')...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  relevant_name <- ifelse(y=='Original', 
                          relevant_row$orig_name,
                          relevant_row$irnt_name)
  phenos_with_stats <- phenos_random_stats %>% 
    subset(PHENOTYPE==relevant_name) %>% 
    as.data.frame()
  colnames(phenos_with_stats) <- c('Phenotype',
                                   'Cosine Similarity Evenness',
                                   'Pearson r Evenness',
                                   'Max Prevalence-Percentile Slope',
                                   'Slope p-value',
                                   'Mean Incremental R\u00B2',
                                   'Incremental R\u00B2 p-value')
  phenos_with_stats[['Mean Incremental R\u00B2']] <- format(
    phenos_with_stats[['Mean Incremental R\u00B2']],
    scientific=TRUE)
  phenos_with_stats[['Incremental R\u00B2 p-value']] <- format(
    phenos_with_stats[['Incremental R\u00B2 p-value']],
    scientific=TRUE)
  return(list(TABLE = phenos_with_stats[,c(2:3,6:7)],
              SENTENCE = paste0('3. Other Statistics for ', x, ' (',y,').')))
}

summarizePheno <- function(name) {
  message(date(), ': Generating summary list for ', name, '...')
  target_row <- ukbb_table %>% subset(english_name==name)
  return(list(No_Avail_PGS=target_row$no_avail_pgs,
              Zero_PGS_Sentence=paste0("No PGS data for ", name,"."),
              Category=target_row$category,
              Field=target_row$ukbb_field,
              URL=paste0("https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=",target_row$ukbb_field)))
}

#' Get PGS stratification summary statistics
#' x = phenotype
#' y = PGS type
getPGSStratStats <- function(x, y) {
  # no. variants in PGS (summarize from directory of PGS files)
  # distribution across autosomes
  # PC-related stratification metrics x 11
  # performance metrics 
  message(date(), ': Generating ', y,' stratification statistics for ', x,'...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  pheno_name <- relevant_row[['pheno_name']]
  if (is.null(y)) {
    pgs_df <- NULL
    pgs_strat_df <- NULL
  } else if (y == "Clumping and thresholding (lenient)") {
    pgs_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-5_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-5_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else if (y == "Clumping and thresholding (stringent)") {
    pgs_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-8_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-8_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else {
    stop("y must be NULL or Clumping and thresholding.")
  }
  
  if (!is.null(pgs_df)) {
    chrom_tab <- table(pgs_df$chrom)
    pie_chart_df <- data.frame(Chromosome=sapply(names(chrom_tab),removeLeadingZero),
                               Count=as.numeric(chrom_tab),
                               row.names=NULL)  
    
    # collate objects to return
    df_to_return <- data.frame(`Unsigned Metric` = c("PC1 Cosine Similarity",
                                                     "PC1 Pearson r",
                                                     "PC1 Spearman \u03C1"),
                               `Value` = c(pgs_strat_df$PC1_ABS_COSSIM,
                                           pgs_strat_df$PC1_ABS_PEARSON,
                                           pgs_strat_df$PC1_ABS_SPEARMAN),
                               `Evenness Metric` = c("Cosine Similarity",
                                                     "Pearson r",
                                                     "Spearman \u03C1"),
                               `Value` = c(pgs_strat_df$EVENNESS_COSSIM,
                                           pgs_strat_df$EVENNESS_PEARSON,
                                           pgs_strat_df$EVENNESS_SPEARMAN),
                               `Linear Model Metric` = c("R\u00B2",
                                                         "Adjusted R\u00B2",
                                                         "No. Significant Features"),
                               `Value` = c(format(pgs_strat_df$R2_PHENO_PCS,scientific=TRUE),
                                           format(pgs_strat_df$ADJ_R2_PHENO_PCS,scientific=TRUE),
                                           pgs_strat_df$LinearModel_NUM_SIG_VARS),
                               check.names=FALSE)
    # plot_to_return <- ggplot(pie_chart_df, aes(x="", y=Count, fill=Chromosome)) +
    #   geom_col(width = 1, color = 1) +
    #   geom_bar(stat="identity", width=1) +
    #   coord_polar("y", start=0) +
    #   theme_void() 
    plot_to_return <- ggplot(pie_chart_df, aes(x=Chromosome, y=Count, fill=Chromosome)) +
      geom_bar(stat="identity", width=1,colour="#636363") +
      theme_bw() +
      theme(legend.position="none",
            plot.title = element_text(face='bold',size=14.5,hjust=0.5),
            strip.text.x = element_text(size = 15),
            axis.title = element_text(size=14),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            text=element_text(family='DejaVu Sans'))
  } else {
    df_to_return <-NULL
    plot_to_return<-NULL
  }
  
  # return
  return(list(SENTENCE=ifelse(is.null(pgs_df),
                              "",
                              paste0('There are ', nrow(pgs_df), 
                                     ' variants in ',y,' PGS for ', x,'.')),
              TABLE=df_to_return,
              PLOT=plot_to_return,
              CHROM_DIST_SENTENCE=ifelse(is.null(pgs_df),"","1. Variant Distribution by Chromosome"),
              PGS_STRAT_SENTENCE=ifelse(is.null(pgs_df),"","2. PC Stratification")))
}

#' Get PGS summary statistics 1 -- perturbed-fixed architecture
#' x = phenotype
#' y = PGS type
#' z = cutoff
getPGSStats1 <- function(x, y, z) {
  # sensitivity metrics (matched to performance metrics)
  # fixed-perturbed architecture
  message(date(), ': Generating ', y,' perturbed-fixed architecture using cutoff = ', z,'...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  pheno_name <- relevant_row[['pheno_name']]

  if (y == "Clumping and thresholding (lenient)") {
    perturb_fix_df <- readr::read_csv(paste0(
      "tables/perturb-fixed-architecture/gwas1e-5_cutoff",z,".csv")) %>%
      subset(PHENOTYPE==pheno_name)
  } else if (y == "Clumping and thresholding (stringent)") {
    perturb_fix_df <- readr::read_csv(paste0(
      "tables/perturb-fixed-architecture/gwas1e-8_cutoff",z,".csv")) %>%
      subset(PHENOTYPE==pheno_name)
  } else {
    perturb_fix_df <- NULL
  }
  
  df_to_return <- data.frame(`Fixed Variant Quantity` = c("Maximum |\u03B2\u2C7C|",
                                                          "Median |\u03B2\u2C7C|",
                                                          "Mean |\u03B2\u2C7C|",
                                                          "Sum of |\u03B2\u2C7C|",
                                                          "Variance of \u03B2\u2C7C"),
                             `Value` = c(perturb_fix_df$MAX_MAG_BKGRD,
                                         perturb_fix_df$MEDIAN_MAG_BKGRD,
                                         perturb_fix_df$AVE_MAG_BKGRD,
                                         perturb_fix_df$SUM_BKGRD,
                                         perturb_fix_df$VAR_BKGRD),
                             `Perturbed Variant Quantity` = c("Maximum |\u03B2\u2C7C|",
                                                              "Median |\u03B2\u2C7C|",
                                                              "Mean |\u03B2\u2C7C|",
                                                              "Sum of |\u03B2\u2C7C|",
                                                              "Variance of \u03B2\u2C7C"),
                             `Value` = c(perturb_fix_df$MAX_MAG_TARGET,
                                         perturb_fix_df$MEDIAN_MAG_TARGET,
                                         perturb_fix_df$AVE_MAG_TARGET,
                                         perturb_fix_df$SUM_TARGET,
                                         perturb_fix_df$VAR_TARGET),
                             `Perturbed vs Fixed` = c("Wilcoxon(Perturbed |\u03B2\u2C7C|,Fixed |\u03B2\u2C7C|) p-value",
                                                      "Sum(Perturbed |\u03B2\u2C7C|)/Sum(Fixed |\u03B2\u2C7C|)",
                                                      "","",""),
                             `Value` = c(format(perturb_fix_df$BETA_WILCOX_PVAL,scientific=TRUE),
                                         round(perturb_fix_df$SUM_RATIO_TARGET_BKGRD,digits=4),
                                         NA,NA,NA),
                             check.names=FALSE) 
  
  return(list(SENTENCE=ifelse(is.null(perturb_fix_df),"",paste0('There are ', 
                                                                perturb_fix_df$N_TARGET_SNPS, 
                                                                ' variants with GWAS p-value at least ',
                                                                z, '. These variants are Perturbed, while the remaining ',
                                                                perturb_fix_df$N_SNPS-perturb_fix_df$N_TARGET_SNPS,' are Fixed.')),
              TABLE=df_to_return,
              PERTURB_FIX_HEADER=ifelse(is.null(perturb_fix_df),"","Perturbed-Fixed Architecture")))
}

#' Get PGS summary statistics 2 -- specific performance metric and sensitivity
#' Show both sPGS and permutation pPGS
#' x = phenotype
#' y = PGS type
#' z = cutoff
#' w = metric
getPGSStats2 <- function(x,y,z,w) {
  # sensitivity metrics (matched to performance metrics)
  message(date(), ': Generating ', 
          y,' sensitivity plots and data using cutoff = ', 
          z,' and metric = ',w,'...')
  relevant_row <- ukbb_table %>% subset(english_name==x)
  pheno_name <- relevant_row[['pheno_name']]
  metric_name <- (metric_names_df %>% subset(english_name==w))[["col_name"]]
  if (y == "Clumping and thresholding (lenient)") {
    sensitivity_df <- readr::read_csv(paste0(
      "tables/raw_perf_vs_perturb/",pheno_name,
      "_gwas_1e-5_cutoff",gsub("e-0", "e-", z),"_metrics_df.csv")) %>% 
      select(c("TYPE",metric_name))
  } else if (y == "Clumping and thresholding (stringent)"){
    sensitivity_df <- readr::read_csv(paste0(
      "tables/raw_perf_vs_perturb/",pheno_name,
      "_gwas_1e-8_cutoff",gsub("e-0", "e-", z),"_metrics_df.csv")) %>% 
      select(c("TYPE",metric_name))
  } else {
    sensitivity_df <- NULL
  }
  colnames(sensitivity_df) <- c("Type","value")
  
  shuffle.and.orig.melted.df <- reshape2::melt(sensitivity_df, id.vars = 'Type')
  
  # Generate plot 
  plot_to_return <- ggplot(sensitivity_df, aes(x=value,fill=factor(Type))) +
    geom_histogram(position = "dodge") +
    theme_bw() +
    ylab('Count') +
    xlab('Value') +
    guides(fill=guide_legend("PGS Type")) +
    ggtitle('Relative Performance of PGS') +
    theme(legend.position = 'right',
          plot.title = element_text(size=16,hjust=0.5),
          strip.text.x = element_text(size = 14.5),
          axis.title = element_text(size=15),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          text=element_text(family='DejaVu Sans')) +
    scale_fill_manual(values = c("#b30000","#6baed6","#74c476"),
                      breaks = c("original", "shuffle", "signflip"),
                      labels=c("original","pPGS","sPGS")) 
  
  # Generate sentence
  orig_perf <- (sensitivity_df %>% subset(Type=="original"))[['value']]
  sPGS_p_val <- format(mean(orig_perf <= 
                              c((sensitivity_df %>% subset(Type=="signflip"))[['value']],orig_perf)),
                       scientific=TRUE)
  pPGS_p_val <- format(mean(orig_perf <= 
                              c((sensitivity_df %>% subset(Type=="shuffle"))[['value']],orig_perf)),
                       scientific=TRUE)
  
  df_to_return <- data.frame(Quantity = c("Original PGS Performance Metric",
                                          "Permuted PGS Relative Performance p-val",
                                          "Sign-flipped PGS Relative Performance p-val"),
                             Value = c(round(orig_perf,digits=4),
                                       pPGS_p_val,
                                       sPGS_p_val))
  # Return
  return(list(PLOT=plot_to_return,
              TABLE=df_to_return,
              HEADER=ifelse(is.null(sensitivity_df),"","PGS Performance and Sensitivity"),
              SENTENCE=ifelse(is.null(sensitivity_df),
                              "",
                              "Performance of the PGS based on the selected metric is reported below, along with its relative performance against perturbed PGSs. There are two types of perturbed PGSs: permuted PGSs (pPGSs) and sign-flipped PGSs (sPGSs). See <b>About</b> tab for details.")))
}

# x <- covar.ids[1]
# plots <- getPCPlots(x)
# getwd()
# readRDS('/deep_learning/aaw/021423/plots/Alcohol_intake_frequency_-3.0_plots.rds')
# List the library paths
# The issue is likely to be in the first directory
# paths = .libPaths()
# 
# ## Try and detect bad files
# list.files(paths, 
#            pattern = "^00LOCK*|*\\.rds$|*\\.RDS$",
#            full.names = TRUE)
# 
# ## List files of size 0
# l = list.files(paths, full.names = TRUE)
# l[sapply(l, file.size) == 0]
