############################################
########### Auxiliary functions ############
############################################
covar_corrs <- readr::read_csv('covar_corrs.csv')
covar_cossim <- readr::read_csv('covar_cossim.csv')
phenos_corrs <- readr::read_csv('phenos_corrs.csv')
phenos_cossim <- readr::read_csv('phenos_cossim.csv')
phenos_random_stats <- readr::read_csv('phenos_randomization_stats.csv')
ukbb_table <- readr::read_csv('ukbb_table.csv')

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
          legend.position = 'none') +
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
          legend.position = 'none') +
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
          legend.position = 'none') +
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
          legend.position = 'none') +
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
  if (y == "Clumping and thresholding (lenient)") {
    pgs_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-5_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-5_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else if (y == "Clumping and thresholding (stringent)") {
    pgs_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-8_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-8_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else {
    pgs_df <- NULL
    pgs_strat_df <- NULL
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
      theme(legend.position="none")
  } else {
    df_to_return <-NULL
    plot_to_return<-NULL
  }
  
  # return
  return(list(SENTENCE=ifelse(is.null(pgs_df),"",paste0('There are ', nrow(pgs_df), ' variants in ',y,'.')),
              TABLE=df_to_return,
              PLOT=plot_to_return,
              CHROM_DIST_SENTENCE=ifelse(is.null(pgs_df),"","1. Variant Distribution"),
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
    perturb_fix_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-5_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-5_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else if (y == "Clumping and thresholding (stringent)") {
    pgs_df <- readr::read_delim(paste0("tables/CnT_PGS_files/",pheno_name,".loci.common.1e-8_no_dups.txt"))
    pgs_strat_df <- readr::read_csv(paste0("tables/train_n288728_prs_gwas_1e-8_gPC_metrics.csv")) %>%
      subset(PHENO==pheno_name)
  } else {
    pgs_df <- NULL
    pgs_strat_df <- NULL
  }
  
  return(list(PERTURB_FIX_SENTENCE=ifelse(is.null(pgs_df),"","3. Perturbed-Fixed Architecture"),
              SENSITIVITY_SENTENCE=ifelse(is.null(pgs_df),"","4. Performance and Sensitivity")))
}

#' Get PGS summary statistics 2 -- specific performance metric and sensitivity
#' x = PGS method
#' y = cutoff
#' z = metric
getPGSStats2 <- function(x,y,z) {
  
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
