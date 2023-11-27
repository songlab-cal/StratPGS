############################################
########### Auxiliary functions ############
############################################
covar_corrs <- readr::read_csv('covar_corrs.csv')
covar_cossim <- readr::read_csv('covar_cossim.csv')
phenos_corrs <- readr::read_csv('phenos_corrs.csv')
phenos_cossim <- readr::read_csv('phenos_cossim.csv')
phenos_random_stats <- readr::read_csv('phenos_randomization_stats.csv')
ukbb_table <- readr::read_csv('ukbb_table.csv')
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
