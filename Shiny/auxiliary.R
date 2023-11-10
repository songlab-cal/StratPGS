############################################
########### Auxiliary functions ############
############################################
covar_corrs <- readr::read_csv('covar_corrs.csv')
covar_cossim <- readr::read_csv('covar_cossim.csv')
phenos_corrs <- readr::read_csv('phenos_corrs.csv')
phenos_cossim <- readr::read_csv('phenos_cossim.csv')
phenos_random_stats <- readr::read_csv('phenos_randomization_stats.csv')

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

getPhenoCorrPlot <- function(x) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating corrplot for ', x, '...')
  corr.with.gpcs <- phenos_corrs %>% subset(PHENOTYPE == x)
  melted.corr.with.gpcs <- reshape2::melt(corr.with.gpcs,id.vars = 'PHENOTYPE')
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

getPhenoCossimPlot <- function(x) {
  #covar.plot.list <- readRDS(file=paste0('plots/',x,'_plots.rds'))
  message(date(), ': Generating cosine similarity plot for ', x, '...')
  cossim.with.gpcs <- phenos_cossim %>% subset(PHENOTYPE == x)
  melted.cossim.with.gpcs <- reshape2::melt(cossim.with.gpcs,id.vars = 'PHENOTYPE')
  melted.cossim.with.gpcs$value <- as.numeric(melted.cossim.with.gpcs$value)
  melted.cossim.with.gpcs$max_mag <- abs(melted.cossim.with.gpcs$value) == max(abs(melted.cossim.with.gpcs$value))
  
  cossim.plot<- ggplot(melted.cossim.with.gpcs, aes(x = variable, y = value)) +
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

getPhenoStats <- function(x) {
  message(date(), ': Generating statistics for ', x, '...')
  phenos.with.stats <- phenos_random_stats %>% 
    subset(PHENOTYPE == x) %>% 
    as.data.frame()
  colnames(phenos.with.stats) <- c('Phenotype',
                                   'Cosine Similarity Evenness',
                                   'Pearson r Evenness',
                                   'Max Prevalence-Percentile Slope',
                                   'Slope p-value',
                                   'Mean Incremental R2',
                                   'Incremental R2 p-value')
  return(list(TABLE = phenos.with.stats[,-1],
              SENTENCE = paste0('Below are some statistics for ', x,'.')))
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
