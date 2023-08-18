## Mining PGS Catalogue metadata -----------------------------------------------
# First is PGS000001
# Last one is PGS003759

library(jsonlite, 
        lib.loc = "/global/home/groups/consultsw/sl-7.x86_64/modules/r-packages/2022-10-14-r4.2")

# Test 
webpage <- read_json('https://www.pgscatalog.org/rest/score/all')
res_nodes <- html_nodes(webpage, "results")
ids <- unlist(lapply(webpage$results, function(x){return(x[['id']])}))
prefix_url <- 'https://www.pgscatalog.org/rest/score/all'
pgs_string <- paste0(paste0('PGS0000',51:99,collapse=','),',PGS000100')
scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
webpage_2 <- read_json(scrap_input)
ids_2 <- unlist2(lapply(webpage_2$results, function(x){return(x[['id']])}))
sampsize_2 <- unlist2(lapply(webpage_2$results, function(x) {return(getSampSize(x))}))
method_name_2 <- unlist2(lapply(webpage_2$results, function(x){return(x[['method_name']])}))
no_variants_2 <- unlist2(lapply(webpage_2$results, function(x){return(x[['variants_number']])}))
params_2 <- unlist2(lapply(webpage_2$results, function(x){return(x[['method_params']])}))
trait_reported_2 <- unlist2(lapply(webpage_2$results, function(x){return(x[['trait_reported']])}))

## Helper functions ------------------------------------------------------------
unlist2 <- function(my_list) {
  null_indexes <- sapply(my_list, is.null)
  my_list[null_indexes] <- NA
  return(unlist(my_list))
}

getScrapInput <- function(i) {
  # This prefixes all cases
  prefix_url <- 'https://www.pgscatalog.org/rest/score/all'
  if (i==1) {
    return(prefix_url)
  } else if (i==2) {
    pgs_string <- paste0(paste0('PGS0000',51:99,collapse=','),',PGS000100')
    scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
    return(scrap_input)
  } else if (i<20) {
    pgs_string <- paste0('PGS000',(50*(i-1)+1):(50*i),collapse=',')
    scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
    return(scrap_input)
  } else if (i==20) {
    pgs_string <- paste0(paste0('PGS000',951:999,collapse=','),',PGS001000')
    scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
    return(scrap_input)
  } else if (i<=75){
    pgs_string <- paste0('PGS00',(50*(i-1)+1):(50*i),collapse=',')
    scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
    return(scrap_input)
  } else if (i==76) {
    pgs_string <- paste0('PGS00375',1:9,collapse=',')
    scrap_input<-paste0(prefix_url,'?filter_ids=',pgs_string)
    return(scrap_input)
  } else {
    warning("Index i needs to be between 1 and 76, inclusive.")
    return(NA)
  }
}

getSampSize <- function(x) {
  if (length(x[['samples_variants']])<1) {
    #search samples_training first
    if (length(x[['samples_training']])>=1) {
      return(unlist(lapply(x[['samples_training']], 
                           function(y) {return(y[['sample_number']])})) %>% sum())
    } else {
      # if not, return NA
      return(NA)
    }
  } else {
    #lapply to extract all sample_number subnode
    #sum
    return(unlist(lapply(x[['samples_variants']], 
                         function(y) {return(y[['sample_number']])})) %>% sum())
  }
}

## Main Body -------------------------------------------------------------------
library(logr)
out_dir <- "/global/scratch/projects/fc_songlab/alan/prs_and_structure/"
tmp <- file.path(paste0(out_dir,
                        "logs/mine_PGS_Catalogue_",
                        format(Sys.Date(), "%m%d%y.log")))
lf <- log_open(tmp)

PRS_metadata_df <- data.frame(ID=character(),
                              PHENOTYPE=character(),
                              SAMPLE_SIZE=numeric(),
                              METHODOLOGY=character(),
                              PARAMETERS=character(),
                              NO_VARIANTS=numeric())
for (i in 1:76) {
  log_print(paste0(date(), ": Scrapping next 50 PRSs starting from ", 50*(i-1)+1))
  scrap_input <- getScrapInput(i)
  webpage <- read_json(scrap_input)
  ids <- unlist2(lapply(webpage$results, function(x){return(x[['id']])}))
  sampsize <- unlist2(lapply(webpage$results, function(x) {return(getSampSize(x))}))
  method_name <- unlist2(lapply(webpage$results, function(x){return(x[['method_name']])}))
  no_variants <- unlist2(lapply(webpage$results, function(x){return(x[['variants_number']])}))
  params <- unlist2(lapply(webpage$results, function(x){return(x[['method_params']])}))
  trait_reported <- unlist2(lapply(webpage$results, function(x){return(x[['trait_reported']])}))
  PRS_metadata_df <- rbind(PRS_metadata_df,
                           data.frame(ID=ids,
                                      PHENOTYPE=trait_reported,
                                      SAMPLE_SIZE=sampsize,
                                      METHODOLOGY=method_name,
                                      PARAMETERS=params,
                                      NO_VARIANTS=no_variants))
}

# Save
readr::write_csv(PRS_metadata_df, 
                 file=paste0(out_dir,"PGS_catalogue_metadata.csv"))

# Close stream and write to log file
log_close()
writeLines(readLines(lf))

## Plotting
PRS_metadata_annotated_df <- readr::read_csv(paste0(out_dir,"PGS_catalogue_metadata_annotated.csv"))
View(PRS_metadata_df %>% 
       group_by(PHENOTYPE) %>% 
       summarise(MEDIAN=median(NO_VARIANTS),MAX=max(NO_VARIANTS)))
max_median_by_trait_df <- PRS_metadata_annotated_df %>% 
  group_by(PHENOTYPE) %>% 
  summarise(MEDIAN=median(NO_VARIANTS),MAX=max(NO_VARIANTS))
readr::write_csv(max_median_by_trait_df, 
                 file=paste0(out_dir,"PGS_catalogue_max_median_by_trait.csv"))
library(ggplot2)
plot1 <- ggplot(data=PRS_metadata_annotated_df) +
  geom_histogram(aes(x=log10(NO_NONZERO_VARIANTS)),
                 fill='grey',
                 color='black',
                 alpha=0.7) +
  theme_bw() +
  xlab(expression(paste(log[10],'(No. Variants Included in PRS)'))) +
  ylab('Count') +
  ylim(c(0,600)) +
  ggtitle('A. Dimensionality of 3,683 PRSs from the PGS Catalogue') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16))

ggsave(plot1,
       file = "/global/scratch/projects/fc_songlab/alan/prs_and_structure/PRS_variant_counts_logplot.jpg",
       width = 6.5, height = 4,
       dpi = 300)

## Creation of new Figure 1  ---------------------------------------------------
dir_prefix <- "/global/scratch/projects/fc_songlab/alan/prs_and_structure/"
prs_dist_summaries_cutoff1e6 <- readr::read_csv(paste0(dir_prefix,
                                                       "Github/MCH_PRS_dist_summaries/prs_dist_summaries_cutoff1e-6.csv"))
combined_cutoff1e6_shuffle_metrics <- readr::read_csv(paste0(dir_prefix,
                                                             "Github/MCH_raw_sensitivity_scores/combined_cutoff1e-6_shuffle_metrics_df.csv"))
prs_max_shuffle_sensitivity_min <- readr::read_csv(paste0(dir_prefix,
                                                             "Github/MCH_min_and_sum_sensitivity_scores/prs_max_shuffle_sensitivity_min.csv"))
orig_prs_ranked_df <- readr::read_csv(paste0(dir_prefix,
                                             "Github/MCH_original_PRS_performance/orig_prs_ranked_performance.csv"),
                                      show_col_types = FALSE)
orig_prs_performance <- readr::read_csv(paste0(dir_prefix,
                                             "Github/MCH_original_PRS_performance/orig_prs_performance.csv"),
                                      show_col_types = FALSE)
# Sensitivity plotting
library(ggrepel)
left_df <- prs_dist_summaries_cutoff1e6 %>% select(c("PHENOTYPE","N_SNPS")); colnames(left_df)[1] <- "PRS_ID"
right_df <- prs_max_shuffle_sensitivity_min %>% select(c("PRS_ID","CUTOFF_1e-6"));colnames(right_df)[2] <- "AGGREGATED"
combined_df_sensitivity <- left_join(left_df,right_df,by="PRS_ID")

new_plot_C <- ggplot(combined_df_sensitivity %>% 
                       subset(AGGREGATED >=0.75),
                     aes(y=AGGREGATED,
                         x=log10(N_SNPS),
                         label=PRS_ID)) +
  geom_point() +
  #geom_text(hjust=0,vjust=0) +
  geom_label_repel(aes(label = PRS_ID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   alpha=0.5,
                   segment.color = 'grey50',
                   max.overlaps=Inf) +
  theme_bw() +
  xlab(expression(paste(log[10],'(No. Variants)'))) +
  xlim(c(2,7)) +
  ylim(c(0.75,1)) +
  #ylab('Minimum Sensitivity Degree Across Metrics') +
  ylab('Sensitivity Score') +
  ggtitle('C. Sensitivity of MCH PRSs') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16))

# Original performance plotting
ranks_mat <- apply(orig_prs_performance[-c(16,17),-1],2,function(x){rank(x)})
agg_ranks <- rowSums(ranks_mat)
PRS_IDs <- orig_prs_performance$PRS_ID[-c(16,17)]
agg_rank_df <- data.frame(PRS_ID=PRS_IDs,AGGREGATED=agg_ranks)
combined_df_orig <- left_join(left_df,agg_rank_df,by="PRS_ID")

two_metrics_df <- orig_prs_performance[-c(16,17),c("PRS_ID","SPEARMAN_RHO","TOP10PCT_OR")]
two_metrics_df <- left_join(left_df,two_metrics_df,by="PRS_ID")

new_plot_B <- ggplot(two_metrics_df,aes(y=SPEARMAN_RHO,
                          x=TOP10PCT_OR,
                          label=PRS_ID)) +
  geom_point(aes(color=log10(N_SNPS))) +
  #geom_text(hjust=0,vjust=0) +
  geom_label_repel(aes(label = PRS_ID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   alpha=0.5,
                   segment.color = 'grey50',
                   max.overlaps=Inf) +
  theme_bw() +
  xlab("OR of Top 10% Scoring Individuals") +
  ylab('Rank Correlation with Phenotype') +
  theme(plot.title = element_text(face='bold',size=17),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = c(.6,.12), 
        legend.direction = "horizontal",
        legend.title=element_text(size=15,hjust=-0.5)) +
  scale_colour_continuous(name=expression(paste(log[10],"(No. Variants)     "))) +
  ggtitle('B. Performance of MCH PRSs')

# Combining all plots
library(gridExtra)
plots_B_and_C <- gridExtra::grid.arrange(new_plot_B,new_plot_C,ncol=2)
all_three_plots <- gridExtra::grid.arrange(plot1, plots_B_and_C,nrow=2)

ggsave(all_three_plots,
       filename = paste0(dir_prefix,"figures/MainText_080523_Fig1_1e-6.jpg"),
       width = 12.5, height = 12.5,
       dpi = 300)
