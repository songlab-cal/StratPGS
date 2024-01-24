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

## Plotting ------------------------------------------------------------
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
  xlab(expression(paste(log[10],'(No. Variants Included in PGS)'))) +
  ylab('Count') +
  ylim(c(0,600)) +
  geom_vline(xintercept=log10(median(PRS_metadata_annotated_df$NO_NONZERO_VARIANTS)),
             lty='dashed',
             colour="#0570b0") +
  annotate(geom="text", 
           x=1.6, y=500, 
           label="Median No. Variants\n= 5,578",
           size=5,
           color = "#0570b0") +
  #ggtitle('A. Dimensionality of 3,683 PRSs from the PGS Catalogue') +
  ggtitle('Dimensionality of 3,683 PGSs\nfrom the PGS Catalogue') +
  theme(plot.title = element_text(face='bold',size=17,hjust=0.5),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        text=element_text(family='DejaVu Sans'))

# Optional: Save
ggsave(plot1,
       file = "PRS_variant_counts_logplot_100423.jpg",
       width = 5, height = 5,
       #width = 6.5, height = 4,
       dpi = 1000)

