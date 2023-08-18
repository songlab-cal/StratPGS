## Mining PGS Catalogue metadata -----------------------------------------------
# Date: 8/4/23
# First is PGS000001
# Last one is PGS003759
# This version only considers the number of nonzero variant effects
library(jsonlite, 
        lib.loc = "/global/home/groups/consultsw/sl-7.x86_64/modules/r-packages/2022-10-14-r4.2")

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
                        "logs/mine_PGS_Catalogue_annotated_",
                        format(Sys.Date(), "%m%d%y.log")))
lf <- log_open(tmp)

# PRS_metadata_df <- data.frame(ID=character(),
#                               PHENOTYPE=character(),
#                               SAMPLE_SIZE=numeric(),
#                               METHODOLOGY=character(),
#                               PARAMETERS=character(),
#                               NO_VARIANTS=numeric(),
#                               FTP_WEBPAGE=character())
# for (i in 1:76) {
#   log_print(paste0(date(), ": Scrapping next 50 PRSs starting from ", 50*(i-1)+1))
#   scrap_input <- getScrapInput(i)
#   webpage <- read_json(scrap_input)
#   ids <- unlist2(lapply(webpage$results, function(x){return(x[['id']])}))
#   sampsize <- unlist2(lapply(webpage$results, function(x) {return(getSampSize(x))}))
#   method_name <- unlist2(lapply(webpage$results, function(x){return(x[['method_name']])}))
#   no_variants <- unlist2(lapply(webpage$results, function(x){return(x[['variants_number']])}))
#   params <- unlist2(lapply(webpage$results, function(x){return(x[['method_params']])}))
#   trait_reported <- unlist2(lapply(webpage$results, function(x){return(x[['trait_reported']])}))
#   ftp_webpage <- unlist2(lapply(webpage$results, function(x){return(x[['ftp_scoring_file']])}))
#   PRS_metadata_df <- rbind(PRS_metadata_df,
#                            data.frame(ID=ids,
#                                       PHENOTYPE=trait_reported,
#                                       SAMPLE_SIZE=sampsize,
#                                       METHODOLOGY=method_name,
#                                       PARAMETERS=params,
#                                       NO_VARIANTS=no_variants,
#                                       FTP_WEBPAGE=ftp_webpage))
# }

PRS_metadata_df <- readr::read_csv(paste0(out_dir,
                                          "PGS_catalogue_metadata.csv"))
# Add nonzero variant count 
no_nonzero_vars <- c()
for (j in 1:nrow(PRS_metadata_df)) {
  log_print(paste0(date(), ": Counting non-zero variants for ", 
                   j,"th Polygenic Score: ", 
                   PRS_metadata_df$ID[j]))
  prs_file <- data.table::fread(PRS_metadata_df$FTP_WEBPAGE[j],
                                skip=14)
  no_nonzero_vars <- c(no_nonzero_vars, 
                       sum(prs_file$effect_weight!=0))
  log_print(paste0(date(), ": No. non-zero variants for ", 
                   j,"th Polygenic Score = ", 
                   sum(prs_file$effect_weight!=0)))
  rm(prs_file)
}
PRS_metadata_df$NO_NONZERO_VARIANTS <- no_nonzero_vars

# Save
readr::write_csv(PRS_metadata_df, 
                 file=paste0(out_dir,"PGS_catalogue_metadata_annotated.csv"))

# Close stream and write to log file
log_close()
writeLines(readLines(lf))

## Analyze ---------------------------------------------------------------------
PRS_metadata_annotated_df <- readr::read_csv(paste0(out_dir,
                                                    "PGS_catalogue_metadata_annotated.csv"))
PRS_metadata_df <- readr::read_csv(paste0(out_dir,
                                          "PGS_catalogue_metadata.csv"))
PRS_metadata_annotated_df$NONZERO_FRACTION <- PRS_metadata_annotated_df$NO_NONZERO_VARIANTS/
  PRS_metadata_annotated_df$NO_VARIANTS
summary(PRS_metadata_annotated_df$NONZERO_FRACTION)

# Fix the 0 fraction rows 
# (PRS file has one fewer metadata row hence skipped column names when counting)
PRS_metadata_annotated_df %>% subset(NO_NONZERO_VARIANTS == 0)
which(PRS_metadata_annotated_df$NO_NONZERO_VARIANTS ==  0)
# Manual checking shows that both files have nonzero fraction = 1
# https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002295/ScoringFiles/
# prs_file <- data.table::fread(PRS_metadata_df$FTP_WEBPAGE[1707],skip=13) 
# sum(prs_file$effect_weight!=0) # 21
# prs_file <- data.table::fread(PRS_metadata_df$FTP_WEBPAGE[2223],skip=13)
# sum(prs_file$effect_weight!=0) # 413
# https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS001779/ScoringFiles/ 
PRS_metadata_annotated_df$NO_NONZERO_VARIANTS[c(1707,2223)] <- c(21,413)
PRS_metadata_annotated_df$NONZERO_FRACTION <- PRS_metadata_annotated_df$NO_NONZERO_VARIANTS/
  PRS_metadata_annotated_df$NO_VARIANTS
summary(PRS_metadata_annotated_df$NONZERO_FRACTION)
summary(PRS_metadata_annotated_df$NO_NONZERO_VARIANTS) 
summary(PRS_metadata_annotated_df$NO_VARIANTS) 

readr::write_csv(PRS_metadata_annotated_df, 
                 file=paste0(out_dir,"PGS_catalogue_metadata_annotated.csv"))
