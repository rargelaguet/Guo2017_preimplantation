################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo2017_preimplantation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/Guo2017_preimplantation/settings.R")
} else {
  stop("Computer not recognised")
}

# Output directory
io$outdir <- paste0(io$basedir,"/metacc/stats")

####################
## Define options ##
####################

# Define lineages
# opts$lineages <- c(
#   # "Zygote", 
#   # "2cell", 
#   # "4cell", 
#   "8cell", 
#   "Morula", 
#   "ICM", 
#   "TE"
#   # "hESC"
# )

# Define genomic contexts for methylation
opts$annos <- NULL
# opts$annos <- c("CGI","genebody","prom_2000_2000_noncgi","L1")
if (is.null(opts$annos)) {
  # opts$met.annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)# %>% head(n=2)
}

#####################
## Update metadata ##
#####################

# sample_metadata <- sample_metadata %>% 
#   .[lineage%in%opts$lineage] %>%
#   # .[sex%in%c("Female","Male")] %>% 
#   droplevels 
