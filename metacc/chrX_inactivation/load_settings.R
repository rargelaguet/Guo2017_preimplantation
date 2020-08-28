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

sample_metadata[sex%in%c("Male","Female"),length(unique(embryo)),by=c("sex","stage")]

# Output directory
io$outdir <- paste0(io$basedir,"/metacc/chrX_inactivation")

####################
## Define options ##
####################

opts$stages <- c(
  # "Zygote",
  # "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

# Define genomic contexts for methylation
# opts$annos <- NULL
opts$annos <- c("CGI","genebody","L1")
if (is.null(opts$annos)) {
  # opts$met.annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)
}

#####################
## Update metadata ##
#####################

sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stage] %>%
  .[sex%in%c("Female","Male")] %>% 
  droplevels 

table(sample_metadata$stage)
table(sample_metadata$sex)
