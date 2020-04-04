################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
} else {
  stop("Computer not recognised")
}

# Folders with the differential analysis results
io$diff.met <- paste0(io$basedir,"/met/differential/feature_level")

# Output directory
io$outdir <- paste0(io$basedir,"/met/boxplots")

####################
## Define options ##
####################

# Define stages
opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

# Define genomic contexts for methylation
opts$met.annos <- c(
  "prom_2000_2000" = "Promoters"
)

# Options for selecting differential hits
opts$min.fdr <- 0.10
opts$min.met.diff <- 25

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[stage%in%opts$stages]
opts$met_cells <- sample_metadata %>% .[,id_met]
opts$acc_cells <- sample_metadata %>% .[,id_acc]
