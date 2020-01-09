library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Guo_2017/settings.R")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/met/dimensionality_reduction")

####################
## Define options ##
####################


# Define which annotations to look at
opts$annos <- c(
  "prom_2000_2000" = "Promoters"
)

opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)


# Define which cells to  use
opts$cells <- sample_metadata %>% 
  .[stage%in%opts$stages,id_met] %>%
  as.character

# Filtering options
opts$min.CpGs <- 1          # minimum number of CpG sites per feature and cell
opts$min.coverage <- 0.20   # minimum coverage (fraction of cells with at least min.CpG measurements)
opts$nfeatures <- 2500     # number of features per view (filter based on variance)

# Output file
io$outfile = sprintf("%s/hdf5/model_%s_%s.hdf5",io$outdir,paste(names(opts$annos), collapse="_"), paste(opts$stages, collapse="_"))

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>% 
  .[id_met%in%opts$cells] %>%
  .[,c("id_met","stage")] %>%
  merge(fread(io$met.stats)[,c("id_met","mean")], by="id_met")

table(sample_metadata$stage)

