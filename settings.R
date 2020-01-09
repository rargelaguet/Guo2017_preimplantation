
#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/Guo_2017"
} else if (grepl("yoda",Sys.info()['nodename'])) {
  # io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/Guo_2017"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$gene_metadata <- paste0(io$basedir,"/features/genes/Mmusculus_genes_BioMart.87.txt")


io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met.stats <- paste0(io$basedir,"/met/stats/sample_stats.txt")

io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc.stats <- paste0(io$basedir,"/acc/stats/sample_stats.txt")

io$features.dir <- paste0(io$basedir,"/features/filt")

#############
## Options ##
#############

opts <- list()

opts$stages <- c(
  "Zygote",
  "2-cell",
  "4-cell",
  "8-cell",
  "Morula",
  "ICM",
  "TE"
)

##########################
## Load sample metadata ##
##########################

factor.cols <- c("sample","id_met","id_acc","stage")

sample_metadata <- fread(io$metadata) %>% 
  .[stage%in%opts$stages] %>%
  .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>%
  .[,stage:=factor(stage,levels=opts$stages)] %>%
  droplevels
  
  