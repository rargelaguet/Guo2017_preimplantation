
################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/Guo_2017"
} else {
  stop()
}

# Sample metadata
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

# Genomic contects
io$annos_dir <- paste0(io$basedir,"/features/filt")

# DNA methylation and chromatin accessibility data
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/stats/stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/stats.txt")

# Output directory
io$outdir <- "/Users/ricard/data/Guo_2017/metacc/profiles_enhancers"

####################
## Define options ##
####################

opts <- list()

# Define stages
opts$stage <- c(
  "2-cell",
  "4-cell",
  "8-cell",
  "16-cell",
  "Morula",
  "ICM"
)

# Define genomic contexts for methylation
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers"
)

# Define colors
opts$colors <- c(
  "2-cell" = "grey70",
  "4-cell" = "grey70",
  "8-cell" = "grey70",
  "16-cell" = "grey70",
  "Morula" = "grey70",
  "ICM" = "grey70"
)


# Output figure settings
opts$width = 10
opts$height = 6

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% .[stage%in%opts$stage] 
opts$met_cells <- tmp %>% .[,id_met]
opts$acc_cells <- tmp %>% .[,id_acc]
