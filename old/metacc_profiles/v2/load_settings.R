
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
io$met.dir <- paste0(io$basedir,"/met/cpg_level")
io$acc.dir <- paste0(io$basedir,"/acc/gpc_level")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/stats/stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/stats.txt")

# Pluripotency vs Midbrain H3K27ac data
io$esc_vs_brain <- paste0(io$basedir,"/H3K27ac/E7.5_enhancers_Ect_ESC_brain_top250.txt")

# Output directory
io$pdfdir <- paste0(io$basedir,"/metacc/profiles_enhancers")


####################
## Define options ##
####################

opts <- list()

# Define stages
opts$stage <- c(
  # "2-cell",
  # "4-cell",
  "8-cell",
  "16-cell",
  "Morula",
  "ICM"
)

# Define genomic contexts
opts$annos <- c(
  # "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)

# Define window positions and characteristics
opts$positions <- c(
  # "H3K27ac_distal_E7.5_Mes_intersect12"="center",
  "H3K27ac_distal_E7.5_End_intersect12"="center",
  "H3K27ac_distal_E7.5_Ect_intersect12"="center"
)
opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150


# Define which cells to use
tmp <- fread(io$sample.metadata) %>% .[stage%in%opts$stage] 
opts$met.cells <- tmp[,id_met]
opts$acc.cells <- tmp[,id_acc]
