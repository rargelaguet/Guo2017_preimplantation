library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/Guo_2017"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/Guo_2017/mofa"
} else {
  # io$basedir <- "/hps/nobackup/stegle/users/ricard/Guo_2017"
  # io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  # io$outdir <- "/homes/ricard/gastrulation/metaccrna/mofa/E7.5/out"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$annos_dir  <- paste0(io$basedir1, "/features/filt")

####################
## Define options ##
####################

opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  "prom_2000_2000",
  # "genebody",
  # "ESC_DHS",
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500"
)

opts$acc.annos <- c(
  "prom_200_200",
  # "genebody",
  # "ESC_DHS",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
)

opts$stage_lineage <- c(
  # "ESC",
  "TE",
  "ICM",
  "Morula",
  "16-cell"
  # "8-cell",
  # "4-cell",
  # "2-cell",
  # "Zygote"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature
opts$met_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5       # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature
opts$acc_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define which cells to use
tmp <- fread(io$sample.metadata)
opts$met_cells <- tmp %>% .[stage%in%opts$stage_lineage,id_met]
opts$acc_cells <- tmp %>% .[stage%in%opts$stage_lineage,id_acc]


###################
## Load metadata ##
###################

sample_metadata <- fread(io$sample.metadata) %>%
  .[id_met%in%opts$met_cells | id_acc %in% opts$acc_cells ]
